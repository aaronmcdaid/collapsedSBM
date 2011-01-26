#include "sbm_state.hpp"
#include "aaron_utils.hpp"
using namespace shmGraphRaw;
namespace sbm {
	State::State(const GraphType * const g) : _g(g), _N(g->numNodes()) {
		// initialize it with every node in one giant cluster
		this->_k = 1;
		this->clusters.reserve(100 + 2*_N); // this is important. We mustn't allow the clusters vector to be moving about in RAM.
		this->clusters.push_back(Cluster());
		assert(this->clusters.size()==1);

		assert(this->cluster_id.size()==0);
		assert(this->its.size()==0);
		for(int i=0; i<this->_N; i++) {
			this->cluster_id.push_back(0);
			Cluster &cl = this->clusters.back();
			this->its.push_back( cl.newMember(i) );
			assert(*this->its.at(i) == i);
		}

		assert((int)this->cluster_id.size()==this->_N);
		assert((int)this->its.size()==this->_N);
		assert((int)this->clusters.back().members.size()==this->_N);

		// inform EdgeCounts of all the edges
		for(int relId = 0; relId < this->_g->numRels(); relId++) {
			const std::pair<int, int> & eps = this->_g->EndPoints(relId);
			const int cl1 = this->cluster_id.at(eps.first);
			const int cl2 = this->cluster_id.at(eps.second);
			_edgeCounts.inform(cl1,cl2);
		}
	}
	const int State::Cluster::order() const {
		return this->members.size();
	}
	std::list<int>::iterator State::Cluster::newMember(const int n) {
		this->members.push_front(n);
		return this->members.begin();
	}

	int State::appendEmptyCluster() {
		const int newClusterID = this->_k;
		this->_k ++;
		assert(this->_k < (int)this->clusters.capacity());
		this->clusters.push_back(Cluster());
		assert((int)this->clusters.back().members.size()==0);
		return newClusterID;
	}
	int State::isolateNode(const int n) { // create a new (probably temporary) cluster to hold this one node
		assert(n>=0 && n<this->_N);
		const int newClusterID = this->appendEmptyCluster();
		const int oldClusterID = this->cluster_id.at(n);
		const int oldClusterSize = this->clusters.at(oldClusterID).order();
		const list<int>::iterator it = this->its.at(n);
		assert(*it == n);
		this->clusters.at(oldClusterID).members.erase(it);
		assert(oldClusterSize-1 == this->clusters.at(oldClusterID).order());

		this->cluster_id.at(n) = newClusterID;
		Cluster &cl = this->clusters.at(newClusterID);
		assert(cl.order()==0);
		const list<int>::iterator newit = cl.newMember(n);
		this->its.at(n) = newit;

		this->informNodeMove(n, oldClusterID, newClusterID);

		return newClusterID;
	}

	void State::internalCheck() const {
		assert(this->_k>0);
		assert(this->_k == (int)this->clusters.size());
		assert(this->_N == (int)this->cluster_id.size());
		assert(this->_N == (int)this->its.size());
		boost::unordered_set<int> alreadyConsidered;
		for(int CL = 0; CL < this->_k; CL++) {
			const Cluster &cl = this->clusters.at(CL);
			for(list<int>::const_iterator i = cl.members.begin(); i!=cl.members.end(); i++) {
				const int n = *i;
				assert(n>=0 && n<this->_N);
				bool wasAccepted = alreadyConsidered.insert(n).second;
				assert(wasAccepted);
				assert(n == *this->its.at(n));
				assert(CL == this->cluster_id.at(n));
				assert(i == this->its.at(n));
			}
		}
		assert((int)alreadyConsidered.size() == this->_N);

		typedef boost::unordered_map< int , boost::unordered_map<int,int> >::value_type outer_value_type;
		typedef                             boost::unordered_map<int,int>  ::value_type inner_value_type;
		int numEdgeEnds = 0; // should end up at twice g->numRels() // assuming no self-loops of course
		forEach(const outer_value_type & outer, amd::mk_range(this->_edgeCounts.counts)) {
			assert(outer.first >= 0 && outer.first < this->_k);
			forEach(const inner_value_type & inner, amd::mk_range(outer.second)) {
				cout << outer.first << ' ' << inner.first << '\t' << inner.second << endl;
				assert(inner.first >= 0 && inner.first < this->_k);
				assert(inner.second >= 0); // maybe > 0 in future ? TODO SPEED ?
				numEdgeEnds += inner.second;
				if(outer.first == inner.first) // we should double count those inside a cluster during this verification process
					numEdgeEnds += inner.second;
			}
		}
		assert(numEdgeEnds == 2*this->_g->numRels());
	}

	void State::shortSummary() const {
		cout << " == Summary: ==" << endl;
		PP(this->_k);
		for(int n=0; n<this->_N;n++) {
			const int id_of_cluster = this->cluster_id.at(n);
			if(id_of_cluster<10)
				cout << (char)('0' + id_of_cluster);
			else
				cout << (char)('a' + id_of_cluster - 10);
		}
		cout << endl;
	}


	void State::EdgeCounts::inform(const int cl1, const int cl2) { // inform us of an edge between cl1 and cl2
		assert(cl1 >= 0); // && cl1 < this->_k);
		assert(cl2 >= 0); // && cl2 < this->_k);
		this->counts[cl1][cl2]++;
		if(cl1 != cl2)
			this->counts[cl2][cl1]++;
	}
	void State::EdgeCounts::uninform(const int cl1, const int cl2) { // inform us of an edge between cl1 and cl2
		assert(cl1 >= 0); // && cl1 < this->_k);
		assert(cl2 >= 0); // && cl2 < this->_k);
		this->counts[cl1][cl2]--;
		if(cl1 != cl2)
			this->counts[cl2][cl1]--;
	}
	void State::informNodeMove(const int n, const int oldcl, const int newcl) { // a node has just moved from one cluster to another. We must consider it's neighbours for _edgeCounts
		const PlainMem::mmap_uset_of_ints & rels = this->_g->myRels(n);
		forEach(const int relId, amd::mk_range(rels)) {
			const int otherEnd = this->_g->oppositeEndPoint(relId, n);
			const int otherCl = this->cluster_id.at(otherEnd);
			this->_edgeCounts.uninform(otherCl, oldcl);
			this->_edgeCounts.  inform(otherCl, newcl);
		}
	}
} // namespace sbm
