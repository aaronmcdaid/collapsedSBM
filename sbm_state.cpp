#include <gsl/gsl_sf.h>
#include "sbm_state.hpp"
#include "aaron_utils.hpp"
using namespace shmGraphRaw;
namespace sbm {
	struct SelfLoopsNotSupported : public std::exception {
	};
	State::State(const GraphType * const g) : _g(g), _N(g->numNodes()) {
		// initialize it with every node in one giant cluster
		this->_k = 1;
		this->clusters.push_back(new Cluster());
		assert(this->clusters.size()==1);

		assert(this->cluster_id.size()==0);
		assert(this->its.size()==0);
		for(int i=0; i<this->_N; i++) {
			this->cluster_id.push_back(0);
			Cluster *cl = this->clusters.back();
			assert(cl);
			this->its.push_back( cl->newMember(i) );
			assert(*this->its.at(i) == i);
		}

		assert((int)this->cluster_id.size()==this->_N);
		assert((int)this->its.size()==this->_N);
		assert((int)this->clusters.back()->members.size()==this->_N);

		// inform EdgeCounts of all the edges
		for(int relId = 0; relId < this->_g->numRels(); relId++) {
			const std::pair<int, int> & eps = this->_g->EndPoints(relId);
			if(eps.first == eps.second)
				throw SelfLoopsNotSupported();
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

	int State:: appendEmptyCluster() {
		const int newClusterID = this->_k;
		this->_k ++;
		Cluster * newCluster = new Cluster();
		this->clusters.push_back(newCluster);
		assert(newCluster->members.size()==0);
		return newClusterID;
	}
	void State:: deleteClusterFromTheEnd() {
		assert(this->_k >= 1);
		this->_k --;
		Cluster * clusterToDelete = this->clusters.at(this->_k);
		assert(clusterToDelete->order() == 0);
		this->clusters.pop_back();
		delete clusterToDelete;
	}
	void State::moveNode(const int n, const int newClusterID) {
		const int oldClusterID = this->cluster_id.at(n);
		const int oldClusterSize = this->clusters.at(oldClusterID)->order();
		assert(newClusterID >= 0 && newClusterID < this->_k);
		assert(oldClusterID >= 0 && oldClusterID < this->_k);
		assert(newClusterID != oldClusterID);

		Cluster *cl = this->clusters.at(newClusterID);
		assert(cl);

		const list<int>::iterator it = this->its.at(n);
		assert(*it == n);
		this->clusters.at(oldClusterID)->members.erase(it); // fix up this->clusters.members
		const list<int>::iterator newit = cl->newMember(n); // fix up this->clusters.members
		this->its.at(n) = newit; // fix up this->its

		this->cluster_id.at(n) = newClusterID; // fix up this->cluster_id

		assert(oldClusterSize-1 == this->clusters.at(oldClusterID)->order());
	}
	int State::isolateNode(const int n) { // create a new (probably temporary) cluster to hold this one node
		assert(n>=0 && n<this->_N);
		const int newClusterID = this->appendEmptyCluster();
		const int oldClusterID = this->cluster_id.at(n);

		Cluster *cl = this->clusters.at(newClusterID);
		assert(cl->order()==0);

		this->moveNode(n, newClusterID);

		this->informNodeMove(n, oldClusterID, newClusterID);

		return newClusterID;
	}
	void State::unIsolateTempNode(const int n, const int newClusterID) { // move a node from its 'temporary' cluster to an existing cluster
		const int oldClusterID = this->cluster_id.at(n);
		this->moveNode(n, newClusterID);
		assert(oldClusterID+1 == this->_k);
		this->deleteClusterFromTheEnd();

		this->informNodeMove(n, oldClusterID, newClusterID);
	}

	void State:: summarizeEdgeCounts() const {
		forEach(const EdgeCounts::outer_value_type & outer, amd::mk_range(this->_edgeCounts.counts)) {
			assert(outer.first >= 0 && outer.first < this->_k);
			forEach(const EdgeCounts::inner_value_type & inner, amd::mk_range(outer.second)) {
				cout << " edges " << outer.first << ',' << inner.first << '\t' << inner.second << endl;
				assert(inner.first >= 0 && inner.first < this->_k);
				assert(inner.second >= 0); // maybe > 0 in future ? TODO SPEED ?
			}
		}
	}
	void State:: blockDetail() const {
			cout << "    | ";
			for(int j=0; j<this->_k; j++) {
				const Cluster *J = this->clusters.at(j);
				assert(J);
				const int nj = J->order();
				cout << printfstring("      %5d     ", nj) << "   ";
			}
			cout << endl;
		for(int i=0; i<this->_k; i++) {
			const Cluster *I = this->clusters.at(i);
			assert(I);
			const int ni = I->order();
			cout << printfstring("%3d", ni) << " | ";
			for(int j=0; j<this->_k; j++) {
				const Cluster *J = this->clusters.at(j);
				assert(J);
				const int nj = J->order();
				const int edges = this->_edgeCounts.get(i,j);
				const int pairs = i==j ? (ni * (nj-1) / 2) : (ni*nj);
				assert(edges <= pairs);
				cout << printfstring("%10s %-#5.2f", printfstring("%d/%d", edges, pairs).c_str(), double(edges)/double(pairs)) << " | ";
			}
			cout << endl;
		}
	}

	void State::internalCheck() const {
		assert(this->_k>0);
		assert(this->_k == (int)this->clusters.size());
		assert(this->_N == (int)this->cluster_id.size());
		assert(this->_N == (int)this->its.size());
		boost::unordered_set<int> alreadyConsidered;
		for(int CL = 0; CL < this->_k; CL++) {
			const Cluster *cl = this->clusters.at(CL);
			assert(cl);
			for(list<int>::const_iterator i = cl->members.begin(); i!=cl->members.end(); i++) {
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

		EdgeCounts edgeCountsVerification;
		for(int relId = 0; relId < this->_g->numRels(); relId++) {
			const std::pair<int, int> & eps = this->_g->EndPoints(relId);
			const int cl1 = this->cluster_id.at(eps.first);
			const int cl2 = this->cluster_id.at(eps.second);
			edgeCountsVerification.inform(cl1,cl2);
		}
		DYINGWORDS(edgeCountsVerification.counts.size() == this->_edgeCounts.counts.size()) {
			PP(edgeCountsVerification.counts.size());
			PP(this->_edgeCounts.counts.size());
		}
		assert(edgeCountsVerification.counts.size() == this->_edgeCounts.counts.size());

		int numEdgeEnds = 0; // should end up at twice g->numRels() // assuming no self-loops of course
		forEach(const EdgeCounts::outer_value_type & outer, amd::mk_range(this->_edgeCounts.counts)) {
			assert(outer.first >= 0 && outer.first < this->_k);
			const EdgeCounts::outer_value_type::second_type & outerVerification = edgeCountsVerification.counts.at(outer.first);
			assert(outerVerification.size() == outer.second.size());
			forEach(const EdgeCounts::inner_value_type & inner, amd::mk_range(outer.second)) {
				assert(inner.first >= 0 && inner.first < this->_k);
				assert(inner.second > 0); // maybe > 0 in future ? TODO SPEED ?
				assert(inner.second == outerVerification.at(inner.first));
				numEdgeEnds += inner.second;
				if(outer.first == inner.first) // we should double count those inside a cluster during this verification process
					numEdgeEnds += inner.second;
			}
		}
		assert(numEdgeEnds == 2*this->_g->numRels());


	}

	void State::shortSummary() const {
		cout << endl << " == Summary: ==" << endl;
		PP(this->_k);
		PP(this->pmf());
		for(int n=0; n<this->_N;n++) {
			const int id_of_cluster = this->cluster_id.at(n);
			if(id_of_cluster<10)
				cout << (char)('0' + id_of_cluster);
			else if(id_of_cluster < 36)
				cout << (char)('a' + id_of_cluster - 10);
			else
				cout << '<' << id_of_cluster << '>';
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
	void State::EdgeCounts::partialUnInform(const int clA, const int clB) { // private to EdgeCounts
		EdgeCounts::outer_value_type::second_type & edgeMapForOneCluster = this->counts.at(clA);
		EdgeCounts::outer_value_type::second_type::iterator aSingleEntry = edgeMapForOneCluster.find(clB);
		assert(edgeMapForOneCluster.at(clB) == aSingleEntry->second);
		aSingleEntry->second --;
		assert(edgeMapForOneCluster.at(clB) == aSingleEntry->second);
		if(aSingleEntry->second == 0)
			edgeMapForOneCluster.erase(aSingleEntry);
		if(edgeMapForOneCluster.size()==0)
			this->counts.erase(clA);
	}
	void State::EdgeCounts::uninform(const int cl1, const int cl2) { // inform us of an edge between cl1 and cl2
		assert(cl1 >= 0); // && cl1 < this->_k);
		assert(cl2 >= 0); // && cl2 < this->_k);
		this->partialUnInform(cl1,cl2);
		if(cl1 != cl2)
			this->partialUnInform(cl2,cl1);
	}
	int  State::EdgeCounts:: get(const int cl1, const int cl2) const throw() {
		try {
			return this->counts.at(cl1).at(cl2);
		} catch (const std::out_of_range &e) {
			return 0;
		}
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
	long double assertNonPositiveFinite(const long double x) {
		assert(isfinite(x));
		assert(x<=0.0L);
		return x;
	}

#define LOG2GAMMA(x) (M_LOG2E * gsl_sf_lngamma(x))
#define LOG2FACT(x)  (M_LOG2E * gsl_sf_lnfact(x))
	long double State:: P_z() const { // given our current this->_k, what's P(z | k)
		// const long double priorOnK = -this->_k; // Exponential prior on K
		const long double priorOnK = -LOG2FACT(this->_k); // Poisson(1) prior on K
		const long double K_dependant_bits = priorOnK + LOG2GAMMA(this->_k) - LOG2GAMMA(this->_k + this->_N);
		long double perCluster_bits = 0.0L;
		for(int CL=0; CL < this->_k; CL++) {
			const Cluster *cl = this->clusters.at(CL);
			assert(cl);
			perCluster_bits += LOG2FACT(cl->order());
		}
		// PP(K_dependant_bits + perCluster_bits);
		return assertNonPositiveFinite(K_dependant_bits + perCluster_bits);
	}
	long double State:: P_edges_given_z_slow() const {
		long double edges_bits = 0.0L;
		for(int i=0; i<this->_k; i++) {
			const Cluster *I = this->clusters.at(i);
			assert(I);
			for(int j=0; j<=i; j++) {
				const Cluster *J = this->clusters.at(j);
				assert(J);
				const int ni = I->order();
				const int nj = J->order();
				const int edges = this->_edgeCounts.get(i,j);
				/*
				cout << ni << " x " << nj;
				if(i==j)
					cout << "-1 x1/2";
				cout	<< " pairs. "
					<< edges << " edges."
					<< endl;
					*/
				const int pairs = i==j ? (ni * (nj-1) / 2) : (ni*nj);
				assert(edges <= pairs);
				// PP2(pairs,edges);
				if(pairs > 0) {
					edges_bits -= log2(pairs);
					edges_bits -= M_LOG2E * gsl_sf_lnchoose(pairs, edges);
				}
				// PP(edges_bits);
			}
		}
		// PP(edges_bits);
		return assertNonPositiveFinite(edges_bits);
	}
	long double State:: P_edges_given_z() const { // this function might be used to try faster ways to calculate the same data. Especially where there are lots of clusters in a small graph.
		const long double slow = this->P_edges_given_z_slow();
		// TODO: Use a faster implementation
		return slow;
	}
	long double State:: pmf_slow() const {
		return this->P_edges_given_z_slow() + this->P_z();
	}
	long double State:: pmf() const {
		return this->P_edges_given_z() + this->P_z();
	}
} // namespace sbm
