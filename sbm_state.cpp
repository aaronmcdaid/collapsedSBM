#include "sbm_state.hpp"
#include "aaron_utils.hpp"
using namespace shmGraphRaw;
namespace sbm {
	State::State(const ReadableShmGraphBase * const g) : _g(g), _N(g->numNodes()) {
		// initialize it with every node in one giant cluster
		this->_k = 1;
		this->clusters.reserve(100 + 2*_N); // this is important. We mustn't allow the clusters vector to be moving about in RAM.
		this->clusters.push_back(Cluster());
		assert(this->clusters.size()==1);

		assert(this->cluster_id.size()==0);
		assert(this->its.size()==0);
		for(int i=0; i<this->_N; i++) {
			this->cluster_id.push_back(0);
			this->clusters.back().members.push_front(i);
			this->its.push_back( this->clusters.back().members.begin() );
			assert(*this->its.at(i) == i);
		}

		assert((int)this->cluster_id.size()==this->_N);
		assert((int)this->its.size()==this->_N);
		assert((int)this->clusters.back().members.size()==this->_N);
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
	}
} // namespace sbm
