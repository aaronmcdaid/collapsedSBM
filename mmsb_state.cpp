#include "aaron_utils.hpp"
#include "mmsb_state.hpp"
namespace sbm {
	int & PairCounts:: set(const int i, const int j) {
			return (*this)[std::make_pair(i,j)];
	}
	int   PairCounts:: get(const int i, const int j) const {
			if (this->count(std::make_pair(i,j))==0)
				return 0;
			return this-> at(std::make_pair(i,j));
	}
	MMSBstate:: MMSBstate (const GraphType * const g): _N(g->numNodes()), _k(1) {
		for(int i=0; i< this->_N; i++) {
			ls.push_back(new Labelling(_N));
		}
		for(int i=0; i< this->_N; i++) {
			assert((int)ls.at(i)->clusters.size()==this->_k);
		}
		const int totalNumberOfPairs = this->_N * (this->_N-1) / 2;
		this->numPairs.set(0,0) = totalNumberOfPairs;
		PP(this->numPairs.get(0,0));
		this->numEdges.set(0,0) = g->numRels();
		PP(this->numEdges.get(0,0));
	};
	void MMSBstate:: appendEmptyCluster() {
		for(int i=0; i< this->_N; i++) {
			Labelling * l = ls.at(i);
			assert(l);
			assert((int)l->clusters.size() == this->_k);
			l->appendEmptyCluster();
		}
		this->_k ++;
	}
	long double MMSBstate:: P_z_K() const { // 1 and 2
		// const long double priorOnK = -this->_k; // Exponential prior on K
		const long double priorOnK = -LOG2FACT(this->_k); // Poisson(1) prior on K
		const long double K_dependant_bits = priorOnK + LOG2GAMMA(this->_k) - LOG2GAMMA(this->_k + this->_N);
		return sbm:: assertNonPositiveFinite(K_dependant_bits);
	}
	void MMSBstate:: P_zs_given_K() const {
		for(int i=0; i< this->_N; i++) {
			Labelling * l = ls.at(i);
			assert(l);
			assert((int)l->clusters.size() == this->_k);
			for(int j=0; j< this->_N; j++) {
				const int clid = l->cluster_id.at(j);
				cout << ' ' << clid;
			}
			cout << '\t' << this->P_z_K() + l->SumOfLog2LOrders;
			cout << endl;
		}
	}
} // namespace sbm
