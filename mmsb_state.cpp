#include "aaron_utils.hpp"
#include "mmsb_state.hpp"
namespace sbm {
	int & PairCounts:: set(const int i, const int j) {
		return this-> operator [] (std::make_pair(i,j));
	}
	int   PairCounts:: get(const int i, const int j) const {
			if (this->count(std::make_pair(i,j))==0)
				return 0;
			return this-> at(std::make_pair(i,j));
	}
	MMSBstate:: MMSBstate (const GraphType * const g): _g(g), _N(g->numNodes()), _k(1) {
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
	long double P_beta_binomial(const int n, const int k) {
		assert(k >= 0);
		assert(n >= k);
		if(n==0)
			return 0.0L;
		long double x = - log2l(n+1);
		PP(x);
		x -= M_LOG2E * gsl_sf_lnchoose(n, k);
		PP(x);
		return sbm:: assertNonPositiveFinite(x);
	}
	void MMSBstate:: moveOnePair(const int w,const int v,const int clid) { // node w, when interacting with v, should take identity clid
		assert(w!=v);
		assert(clid >= 0 && clid < this->_k);
		const bool areConnected = this->_g->are_connected(w,v);
		PP2(w,v);
		PP(areConnected);
		Labelling * l = ls.at(w);
		assert(l);
		const int old_clid = l->cluster_id.at(v);
		assert(clid != old_clid);
		PP2(w,v);
		// PP2(old_clid, clid);
		// moving the node is easy enough. It's updating the pair- and edge- counts that's interesting
		
		// Let's check on the other side of this pair
		Labelling * l2 = ls.at(v);
		assert(l2);
		const int fixed_id = l2->cluster_id.at(w);
		PP2(fixed_id, old_clid);
		const int OrigWasPairs = this->numPairs.get(fixed_id, old_clid);
		const int OrigWasEdges = this->numEdges.get(fixed_id, old_clid);
		PP2(OrigWasPairs, OrigWasEdges);
		const long double origWasBits = P_beta_binomial(OrigWasPairs, OrigWasEdges);
		const int OrigWillPairs = OrigWasPairs - 1;
		const int OrigWillEdges = OrigWasEdges - (areConnected ? 1 : 0);
		PP2(OrigWillPairs, OrigWillEdges);
		const long double origWillBits = P_beta_binomial(OrigWillPairs, OrigWillEdges);
		PP (origWillBits - origWasBits);

		PP2(fixed_id, clid);
		const int NewWasPairs = this->numPairs.get(fixed_id, clid);
		const int NewWasEdges = this->numEdges.get(fixed_id, clid);
		PP2(NewWasPairs, NewWasEdges);
		const long double newWasBits = P_beta_binomial(NewWasPairs, NewWasEdges);
		const int NewWillPairs = NewWasPairs + 1;
		const int NewWillEdges = NewWasEdges + (areConnected ? 1 : 0);
		PP2(NewWillPairs, NewWillEdges);
		const long double newWillBits = P_beta_binomial(NewWillPairs, NewWillEdges);
		PP(newWillBits - newWasBits);
		const long double delta_g_given_z = origWillBits - origWasBits + newWillBits - newWasBits;
		assert(isfinite(delta_g_given_z));
		PP(delta_g_given_z);

		if(areConnected) {
			assert(delta_g_given_z > 0.0L); // this isn't a serious assert, but it does work in this initialization
		} else {
			assert(delta_g_given_z < 0.0L);
		}

	}
} // namespace sbm
