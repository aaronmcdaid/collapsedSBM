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

		for(int i=0; i<this->_N; i++) {
			const bool wasInserted = this->nodeNames.insert(atoi(this->_g->NodeAsString(i))).second;
			assert(wasInserted);
		}
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
		forEach(int nodeName, amd::mk_range(this->nodeNames))
		// for(int i=0; i< this->_N; i++)
		{
			const int i = this->_g->StringToNodeId(printfstring("%d", nodeName).c_str());
			cout << this->_g->NodeAsString(i) << '\t';
			Labelling * l = ls.at(i);
			assert(l);
			assert((int)l->clusters.size() == this->_k);
			for(int j=0; j< this->_N; j++) {
				const int clid = l->cluster_id.at(j);
				cout << ' ' << clid;
			}
			cout << '\t' << this->P_z_K() + l->SumOfLog2Facts;
			const int sz0 = l->clusters.at(0)->order();
			cout << '\t';
			if(sz0*2 > this->_N)
				cout << "1 ";
			else
				cout << " 2";
			cout << "\t" << this->_g->degree(i);
			cout << "\t(";
			for(int k=0; k<this->_k; k++)
				cout << ',' << l->clusters.at(k)->order();
			cout << ")";
			cout << endl;
		}
	}
	long double P_beta_binomial(const int n, const int k) {
		assert(k >= 0);
		assert(n >= k);
		if(n==0)
			return 0.0L;
		long double x = - log2l(n+1);
		// PP(x);
		x -= M_LOG2E * gsl_sf_lnchoose(n, k);
		// PP(x);
		return sbm:: assertNonPositiveFinite(x);
	}
	void MMSBstate:: moveOnePair(const int w,const int v,const int clid) { // node w, when interacting with v, should take identity clid
		const bool verbose = false;
		assert(w!=v);
		assert(clid >= 0 && clid < this->_k);
		Labelling * l = ls.at(w);
		const int old_clid = l->cluster_id.at(v);
		if(clid == old_clid)
			return;
		if(verbose) cout << endl << "   moveOnePair()" << endl;
		assert(clid != old_clid);
		const bool areConnected = this->_g->are_connected(w,v);
		if(verbose) PP2(w,v);
		if(verbose) PP(areConnected);
		assert(l);
		if(verbose) PP2(w,v);
		// PP2(old_clid, clid);
		// moving the node is easy enough. It's updating the pair- and edge- counts that's interesting
		
		// Let's check on the other side of this pair
		Labelling * l2 = ls.at(v);
		assert(l2);
		const int fixed_id = l2->cluster_id.at(w);
		if(verbose) PP2(fixed_id, old_clid);
		const int OrigWasPairs = this->numPairs.get(fixed_id, old_clid);
		const int OrigWasEdges = this->numEdges.get(fixed_id, old_clid);
		if(verbose) PP2(OrigWasPairs, OrigWasEdges);
		const long double origWasBits = P_beta_binomial(OrigWasPairs, OrigWasEdges);
		const int OrigWillPairs = OrigWasPairs - 1;
		const int OrigWillEdges = OrigWasEdges - (areConnected ? 1 : 0);
		if(verbose) PP2(OrigWillPairs, OrigWillEdges);
		const long double origWillBits = P_beta_binomial(OrigWillPairs, OrigWillEdges);
		if(verbose) PP (origWillBits - origWasBits);

		if(verbose) PP2(fixed_id, clid);
		const int NewWasPairs = this->numPairs.get(fixed_id, clid);
		const int NewWasEdges = this->numEdges.get(fixed_id, clid);
		if(verbose) PP2(NewWasPairs, NewWasEdges);
		const long double newWasBits = P_beta_binomial(NewWasPairs, NewWasEdges);
		const int NewWillPairs = NewWasPairs + 1;
		const int NewWillEdges = NewWasEdges + (areConnected ? 1 : 0);
		if(verbose) PP2(NewWillPairs, NewWillEdges);
		const long double newWillBits = P_beta_binomial(NewWillPairs, NewWillEdges);
		if(verbose) PP(newWillBits - newWasBits);
		const long double delta_g_given_z = origWillBits - origWasBits + newWillBits - newWasBits;
		assert(isfinite(delta_g_given_z));
		if(verbose) PP(delta_g_given_z);

		const long double wasOrdersBits = l->SumOfLog2Facts;
		l->moveNode(v, clid);
		const long double willOrdersBits = l->SumOfLog2Facts;
		l->moveNode(v, old_clid);
		const long double verifywasOrdersBits = l->SumOfLog2Facts;
		assert(verifywasOrdersBits == wasOrdersBits);
		if(verbose) PP2(wasOrdersBits , willOrdersBits);
		const long delta_P_z = willOrdersBits-wasOrdersBits;
		if(verbose) PP (delta_P_z);
		assert(isfinite(delta_P_z));
		const long double delta = delta_P_z + delta_g_given_z;
		if(verbose) PP(delta);

		if(log2(drand48()) < delta) { // accept
			l->moveNode(v, clid);
			const long double verifywillOrdersBits = l->SumOfLog2Facts;
			assert(verifywillOrdersBits == willOrdersBits);
			this->numPairs.set(fixed_id, old_clid) = OrigWillPairs;
			this->numEdges.set(fixed_id, old_clid) = OrigWillEdges;
			this->numPairs.set(old_clid, fixed_id) = OrigWillPairs;
			this->numEdges.set(old_clid, fixed_id) = OrigWillEdges;
			this->numPairs.set(fixed_id, clid) = NewWillPairs;
			this->numEdges.set(fixed_id, clid) = NewWillEdges;
			this->numPairs.set(clid, fixed_id) = NewWillPairs;
			this->numEdges.set(clid, fixed_id) = NewWillEdges;
		} else {
		}

	}
} // namespace sbm
