#include "aaron_utils.hpp"
#include "mmsb_state.hpp"

#define SIMPLER_Z 2
#define SIMPLER_EDGES 2
namespace sbm {
	void PairCounts :: set(const int i, const int j, const int to) {
		assert(to >= 0);
		this-> operator [] (std :: make_pair(i,j)) = to;
		this-> operator [] (std :: make_pair(j,i)) = to;
		if(this-> at(std :: make_pair(i,j)) == 0) {
			this-> erase(std :: make_pair(i,j));
			if(i!=j)
				this-> erase(std :: make_pair(j,i));
			cout << __LINE__ << endl;
			this-> at(std :: make_pair(i,j));
			cout << __LINE__ << endl;
		}
	}
	void PairCounts :: inc(const int i, const int j) {
		this-> operator [] (std :: make_pair(i,j)) ++;
		if(i!=j)
			this-> operator [] (std :: make_pair(j,i)) ++;
	}
	void PairCounts :: dec(const int i, const int j) {
		assert(this-> count (std :: make_pair(i,j)) > 0);
		assert(this-> count (std :: make_pair(j,i)) > 0);
		assert(this-> at    (std :: make_pair(i,j)) > 0);
		assert(this-> at    (std :: make_pair(j,i)) > 0);
		this-> operator [] (std :: make_pair(i,j)) --;
		if(i!=j)
			this-> operator [] (std :: make_pair(j,i)) --;
		if(this-> at(std :: make_pair(i,j)) == 0) {
			this-> erase(std :: make_pair(i,j));
		}
		if(i != j && this-> at(std :: make_pair(j,i)) == 0) {
			this-> erase(std :: make_pair(j,i));
		}
	}
	int   PairCounts :: get(const int i, const int j) const {
			if (this->count(std :: make_pair(i,j))==0)
				return 0;
			return this-> at(std :: make_pair(i,j));
	}
	MMSBstate :: MMSBstate (const GraphType * const g): _g(g), _N(g->numNodes()), _k(1) {
		for(int i=0; i< this->_N; i++) {
			ls.push_back(new Labelling(_N, 1.0L));
		}
		for(int i=0; i< this->_N; i++) {
			assert((int)ls.at(i)->clusters.size()==this->_k);
		}
		const int totalNumberOfPairs = this->_N * (this->_N-1) / 2;
		this->numPairs.set(0,0,totalNumberOfPairs);
		PP(this->numPairs.get(0,0));
		this->numEdges.set(0,0,g->numRels());
		PP(this->numEdges.get(0,0));

		for(int i=0; i<this->_N; i++) {
			const bool wasInserted = this->nodeNames.insert(atoi(this->_g->NodeAsString(i))).second;
			assert(wasInserted);
		}
	};
	void MMSBstate :: appendEmptyCluster() {
		for(int i=0; i< this->_N; i++) {
			Labelling * l = ls.at(i);
			assert(l);
			assert((int)l->clusters.size() == this->_k);
			l->appendEmptyCluster();
		}
		this->_k ++;
	}
	long double MMSBstate :: P_z_K() const { // 1 and 2
		// const long double priorOnK = -this->_k; // Exponential prior on K
		const long double priorOnK = -LOG2FACT(this->_k); // Poisson(1) prior on K
		const long double K_dependant_bits = priorOnK + LOG2GAMMA(this->_k) - LOG2GAMMA(this->_k + this->_N);
		return sbm :: assertNonPositiveFinite_line(K_dependant_bits, __LINE__);
	}
	void MMSBstate :: P_zs_given_K() const {
		cout << endl << "  P_zs_given_K()" << endl;
		forEach(int nodeName, amd :: mk_range(this->nodeNames))
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
			cout << '\t' << /*this->P_z_K() + */l->SumOfLog2Facts; // TODO: Make this depend on K
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
		PP(this->pmf_slow());
		for(int k=0; k<this->_k; k++) {
			for(int l=0; l<this->_k; l++) {
				const int pairs = this->numPairs.get(k,l);
				const int edges = this->numEdges.get(k,l);
				cout << k << ',' << l
					<< ":\t" << edges << "/" << pairs
					<< "\t" << double(edges)/double(pairs)
					<< endl;
			}
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
		return sbm :: assertNonPositiveFinite_line(x, __LINE__);
	}
	long double MMSBstate :: pmf_slow() const {
		// this will check that numEdges and numPairs as well as calculate everything
		long double log2facts = 0.0L;
		for(int i=0; i<this->_N; i++) {
			const Labelling * l = ls.at(i);
			{
				long double verify = 0.0L;
				for(int k=0; k<this->_k;k++) {
					const Cluster *cl = l->clusters.at(k);
					const int o = cl->order();
					verify += LOG2FACT(o);
				}
				assert(verify == l->SumOfLog2Facts);
			}
			log2facts += l->SumOfLog2Facts;
		}
		long double edges_bits = 0.0L;
		{
			PairCounts numPairsVerify;
			PairCounts numEdgesVerify;
			for(int w = 0; w< this->_N; w++) {
				for(int v = 0; v< w; v++) {
					const bool areConnected = this->_g->are_connected(w,v);
					// PP2(w,v);
					const int clid1 = this->ls.at(w)->cluster_id.at(v);
					const int clid2 = this->ls.at(v)->cluster_id.at(w);
					assert(clid1 >= 0 && clid2 >= 0);
					assert(clid1 < this->_k && clid2 < this->_k);
					// PP2(clid1,clid2);
					if(clid1 == clid2) {
						numPairsVerify.inc(clid1,clid2);
						if(areConnected)
							numEdgesVerify.inc(clid1,clid2);
					} else {
						numPairsVerify.inc(clid2,clid1);
						if(areConnected) {
							numEdgesVerify.inc(clid1,clid2);
						}
					}
				}
			}
			DYINGWORDS(this->numEdges.size() == numEdgesVerify.size()) {
				PP2(this->numEdges.size() , numEdgesVerify.size());
				cout << " numEdgesVerify" << endl;
				forEach(const typeof(pair< pair<int,int> ,int>) &x, amd :: mk_range(numEdgesVerify)) {
					PP3(x.first.first,x.first.second,x.second);
				}
				cout << " this->numEdges" << endl;
				forEach(const typeof(pair< pair<int,int> ,int>) &x, amd :: mk_range(this->numEdges)) {
					PP3(x.first.first,x.first.second,x.second);
				}
			}
			assert(this->numEdges.size() == numEdgesVerify.size());
			assert(this->numEdges == numEdgesVerify);
			DYINGWORDS(this->numPairs.size() == numPairsVerify.size()) {
				PP2(this->numPairs.size() , numPairsVerify.size());
				cout << " numPairsVerify" << endl;
				forEach(const typeof(pair< pair<int,int> ,int>) &x, amd :: mk_range(numPairsVerify)) {
					PP3(x.first.first,x.first.second,x.second);
				}
				cout << " this->numPairs" << endl;
				forEach(const typeof(pair< pair<int,int> ,int>) &x, amd :: mk_range(this->numPairs)) {
					PP3(x.first.first,x.first.second,x.second);
				}
			}
			assert(this->numPairs.size() == numPairsVerify.size());
			assert(this->numPairs == numPairsVerify);
		}
		for(int k=0; k<this->_k; k++) {
			for(int l=0; l<=k; l++) {
				const int pairs = this->numPairs.get(k,l);
				const int edges = this->numEdges.get(k,l);
				cout << k << ',' << l
					<< ":\t" << edges << "/" << pairs
					<< "\t" << double(edges)/double(pairs)
					<< endl;
				edges_bits += P_beta_binomial(pairs,edges);
			}
		}
		PP3(          log2facts , edges_bits, log2facts+edges_bits);
		PP3(SIMPLER_Z*log2facts , edges_bits, SIMPLER_Z*log2facts+SIMPLER_EDGES*edges_bits);
		return SIMPLER_Z*log2facts + SIMPLER_EDGES*edges_bits;
	}
	long double MMSBstate :: moveOnePair(const int w,const int v,const int clid) const { // node w, when interacting with v, should take identity clid
		const bool verbose = false;
		assert(w!=v);
		assert(clid >= 0 && clid < this->_k);
		Labelling * l = ls.at(w);
		const int old_clid = l->cluster_id.at(v);
		if(clid == old_clid)
			return 0.0L;
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
		// if(verbose)
			PP(delta_g_given_z);

		const long double wasOrdersBits = l->SumOfLog2Facts;
		cout << __LINE__ << endl;
		int old_clid_verify = l->moveNode(v, clid);
		assert(old_clid_verify == old_clid);
		cout << __LINE__ << endl;
		const long double willOrdersBits = l->SumOfLog2Facts;
		int pretend_clid_verify = l->moveNode(v, old_clid);
		assert(pretend_clid_verify == clid);
		cout << __LINE__ << endl;
		const long double verifywasOrdersBits = l->SumOfLog2Facts;
		assert(verifywasOrdersBits == wasOrdersBits);
		if(verbose) PP2(wasOrdersBits , willOrdersBits);
		const long double delta_P_z = willOrdersBits-wasOrdersBits;
		// if(verbose)
			PP (delta_P_z);
		assert(isfinite(delta_P_z));
		const long double delta = SIMPLER_Z*delta_P_z + SIMPLER_EDGES*delta_g_given_z;
		if(verbose) PP(delta);

		return delta;
	}
	int MMSBstate :: performMoveAndUpdateEdges(const int w, const int v, const int clid) {
		assert(w!=v); // TODO: need to move efficiently store the Labelling, given that self loops are meaningless
		const bool areConnected = this->_g->are_connected(w,v);
		Labelling * l = ls.at(w);
		const int old_clid = l->moveNode(v, clid);
		Labelling * l2 = ls.at(v);
		assert(l2);
		const int fixed_id = l2->cluster_id.at(w);
		this->numPairs.dec(fixed_id, old_clid);
		if(areConnected)
			this->numEdges.dec(fixed_id, old_clid);
		this->numPairs.inc(fixed_id, clid);
			if(areConnected)
				this->numEdges.inc(fixed_id, clid);
		return old_clid;
	}
	long double MMSBstate :: MetropolisMoveOnePair(const int w,const int v,const int clid) { // node w, when interacting with v, should take identity clid
		const Labelling * l = ls.at(w);
		const int old_clid = l->cluster_id.at(v);
		if(old_clid == clid) {
			cout << "abandon" << endl;
			return 0.0L;
		}

		const long double delta = this->moveOnePair(w,v,clid);
		const long double delta2= this->moveOnePair(w,v,clid);
		PP3(w,v,clid);
		PP2(delta , delta2);
		assert(delta == delta2);
		if(log2(drand48()) < delta) { // accept
			cout << " ACCEPT" << endl;
			cout << __LINE__ << endl;
			int old_clid_verify = this-> performMoveAndUpdateEdges(w,v,clid);
			assert(old_clid_verify == old_clid);
			// const long double verifywillOrdersBits = l->SumOfLog2Facts;
			// assert(verifywillOrdersBits == willOrdersBits);
			{
				Labelling * l = ls.at(w);
				const int new_clid = l->cluster_id.at(v);
				assert(new_clid == clid);
				const long double reverse_delta = this->moveOnePair(w,v,old_clid);
				PP2(reverse_delta, delta);
				assert(VERYCLOSE(0.0L, reverse_delta+delta));
			}
			return delta;
		} else {
			return 0.0L;
		}
	}
} // namespace sbm
