#include "sbm_state.hpp"
#include "aaron_utils.hpp"
using namespace shmGraphRaw;
namespace sbm {
#define assertNonPositiveFinite(x) assertNonPositiveFinite_line(x, __LINE__)
	long double assertNonPositiveFinite_line(const long double x, const int lineno) {
		DYINGWORDS(isfinite(x) && x<=0.0L) {
			PP2(x,lineno);
		}
		assert(isfinite(x));
		assert(x<=0.0L);
		return x;
	}
	Labelling::Labelling(const int _N, const long double alpha) : _N(_N), _k(1), _alpha(alpha) {
		this->SumOfLog2Facts   = LOG2GAMMA(this->_N + this->_alpha);
		this->NonEmptyClusters = 1;
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
	}
	State::State(const GraphType * const g, const shmGraphRaw:: EdgeDetailsInterface * edge_details) : _g(g), _edge_details(edge_details), _N(g->numNodes()), _alpha(1.0L), labelling(this->_N, this->_alpha), _edgeCounts(edge_details) {
		// initialize it with every node in one giant cluster
		this->_k = 1;

		// inform EdgeCounts of all the edges
		this->total_edge_weight = 0.0L;
		for(int relId = 0; relId < this->_g->numRels(); relId++) {
			const std::pair<int, int> & eps = this->_g->EndPoints(relId);
			// if(eps.first == eps.second) throw SelfLoopsNotSupported(); // We will assume that self loops have been dealt with appropriately elsewhere.
			const int cl1 = this->labelling.cluster_id.at(eps.first);
			const int cl2 = this->labelling.cluster_id.at(eps.second);
			this->_edgeCounts.inform(cl1,cl2, relId);
			this->total_edge_weight += edge_details->getl2h(relId) + edge_details->geth2l(relId);
		}

		// to ensure the nodes are dealt with in the order of their integer node name
		for(int n=0; n<this->_N;n++) {
			nodeNamesInOrder.insert(atoi(this->_g->NodeAsString(n)));
		}
		assert( (int)nodeNamesInOrder.size() == this->_N);
	}
	const int Cluster::order() const {
		return this->members.size();
	}
	std::list<int>::iterator Cluster::newMember(const int n) {
		this->members.push_front(n);
		return this->members.begin();
	}

	void Labelling:: appendEmptyCluster() {
		Cluster * newCluster = new Cluster();
		this->clusters.push_back(newCluster);
		this->_k ++;
		this->SumOfLog2Facts += LOG2GAMMA(this->_alpha);
	}
	int State:: appendEmptyCluster() {
		this->labelling.appendEmptyCluster();
		const int newClusterID = this->_k;
		this->_k ++;
		return newClusterID;
	}
	void Labelling:: deleteClusterFromTheEnd() {
		assert(this->_k >= 1);
		Cluster * clusterToDelete = this->clusters.at(this->_k - 1);
		assert(clusterToDelete->order() == 0);
		this->clusters.pop_back();
		delete clusterToDelete;
		this->_k --;
		this->SumOfLog2Facts -= LOG2GAMMA(this->_alpha);
	}
	void State:: deleteClusterFromTheEnd() {
		assert(this->_k >= 1);
		this->_k --;
		this->labelling.deleteClusterFromTheEnd();
	}
	int Labelling:: moveNode(const int n, const int newClusterID) {
		const int _k = this->clusters.size();
		const int oldClusterID = this->cluster_id.at(n);
		const int oldClusterSize = this->clusters.at(oldClusterID)->order();
		assert(newClusterID >= 0 && newClusterID < _k);
		assert(oldClusterID >= 0 && oldClusterID < _k);
		assert(newClusterID != oldClusterID);

		Cluster *cl = this->clusters.at(newClusterID);
		Cluster *oldcl = this->clusters.at(oldClusterID);
		assert(cl);

		const list<int>::iterator it = this->its.at(n);
		assert(*it == n);
		oldcl->members.erase(it); // fix up this->clusters.members
		const list<int>::iterator newit = cl->newMember(n); // fix up this->clusters.members
		this->its.at(n) = newit; // fix up this->its

		this->cluster_id.at(n) = newClusterID; // fix up this->cluster_id

		assert(oldClusterSize-1 == oldcl->order());
		if(oldcl->order() == 0) {
			this->NonEmptyClusters--;
		}
		if(cl->order() == 1) {
			this->NonEmptyClusters++;
		}
		const int   to_order =    cl->order();
		assert(to_order > 0);
		this->SumOfLog2Facts +=
					+ (oldcl->order()<0?0.0L:LOG2GAMMA(this->_alpha+oldcl->order()))
					- (oldcl->order()<0?0.0L:LOG2GAMMA(this->_alpha+oldcl->order()+1))
					;
		this->SumOfLog2Facts +=
					+ (cl->order()<0   ?0.0L:LOG2GAMMA(this->_alpha+cl->order()))
					- (cl->order()<0   ?0.0L:LOG2GAMMA(this->_alpha+cl->order()-1))
					;
		assert(isfinite(this->SumOfLog2Facts));
		return oldClusterID;
	}
	void State::moveNode(const int n, const int newClusterID) {
		this->labelling.moveNode(n, newClusterID);
	}
	int State::isolateNode(const int n) { // create a new (probably temporary) cluster to hold this one node
		assert(n>=0 && n<this->_N);
		const int newClusterID = this->appendEmptyCluster();
		const int oldClusterID = this->labelling.cluster_id.at(n);

		Cluster *cl = this->labelling.clusters.at(newClusterID);
		assert(cl->order()==0);

		this->moveNode(n, newClusterID);

		this->informNodeMove(n, oldClusterID, newClusterID);

		return newClusterID;
	}
	void State::unIsolateTempNode(const int n, const int newClusterID) { // move a node from its 'temporary' cluster to an existing cluster
		const int oldClusterID = this->labelling.cluster_id.at(n);
		this->moveNode(n, newClusterID);
		assert(oldClusterID+1 == this->_k);
		this->deleteClusterFromTheEnd();

		this->informNodeMove(n, oldClusterID, newClusterID);
	}

	void State:: summarizeEdgeCounts() const {
		forEach(const State :: EdgeCounts :: map_type :: value_type  &x, amd::mk_range(this->_edgeCounts.counts)) {
			if(x.second != 0.0L) {
				cout << x.first.first
					<< ',' << x.first.second
					<< '\t' << x.second
					<< endl ;
			}
		}
	}
	void State:: blockDetail(const ObjectiveFunction *obj) const {
		int pairsEncountered = 0;
		long double total_edge_weight_verification = 0.0L;

		cout << "    | ";
		for(int j=0; j<this->_k; j++) {
			const Cluster *J = this->labelling.clusters.at(j);
			assert(J);
			const int nj = J->order();
			cout << printfstring("      %5d     ", nj) << "   ";
		}
		cout << endl;

		for(int i=0; i<this->_k; i++) {
			const Cluster *I = this->labelling.clusters.at(i);
			assert(I);
			const int ni = I->order();
			cout << printfstring("%3d", ni) << " | ";
			for(int j=0; j<this->_k; j++) {
				const Cluster *J = this->labelling.clusters.at(j);
				assert(J);
				const int nj = J->order();
				const long double edges = this->_edgeCounts.get(i,j) + ( obj->isTwoSided(i,j) ? this->_edgeCounts.get(j,i) : 0);
				int pairs = ni*nj; // if i==j, then ni==nj
				if (i==j) {
					assert(ni==nj);
					if(obj->directed)
						pairs = ni * (nj - 1); // both directions included
					else
						pairs = ni * (nj - 1) / 2; // just one direction included
				}
				if(i==j && obj->selfloops)
					pairs += ni;
				if(obj->directed || j<=i) {
					pairsEncountered += pairs;
					total_edge_weight_verification += edges;
				}
				cout << printfstring("%10s %-#5.2Lf", printfstring("%Lg/%d", edges, pairs).c_str(), edges/double(pairs)) << " | ";
			}
			cout << endl;
		}
		if(obj->directed) {
			assert(pairsEncountered == this->_N * (this->_N + (obj->selfloops?0:-1) ));
		} else {
			assert(pairsEncountered == this->_N * (this->_N + (obj->selfloops?1:-1) ) / 2);
		}
		assert(total_edge_weight_verification == this->total_edge_weight);
	}

	void State::internalCheck() const {
		assert(this->_k>0);
		assert(this->_k == (int)this->labelling.clusters.size());
		assert(this->_N == (int)this->labelling.cluster_id.size());
		assert(this->_N == (int)this->labelling.its.size());
		boost::unordered_set<int> alreadyConsidered;
		long double sumVerify = 0.0L;
		long double sumVerifyFacts = 0.0L;
		long double sumVerifyInternal = 0.0L;
		int NonEmptyVerify = 0;
		for(int CL = 0; CL < this->_k; CL++) {
			const Cluster *cl = this->labelling.clusters.at(CL);
			assert(cl);
			for(list<int>::const_iterator i = cl->members.begin(); i!=cl->members.end(); i++) {
				const int n = *i;
				assert(n>=0 && n<this->_N);
				bool wasAccepted = alreadyConsidered.insert(n).second;
				assert(wasAccepted);
				assert(n == *this->labelling.its.at(n));
				assert(CL == this->labelling.cluster_id.at(n));
				assert(i == this->labelling.its.at(n));
			}
			if(cl->order() >0)
				NonEmptyVerify++;
			sumVerifyFacts += LOG2GAMMA(cl->order() + this->_alpha);
			if(cl->order() >=2) {
				sumVerify += log2l(cl->order());
				sumVerifyInternal += log2l(cl->order()*(cl->order()-1)/2);
			}
		}
		assert(NonEmptyVerify == this->labelling.NonEmptyClusters);
		assert(VERYCLOSE(sumVerifyFacts , this->labelling.SumOfLog2Facts));
		this->labelling.SumOfLog2Facts = sumVerifyFacts;
		assert((int)alreadyConsidered.size() == this->_N);

		EdgeCounts edgeCountsVerification(this->_edge_details);
		for(int relId = 0; relId < this->_g->numRels(); relId++) {
			const std::pair<int, int> & eps = this->_g->EndPoints(relId);
			const int cl1 = this->labelling.cluster_id.at(eps.first);
			const int cl2 = this->labelling.cluster_id.at(eps.second);
			edgeCountsVerification.inform(cl1,cl2,relId);
		}
		forEach(const State :: EdgeCounts :: map_type :: value_type  &x, amd::mk_range(this->_edgeCounts.counts)) {
			if(x.second == 0) {
				// we can ignore this
			} else {
				assert(x.second > 0);
				assert(edgeCountsVerification.counts.count(x.first)==1);
				assert(edgeCountsVerification.counts.at(x.first)==x.second);
				edgeCountsVerification.counts.erase(x.first);
			}
		}
		forEach(const State :: EdgeCounts :: map_type :: value_type &x, amd::mk_range(edgeCountsVerification.counts)) {
			assert(x.second == 0.0L);
		}
		//assert(edgeCountsVerification.counts.size()==0);

		// assert(VERYCLOSE(this->pmf(obj) , this->pmf_slow(obj)));

	}

	void State:: shortSummary(const ObjectiveFunction *obj) const {
		cout << endl << " == Summary: ==" << endl;
		PP(this->_k);
		int nonEmptyK = 0;
		for(int k=0; k<this->_k; k++) {
			if(this->labelling.clusters.at(k)->order() > 0)
				nonEmptyK++;
		}
		PP(nonEmptyK);
		PP(this->pmf(obj));
		forEach(int node_name, amd::mk_range(this->nodeNamesInOrder))
		// for(int n=0; n<this->_N;n++)
		{
			const int n = this->_g->StringToNodeId(printfstring("%d", node_name).c_str());
			// PP2(n, this->_g->NodeAsString(n));
			const int id_of_cluster = this->labelling.cluster_id.at(n);
			if(id_of_cluster<10)
				cout << (char)('0' + id_of_cluster);
			else if(id_of_cluster < 36)
				cout << (char)('a' + id_of_cluster - 10);
			else
				cout << '<' << id_of_cluster << '>';
		}
		cout << endl;
	}


	State:: EdgeCounts:: EdgeCounts(const EdgeDetailsInterface *edge_details) : _edge_details(edge_details) {
	}
	void State::EdgeCounts::inform(const int cl1, const int cl2, int relId) { // inform us of an edge between cl1 and cl2
		assert(cl1 >= 0); // && cl1 < this->_k);
		assert(cl2 >= 0); // && cl2 < this->_k);
		long double l2h = this->_edge_details->getl2h(relId);
		long double h2l = this->_edge_details->geth2l(relId);
		this->counts[make_pair(cl1,cl2)]+=l2h;
		this->counts[make_pair(cl2,cl1)]+=h2l;
		// this->counts[make_pair(cl1,cl2)]++;
	}
	void State::EdgeCounts::uninform(const int cl1, const int cl2, int relId) { // inform us of an edge between cl1 and cl2
		assert(cl1 >= 0); // && cl1 < this->_k);
		assert(cl2 >= 0); // && cl2 < this->_k);
		long double l2h = this->_edge_details->getl2h(relId);
		long double h2l = this->_edge_details->geth2l(relId);
		this->counts[make_pair(cl1,cl2)]-=l2h;
		this->counts[make_pair(cl2,cl1)]-=h2l;
		// this->counts[make_pair(cl1,cl2)]--;
		DYINGWORDS(this->counts[make_pair(cl1,cl2)] >= 0) {
			PP2(cl1,cl2);
		}
	}
	long double  State::EdgeCounts:: get(const int cl1, const int cl2) const throw() {
		long double lowToHigh;
		try { lowToHigh = this->counts.at(make_pair(cl1,cl2)); } catch (const std::out_of_range &e) { lowToHigh = 0.0L; }
		return lowToHigh;
	}
	void State::informNodeMove(const int n, const int oldcl, const int newcl) { // a node has just moved from one cluster to another. We must consider it's neighbours for _edgeCounts
		assert(oldcl != newcl);
		const PlainMem::mmap_uset_of_ints & rels = this->_g->myRels(n);
		forEach(const int relId, amd::mk_range(rels)) {
			const pair<int,int> eps = this->_g->EndPoints(relId);
			const int fstClusterNew = this->labelling.cluster_id.at(eps.first);
			const int sndClusterNew = this->labelling.cluster_id.at(eps.second);
			int fstClusterOld = fstClusterNew;
			if(n == eps.first) {
				assert(fstClusterOld == newcl);
				fstClusterOld = oldcl;
			}
			int sndClusterOld = sndClusterNew;
			if(n == eps.second) {
				assert(sndClusterOld == newcl);
				sndClusterOld = oldcl;
			}
			this->_edgeCounts.uninform(fstClusterOld, sndClusterOld, relId);
			this->_edgeCounts.  inform(fstClusterNew, sndClusterNew, relId);
		}
	}
	void State:: moveNodeAndInformOfEdges(const int n, const int newcl) {
		const int oldClusterID = this->labelling.cluster_id.at(n);
		this->moveNode(n, newcl);
		this->informNodeMove(n, oldClusterID, newcl);
	}
	void Labelling :: swapClusters(const int cl1, const int cl2) {
		assert(cl1 != cl2);
		// PP2(cl1,cl2);
		const Cluster * CL1 = this->clusters.at(cl1);
		const Cluster * CL2 = this->clusters.at(cl2);
		forEach(int n, amd::mk_range(CL1->members)) {
			assert(this->cluster_id.at(n) == cl1);
			       this->cluster_id.at(n) =  cl2 ;
		}
		forEach(int n, amd::mk_range(CL2->members)) {
			assert(this->cluster_id.at(n) == cl2);
			       this->cluster_id.at(n) =  cl1 ;
		}
		swap(this->clusters.at(cl1), this->clusters.at(cl2));
	}
	void State :: swapClusters(const int cl1, const int cl2) {
		this->labelling.swapClusters(cl1,cl2);
		// I must also swap the _edgeCounts structure
		const size_t preSize = this->_edgeCounts.counts.size();
		// PP(this->_edgeCounts.counts.size());
		// PP(this->_k);
		if(0)
		for(State :: EdgeCounts :: map_type :: const_iterator i = this->_edgeCounts.counts.begin(); i != this->_edgeCounts.counts.end(); i++ ) {
			PP3(i->first.first, i->first.second, i->second);
		}
		// cout << "swapping" << endl;
		State :: EdgeCounts :: map_type merge_these_back_in_again;
		for(State :: EdgeCounts :: map_type ::       iterator i = this->_edgeCounts.counts.begin(); i != this->_edgeCounts.counts.end();  ) {
			if( i->first.first  == cl1
			 || i->first.first  == cl2
			 || i->first.second == cl1
			 || i->first.second == cl2
			 ) {
				// flip them
				pair< pair<int,int> , int> tmp = *i;
				// PP3(tmp.first.first, tmp.first.second, tmp.second);
				if(tmp.first.first == cl1) tmp.first.first = cl2; else if(tmp.first.first == cl2) tmp.first.first = cl1;
				if(tmp.first.second == cl1) tmp.first.second = cl2; else if(tmp.first.second == cl2) tmp.first.second = cl1;
				// remove it
				// readd them
				// PP3(tmp.first.first, tmp.first.second, tmp.second);
				i = this->_edgeCounts.counts.erase(i); // erase it and leave the iterator pointing at the next item
				merge_these_back_in_again.insert(tmp);
			} else
				++i;
		}
		// PP(merge_these_back_in_again.size());
		for(State :: EdgeCounts :: map_type ::       iterator i = merge_these_back_in_again.begin(); i != merge_these_back_in_again.end(); i++) {
			this->_edgeCounts.counts.insert(*i);
		}
		// PP2(preSize, this->_edgeCounts.counts.size());
		assert(preSize == this->_edgeCounts.counts.size());
	}

	long double State:: P_z_K() const { // 1 and 2
		// const long double priorOnK = -this->_k; // Geometric(0.5) prior on K
		const long double priorOnK = -LOG2FACT(this->_k); // Poisson(1) prior on K
		return assertNonPositiveFinite(priorOnK);
	}
	long double State:: P_z_orders() const { // given our current this->_k, what's P(z | k)
		return                        LOG2GAMMA(this->_k * this->_alpha) - LOG2GAMMA(this->_k * this->_alpha + this->_N) - this->_k*LOG2GAMMA(this->_alpha) + this->labelling.SumOfLog2Facts;
	}

	long double State:: P_z_slow() const { // given our current this->_k, what's P(z | k)
		const long double K_prior = this->P_z_K();
		long double perCluster_bits = LOG2GAMMA(this->_k * this->_alpha) - LOG2GAMMA(this->_k * this->_alpha + this->_N) - this->_k*LOG2GAMMA(this->_alpha);
		for(int CL=0; CL < this->_k; CL++) {
			const Cluster *cl = this->labelling.clusters.at(CL);
			assert(cl);
			perCluster_bits += LOG2GAMMA(cl->order() + this->_alpha);
		}
		assert(perCluster_bits == P_z_orders());
		if(VERYCLOSE(perCluster_bits, 0.0L)) {
			perCluster_bits = 0.0L;
		}
		return assertNonPositiveFinite(K_prior + perCluster_bits);
	}

	
	long double State:: P_edges_given_z_slow(const ObjectiveFunction *obj) const {
		long double edges_bits = 0.0L;
		int pairsEncountered = 0;
		long double total_edge_weight_verification = 0.0L;
		int blocksEncountered = 0; // should be K*K, or 1/2 * K * (K+1)
		for(int i=0; i<this->_k; i++) {
			const Cluster *I = this->labelling.clusters.at(i);
			assert(I);
			for(int j=0; j<this->_k; j++) {
				if(!obj->isValidBlock(i,j))
					break;
				assert(obj->directed || j <= i);
				++ blocksEncountered;
				const Cluster *J = this->labelling.clusters.at(j);
				assert(J);
				const int ni = I->order();
				const int nj = J->order();
				const long double edges = this->_edgeCounts.get(i,j) + ( obj->isTwoSided(i,j) ? this->_edgeCounts.get(j,i) : 0);
				{
					const long double edges_= this->_edgeCounts.get(i,j) + ( (obj->directed || j==i) ? 0 : this->_edgeCounts.get(j,i));
					assert(edges == edges_);
				}
				total_edge_weight_verification += edges;
				/*
				cout << ni << " x " << nj;
				if(i==j)
					cout << "-1 x1/2";
				cout	<< " pairs. "
					<< edges << " edges."
					<< endl;
					*/
				int pairs = ni*nj; // if i==j, then ni==nj
				if (i==j) {
					assert(ni==nj);
					if(obj->directed)
						pairs = ni * (nj - 1); // both directions included
					else
						pairs = ni * (nj - 1) / 2; // just one direction included
				}
				if(i==j && obj->selfloops)
					pairs += ni;
				// PP2(i,j);
				// PP2(ni,nj);
				// PP2(edges, pairs);
				if(pairs > 0) {
					pairsEncountered += pairs;
					edges_bits += obj->log2OneBlock(edges, pairs, i==j);
				}
				// PP(edges_bits);
			}
		}
		if(obj->directed) {
			assert(blocksEncountered == this->_k * this->_k);
			assert(pairsEncountered == this->_N * (this->_N + (obj->selfloops?0:-1) ));
		} else {
			assert(blocksEncountered == this->_k * (this->_k+1) / 2);
			assert(pairsEncountered == this->_N * (this->_N + (obj->selfloops?1:-1) ) / 2);
		}
		DYINGWORDS(total_edge_weight_verification == this->total_edge_weight) {
			PP2(total_edge_weight_verification , this->total_edge_weight);
			PP (total_edge_weight_verification - this->total_edge_weight);
		}
		/*
		DYINGWORDS(VERYCLOSE(edges_bits, this->P_edges_given_z_baseline() + this->P_edges_given_z_correction())) {
			cout << endl << "DYINGWORDS:" << __LINE__ << endl;
			PP(this->P_edges_given_z_baseline() + this->P_edges_given_z_correction() - edges_bits);
			PP(this->P_edges_given_z_baseline() + this->P_edges_given_z_correction());
			PP(this->P_edges_given_z_baseline());
			PP(edges_bits_no_edges);
			PP(this->P_edges_given_z_correction());
			PP(- edges_bits);
		}
		*/
		return assertNonPositiveFinite(edges_bits);
	}
	long double State:: P_edges_given_z(const ObjectiveFunction *obj) const { // this function might be used to try faster ways to calculate the same data. Especially where there are lots of clusters in a small graph.
		const long double slow = this->P_edges_given_z_slow(obj);
		return assertNonPositiveFinite(slow);
	}
	struct BaseLineNotCorrectException {
	};
	long double State:: pmf_slow(const ObjectiveFunction *obj) const {
		return assertNonPositiveFinite( this->P_edges_given_z_slow(obj) + this->P_z_slow() );
	}
	long double State:: pmf(const ObjectiveFunction *obj) const {
		return this->pmf_slow(obj);
		/*
		const long double fast = P_z_K() + P_z_orders() + this->P_edges_given_z_baseline() + this->P_edges_given_z_correction();
		// const long double slow = this->P_z() + this->P_edges_given_z();
		// assert(VERYCLOSE(fast, slow));
		return assertNonPositiveFinite(fast);
		*/
	}

	ObjectiveFunction :: ObjectiveFunction(const bool s, const bool d, const bool w) : selfloops(s), directed(d), weighted(w) {}
	bool ObjectiveFunction:: isValidBlock(const int i, const int j) const {
		if(this->directed)
			return true;
		else
			return j <=  i;
	}
	bool ObjectiveFunction:: isTwoSided(const int i, const int j) const{ // should the other direction be included when counting the edges?
		if(j==i) // get(i,j) == get(j,i) already
			return false;
		if(this->directed)
			return false;
		return true; // it's undirected and non-diagonal, hence the reverse direction should be included
	}
	ObjectiveFunction_Bernoulli :: ObjectiveFunction_Bernoulli(const bool s, const bool d, const bool w) : ObjectiveFunction(s,d,w) {}
	long double ObjectiveFunction_Bernoulli :: log2OneBlock(const long double edge_total, const int pairs, bool isDiagonal) const { // virtual
		assert(pairs > 0);
		const long double beta_1 = isDiagonal ? 1.0L : 1.0L; // prior
		const long double beta_2 = isDiagonal ? 1.0L : 1.0L; // prior. number of NON-edges in the urn. A large number here for !isDiagonal will give sparse offDiagonal (i.e. community finding)
		assert(isfinite(edge_total));
		assert(edge_total == floor(edge_total));
		assert(edge_total >= 0);
		assert(edge_total <= pairs);
		// edges_bits += - log2l(pairs + 1) - M_LOG2E * gsl_sf_lnchoose(pairs, edges);
		long double p2 = M_LOG2E * gsl_sf_lnbeta(edge_total + beta_1, pairs - edge_total + beta_2);
		// long double p = - log2l(pairs + 1) - M_LOG2E * gsl_sf_lnchoose(pairs, edge_total);
		// assert(VERYCLOSE(p,p2));
		// assert(0.0 ==         gsl_sf_lnbeta(beta_1, beta_2)); // this'll be zero while \beta_1 = \beta_2 = 1
		const long double p3 = p2 - M_LOG2E * gsl_sf_lnbeta(beta_1, beta_2);
		assertNonPositiveFinite(p3);
		return p3;
	}
	ObjectiveFunction_Poisson :: ObjectiveFunction_Poisson(const bool s, const bool d, const bool w) : ObjectiveFunction(s,d,w) {}
	long double ObjectiveFunction_Poisson :: log2OneBlock(const long double y_b, const int p_b, bool isDiagonal) const { // virtual
		assert(p_b > 0);
		// y_b should be a non-negative integer
		assert(isfinite(y_b));
		// assert(y_b == floor(y_b));
		assert(y_b >= 0);
		const long double s = 1.0L;
		const long double theta = 2.0L;
		const long double log2p = LOG2GAMMA(s + y_b) + ( s+y_b) * -log2(p_b + 1.0L/theta) 
		       	-   LOG2GAMMA(s)  - s*log2(theta) // this denominator is important because it depends on the number of blocks.
			;
		assertNonPositiveFinite(log2p);
		return log2p;
	}
} // namespace sbm
