#include "sbm_state.hpp"

#include "macros.hpp"
#include "Range.hpp"
#include <tr1/unordered_set>
#include <vector>
#include <map>
#include <cmath>

using namespace std;
#define printfstring(...) ({ char str[1000]; sprintf(str, __VA_ARGS__) ; (std :: string (str)); })

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
	Labelling :: Labelling(const int _N, const long double alpha) : _N(_N), _k(1), _alpha(alpha) ,log2GammaAlphaPlus(_N+1) {
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
		assert((int)this->clusters.back()->order()==this->_N);

		for(int i=0; i < this->_N+1; i++) {
			log2GammaAlphaPlus.at(i) = LOG2GAMMA(this->_alpha+i);
		}
	}
	State :: State(const GraphType * const g, const bool mega /* = false */, const double alpha /*=1.0*/) : _g(g), vsg(g->get_plain_graph()), _edge_details(g->get_edge_weights()), _N(g->numNodes()), _alpha(alpha), _mega(mega), labelling(this->_N, this->_alpha), _edgeCounts(g->get_edge_weights()) {
		// initialize it with every node in one giant cluster
		this->_k = 1;

		// inform EdgeCounts of all the edges
		this->total_edge_weight = 0.0L;
		for(int relId = 0; relId < this->_g->numRels(); relId++) {
			const std :: pair<int, int> & eps = this->vsg->EndPoints(relId);
			// if(eps.first == eps.second) throw SelfLoopsNotSupported(); // We will assume that self loops have been dealt with appropriately elsewhere.
			const int cl1 = this->labelling.cluster_id.at(eps.first);
			const int cl2 = this->labelling.cluster_id.at(eps.second);
			this->_edgeCounts.inform(cl1,cl2, relId);
			this->total_edge_weight += this->_edge_details->getl2h(relId) + this->_edge_details->geth2l(relId);
		}

		assert(this->sum_weights_BOTH_directions__.empty());
		long double w_verify = 0.0L;
		for(int r = 0; r<this->vsg->numRels(); ++r) {
			pair<int32_t,int32_t> eps = this->vsg->EndPoints(r);
			const long double w = this->_edge_details->getl2h(r) + this->_edge_details->geth2l(r);
			this->sum_weights_BOTH_directions__[ eps ] += w;
			w_verify += w;
			if(eps.first<eps.second) {
				swap(eps.first, eps.second);
				this->sum_weights_BOTH_directions__[ eps ] += w;
			}
		}
		assert(w_verify == this->total_edge_weight);
	}
	Cluster :: Cluster() : _order(0) {
	}
	int Cluster :: order() const {
		// assert( this->_order == (int)this->members.size() ) ;
		return this->_order;
	}
	std :: list<int>:: iterator Cluster :: newMember(const int n) {
		this->members.push_front(n);
		++ this->_order;
		return this->members.begin();
	}
	void Cluster :: eraseMember(const std :: list<int>:: iterator it) {
		this->members.erase(it); // fix up this->clusters.members
		this->_order --;
	}
	const std :: list<int> & Cluster :: get_members() const {
		return members;
	}

	int Labelling :: appendEmptyCluster() {
		const int emptyClusterId = this->_k;
		Cluster * newCluster = new Cluster();
		this->clusters.push_back(newCluster);
		this->_k ++;
		this->SumOfLog2Facts += LOG2GAMMA(this->_alpha);
		return emptyClusterId;
	}
	int State :: appendEmptyCluster() {
		this->labelling.appendEmptyCluster();
		const int newClusterID = this->_k;
		this->_k ++;
		return newClusterID;
	}
	void Labelling :: deleteClusterFromTheEnd() {
		assert(this->_k >= 1);
		Cluster * clusterToDelete = this->clusters.at(this->_k - 1);
		assert(clusterToDelete->order() == 0);
		this->clusters.pop_back();
		delete clusterToDelete;
		this->_k --;
		this->SumOfLog2Facts -= LOG2GAMMA(this->_alpha);
	}
	void State :: deleteClusterFromTheEnd() {
		assert(this->_k >= 1);
		this->_k --;
		this->labelling.deleteClusterFromTheEnd();
	}
	void Labelling :: deltaSumOfLog2Facts(const int oldClOrder, const int clOrder) {
		if(oldClOrder == 0) {
			this->NonEmptyClusters--;
		}
		if(clOrder == 1) {
			this->NonEmptyClusters++;
		}
		this->SumOfLog2Facts +=
					+ (this->log2GammaAlphaPlus.at(oldClOrder   ))
					- (this->log2GammaAlphaPlus.at(oldClOrder+1 ))
					+ (this->log2GammaAlphaPlus.at(clOrder      ))
					- (this->log2GammaAlphaPlus.at(clOrder-1    ))
					;
	}
	void Labelling :: fixUpIterators(const int n, Cluster *cl, Cluster *oldcl) {
		const list<int>:: iterator it = this->its.at(n);
		assert(*it == n);
		oldcl->eraseMember(it);
		const list<int>:: iterator newit = cl->newMember(n); // fix up this->clusters.members
		this->its.at(n) = newit; // fix up this->its
	}
	int Labelling :: moveNode(const int n, const int newClusterID) {
		const int oldClusterID = this->cluster_id.at(n);
		const int oldClusterSize = this->clusters.at(oldClusterID)->order();
		const int _k = static_cast<int>(this->clusters.size());
		assert(newClusterID >= 0 && newClusterID < _k);
		assert(oldClusterID >= 0 && oldClusterID < _k);
		assert(newClusterID != oldClusterID);

		Cluster *cl = this->clusters.at(newClusterID);
		Cluster *oldcl = this->clusters.at(oldClusterID);
		assert(cl);

		this->fixUpIterators(n, cl, oldcl);

		this->cluster_id.at(n) = newClusterID; // fix up this->cluster_id

		assert(oldClusterSize-1 == oldcl->order());
		const int   clOrder =    cl->order();
		const int   oldClOrder =  oldcl->order();
		assert(clOrder > 0);
		this->deltaSumOfLog2Facts(oldClOrder, clOrder);
		// assert(isfinite(this->SumOfLog2Facts));
		return oldClusterID;
	}
	int State :: moveNode(const int n, const int newClusterID) {
		return this->labelling.moveNode(n, newClusterID);
	}
	int State :: isolateNode(const int n) { // create a new (probably temporary) cluster to hold this one node
		assert(n>=0 && n<this->_N);
		const int newClusterID = this->appendEmptyCluster();
		const int oldClusterID = this->labelling.cluster_id.at(n);

		Cluster *cl = this->labelling.clusters.at(newClusterID);
		assert(cl->order()==0);

		this->moveNode(n, newClusterID);

		this->informNodeMove(n, oldClusterID, newClusterID);

		return newClusterID;
	}
	void State :: unIsolateTempNode(const int n, const int newClusterID) { // move a node from its 'temporary' cluster to an existing cluster
		const int oldClusterID = this->labelling.cluster_id.at(n);
		this->moveNode(n, newClusterID);
		assert(oldClusterID+1 == this->_k);
		this->deleteClusterFromTheEnd();

		this->informNodeMove(n, oldClusterID, newClusterID);
	}

	void State :: summarizeEdgeCounts() const {
		for (int i=0; i<this->_k; i++)
		for (int j=0; j<this->_k; j++)
		{
			const long double x = this->_edgeCounts.read(i,j);
			if(x != 0.0L) {
				cout << i
					<< ',' << j
					<< '\t' << x
					<< endl ;
			}
		}
	}
	void State :: blockDetail(const ObjectiveFunction *obj) const {
		long int pairsEncountered = 0;
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
			const long int ni = I->order();
			cout << printfstring("%3ld", ni) << " | ";
			for(int j=0; j<this->_k; j++) {
				const Cluster *J = this->labelling.clusters.at(j);
				assert(J);
				const long int nj = J->order();
				const long double edges = obj->relevantWeight(i,j, &this->_edgeCounts);
				{
					const long double edges_ = this->_edgeCounts.read(i,j) + ( obj->isTwoSided(i,j) ? this->_edgeCounts.read(j,i) : 0);
					assert(edges==edges_);
				}
				long int pairs = ni*nj; // if i==j, then ni==nj
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
				cout << printfstring("%10s %-#5.2Lf", printfstring("%Lg/%ld", edges, pairs).c_str(), edges/static_cast<double>(pairs)) << " | ";
			}
			cout << endl;
		}
		if(obj->directed) {
			assert(pairsEncountered == long(this->_N) * long(this->_N + (obj->selfloops?0:-1) ));
		} else {
			assert(pairsEncountered == long(this->_N) * long(this->_N + (obj->selfloops?1:-1) ) / 2);
		}
		assert(total_edge_weight_verification == this->total_edge_weight);
	}

	void State :: internalCheck() const {
		assert(this->_k>0);
		assert(this->_k == (int)this->labelling.clusters.size());
		assert(this->_N == (int)this->labelling.cluster_id.size());
		assert(this->_N == (int)this->labelling.its.size());
		std :: tr1 :: unordered_set<int> alreadyConsidered;
		long double sumVerify = 0.0L;
		long double sumVerifyFacts = 0.0L;
		long double sumVerifyInternal = 0.0L;
		int NonEmptyVerify = 0;
		for(int CL = 0; CL < this->_k; CL++) {
			const Cluster *cl = this->labelling.clusters.at(CL);
			assert(cl);
			for(list<int>:: const_iterator i = cl->get_members().begin(); i!=cl->get_members().end(); i++) {
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

		{ // this->total_edge_weight
			long double total_edge_weight_verification = 0.0L;
			long double total_edge_weight_Inside_verification = 0.0L; // just the edges that are inside a cluster
			long double total_edge_weight_Outside_verification = 0.0L; // just the edges that are from one cluster to another
			// forEach(const State :: EdgeCounts :: map_type :: value_type  &x, amd :: mk_range(this->_edgeCounts.counts))
			for(int i=0; i<this->_k; i++)
			for(int j=0; j<this->_k; j++)
			{
				pair< pair<int,int>, long double> x(make_pair(i,j),this->_edgeCounts.read(i,j));
				if(x.second == 0) {
					// we can ignore this
				} else {
					assert(x.second > 0);
					assert(x.first.first < this->_k);
					assert(x.first.second < this->_k);
					assert(x.first.first >= 0);
					assert(x.first.second >= 0);

					total_edge_weight_verification += x.second;
					if(x.first.first == x.first.second)
						total_edge_weight_Inside_verification += x.second;
					else
						total_edge_weight_Outside_verification += x.second;
				}
			}
			assert(this->total_edge_weight == total_edge_weight_verification);
			assert(this->total_edge_weight == total_edge_weight_Outside_verification + total_edge_weight_Inside_verification);
			assert(this->_edgeCounts.externalEdgeWeight == total_edge_weight_Outside_verification);
		}
		EdgeCounts edgeCountsVerification(this->_edge_details);
		for(int relId = 0; relId < this->_g->numRels(); relId++) {
			const std :: pair<int, int> & eps = this->vsg->EndPoints(relId);
			const int cl1 = this->labelling.cluster_id.at(eps.first);
			const int cl2 = this->labelling.cluster_id.at(eps.second);
			edgeCountsVerification.inform(cl1,cl2,relId);
		}
		assert(edgeCountsVerification.externalEdgeWeight == this->_edgeCounts.externalEdgeWeight);
		for(int i=0; i<this->_k; i++) {
			edgeCountsVerification.at(i,this->_k-1);
			const std :: vector<long double> & xx = edgeCountsVerification.counts.at(i);

			for(int j=0; j<this->_k; j++) {
				assert(xx.at(j) == this->_edgeCounts.read(i,j));
			}
		}
		/*
		// forEach(const State :: EdgeCounts :: map_type :: value_type  &x, amd :: mk_range(this->_edgeCounts.counts))
		for(int i=0; i<this->_k; i++)
		for(int j=0; j<this->_k; j++)
		{
			if(x.second == 0) {
				// we can ignore this
			} else {
				assert(x.second > 0);
				assert(edgeCountsVerification.counts.count(x.first)==1);
				assert(edgeCountsVerification.counts.at(x.first)==x.second);
				edgeCountsVerification.counts.erase(x.first);
			}
		}
		forEach(const State :: EdgeCounts :: map_type :: value_type &x, amd :: mk_range(edgeCountsVerification.counts)) {
			assert(x.second == 0.0L);
		}
		*/
		//assert(edgeCountsVerification.counts.size()==0);

		// assert(VERYCLOSE(this->pmf(obj) , this->pmf_slow(obj)));

	}

	/*static*/ double State :: NMI(const vector<int> &left, const vector<int> &top) {
		const int N = static_cast<int>(left.size());
		assert(left.size() == top.size());
		map< pair<int,int> , int> joint;
		map< int , int> left_marginal;
		map< int , int> top_marginal;
		for(int n=0; n<N; n++) {
			joint[make_pair(left.at(n), top.at(n))] ++;
			left_marginal[left.at(n)] ++;
			top_marginal[top.at(n)] ++;
		}
		double mutual_information = 0.0;
		double H_left = 0.0;
		forEach( typeof(pair<const int,int>) & l, amd :: mk_range(left_marginal) ) {
			const double P_left  = static_cast<double>(left_marginal[l.first]) / N;
			if(P_left > 0.0)
				H_left += - P_left * log2(P_left);
		}
		assert(isfinite(H_left));
		double H_top = 0.0;
		forEach( typeof(pair<const int,int>) & t, amd :: mk_range(top_marginal) ) {
			const double P_top  = static_cast<double>(top_marginal[t.first]) / N;
			if(P_top > 0.0)
				H_top += - P_top * log2(P_top);
		}
		assert(isfinite(H_top));
		forEach( typeof(pair<const int,int>) & l, amd :: mk_range(left_marginal) ) {
			forEach( typeof(pair<const int,int>) & t, amd :: mk_range(top_marginal) ) {
				const double P_joint = static_cast<double>(joint[make_pair(l.first,t.first)]) / N;
				const double P_left  = static_cast<double>(left_marginal[l.first]) / N;
				const double P_top   = static_cast<double>( top_marginal[t.first]) / N;
				if(P_joint > 0.0)
					mutual_information += P_joint * (log2(P_joint) - log2(P_left) - log2(P_top));
				else
					assert(P_joint == 0.0);
			}
		}
		assert(isfinite(mutual_information));
		const double normalization = std :: max(H_left, H_top);
		if(VERYCLOSE(mutual_information, 0.0L)) {
			mutual_information = 0.0;
		}
		assert(mutual_information >= 0.0);
		assert(mutual_information <= normalization);

		if(normalization > 0.0) {
			const double nmi = mutual_information / normalization;
			assert(isfinite(nmi));
			return nmi;
		} else {
			return 0;
		}
	}

	static string cluster_id_to_string(const int cluster_id) {
		ostringstream oss;
		if(cluster_id<10)
			cout << (char)('0' + cluster_id);
		else if(cluster_id < 36)
			cout << (char)('a' + cluster_id - 10);
		else
			cout << '<' << cluster_id << '>';
		return oss.str();
	}
	static void print_positions(const State *s) {
		if(s->cluster_to_points_map.empty())
			return;
		for(int n=0; n < s->_N; n++) {
			const int k = s->labelling.cluster_id.at(n);
			sbm :: State :: point_type current_position = s->cluster_to_points_map.at(k).at(n);
			cout << "\"" << s->_g->node_name_as_string(n) << "\"," << k << current_position << endl;
		}
	}
	void State :: shortSummary(const ObjectiveFunction *obj, const vector<int> *groundTruth) const {
		cout << endl << " == Summary: ==" << endl;
		PP(this->_k);
		int nonEmptyK = 0;
		for(int k=0; k<this->_k; k++) {
			if(this->labelling.clusters.at(k)->order() > 0)
				nonEmptyK++;
		}
		PP(nonEmptyK);
		PP(this->pmf(obj));
		if (!this->_mega) {
			for (int n=0; n < this->_N; n++) {
				// PP2(n, this->_g->NodeAsString(n));
				const int id_of_cluster = this->labelling.cluster_id.at(n);
				cout << cluster_id_to_string(id_of_cluster);
			}
			cout << endl;
		}
		if(groundTruth) {
			assert(!groundTruth->empty());
			vector<int> z_vector(this->_N);
			for (int n=0; n < this->_N; n++) {
				cout << cluster_id_to_string(groundTruth->at(n));
				const int id_of_cluster = this->labelling.cluster_id.at(n);
				z_vector.at(n) = id_of_cluster;
			}
			assert(z_vector.size() == groundTruth->size());
			cout << " ground truth. NMI=";
			cout
				<< NMI(z_vector, *groundTruth) * 100
				<< " %"
				<< endl;
		}
		this->KandClusterSizes();
		print_positions(this);
	}
	void State :: KandClusterSizes() const {
		vector<int> cluster_sizes;
		for(int k=0; k<this->_k; k++) {
			cluster_sizes.push_back(this->labelling.clusters.at(k)->order());
		}
		sort(cluster_sizes.begin(), cluster_sizes.end(), greater<int>());
		cout << "K and cluster sizes:" << this->_k;
		forEach(int size, amd :: mk_range(cluster_sizes)) {
			cout << '\t' << size;
		}
		cout << endl;
	}


	long double & State :: EdgeCounts :: at(int i,int j) {
			assert(i>=0);
				assert(j>=0);
				if((int) this->counts.size() <= i)
					this->counts.resize(i+1);
				std :: vector<long double> &y = this->counts.at(i);
				if((int) y.size() <= j)
					y.resize(j+1);
				return y.at(j);
	}
	long double State :: EdgeCounts :: read(int i, int j) const {
				assert(i>=0);
				assert(j>=0);
				if((int) this->counts.size() <= i)
					this->counts.resize(i+1);
				std :: vector<long double> &y = this->counts.at(i);
				if((int) y.size() <= j)
					y.resize(j+1);
				return y.at(j);
	}
	State :: EdgeCounts :: EdgeCounts(const graph :: weights :: EdgeDetailsInterface *edge_details) : _edge_details(edge_details), externalEdgeWeight(0.0L) {
	}
	void State :: EdgeCounts :: inform(const int cl1, const int cl2, int relId) { // inform us of an edge between cl1 and cl2
		assert(cl1 >= 0); // && cl1 < this->_k);
		assert(cl2 >= 0); // && cl2 < this->_k);
		long double l2h = this->_edge_details->getl2h(relId);
		long double h2l = this->_edge_details->geth2l(relId);
		this->at(cl1,cl2)+=l2h;
		this->at(cl2,cl1)+=h2l;
		// this->counts[make_pair(cl1,cl2)]++;
		if(cl1 != cl2) {
			this->externalEdgeWeight += l2h+h2l;
		}
	}
	void State :: EdgeCounts :: uninform(const int cl1, const int cl2, int relId) { // inform us of an edge between cl1 and cl2
		assert(cl1 >= 0); // && cl1 < this->_k);
		assert(cl2 >= 0); // && cl2 < this->_k);
		long double l2h = this->_edge_details->getl2h(relId);
		long double h2l = this->_edge_details->geth2l(relId);
		this->at(cl1,cl2)-=l2h;
		this->at(cl2,cl1)-=h2l;
		// this->counts[make_pair(cl1,cl2)]--;
		if(cl1 != cl2) {
			this->externalEdgeWeight -= l2h+h2l;
		}
	}
	void State :: informNodeMove(const int n, const int oldcl, const int newcl) { // a node has just moved from one cluster to another. We must consider it's neighbours for _edgeCounts
		assert(oldcl != newcl);
		const std :: vector<int32_t> & rels = this->vsg->neighbouring_rels_in_order(n);
		forEach(const int relId, amd :: mk_range(rels)) {
			const pair<int,int> eps = this->vsg->EndPoints(relId);
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
	int State :: moveNodeAndInformOfEdges(const int n, const int newcl) {
		const int oldClusterID = this->labelling.cluster_id.at(n);
		this->moveNode(n, newcl);
		this->informNodeMove(n, oldClusterID, newcl);
		return oldClusterID;
	}
	int State :: moveNodeAndInformOfEdges2(const int n, const int newcl) {
		const int oldClusterID = this->labelling.cluster_id.at(n);
		if(oldClusterID != newcl)
			return this->moveNodeAndInformOfEdges(n, newcl);
		return oldClusterID;
	}
	void Labelling :: swapClusters(const int cl1, const int cl2) {
		assert(cl1 != cl2);
		// PP2(cl1,cl2);
		const Cluster * CL1 = this->clusters.at(cl1);
		const Cluster * CL2 = this->clusters.at(cl2);
		forEach(int n, amd :: mk_range(CL1->get_members())) {
			assert(this->cluster_id.at(n) == cl1);
			       this->cluster_id.at(n) =  cl2 ;
		}
		forEach(int n, amd :: mk_range(CL2->get_members())) {
			assert(this->cluster_id.at(n) == cl2);
			       this->cluster_id.at(n) =  cl1 ;
		}
		swap(this->clusters.at(cl1), this->clusters.at(cl2));
	}
	void State :: swapClusters(const int cl1, const int cl2) {
		assert(cl1 != cl2);
		this->labelling.swapClusters(cl1,cl2);
		// I must also swap the _edgeCounts structure
		// PP(this->_edgeCounts.counts.size());
		// PP(this->_k);
		for(int i=0; i<this->_k; i++) {
			this->_edgeCounts.read(i,this->_k - 1); // just to ensure counts :: x is big enough
			std :: vector<long double> &row = this->_edgeCounts.counts.at(i);
			// PP(row.size());
			// PP2(cl1,cl2);
			// PP2(i,this->_k);
			swap(row.at(cl1),row.at(cl2));
		}
		swap(this->_edgeCounts.counts.at(cl1), this->_edgeCounts.counts.at(cl2));
	}

	long double State :: P_z_K() const { // 1 and 2
		// const long double priorOnK = -this->_k; // Geometric(0.5) prior on K
		const long double priorOnK = -LOG2FACT(this->_k); // Poisson(1) prior on K
		return assertNonPositiveFinite(priorOnK);
	}
	long double State :: P_z_orders() const { // given our current this->_k, what's P(z | k)
		return                        LOG2GAMMA(this->_k * this->_alpha) - LOG2GAMMA(this->_k * this->_alpha + this->_N) - this->_k*LOG2GAMMA(this->_alpha) + this->labelling.SumOfLog2Facts;
	}

	long double State :: P_z_slow() const { // given our current this->_k, what's P(z | k)
		const long double K_prior = this->P_z_K();
		long double perCluster_bits = LOG2GAMMA(this->_k * this->_alpha) - LOG2GAMMA(this->_k * this->_alpha + this->_N) - this->_k*LOG2GAMMA(this->_alpha);
		for(int CL=0; CL < this->_k; CL++) {
			const Cluster *cl = this->labelling.clusters.at(CL);
			assert(cl);
			perCluster_bits += LOG2GAMMA(cl->order() + this->_alpha);
		}
		if(VERYCLOSE(perCluster_bits, this->P_z_orders())) {
			perCluster_bits = this->P_z_orders();
		}
		assert(perCluster_bits == this->P_z_orders());
		if(VERYCLOSE(perCluster_bits, 0.0L)) {
			perCluster_bits = 0.0L;
		}
		return assertNonPositiveFinite(K_prior + perCluster_bits);
	}
	long double State :: P_z() const { // given our current this->_k, what's P(z | k)
		const long double K_prior = this->P_z_K();
		return assertNonPositiveFinite(K_prior + P_z_orders());
	}

	
	long double State :: P_edges_given_z_slow(const ObjectiveFunction *obj) const {
		long double edges_bits = 0.0L;
		long int pairsEncountered = 0;
		long double total_edge_weight_verification = 0.0L;
		int blocksEncountered = 0; // should be K*K, or 1/2 * K * (K+1)
		for(int i=0; i<this->_k; i++) {
			for(int j=0; j<this->_k; j++) {
				if(!obj->isValidBlock(i,j))
					break;
				assert(obj->directed || j <= i);
				++ blocksEncountered;


				const long double edges = obj->relevantWeight(i,j, &this->_edgeCounts);
				{
					const long double edges_= this->_edgeCounts.read(i,j) + ( (obj->directed || j==i) ? 0 : this->_edgeCounts.read(j,i));
					assert(edges == edges_);
					const long double edges__ = this->_edgeCounts.read(i,j) + ( obj->isTwoSided(i,j) ? this->_edgeCounts.read(j,i) : 0);
					assert(edges == edges__);
				}
				total_edge_weight_verification += edges;

				const long int pairs_ = obj->numberOfPairsInBlock(i,j, &this->labelling);
				{
					long int pairs = 0;
					const Cluster *I = this->labelling.clusters.at(i);
					assert(I);
					const Cluster *J = this->labelling.clusters.at(j);
					assert(J);
					const long int ni = I->order();
					const long int nj = J->order();
					pairs += ni*nj; // if i==j, then ni==nj
					if (i==j) {
						assert(ni==nj);
						if(obj->directed)
							pairs = ni * (nj - 1); // both directions included
						else
							pairs = ni * (nj - 1) / 2; // just one direction included
					}
					if(i==j && obj->selfloops)
						pairs += ni;
					assert(pairs == pairs_);
				}
				pairsEncountered += pairs_;

				if(!this->cluster_to_points_map.empty() && i==j) { // latentspace on the diagonal, pick them up later
				} else {
					edges_bits += obj->log2OneBlock(edges, pairs_, i==j);
				}
			}
		}
		if(!this->cluster_to_points_map.empty()) { // latentspace on the diagonal,
			assert(!obj->weighted);
			// first, assume every pairs of nodes inside a cluster is *disconnected*, then correct that later.
			long double ls_bits = 0.0L;
			for(int k=0; k < this->_k; ++k) {
				const Cluster *I = this->labelling.clusters.at(k);
				for(auto n : I->get_members()) {
					for(auto m : I->get_members()) {
						if(n==m && !obj->selfloops)
							continue;
						if(!obj->directed && n>m)
							continue;
						if(!obj->directed)
							assert(n<=m);
						// PP2(__LINE__, n,m);
						sbm :: State :: point_type pn = this->cluster_to_points_map.at(k).at(n);
						sbm :: State :: point_type pm = this->cluster_to_points_map.at(k).at(m);
						/*
						 * P(connect) = exp(ls_alpha_k - dist_2)              / (1 + exp(ls_alpha_k - dist_2));
						 * P(disc   ) = (1 + exp(ls_alpha_k - dist_2))/(1 + exp(ls_alpha_k - dist_2)) - exp(ls_alpha_k - dist_2)              / (1 + exp(ls_alpha_k - dist_2));
						 * P(disc   ) = (1+exp() - exp(ls_alpha_k - dist_2))  / (1 + exp(ls_alpha_k - dist_2));
						 * P(disc   ) =  1                                    / (1 + exp(ls_alpha_k - dist_2));
						 */
						const double l2_p_DISconnecting = l2_likelihood( pn, pm, false);
						assert(isfinite(l2_p_DISconnecting));
						ls_bits += l2_p_DISconnecting;
					}
				}
			}
			for(int relId = 0; relId<this->_g->numRels(); ++relId) {
					pair<int32_t, int32_t> eps = this->vsg->EndPoints(relId);
					assert(eps.first <= eps.second);
					if(eps.first == eps.second) {
						assert(obj->selfloops);
					}
					const int n = eps.first;
					const int m = eps.second;
					const int k  = this->labelling.cluster_id.at(n);
					const int km = this->labelling.cluster_id.at(m);

					if(k !=km)
						continue;

					sbm :: State :: point_type pn = this->cluster_to_points_map.at(k).at(n);
					sbm :: State :: point_type pm = this->cluster_to_points_map.at(k).at(m);
					const double l2_p_DISconnecting = l2_likelihood( pn, pm, false);
					const double l2_p_connecting    = l2_likelihood( pn, pm, true);

					// for each edge, we must undo the l2_p_DISconnecting and replace it with a l2_p_connecting
					if(this->_edge_details->getl2h(relId)) {
						ls_bits -= l2_p_DISconnecting;
						ls_bits += l2_p_connecting;
					}
					if(this->_edge_details->geth2l(relId)) {
						assert(obj->directed);
						ls_bits -= l2_p_DISconnecting;
						ls_bits += l2_p_connecting;
					}
			}
			if(1){ // prior on the positions
				for(int n=0; n < this->_N; ++n) {
					const int z_n = this->labelling.cluster_id.at(n);
					sbm :: State :: point_type p = this->cluster_to_points_map.at(z_n).at(n);
					sbm :: State :: point_type pzero;
					pzero.zero();
					const double dist_2 = pzero.dist_2(p);
					// PP2(__LINE__, dist_2);
					const double ln_prior_position = - dist_2 / (2*ls_prior_sigma_2);
					const double l2_prior_position = M_LOG2E * ln_prior_position;
					ls_bits += l2_prior_position;
				}
			}
			assert(isfinite(ls_bits));
			edges_bits += ls_bits;
		}
		if(obj->directed) {
			assert(blocksEncountered == this->_k * this->_k);
			assert(pairsEncountered == long(this->_N) * long(this->_N + (obj->selfloops?0:-1) ));
		} else {
			assert(blocksEncountered == this->_k * (this->_k+1) / 2);
			assert(pairsEncountered == long(this->_N) * long(this->_N + (obj->selfloops?1:-1) ) / 2);
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
		assert(isfinite(edges_bits));
		return edges_bits;
		// return assertNonPositiveFinite(edges_bits); // it doesn't *have* to be finite, just as in our Poisson model (due to the ignoring of x!)
	}
	long double State :: P_edges_given_z(const ObjectiveFunction *obj) const { // this function might be used to try faster ways to calculate the same data. Especially where there are lots of clusters in a small graph.
		const long double slow = this->P_edges_given_z_slow(obj);
		return assertNonPositiveFinite(slow);
	}
	struct BaseLineNotCorrectException {
	};
	long double State :: pmf_slow(const ObjectiveFunction *obj) const {
		const long double t = this->P_edges_given_z_slow(obj) + this->P_z_slow() ;
		assert(isfinite(t));
		return t;
	}
	long double State :: pmf(const ObjectiveFunction *obj) const {
		return this->pmf_slow(obj);
		/*
		const long double fast = P_z_K() + P_z_orders() + this->P_edges_given_z_baseline() + this->P_edges_given_z_correction();
		// const long double slow = this->P_z() + this->P_edges_given_z();
		// assert(VERYCLOSE(fast, slow));
		return assertNonPositiveFinite(fast);
		*/
	}

	ObjectiveFunction :: ObjectiveFunction(const bool s, const bool d, const bool w) : selfloops(s), directed(d), weighted(w) {}
	bool ObjectiveFunction :: isValidBlock(const int i, const int j) const {
		if(this->directed)
			return true;
		else
			return j <=  i;
	}
	bool ObjectiveFunction :: isTwoSided(const int i, const int j) const{ // should the other direction be included when counting the edges?
		if(j==i) // get(i,j) == get(j,i) already
			return false;
		if(this->directed)
			return false;
		return true; // it's undirected and non-diagonal, hence the reverse direction should be included
	}
	long double ObjectiveFunction :: relevantWeight(const int i, const int j, const State :: EdgeCounts *edge_counts) const { // this might include the other direction, if it's undirected
		const long double edges = edge_counts->read(i,j) + ( this->isTwoSided(i,j) ? edge_counts->read(j,i) : 0);
		return edges;
	}
	long int ObjectiveFunction :: numberOfPairsInBlock(const int i, const int j, const Labelling *labelling) const {
		long int pairs=0;
					const Cluster *I = labelling->clusters.at(i);
					assert(I);
					const Cluster *J = labelling->clusters.at(j);
					assert(J);
					const long int ni = I->order();
					const long int nj = J->order();
					pairs += ni*nj; // if i==j, then ni==nj
					if (i==j) {
						assert(ni==nj);
						if(this->directed)
							pairs = ni * (nj - 1); // both directions included
						else
							pairs = ni * (nj - 1) / 2; // just one direction included
						if(this->selfloops)
							pairs += ni;
					}
		return pairs;
	}
	ObjectiveFunction_Bernoulli :: ObjectiveFunction_Bernoulli(const bool s, const bool d, const bool w) : ObjectiveFunction(s,d,w) {}
	long double ObjectiveFunction_Bernoulli :: log2OneBlock(const long double edge_total, const long int pairs, bool isDiagonal) const { // virtual
		if(pairs==0)
			return 0.0L;
		assert(pairs > 0);
		const long double beta_1 = isDiagonal ? 1.0L : 1.0L; // prior
		const long double beta_2 = isDiagonal ? 1.0L : 1.0L; // prior. number of NON-edges in the urn. A large number here for !isDiagonal will give sparse offDiagonal (i.e. community finding)
		assert(isfinite(edge_total));
		assert(edge_total == floor(edge_total));
		assert(edge_total >= 0);
		assert(edge_total <= pairs);
		// edges_bits += - log2l(pairs + 1) - M_LOG2E * gsl_sf_lnchoose(pairs, edges);
		long double p2 = M_LOG2E * gsl_sf_lnbeta(static_cast<double>(edge_total + beta_1), static_cast<double>(pairs - edge_total + beta_2));
		// long double p = - log2l(pairs + 1) - M_LOG2E * gsl_sf_lnchoose(pairs, edge_total);
		// assert(VERYCLOSE(p,p2));
		// assert(0.0 ==         gsl_sf_lnbeta(beta_1, beta_2)); // this'll be zero while \beta_1 = \beta_2 = 1
		const long double p3 = p2 - M_LOG2E * gsl_sf_lnbeta(beta_1, beta_2);
		assertNonPositiveFinite(p3);
		return p3;
	}
	ObjectiveFunction_Poisson :: ObjectiveFunction_Poisson(const bool s, const bool d, const bool w) : ObjectiveFunction(s,d,w) {}
	/*static*/ long double ObjectiveFunction_Poisson :: s = 1.0L;
	/*static*/ long double ObjectiveFunction_Poisson :: theta = 1.0L;
	long double ObjectiveFunction_Poisson :: log2OneBlock(const long double y_b, const long int p_b, bool /*isDiagonal*/) const { // virtual
		if(p_b==0)
			return 0.0L;
		assert(p_b > 0);
		// y_b should be a non-negative integer
		assert(isfinite(y_b));
		// assert(y_b == floor(y_b));
		assert(y_b >= 0);
		const long double s = ObjectiveFunction_Poisson :: s;
		const long double theta = ObjectiveFunction_Poisson :: theta;
		const long double log2p = LOG2GAMMA(s + y_b) + ( s+y_b) * -log2l(p_b + 1.0L/theta) 
		       	-   LOG2GAMMA(s)  - s*log2(theta) // this denominator is important because it depends on the number of blocks.
			;
		// PP2(y_b, p_b);
		// PP(log2p);
		// assertNonPositiveFinite(log2p); // this needn't hold after all. We've ignored \prod x_i, hence this might be positive sometimes.
		return log2p;
	}

double l2_likelihood( sbm :: State :: point_type near, sbm :: State :: point_type far, bool connected) {
	const double dist_2 = near.dist_2(far);
	const double ln_p_connecting       =                            - log(1 + exp( dist_2 - sbm :: ls_alpha_k ));
	const double ln_p_DISconnecting    =                            - log(1 + exp( sbm :: ls_alpha_k - dist_2 ));
	assert(VERYCLOSE(exp(ln_p_connecting) + exp(ln_p_DISconnecting) , 1.0));
	if(connected) {
		const double l2_p_connecting       = M_LOG2E * ln_p_connecting;
		return l2_p_connecting;
	} else {
		//const double ln_p_DISconnecting    = sbm :: ls_alpha_k - dist_2 - log(1 + exp(sbm :: ls_alpha_k - dist_2));
		const double l2_p_DISconnecting    = M_LOG2E * ln_p_DISconnecting;
		return l2_p_DISconnecting;
	}
}
double ls_alpha_k = 0; // TODO: fixed, or put a prior on it? 
} // namespace sbm
