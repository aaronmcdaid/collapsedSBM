#include "graph/network.hpp"
#include "graph/loading.hpp"
using namespace std;
#include <algorithm>
#include <vector>
#include <fstream>
#include <map>

#include <getopt.h>
#include <unistd.h>
#include <libgen.h>
#include <float.h>
// #include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>

#include "aaron_utils.hpp"
#include "sbm_state.hpp"
#include "scf.hpp"
#include "sbm.hpp"

const char gitstatus[] = 
#include "comment.txt"
#include "gitstatus.txt"
;
#include "cmdline.h"


template<typename T,typename V> T up_cast(V x) { return x; }


struct UsageMessage {
};

static void runSBM(const graph :: NetworkInterfaceConvertedToStringWithWeights *g, const int commandLineK, const sbm :: ObjectiveFunction * const obj, const bool initializeToGT, const vector<int> * const groundTruth, const int iterations, const bool algo_gibbs, const bool algo_m3 , const  gengetopt_args_info &args_info) ;


int main(int argc, char **argv) {
	gengetopt_args_info args_info;
	if (cmdline_parser (argc, argv, &args_info) != 0)
		exit(1) ;
	if(args_info.git_version_flag) {
		PP(gitstatus);
		for (int i=0; i<argc; i++) {
			PP(argv[i]);
		}
	}
	if(args_info.inputs_num != 1) {
		cmdline_parser_print_help();
		exit(1);
	}

	const char * edgeListFileName   = args_info.inputs[0];
	// const char * directoryForOutput = args_info.inputs[1];
	PP(edgeListFileName);
	// PP(directoryForOutput);
	PP(args_info.K_arg);
	PP(args_info.directed_flag);
	PP(args_info.weighted_flag);
	PP(args_info.selfloop_flag);
	PP(args_info.seed_arg);
	PP(args_info.iterations_arg);
	PP(args_info.algo_gibbs_arg);
	PP(args_info.algo_m3_arg);
	PP(args_info.initGT_flag);
	PP(args_info.stringIDs_flag);
	PP(args_info.mega_flag);
	if(args_info.GT_vector_given)
		PP(args_info.GT_vector_arg);
	else
		cout << "args_info.GT_vector_arg:<undefined>";
	PP(args_info.model_scf_flag);

	if(args_info.model_scf_flag) {
		unless(args_info.K_arg == 2) {
			cerr << "Usage error: Currently, the stochastic community finding model (uncollapsed) requires -K 2 as and argument. Exiting." << endl;
			exit(1);
		};
	}
	PP(sbm :: ObjectiveFunction_Poisson :: s);
	PP(sbm :: ObjectiveFunction_Poisson :: theta);

	std :: auto_ptr<sbm :: ObjectiveFunction> obj ( args_info.weighted_flag
		? up_cast<sbm :: ObjectiveFunction*>(new sbm :: ObjectiveFunction_Poisson(args_info.selfloop_flag, args_info.directed_flag, args_info.weighted_flag)    )
		: up_cast<sbm :: ObjectiveFunction*>(new sbm :: ObjectiveFunction_Bernoulli(args_info.selfloop_flag, args_info.directed_flag, args_info.weighted_flag)  )
		);
	std :: auto_ptr<graph :: NetworkInterfaceConvertedToStringWithWeights > network;
	// Try the new graph loader
	if(args_info.stringIDs_flag) {
		network = graph :: loading :: make_Network_from_edge_list_string(edgeListFileName, args_info.directed_flag, args_info.weighted_flag);
	} else {
		network = graph :: loading :: make_Network_from_edge_list_int32(edgeListFileName, args_info.directed_flag, args_info.weighted_flag);
	}

	vector<int> groundTruth;
	if(args_info.GT_vector_given) {
		cerr << " --GT.vector: ground truth temporarily not supported" << endl;
		exit(1);
	}
#if 0
	if(args_info.GT_vector_given) { //populate groundTruth based on the --GT.vector option
		const int N = network->numNodes();
		groundTruth.resize(N);
		std :: ifstream GTvectorStream(args_info.GT_vector_arg);
		set<int> nodeNamesInOrder__;
		for(int v=0; v<N; v++) {
			const bool wasInserted = nodeNamesInOrder__.insert( atoi(network->node_name_as_string(v).c_str())).second;
			assert(wasInserted);
		}
		assert((int)nodeNamesInOrder__.size() == N);
		forEach(const int nodeName, amd :: mk_range(nodeNamesInOrder__)) {
			if(!GTvectorStream.is_open())
				break;
			int z_i;
			GTvectorStream >> z_i;
			if(GTvectorStream.eof())
				break;
			const int i = network->node_name_as_string(printfstring("%d", nodeName).c_str());
			groundTruth.at(i) = z_i;
		}
		if((int) groundTruth.size() != g->numNodes()) {
			cerr << "Error: the GT.vector file \"" << args_info.GT_vector_arg << "\" has " << groundTruth.size() << " lines. "
				<< "I was expecting " << g->numNodes() << "lines."
				<< endl;
			groundTruth.clear();
		}
	}
#endif

	srand48(args_info.seed_arg);
	if(args_info.model_scf_flag) {
		gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);
		gsl_rng_set(r, args_info.seed_arg);
		runSCF(network.get(), args_info.K_arg, args_info.initGT_flag, groundTruth.empty() ? NULL : &groundTruth, args_info.iterations_arg, r);
	} else {
		runSBM(network.get(), args_info.K_arg, obj.get(), args_info.initGT_flag, groundTruth.empty() ? NULL : &groundTruth, args_info.iterations_arg, args_info.algo_gibbs_arg, args_info.algo_m3_arg, args_info);
	}
}

void randomize(sbm :: State &s, const int K) { // randomize the partition and have K clusters in it
	assert(s._k <= K);
	assert(K >= 1);
	while(s._k < K) {
		s.appendEmptyCluster();
	}
	cout << "Randomizing.. ";
	for(int n=0; n<s._N; n++) {
		const int newCluster = drand48() * s._k;
		if(newCluster != s.labelling.cluster_id.at(n))
			s.moveNodeAndInformOfEdges2(n, newCluster);
	}
	PP2(s._k, s.labelling.NonEmptyClusters);
	assert(s._k == K);
	s.internalCheck();
}

bool fiftyfifty() {
	if(drand48() < 0.5)
		return true;
	else
		return false;
}

amd :: FormatFlagStack stack;

	AcceptanceRate :: AcceptanceRate(const char * name) : n(0), a(0), _name(name) {
		mostRecent.insert(make_pair(100,0));
		mostRecent.insert(make_pair(1000,0));
		mostRecent.insert(make_pair(10000,0));
		mostRecent.insert(make_pair(100000,0));
	}
	void AcceptanceRate :: notify(bool accepted) {
		this->acc.push_back(accepted);
		this->n++;
		assert(this->n == (int)this->acc.size());
		if(accepted)
			this->a++;
		forEach(typeof(pair<const int,int>) &x, amd :: mk_range(this->mostRecent)) {
			x.second += accepted ? 1 : 0;
			if(this->n > x.first) {
				x.second -= this->acc.at(this->n - x.first - 1) ? 1 : 0;
			}
			assert(x.second >= 0);
			assert(x.second <= x.first);
		}
	}
	void AcceptanceRate :: dump() const {
		cout << "Acceptance Rate " << '"' << this->_name << "\": ";
		 cout << stack.push << fixed << setw(4) << setprecision(1) << 100.0L * double(this->a)/this->n << " %" << stack.pop; // cout << double(this->a)/this->n;
		cout << "\t" << this->a << " / " << this->n;
		cout << endl;
		forEach(const typeof(pair<const int,int>) &x, amd :: mk_range(this->mostRecent)) {
			if(this->n >= x.first) {
				cout << stack.push
					<< '\t' << setw(7) << x.first << " // " << setw(7) << x.second << "\t "
					<< fixed
					<< setw(4)
					<< setprecision(1)
					<< 100.0L * double(x.second)/x.first << " %" << endl
				<< stack.pop ;
			}
		}
	}

static
bool acceptTest(const long double delta, AcceptanceRate *AR = NULL) {
	bool b;
	if(log2(drand48()) < delta)
		b = true;
	else
		b = false;
	if(AR)
		AR->notify(b);
	return b;
}

long double delta_P_z_x__1RowOfBlocks(sbm :: State &s, const sbm :: ObjectiveFunction *obj, const int pre_k, const int t, const int isolatedClusterId, const long double isolatedNodesSelfLoop);

long double M3(sbm :: State &s, const sbm :: ObjectiveFunction *obj, AcceptanceRate * const AR, AcceptanceRate * const AR_alittleConservative, AcceptanceRate * const AR_veryConservative) {
	// 1. Choose two clusters at random

	if(s._k < 2) {
		AR                    ->notify(false);
		AR_alittleConservative->notify(false);
		AR_veryConservative->notify(false);
		return 0.0L; // should this be recorded as a rejection for the purpose of the acceptance rate?
	}
	assert(s._k >= 2);

// #define M3_debugPrinting
// #define M3_Paranoid
#ifdef M3_Paranoid
	const long double preRandom_P_x_z = s.P_edges_given_z_slow(obj);
	const long double preRandom_P_zK_x = s.pmf_slow(obj);
#endif
	// const long double preRandom_P_z_K = s.P_z_orders();
	assert(AR);
	const int pre_k = s._k;
	const int left = drand48() * s._k;
	int right_; do { right_ = drand48() * s._k; } while (left==right_);
	const int right = right_;
#ifdef M3_debugPrinting
	PP3(s._k, left,right);
#endif

	vector<int> randomizedNodeIDs; // the nodes in both clusters. Will be considered in a random order
	const sbm :: Cluster * const lCluster = s.labelling.clusters.at(left);
	const sbm :: Cluster * const rCluster = s.labelling.clusters.at(right);
	forEach(int n, amd :: mk_range(lCluster->get_members())) { randomizedNodeIDs.push_back(n); }
	forEach(int n, amd :: mk_range(rCluster->get_members())) { randomizedNodeIDs.push_back(n); }

	const int M = randomizedNodeIDs.size();
	assert(M == lCluster->order() + rCluster->order());

	random_shuffle(randomizedNodeIDs.begin(), randomizedNodeIDs.end());

	vector<int> clusterIDs;
	for(int m=0; m<M; m++) {
		const int n = randomizedNodeIDs.at(m);
		clusterIDs.push_back(s.labelling.cluster_id.at(n));
	}
	assert((int) clusterIDs.size() == M);

	long double sumToOrig = 0.0L; // The total delta in P(x|z), but the clusters < pre_k
	long double proposalProbOld_log2 = 0.0L;
	for(int m = 0; m<M; m++) {
		const int nToRemove = randomizedNodeIDs.at(m);
		const int origClusterID = clusterIDs.at(m);
		assert(origClusterID == left || origClusterID == right);
		const int isolatedClusterId = s.appendEmptyCluster();
		const int origClusterID_verify = s.moveNodeAndInformOfEdges(nToRemove, isolatedClusterId);
		assert(origClusterID == origClusterID_verify);
		const long double isolatedNodesSelfLoop = obj->relevantWeight(isolatedClusterId, isolatedClusterId, &s._edgeCounts);

			const long double pre_z_K = s.P_z_orders();
			s.moveNode(nToRemove, left);
			const long double post_z_K_left = s.P_z_orders();
			s.moveNode(nToRemove, isolatedClusterId);
			assert(pre_z_K == s.P_z_orders());
			const long double delta_z_K_left = post_z_K_left - pre_z_K;
			assert(isfinite(delta_z_K_left));
			s.moveNode(nToRemove, right);
			const long double post_z_K_right = s.P_z_orders();
			s.moveNode(nToRemove, isolatedClusterId);
			assert(pre_z_K == s.P_z_orders());
			const long double delta_z_K_right = post_z_K_right - pre_z_K;
			assert(isfinite(delta_z_K_right));
			// PP2(delta_z_K_left,delta_z_K_right);

		const long double toTheLeft = delta_z_K_left  + delta_P_z_x__1RowOfBlocks(s, obj, pre_k, left, isolatedClusterId, isolatedNodesSelfLoop);
		const long double toTheRight= delta_z_K_right + delta_P_z_x__1RowOfBlocks(s, obj, pre_k, right, isolatedClusterId, isolatedNodesSelfLoop);
		const long double toOrig = origClusterID == left ? toTheLeft : toTheRight;
		sumToOrig += toOrig;

		// would would the proposal probability have been in reverse, so to speak?
		const long double toTheLeft__  = toTheLeft  - max(toTheLeft, toTheRight);
		const long double toTheRight__ = toTheRight - max(toTheLeft, toTheRight);
		assert(toTheLeft__ <= 0.0L);
		assert(toTheRight__ <= 0.0L);
		const long double toTheLeft____  = exp2l(toTheLeft__);
		const long double toTheRight____ = exp2l(toTheRight__);
#ifdef M3_debugPrinting
		PP2(nToRemove, origClusterID);
		PP2(toTheLeft, toTheRight);
		PP2(toTheLeft__, toTheRight__);
		PP2(toTheLeft____, toTheRight____);
#endif
		const long double totalProb____ = toTheLeft____ + toTheRight____;
		proposalProbOld_log2 +=          toOrig - max(toTheLeft, toTheRight);
		proposalProbOld_log2 -=          log2l(totalProb____);
	}
	assert(s._k == pre_k + M);

	long double delta_P_z_x_forRandom = 0.0L;
	long double proposalProbNew_log2 = 0.0L;
	// randomly assign them to the two clusters
		for(int m = M-1; m>=0; m--) {
			const int n = randomizedNodeIDs.at(m);
			// cout << "randomly assign " << n << endl;

			const int isolatedClusterId = s.labelling.cluster_id.at(n);
			assert(isolatedClusterId == s._k-1);
			const long double isolatedNodesSelfLoop = obj->relevantWeight(isolatedClusterId, isolatedClusterId, &s._edgeCounts);

			const long double pre_z_K = s.P_z_orders();
			s.moveNode(n, left);
			const long double post_z_K_left = s.P_z_orders();
			s.moveNode(n, isolatedClusterId);
			assert(pre_z_K == s.P_z_orders());
			const long double delta_z_K_left = post_z_K_left - pre_z_K;
			assert(isfinite(delta_z_K_left));
			s.moveNode(n, right);
			const long double post_z_K_right = s.P_z_orders();
			s.moveNode(n, isolatedClusterId);
			assert(pre_z_K == s.P_z_orders());
			const long double delta_z_K_right = post_z_K_right - pre_z_K;
			assert(isfinite(delta_z_K_right));
			// PP2(delta_z_K_left,delta_z_K_right);

			const long double toTheLeft = delta_z_K_left  + delta_P_z_x__1RowOfBlocks(s, obj, pre_k, left, isolatedClusterId, isolatedNodesSelfLoop);
			const long double toTheRight= delta_z_K_right + delta_P_z_x__1RowOfBlocks(s, obj, pre_k, right, isolatedClusterId, isolatedNodesSelfLoop);

			// choose left or right randomly, weighted by toTheLeft and toTheRight
			const long double toTheLeft__  = toTheLeft  - max(toTheLeft, toTheRight);
			const long double toTheRight__ = toTheRight - max(toTheLeft, toTheRight);
			assert(toTheLeft__ <= 0.0L);
			assert(toTheRight__ <= 0.0L);
			const long double toTheLeft____  = exp2l(toTheLeft__);
			const long double toTheRight____ = exp2l(toTheRight__);
			const long double totalProb____ = toTheLeft____ + toTheRight____;
			const long double u = drand48() * totalProb____;
			const int randomClusterId = (u < toTheLeft____) ? left : right;
			proposalProbNew_log2 +=          (u < toTheLeft____) ? toTheLeft__ : toTheRight__;
			proposalProbNew_log2 -=          log2l(totalProb____);
#ifdef M3_debugPrinting
			PP2(toTheLeft, toTheRight);
			PP2(toTheLeft__, toTheRight__);
			PP2(toTheLeft____, toTheRight____);
			PP3(toTheLeft____ , u , toTheRight____);
			PP3(left, randomClusterId, right);
			PP(bool(randomClusterId == right));
#endif
			const int isolatedCluster_verify = s.moveNodeAndInformOfEdges(n, randomClusterId);
			s.deleteClusterFromTheEnd();
			assert(s._k == isolatedCluster_verify);

			delta_P_z_x_forRandom += randomClusterId == left ? toTheLeft : toTheRight;
		}
		assert(s._k == pre_k);

		// must assert that the delta_P_x_z is as was expected
		// const long double postRandom_P_z_K = s.P_z_orders();
#ifdef M3_Paranoid
		const long double postRandom_P_zK_x = s.pmf_slow(obj);
		const long double postRandom_P_x_z = s.P_edges_given_z_slow(obj) + s.P_z_orders();
		PP2(preRandom_P_x_z, postRandom_P_x_z);
		PP2(sumToOrig, delta_P_z_x_forRandom);
		PP2(preRandom_P_x_z, postRandom_P_x_z + sumToOrig - delta_P_z_x_forRandom);
		assert(VERYCLOSE(preRandom_P_zK_x + delta_P_z_x_forRandom - sumToOrig, postRandom_P_zK_x));
		// PP2(preRandom_P_zK_x , postRandom_P_zK_x);
		// PP(delta_P_z_x_forRandom - sumToOrig + postRandom_P_z_K - preRandom_P_z_K);
		// assert(VERYCLOSE( -preRandom_P_zK_x + postRandom_P_zK_x , delta_P_x_z_forRandom - sumToOrig + postRandom_P_z_K - preRandom_P_z_K));

		cout << "Did it randomly find a good one?" << endl;
		PP(delta_P_z_x_forRandom - sumToOrig);
		PP2(proposalProbOld_log2, proposalProbNew_log2);
#endif

	const long double delta_P_zK = delta_P_z_x_forRandom - sumToOrig; // + postRandom_P_z_K - preRandom_P_z_K;
#ifdef M3_Paranoid
	assert(VERYCLOSE(delta_P_zK , postRandom_P_zK_x - preRandom_P_zK_x));
#endif
	const long double acceptanceProbability =
		+ delta_P_zK
		- proposalProbNew_log2 + proposalProbOld_log2
		;

	if(acceptTest(acceptanceProbability, AR)) {
		// how many nodes haven't actually moved?
		int changed = 0;
		for(int m = 0; m<M; m++) {
			const int n = randomizedNodeIDs.at(m);
			const int origClusterID = clusterIDs.at(m);
			const int currentClusterID = s.labelling.cluster_id.at(n);
			assert(currentClusterID == left || currentClusterID == right);
			assert(origClusterID == left || origClusterID == right);
			if(currentClusterID != origClusterID) {
				changed ++;
			}
		}
		assert(changed >= 0 && changed <= M);
		AR_alittleConservative->notify(changed>0);
		AR_veryConservative->notify(changed>0 && changed < M);

		assert(pre_k == s._k);
		return delta_P_zK;
	}
	else
	{
	for(int i = M-1; i>=0; i--) { // regardless of whether the nodes are still isolated, or have already been assigned randomly, this will put everything back the way it was.
		const int nToPutBack = randomizedNodeIDs.at(i);
		const int currentClusterID = s.labelling.cluster_id.at(nToPutBack);
		const int origClusterID = clusterIDs.at(i);
		const int theOtherCluster = origClusterID == left ? right : left;
		if(currentClusterID == origClusterID) {
			continue; // this is OK where it is
		} else if(currentClusterID == theOtherCluster) {
			s.moveNodeAndInformOfEdges(nToPutBack, origClusterID);
		} else { // it's current isolated
			assert(s.labelling.clusters.at(currentClusterID)->order()==1);
			assert(currentClusterID >= pre_k);
			s.moveNodeAndInformOfEdges(nToPutBack, origClusterID);
			s.deleteClusterFromTheEnd();
		}
	}
	assert(pre_k == s._k);
	AR_alittleConservative->notify(false);
	AR_veryConservative->notify(false);
	return 0.0L;
	}
}

// static
long double gibbsOneNode(sbm :: State &s, const sbm :: ObjectiveFunction *obj, AcceptanceRate *AR) {
	if(s._k == 1) {
		AR->notify(false);
		return 0.0L;
	}
	const int pre_k = s._k;
	const int n = drand48() * s._N;
	const int origClusterID = s.labelling.cluster_id.at(n);
	const int isolatedClusterId = s._k;
	assert(pre_k == isolatedClusterId);
	s.appendEmptyCluster();
	s.moveNodeAndInformOfEdges(n, isolatedClusterId);

	// the randomly chosen node is now in a little cluster of is own.
	// We proceed to propose to add it to every existing

	// all the blocks involving isolatedClusterId will be destroyed, no matter where it's merged.
// #define gibbsOneNode_Paranoid
#ifdef gibbsOneNode_Paranoid
	long double P_x_z_forIsolated = 0.0L;
	for(int i=0; i<pre_k+1; i++) { // +1 so as to include iso<=>iso (which will only have an effect with selfloops)
		P_x_z_forIsolated += obj->log2OneBlock(obj->relevantWeight(isolatedClusterId, i, &s._edgeCounts) , obj->numberOfPairsInBlock(isolatedClusterId,i, &s.labelling), i==isolatedClusterId);
		if(i<pre_k && obj->directed) { // don't forget the block in the other direction
			P_x_z_forIsolated += obj->log2OneBlock(obj->relevantWeight(i, isolatedClusterId, &s._edgeCounts) , obj->numberOfPairsInBlock(i, isolatedClusterId, &s.labelling), i==isolatedClusterId);
		}
	}
#endif

	vector<long double> delta_P_x_z_IfIMoveIntoClusterT(pre_k);
	vector<long double> delta_P_z_K_IfIMoveIntoClusterT(pre_k);
	const long double isolatedNodesSelfLoop = obj->relevantWeight(isolatedClusterId, isolatedClusterId, &s._edgeCounts);
	for(int t=0; t<pre_k; t++) {
		long double delta_blocks = delta_P_z_x__1RowOfBlocks(s, obj, pre_k, t, isolatedClusterId, isolatedNodesSelfLoop);

		delta_P_x_z_IfIMoveIntoClusterT.at(t) = delta_blocks;
		{ // this following piece of code should be pretty fast. P_z_orders and moveNode are relatively cheap, as no edges are involved.
				const long double pre_z_K = s.P_z_orders();
				s.moveNode(n, t);
				const long double post_z_K = s.P_z_orders();
				s.moveNode(n, isolatedClusterId);
				assert(pre_z_K == s.P_z_orders());
				delta_P_z_K_IfIMoveIntoClusterT.at(t) = post_z_K - pre_z_K;
		}

#ifdef gibbsOneNode_Paranoid
		{ // paranoid verification that the above calculation was correct
				const long double pre_x_z = s.P_edges_given_z_slow(obj);
				s.moveNodeAndInformOfEdges(n, t);
				const long double post_x_z = s.P_edges_given_z_slow(obj);
				// PP2(post_x_z - pre_x_z, delta_blocks);
				// PP (post_x_z - pre_x_z - delta_blocks);
				// PP2(post_x_z - pre_x_z,  delta_blocks - P_x_z_forIsolated);
				assert(VERYCLOSE(post_x_z - pre_x_z , delta_blocks - P_x_z_forIsolated));
				s.moveNodeAndInformOfEdges(n, isolatedClusterId);
				assert(pre_x_z == s.P_edges_given_z_slow(obj));
		}
#endif
	}


	// now to select cluster with probability proportional to exp(delta_P_x_z_IfIMoveIntoClusterT + delta_P_z_K_IfIMoveIntoClusterT)
	vector<long double> combinedDelta;
	for(int t=0; t<pre_k; t++) {
		combinedDelta.push_back(delta_P_x_z_IfIMoveIntoClusterT.at(t) + delta_P_z_K_IfIMoveIntoClusterT.at(t));
		// PP2(t,combinedDelta.at(t));
	}
	const long double maxDelta = * std :: max_element(combinedDelta.begin(), combinedDelta.end());

	long double sumOfTruncatedExpDeltas = 0.0L;
	for(int t=0; t<pre_k; t++) {
		combinedDelta.at(t) -= maxDelta;
		// PP(exp2l(combinedDelta.at(t)));
		sumOfTruncatedExpDeltas += exp2l(combinedDelta.at(t));
	}

	int newCluster = 0;
	long double unif = sumOfTruncatedExpDeltas * drand48();
	while(unif > exp2l(combinedDelta.at(newCluster))) {
		unif -= exp2l(combinedDelta.at(newCluster));
		newCluster ++;
		assert(newCluster < pre_k);
	}
	assert(isfinite(unif));

#ifdef gibbsOneNode_Paranoid
	{
		// going from isolated to back_home will:
		// SUBTRACT P_x_z_forIsolated
		// ADD delta_P_x_z_IfIMoveIntoClusterT.at(newCluster)
		// ADD delta_P_z_K_IfIMoveIntoClusterT.at(newCluster))

		const long double isolated_pmf = s.pmf(obj);
		s.moveNodeAndInformOfEdges(n, newCluster);
		const long double inNewCluster_pmf = s.pmf(obj);
		// PP2(isolated_pmf, inNewCluster_pmf);
		// PP3(P_x_z_forIsolated, delta_P_x_z_IfIMoveIntoClusterT.at(newCluster), delta_P_z_K_IfIMoveIntoClusterT.at(newCluster));
		// PP2(P_x_z_forIsolated, delta_P_x_z_IfIMoveIntoClusterT.at(newCluster) + delta_P_z_K_IfIMoveIntoClusterT.at(newCluster));
		// PP ( - isolated_pmf      + inNewCluster_pmf);
		// PP ( - P_x_z_forIsolated + delta_P_x_z_IfIMoveIntoClusterT.at(newCluster) + delta_P_z_K_IfIMoveIntoClusterT.at(newCluster));
		assert (VERYCLOSE(  - isolated_pmf      + inNewCluster_pmf
			,   - P_x_z_forIsolated + delta_P_x_z_IfIMoveIntoClusterT.at(newCluster) + delta_P_z_K_IfIMoveIntoClusterT.at(newCluster)
			));
		s.moveNodeAndInformOfEdges(n, isolatedClusterId);
	}
#endif

	s.moveNodeAndInformOfEdges(n, newCluster);

	s.deleteClusterFromTheEnd();
	AR->notify(newCluster != origClusterID);
	return + delta_P_x_z_IfIMoveIntoClusterT.at(newCluster) + delta_P_z_K_IfIMoveIntoClusterT.at(newCluster)
	       - delta_P_x_z_IfIMoveIntoClusterT.at(origClusterID) - delta_P_z_K_IfIMoveIntoClusterT.at(origClusterID)
		;
}

long double delta_P_z_x__1RowOfBlocks(sbm :: State &s, const sbm :: ObjectiveFunction *obj, const int pre_k, const int t, const int isolatedClusterId, const long double isolatedNodesSelfLoop) {
		// moving 1 isolated node into one cluster will affect some blocks.
		// We're only interested in the blocks of clusters whose id is < pre_k
		// i.e. there may have been many nodes isolated as part of the M3 move, but we're not interested in the interactions between them.
		// This function will be called from gibbsOneNode() and from M3()

		assert(t >= 0);
		assert(t < pre_k);
		assert(pre_k < s._k); // we currently have some small "temporary" clusters, as part of gibbsOneNode() or M3()
		assert(isolatedClusterId == s._k-1); // for M3, it'll often be bigger, not just equal. And maybe equal to s._k usually.
		assert(s.labelling.clusters.at(isolatedClusterId)->order()==1);

		long double delta_blocks = 0.0L;
		// pretend we're merging into cluster t

		// The easiest thing to consider first is the block for the inside of t, block t<>t.
		// The edges to and from isolatedClusterId will become internal edges of t
		{
			const long double old_weight = obj->relevantWeight(t,t, &s._edgeCounts);
			const long double new_weight = old_weight + s._edgeCounts.read(isolatedClusterId, t)
								  + s._edgeCounts.read(t, isolatedClusterId)
								+ isolatedNodesSelfLoop
								;
			const long int old_p_b = obj->numberOfPairsInBlock(t,t, &s.labelling);
			const long int new_p_b = old_p_b          + s.labelling.clusters.at(t)->order()
					+ ( (obj->directed) ? (s.labelling.clusters.at(t)->order()) : 0)
					+ ( (obj->selfloops) ? (1) : 0)
					;
			const long double old_log2 = obj -> log2OneBlock(old_weight, old_p_b, true);
			const long double new_log2 = obj -> log2OneBlock(new_weight, new_p_b, true);
			const long double delta_1_block = new_log2 - old_log2;
			delta_blocks += delta_1_block;
		}
		for(int j=0; j<pre_k;j++) {
			if(t==j) continue; // the internal t-block has already been handled in the previous few lines.
			// the blocks from  t->j and j->t will be enlarged
			// PP3(isolatedClusterId, t, j);
			// PP(obj->relevantWeight(isolatedClusterId, j, &s._edgeCounts));
			// PP(obj->relevantWeight(j, isolatedClusterId, &s._edgeCounts));

			// isolatedClusterId is a small cluster distinct from t and j. t might be equal to j, but first we'll pretend they're different.
			// Assuming t!=j, there are up to four types of links that must be considered. iso->j, j->iso, iso->t, t->iso. Be careful with directions

			assert(t!=j);
			const int numDirectionsToConsider = (obj->directed) ? 2 : 1;
			for(int direction = 0; direction < numDirectionsToConsider; direction++) {
				const long double old_weight = direction==0 ? obj->relevantWeight(t,j, &s._edgeCounts) : obj->relevantWeight(j,t, &s._edgeCounts);
				// PP3(isolatedClusterId, t,j);
				// PP(direction);
				// if direction==1, then we have to think about j->t, not t->j
				if(direction==1) {
					assert(obj->directed);
				}
				// PP(obj->relevantWeight(isolatedClusterId, j, &s._edgeCounts));
				// PP(obj->relevantWeight(j, isolatedClusterId, &s._edgeCounts));
				const long double new_weight = old_weight + ( (direction==0) ?  obj->relevantWeight(isolatedClusterId, j, &s._edgeCounts) : obj->relevantWeight(j, isolatedClusterId, &s._edgeCounts)) ;

				const long int old_pairs = obj->numberOfPairsInBlock(t,j, &s.labelling);
				// PP2(old_weight, old_pairs);
				const long int new_pairs = old_pairs + s.labelling.clusters.at(j)->order();

				const long double old_log2 = obj -> log2OneBlock(old_weight, old_pairs, false);
				const long double new_log2 = obj -> log2OneBlock(new_weight, new_pairs, false);
				const long double delta_1_block = new_log2 - old_log2;
				delta_blocks += delta_1_block;
			}
		}
		assert(isfinite(delta_blocks));
		return delta_blocks;
}

// static
long double MoneNode(sbm :: State &s, sbm :: ObjectiveFunction *obj, AcceptanceRate *AR) {
	if(s._k == 1)
	       return 0.0L;	// can't move a node unless there exist other clusters
	assert(s._k > 1); // can't move a node unless there exist other clusters
	const int n = drand48() * s._N;
	const int oldClusterID = s.labelling.cluster_id.at(n);
	int newClusterID;
	do {
		newClusterID = drand48() * s._k;
	} while (newClusterID == oldClusterID);
	assert(newClusterID != oldClusterID);
	// PP(oldClusterID);
	// PP(newClusterID);

#ifdef MoneNodeParanoid
	const long double pre = s.pmf(obj);
	const long double pre_x_z = s.P_edges_given_z_slow(obj);
#endif
	const long double pre_z = s.P_z();

	std :: vector < pair< pair<int,int> , pair<long double, long int> > > blocksBefore; // all the blocks that'll be modified
	for(int i=0; i<s._k; i++) {
		for(int j=0; j<s._k; j++) {
			if(!obj->isValidBlock(i,j))
				break;
			assert(obj->directed || j <= i);
			if(i == newClusterID || j==newClusterID
			|| i == oldClusterID || j==oldClusterID) { // this block is of interest
				const long double w = obj->relevantWeight(i,j, &s._edgeCounts);
				const long int pairs = obj->numberOfPairsInBlock(i,j, &s.labelling);
				blocksBefore.push_back (make_pair(  make_pair(i,j), make_pair(w,       pairs     )));
			}
		}
	}

	s.moveNodeAndInformOfEdges(n, newClusterID);

	const long double post_z = s.P_z();
	long double delta_x_z = 0.0L;
	{ // now to calculate those that have changed
		// cout << "  delta x|z" << endl;
		forEach( typeof(pair< pair<int,int> , pair<long double,long int> >) & block, amd :: mk_range(blocksBefore)) {
			const int i = block.first.first;
			const int j = block.first.second;
			// PP2(i,j);
			const long double old_weight = block.second.first;
			const long int old_pairs = block.second.second;
			const long double new_weight = obj->relevantWeight(i,j, &s._edgeCounts);
			const long int new_pairs = obj->numberOfPairsInBlock(i,j, &s.labelling);
			// PP2(new_weight, old_weight);
			// PP2(new_pairs, old_pairs);
			const long double old_x_z = obj->log2OneBlock(old_weight, old_pairs, i==j);
			const long double new_x_z = obj->log2OneBlock(new_weight, new_pairs, i==j);
			// PP3(old_x_z, new_x_z, new_x_z-old_x_z);
			delta_x_z += new_x_z-old_x_z;
		}
		// PP(delta_x_z);
		// cout << " ~delta x|z" << endl;
	}
	const long double delta_z = post_z - pre_z;
	const long double delta_ = delta_z + delta_x_z;
#ifdef MoneNodeParanoid
	{
		const long double post_x_z = s.P_edges_given_z_slow(obj);
		const long double delta_x_z_ = post_x_z - pre_x_z;
		// PP3(pre_x_z, post_x_z, post_x_z - pre_x_z);
		// PP3(delta_x_z, delta_x_z_, delta_x_z - delta_x_z_);
		assert(VERYCLOSE(delta_x_z , delta_x_z_));

		const long double post = s.pmf(obj);
		const long double delta = post - pre;
		assert(VERYCLOSE(delta, delta_));
		// PP(pre);
		// PP(post);
		// PP(delta);
	}
#endif
	if(acceptTest(delta_, AR)) {
		// cout << " + ";
		return delta_;
	} else {
		// cout << "   ";
		s.moveNodeAndInformOfEdges(n, oldClusterID);
		// assert(VERYCLOSE(s.pmf(obj), pre)); // make sure it has undone it properly
		return 0.0L;
	}
}
#if 0
struct OneChoice {
	const long double delta2;
	const long double delta3;
	const long double delta4;
	OneChoice(long double _d2, long double _d3, long double _d4) : delta2(_d2), delta3(_d3), delta4(_d4) {
	}
	long double deltaSum() const {
		return this->delta2+this->delta3+this->delta4;
	}
};
struct TwoChoices {
	const OneChoice left;
	const OneChoice right;
	const long double left_deltaSum;
	const long double right_deltaSum;
	long double Pleft;
	long double Pright;
	TwoChoices(const OneChoice &_l, const OneChoice &_r) : left(_l), right(_r), left_deltaSum(this->left.deltaSum()), right_deltaSum(this->right.deltaSum()) {
		const long double log2LeftOverRight = left_deltaSum - right_deltaSum;
		const long double LeftOverRight = exp2l(log2LeftOverRight);
		assert(isfinite(LeftOverRight));
		Pright = 1.0L / (LeftOverRight + 1.0L);
		Pleft = 1.0L - Pright;
		// PP2(left_deltaSum, right_deltaSum);
		// PP2(Pleft, Pright);
		if(Pleft == 0.0L)
			Pleft = DBL_MIN;
		if(Pright == 0.0L)
			Pright = DBL_MIN;
		assert(Pleft > 0.0L);
		assert(Pleft <= 1.0L);
		assert(Pright > 0.0L);
		assert(Pright <= 1.0L);
		assert(isfinite(Pleft));
		assert(isfinite(Pright));
	}
};
static OneChoice M3_oneNode(sbm :: State &s, const int n, const int candCluster) {
	// const long double pre = s.pmf(); // TODO remove this, and the corresponding check at the end
	// given an isolated node, and a candidate cluster to join, what's the delta in the fitness?
	
	const int isolatedClusterID = s.labelling.cluster_id.at(n);
	const sbm :: Cluster * isoCl = s.labelling.clusters.at(isolatedClusterID);
	assert(isoCl->order() == 1);
	const sbm :: Cluster * candCl = s.labelling.clusters.at(candCluster);
	assert(isolatedClusterID != candCluster);
	const int old_order = candCl->order();
	const int new_order = old_order+1;
	const int old_NonEmpty = s.labelling.NonEmptyClusters;
	const int new_NonEmpty = old_order==0 ? s.labelling.NonEmptyClusters : (s.labelling.NonEmptyClusters-1);
	const long double preSumOfLog2l = s.labelling.SumOfLog2LOrders;
	const long double postReMergeSumOfLog2l = preSumOfLog2l + log2l(new_order) - ( old_order>1 ? log2l(old_order) : 0.0L);
	// PP(postReMergeSumOfLog2l);


	long double delta2 = 0.0L;
	long double delta3 = 0.0L;
	// long double delta4 = 0.0L;

	if(old_order > 1)
		delta2 -= LOG2FACT(old_order);
	delta2 += LOG2FACT(new_order);

	delta3 += postReMergeSumOfLog2l * (new_NonEmpty-1);
	delta3 -= preSumOfLog2l  * (old_NonEmpty-1);
	if(new_order >= 2)
			delta3 += log2l((new_order * new_order - new_order)/2);
	if(old_order >= 2)
			delta3 -= log2l((old_order * old_order - old_order)/2);
	delta3 = -delta3;

	long double delta4 = 0.0L;
	delta4 -= s.P_edges_given_z_correction_JustOneCluster(candCluster);
	delta4 -= s.P_edges_given_z_correction_JustOneCluster(isolatedClusterID);
	// Any edges between isolatedClusterID and candCluster (i.e. old internal edges) will have been double counted.
	const int doubleCounted_edges = s._edgeCounts.get(isolatedClusterID, candCluster);
	delta4 -= M_LOG2E * gsl_sf_lnchoose(old_order/*doubleCounted_pairs*/, doubleCounted_edges);

	s.moveNodeAndInformOfEdges(n, candCluster); // move the node temporarily, leaving an empty cluster behind, it'll be moved back again in a few lines

	delta4 += s.P_edges_given_z_correction_JustOneCluster(candCluster);

	// PP(delta2);
	// PP(delta3);
	// PP(delta4);
	
	// const long double postTempMove = pre + delta2 + delta3 + delta4;
	// assert(VERYCLOSE(postTempMove , s.pmf()));

	s.moveNodeAndInformOfEdges(n, isolatedClusterID); // move the node back again

	// assert(VERYCLOSE(pre, s.pmf())); // to ensure that we've undone the change

	return OneChoice(delta2 , delta3 , delta4);
}
void M3_old(sbm :: State &s) {
	const bool verbose = false;
	// cout << endl << "     ========== M3 =========" << endl;
	const long double preM3   = s.pmf();
	const long double preM3_1 = s.P_z_K();
	const long double preM3_2 = s.P_z_orders();
	const long double preM3_3 = s.P_edges_given_z_baseline();
	const long double preM3_4 = s.P_edges_given_z_correction();
	const int preM3_k = s._k;
	// as per Nobile & Fearnside's allocation sampler
	// 1. Choose two clusters at random
	// 2. Decide the random order in which the nodes are to be considered.
	// 2. Calculate the proposal probability of the status quo
	// 3. Verify 
	if(s._k < 2)
		return;
	const int cl1 = drand48() * s._k;
	const int cl2 = drand48() * s._k;
	if(cl1 == cl2)
		return;
	if(verbose) cout << endl << "     ========== M3 (found two clusters) =========" << endl;
	const sbm :: Cluster * CL1 = s.labelling.clusters.at(cl1);
	const sbm :: Cluster * CL2 = s.labelling.clusters.at(cl2);
	vector<int> allNodes;
	std :: tr1 :: unordered_map<int, int> statusQuoClustering;
	allNodes.insert(allNodes.end(), CL1->members.begin(), CL1->members.end());
	allNodes.insert(allNodes.end(), CL2->members.begin(), CL2->members.end());
	// PP(CL1->members.size());
	// PP(CL2->members.size());
	// PP(allNodes.size());
	assert(CL1->members.size() + CL2->members.size() == allNodes.size());
	forEach(int x, amd :: mk_range(CL1->members)) { statusQuoClustering[x] = cl1; }
	forEach(int y, amd :: mk_range(CL2->members)) { statusQuoClustering[y] = cl2; }
	// forEach(int x, amd :: mk_range(CL1->members)) { PPt(x); } cout << endl;
	// forEach(int y, amd :: mk_range(CL2->members)) { PPt(y); } cout << endl;
	random_shuffle(allNodes.begin(), allNodes.end());
	// forEach(int z2, amd :: mk_range(allNodes    )) { PPt(z2); } cout << endl;
	long double deltaSumOfTheStatusQuo = 0.0L;
	long double log2ProductOfProposalProbabilitiesForStatusQuo = 0.0;
	for(vector<int>:: const_reverse_iterator remover = allNodes.rbegin(); remover != allNodes.rend(); ++remover) {
		const int node_to_remove = *remover;
		// PP(node_to_remove);
		// const long double pre = s.pmf();
		// const long double pre1 = s.P_z_K();
		const long double pre2 = s.P_z_orders();
		const long double pre3 = s.P_edges_given_z_baseline();
		// const long double pre4 = s.P_edges_given_z_correction();
		// assert(pre == pre1+pre2+pre3+pre4);
		const long double preSumOfLog2l = s.labelling.SumOfLog2LOrders;
		// PP(preSumOfLog2l);
		const long double preNonEmpty = s.labelling.NonEmptyClusters;
		// assert(pre == pre1 + pre2 + pre3 + pre4);

		const int old_clusterID = s.labelling.cluster_id.at(node_to_remove);
		const sbm :: Cluster * old_cluster = s.labelling.clusters.at(old_clusterID);
		const int old_order = old_cluster->order();
		assert(old_order>=1);
		long double delta2 = -LOG2FACT(old_order);
		long double delta4 = - s.P_edges_given_z_correction_JustOneCluster(old_clusterID); // this must be remembered before the move

		const int tempClusterID = s.isolateNode(node_to_remove);

		const long double postSumOfLog2l = s.labelling.SumOfLog2LOrders;
		const long double postNonEmpty = s.labelling.NonEmptyClusters;
		const int new_order = old_order-1;
		if(new_order>1)
			delta2 += LOG2FACT(new_order);
		assert(new_order >= 0);
		const long double post2 = s.P_z_orders();
		assert(VERYCLOSE(post2, pre2 + delta2));

		const long double post3 = s.P_edges_given_z_baseline();
		// PP(post3 - pre3);
		long double delta3 = 0.0;
		// PP2(old_order, new_order);
		delta3 += postSumOfLog2l * (postNonEmpty-1);
		delta3 -= preSumOfLog2l  * (preNonEmpty-1);
		if(new_order >= 2)
			delta3 += log2l((new_order * new_order - new_order)/2);
		if(old_order >= 2)
			delta3 -= log2l((old_order * old_order - old_order)/2);
		delta3 = -delta3;
		// PP(delta3);
		assert(VERYCLOSE(delta3 , post3 - pre3));

		// const long double post4 = s.P_edges_given_z_correction();
		delta4 += s.P_edges_given_z_correction_JustOneCluster(old_clusterID);
		delta4 += s.P_edges_given_z_correction_JustOneCluster(tempClusterID);
		// Any edges between tempClusterID and old_clusterID (i.e. old internal edges) will have been double counted.
		const int doubleCounted_edges = s._edgeCounts.get(tempClusterID, old_clusterID);
		delta4 += M_LOG2E * gsl_sf_lnchoose(new_order/*doubleCounted_pairs*/, doubleCounted_edges);

		// PP2(post4,pre4);
		// PP(post4 - pre4);
		// PP(delta4);
		// assert(post4 - pre4 == delta4);

		// const long double post1 = s.P_z_K();
		// long double delta1 = post1 - pre1;

		// PP(delta1);
		// PP(delta2);
		// PP(delta3);
		// PP(delta4);

		// cout << " ==  M3_oneNode ==" << endl;
		TwoChoices two_choices(M3_oneNode(s, node_to_remove, cl1),M3_oneNode(s, node_to_remove, cl2));
		const long double left  = two_choices.left_deltaSum;
		const long double right = two_choices.right_deltaSum;
		assert(cl1 == old_clusterID || cl2 == old_clusterID);
		const long double statusQuo = cl1 == old_clusterID ? left : right;
		// PP2(cl1, cl2);
		// PP (old_clusterID);
		// PP2(left, right);
		// PP2(two_choices.Pleft, two_choices.Pright);
		// cout << " == ~M3_oneNode ==" << endl;
		assert(VERYCLOSE(-statusQuo , delta2 + delta3 + delta4));
		// const long double post = s.pmf();
		// assert(VERYCLOSE(post, post1+post2+post3+post4));
		// assert(VERYCLOSE(post, pre + delta1 + delta2 + delta3 + delta4));

		// PP2(two_choices.Pleft,two_choices.Pright);
		// PP( (cl1 == old_clusterID) );
		const long double prpsl = cl1 == old_clusterID ? two_choices.Pleft : two_choices.Pright;
		// PP(prpsl);
#define assertFinite(x) assert(isfinite(x))
		assertFinite(prpsl);
		// PP2(prpsl , log2(prpsl));
		// PP2(log2ProductOfProposalProbabilitiesForStatusQuo , log2(prpsl));
		log2ProductOfProposalProbabilitiesForStatusQuo += log2(prpsl);
		// PP(log2ProductOfProposalProbabilitiesForStatusQuo);
		assertFinite(log2ProductOfProposalProbabilitiesForStatusQuo);

		deltaSumOfTheStatusQuo += cl1 == old_clusterID ? two_choices.left_deltaSum : two_choices.right_deltaSum;
	}

	const long double midM3 = s.pmf();
	const long double midM3_1 = s.P_z_K();
	// const long double midM3_2 = s.P_z_orders();
	// const long double midM3_3 = s.P_edges_given_z_baseline();
	// const long double midM3_4 = s.P_edges_given_z_correction();
	// PP2(-deltaSumOfTheStatusQuo, midM3 - midM3_1 - preM3 + preM3_1);
	assert(VERYCLOSE(-deltaSumOfTheStatusQuo , midM3 - midM3_1 - preM3 + preM3_1));

	long double log2ProductOfProposalProbabilitiesForNewProposal = 0.0L;
	long double deltaSumOfTheNewProposal = 0.0L;
	bool IsRandomProposalIdenticalToStatusQuo = true;
	{ // random proposal
		// cout << endl << "  random proposals for M3" << endl << endl;
		for(vector<int>:: const_iterator adder = allNodes.begin(); adder != allNodes.end(); ++adder) {
			// cout << endl << "  random proposal for M3" << endl << endl;
			// const long double preM3OneRandom = s.pmf();
			const int node_to_Add = *adder;
			// PP(node_to_Add);
			const int clID = s.labelling.cluster_id.at(node_to_Add);
			assert(clID + 1 == s._k);
			const sbm :: Cluster * clIsolated = s.labelling.clusters.at(clID);
			assert(clIsolated->order()==1);
			assert(clIsolated->members.front()==node_to_Add);
			// which of the two to add to?
			// PP2(left,right);
			assert(cl1 != clID && cl2 != clID);
			assert(cl1 != cl2);
			const TwoChoices two_choices(M3_oneNode(s, node_to_Add, cl1),M3_oneNode(s, node_to_Add, cl2));
			const long double left = two_choices.left_deltaSum;
			const long double right = two_choices.right_deltaSum;
			assert(VERYCLOSE(left  , M3_oneNode(s, node_to_Add, cl1).deltaSum()));
			assert(VERYCLOSE(right , M3_oneNode(s, node_to_Add, cl2).deltaSum()));
			// PP2(left , two_choices.left.deltaSum());
			// PP2(two_choices.left.deltaSum(),two_choices.right.deltaSum());
			// PP2(two_choices.Pleft,two_choices.Pright);
			long double prpsl;
			enum LeftOrRight { Left, Right };
			LeftOrRight lr = drand48() < two_choices.Pright ? Right : Left;
			// LeftOrRight lr = (statusQuoClustering.at(node_to_Add)==cl2) ? Right : Left; // clone the status quo clustering
			if(lr == Right) {
				// cout << " go right" << endl;
				s.moveNodeAndInformOfEdges(node_to_Add, cl2);
				if(statusQuoClustering.at(node_to_Add)!=cl2)
					IsRandomProposalIdenticalToStatusQuo = false;
				// assert(VERYCLOSE(s.pmf() , preM3OneRandom + two_choices.right_deltaSum));
				// the above assert is (correctly) using an 'out-of-date' value of _k. Hence we don't delete this temporary (now empty) cluster until the next line
				s.deleteClusterFromTheEnd();
				prpsl = two_choices.Pright;
				deltaSumOfTheNewProposal += two_choices.right.deltaSum();
			} else {
				// cout << " go left" << endl;
				s.moveNodeAndInformOfEdges(node_to_Add, cl1);
				if(statusQuoClustering.at(node_to_Add)!=cl1)
					IsRandomProposalIdenticalToStatusQuo = false;
				// assert(VERYCLOSE(s.pmf() , preM3OneRandom + two_choices.left_deltaSum));
				// the above assert is (correctly) using an 'out-of-date' value of _k. Hence we don't delete this temporary (now empty) cluster until the next line
				s.deleteClusterFromTheEnd();
				prpsl = two_choices.Pleft;
				deltaSumOfTheNewProposal += two_choices.left.deltaSum();
			}

			assertFinite(prpsl);
			log2ProductOfProposalProbabilitiesForNewProposal += log2(prpsl);
			assertFinite(log2ProductOfProposalProbabilitiesForNewProposal);

			// cout << endl << " ~random proposal for M3" << endl << endl;
		}
		assert(VERYCLOSE(midM3 + deltaSumOfTheNewProposal + preM3_1 - midM3_1 , s.pmf()));
		// cout << endl << " ~random proposals for M3" << endl << endl;
	}
	const long double pmfOfTheNewProposal = midM3 - midM3_1 + preM3_1 + deltaSumOfTheNewProposal;
	assert(VERYCLOSE(pmfOfTheNewProposal, s.pmf()));

	if(verbose) cout << "Accept or Reject? Not yet implemented." << endl;
	long double acceptanceLog2 = pmfOfTheNewProposal - preM3 - log2ProductOfProposalProbabilitiesForNewProposal + log2ProductOfProposalProbabilitiesForStatusQuo;
	if(VERYCLOSE(acceptanceLog2, 0.0L))
		acceptanceLog2 = 0.0L;
	// PP2(log2ProductOfProposalProbabilitiesForStatusQuo, log2ProductOfProposalProbabilitiesForNewProposal);
	// PP2(preM3, pmfOfTheNewProposal);
	// PP(acceptanceLog2);
	// PP(IsRandomProposalIdenticalToStatusQuo);
	// assert(VERYCLOSE(log2ProductOfProposalProbabilitiesForStatusQuo , log2ProductOfProposalProbabilitiesForNewProposal)); // only true if the proposal is for no change

	assert(preM3_k == s._k);

	// Now we either accept it with probability exp2l(acceptanceLog2), or reject it

	if(log2(drand48()) < acceptanceLog2 /*|| IsRandomProposalIdenticalToStatusQuo*/) {
		// accepting.
		// if(verbose)
			cout << "M3 ACCEPT: " << acceptanceLog2 << endl;
	} else {
		// reject. let's put them all back
		for(vector<int>:: const_iterator reAdder = allNodes.begin(); reAdder != allNodes.end(); ++reAdder) {
			// these should all be isolated nodes, and at the end of the list of clusters
			const int node_to_reAdd = *reAdder;
			// PP(node_to_reAdd);
			// const int clID = s.cluster_id.at(node_to_reAdd);
			s.isolateNode(node_to_reAdd);
			s.moveNodeAndInformOfEdges(node_to_reAdd, statusQuoClustering.at(node_to_reAdd));
			s.deleteClusterFromTheEnd();
		}
		assert(preM3_k == s._k);
		const long double postM3 = s.pmf();
		const long double postM3_1 = s.P_z_K();
		const long double postM3_2 = s.P_z_orders();
		const long double postM3_3 = s.P_edges_given_z_baseline();
		const long double postM3_4 = s.P_edges_given_z_correction();
		// PP2(preM3, postM3);
		assert(preM3_1 == postM3_1);
		assert(preM3_2 == postM3_2);
		assert(VERYCLOSE(preM3_3, postM3_3));
		assert(preM3_4 == postM3_4);
		assert(VERYCLOSE(preM3 , postM3));
	}
	if(verbose) cout << "     ========== ~M3 =========" << endl;
}
#endif
static long double MetropolisOnK(sbm :: State &s, const sbm :: ObjectiveFunction *obj, AcceptanceRate *AR) {
	/// const long double prePMF = s.pmf(obj);
	/// const long double prePMF12 = s.P_z_K();
	const int preK = s._k;
	if(fiftyfifty()) { // propose increase in K
		s.appendEmptyCluster();
		// const long double postPMF12 = s.P_z_K();
		// assert(VERYCLOSE(s.pmf(), prePMF - prePMF12 + postPMF12));
		// const long double postPMF = prePMF - prePMF12 + postPMF12; 
		// assert(VERYCLOSE(s.pmf(), postPMF));
		// assert(postPMF < prePMF);
		// assert(postPMF12 < prePMF12);
		// assert(VERYCLOSE(postPMF - prePMF, postPMF12 - prePMF12));
		const int postK = preK + 1;
		assert(s._k == postK);
		const long double presumed_delta =
			- log2(preK+1) // Poisson(1) prior on K
			// - 1 // Geometric(0.5) prior on K

			+ LOG2GAMMA(s._alpha) // the change to SumOfLog2Facts

			// log2(preK) - log2(s._N+preK)
			+(LOG2GAMMA(postK * s._alpha) - LOG2GAMMA(postK * s._alpha + s._N) - postK*LOG2GAMMA(s._alpha) )
			-(LOG2GAMMA(preK * s._alpha) - LOG2GAMMA(preK  * s._alpha + s._N) - preK*LOG2GAMMA(s._alpha) )
			;
		if(acceptTest(presumed_delta, AR)) {
			// cout << "k: acc inc" << endl;
			assert(s._k>preK);
			// PP(s.pmf(obj) - prePMF);
			// PP2(preK, s._k);
			// PP(s._N);
			// PP( presumed_delta );
			// PP(postPMF12 - prePMF12);
			// assert(VERYCLOSE(presumed_delta, s.pmf(obj) - prePMF));
			return presumed_delta;
		} else {
			// cout << "k: rej inc" << endl;
			s.deleteClusterFromTheEnd();
			// assert(s.pmf(obj)==prePMF);
			// assert(s.P_z_K()==prePMF12);
			assert(s._k==preK);
			return 0.0L;
		}
	} else { // propose decrease
		if(s._k == 1) { // can't propose a decrease
			assert(s._k==preK);
			AR->notify(false);
			return 0.0L;
		}
		assert(s._k > 1);
		const int clusterToProposeDelete = drand48() * s._k;
		if(s.labelling.clusters.at(clusterToProposeDelete)->order()>0) {
			// can't delete this. Time to get out
			assert(s._k==preK);
			AR->notify(false);
			return 0.0L;
		}
		if(clusterToProposeDelete != s._k-1) {
			cout << "swapped during deletion"; PP2(clusterToProposeDelete , s._k-1);
			s.swapClusters(clusterToProposeDelete , s._k-1);
		}
		assert(s.labelling.clusters.back()->order()==0);
		{
			const int postK = preK-1;
			const long double presumed_delta =
				+ log2(preK) // Poisson(1) prior on K
				// + 1 // Geometric(0.5) prior on K

				- LOG2GAMMA(s._alpha) // the change to SumOfLog2Facts

				// - log2(preK-1) + log2(s._N+preK-1)
				+(LOG2GAMMA(postK * s._alpha) - LOG2GAMMA(postK * s._alpha + s._N) - postK*LOG2GAMMA(s._alpha) )
				-(LOG2GAMMA(preK * s._alpha) - LOG2GAMMA(preK  * s._alpha + s._N) - preK*LOG2GAMMA(s._alpha) )
				;
			s.deleteClusterFromTheEnd();
			// const long double postPMF12 = s.P_z_K();
			// assert(VERYCLOSE(s.pmf(obj), prePMF - prePMF12 + postPMF12));
			// const long double postPMF = prePMF - prePMF12 + postPMF12; 
			// assert(VERYCLOSE(s.pmf(obj), postPMF));
			// assert(postPMF > prePMF);
			/// PP2(postPMF12 - prePMF12, presumed_delta);
			if(acceptTest(presumed_delta, AR)) {
				// cout << "k: acc dec" << endl;
				assert(s._k<preK);
				// assert(VERYCLOSE(s.pmf(obj) - prePMF, presumed_delta));
				return presumed_delta;
			} else {
				assert(1==2); // it'll always like decreases, except maybe if the prior on K is an increasing function. // i.e. presumed_delta > 0
				// cout << "k: rej dec" << endl;
				s.appendEmptyCluster();
				// assert(s.pmf(obj)==postPMF);
				// assert(s.P_z_K()==prePMF12);
				assert(s._k==preK);
				return 0.0L;
			}
		}
	}
}

struct CountSharedCluster { // for each *pair* of nodes, how often they share the same cluster
	const int N;
	int denominator;
	vector<int> shared;
	CountSharedCluster(const int _N) : N(_N), denominator(0), shared(_N*_N) {}
	void consume(const vector<int> & cluster_id) {
		assert((int)cluster_id.size() == N);
		++denominator;
		for(int n=0; n<this->N; n++) {
			for(int m=n; m<this->N; m++) {
				if(cluster_id.at(n) == cluster_id.at(m)) {
					++ shared.at(n*this->N+m);
					if(n!=m)
					++ shared.at(m*this->N+n);
				}
			}
		}
	}
	void dump(const sbm :: State & s) {
		if(this->denominator == 0)
			return;
		assert((int)s.nodeNamesInOrder.size() == this->N);
		cout << "  CountSharedCluster(" << this->denominator << ")" << endl;
		//for(int n=0; n<N; n++)
		forEach(const typeof(pair<string,int>) & node_name, amd :: mk_range(s.nodeNamesInOrder))
		{
			const int n = node_name.second;
			//for(int m=0; m<N; m++)
			forEach(const typeof(pair<string,int>) & node_name2, amd :: mk_range(s.nodeNamesInOrder))
			{
				const int m = node_name2.second;
				if(n==m)
					cout << "    ";
				else
					cout << stack.push << setw(4) << int(100.0 * double(shared.at(n*N+m)) / denominator + 0.5) << stack.pop;
			}
			const int my_cluster_id = s.labelling.cluster_id.at(node_name.second);
			cout
				<< "\t" << setw(1 + my_cluster_id) << my_cluster_id
				<< setw(s._k - my_cluster_id) << ""
				<< " \"" << node_name.first << "\"";
			cout << endl;
		}
	}
};

#define CHECK_PMF_TRACKER(track, actual) do { const long double _actual = (actual); long double & _track = (track); if(VERYCLOSE(_track,_actual)) { track = _actual; } assert(_track == _actual); } while(0)

static void runSBM(const graph :: NetworkInterfaceConvertedToStringWithWeights *g, const int commandLineK, const sbm :: ObjectiveFunction * const obj, const bool initializeToGT, const vector<int> * const groundTruth, const int iterations, const bool algo_gibbs, const bool algo_m3 , const  gengetopt_args_info &args_info) {
	PP2(g->numNodes(), g->numRels());
	sbm :: State s(g, !args_info.stringIDs_flag, args_info.mega_flag);

	s.shortSummary(obj, groundTruth); /*s.summarizeEdgeCounts();*/ s.blockDetail(obj);
	s.internalCheck();

	CountSharedCluster *count_shared_cluster;
	if(args_info.mega_flag) // the network is probably too big, don't try to CountSharedCluster
		count_shared_cluster	= 0;
	else
		count_shared_cluster	= new CountSharedCluster(g->numNodes());

	/*
	s.isolateNode(0);
	s.isolateNode(1); // to bring us up to three clusters

	s.shortSummary(); s.summarizeEdgeCounts(); s.blockDetail(obj);
	s.internalCheck();
	PP(s.pmf());

	*/
	if(groundTruth && initializeToGT) {
		for(int v = 0; v < g->numNodes(); v++ ) {
			const int z_i = groundTruth->at(v);
			while(z_i >= s._k) {
				s.appendEmptyCluster();
			}
			s.moveNodeAndInformOfEdges2(v, z_i);
		}
	} else {
		if(commandLineK != -1)
			randomize(s, commandLineK);
		else {
			randomize(s, 2);
			// for(int n=0; n+1<s._N; n++) s.isolateNode(n);
			// assert(s._k == s._N);
		}
	}
	s.shortSummary(obj, groundTruth); /*s.summarizeEdgeCounts();*/ s.blockDetail(obj);
	s.internalCheck();

	long double pmf_track = s.pmf(obj);
	PP(pmf_track);

	AcceptanceRate AR_metroK("metroK");
	AcceptanceRate AR_metro1Node("metro1Node");
	AcceptanceRate AR_gibbs("gibbs");
	AcceptanceRate AR_M3("M3");
	AcceptanceRate AR_M3little("M3lConservative");
	AcceptanceRate AR_M3very  ("M3vConservative");
	for(int i=1; i<=iterations; i++) {
		if(0) {
			if(s._k > 1) // && drand48() < 0.01)
			{
				const int cl1 = s._k * drand48();
				const int cl2 = s._k * drand48();
				if(cl1 != cl2) {
					s.swapClusters(cl1,cl2);
				}
			}
		}
		switch( int(drand48() * 3) ) {
			break; case 0:
				if(commandLineK == -1) {
					pmf_track += MetropolisOnK(s, obj, &AR_metroK);
				} else
					assert(commandLineK == s._k);
			break; case 1:
				if(algo_gibbs)
					pmf_track += gibbsOneNode(s, obj, &AR_gibbs);
			break; case 2:
				if(algo_m3)
					pmf_track += M3(s, obj, &AR_M3, &AR_M3little, &AR_M3very);
		}
		if(i > 30000)
		if(count_shared_cluster)
			count_shared_cluster->consume(s.labelling.cluster_id);

		// PP(i);
		// const long double pre = s.pmf(obj);
		// pmf_track += MoneNode(s, obj, &AR_metro1Node);
		// const long double post = s.pmf(obj);
		// assert(pre + delta == post);
		// if(i%50 == 0)
			// M3(s);
	
		// PP(s.pmf());
		// cout << endl;
		// if(i%1000 == 0)
			// cerr << "i:" << i << endl;
		if(i % args_info.printEveryNIters_arg == 0) {
			cout << endl;
			PP(i);
			s.shortSummary(obj, groundTruth);
			// s.summarizeEdgeCounts();
			AR_metroK.dump();
			AR_metro1Node.dump();
			AR_gibbs.dump();
			AR_M3.dump();
			AR_M3little.dump();
			AR_M3very.dump();
			s.blockDetail(obj);
			if(count_shared_cluster)
				count_shared_cluster->dump(s);
			cout << " end of check at i==" << i << endl;
			CHECK_PMF_TRACKER(pmf_track, s.pmf(obj));
			s.internalCheck();
		}
	}
	s.shortSummary(obj, groundTruth); /*s.summarizeEdgeCounts();*/ s.blockDetail(obj);
	s.internalCheck();
}
