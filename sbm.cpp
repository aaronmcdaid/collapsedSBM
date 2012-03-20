#include "graph/network.hpp"
#include "graph/loading.hpp"
using namespace std;
#include <algorithm>
#include <iomanip>
#include <vector>
#include <fstream>
#include <map>
#include <ctime>

#include <getopt.h>
#include <unistd.h>
#include <libgen.h>
#include <float.h>
// #include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "Range.hpp"
#include "format_flag_stack/format_flag_stack.hpp"
#include "macros.hpp"
#include "sbm_state.hpp"
#include "scf.hpp"
#include "sbm.hpp"

const char gitstatus[] = 
#include "comment.txt"
#include "gitstatus.txt"
;
#include "cmdline.h"


template<typename T,typename V> T up_cast(V x) { return x; }

using namespace std;

struct UsageMessage {
};

static void runSBM(const graph :: NetworkInterfaceConvertedToStringWithWeights *g, const int commandLineK, const sbm :: ObjectiveFunction * const obj, const bool initializeToGT, const vector<int> * const groundTruth, const int iterations, const bool algo_gibbs, const bool algo_m3 , const  gengetopt_args_info &args_info, gsl_rng *r) ;

static void runCEM(const graph :: NetworkInterfaceConvertedToStringWithWeights *g, const int commandLineK
		, const sbm :: ObjectiveFunction * const obj
		, const vector<int> * const groundTruth, const int iterations, const  gengetopt_args_info &args_info, gsl_rng *r) ;
static long double my_likelihood(const int n, sbm :: State &s, const sbm :: ObjectiveFunction *obj, int k = -1);

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
	PP(args_info.algo_metroK_arg);
	PP(args_info.algo_gibbs_arg);
	PP(args_info.algo_1node_arg);
	PP(args_info.algo_m3_arg);
	PP(args_info.algo_ejectabsorb_arg);
	PP(args_info.initGT_flag);
	PP(args_info.stringIDs_flag);
	PP(args_info.mega_flag);
	PP(args_info.alpha_arg);
	if(args_info.GT_vector_given)
		PP(args_info.GT_vector_arg);
	else
		cout << "args_info.GT_vector_arg:<undefined>";
	PP(args_info.model_scf_flag);
	PP(args_info.assume_N_nodes_arg);
	PP(args_info.save_z_arg);
	PP(args_info.algo_sbm_cem_flag);
	PP(args_info.latentspace_flag);
	PP(args_info.lsalpha_arg);
	PP(args_info.algo_lspos_arg);
	PP(args_info.algo_lsm3_arg);
	//PP(args_info.gamma_s_arg);
	//PP(args_info.gamma_phi_arg);
	sbm :: ObjectiveFunction_Poisson :: s     = args_info.gamma_s_arg;
	sbm :: ObjectiveFunction_Poisson :: theta = args_info.gamma_phi_arg;
	sbm :: ls_alpha_k = args_info.lsalpha_arg;


	if(args_info.model_scf_flag) {
		unless(args_info.K_arg == 2) {
			cerr << "Usage error: Currently, the stochastic community finding model (uncollapsed) requires -K 2 as an argument. Exiting." << endl;
			exit(1);
		};
	}
	if(args_info.algo_sbm_cem_flag) {
		unless(args_info.K_arg >= 2) {
			cerr << "Usage error: To use the CEM algorithm (--algo.sbm.cem), you must specify a number of communities (-K) as an argument. Exiting." << endl;
			exit(1);
		};
	}
	if(args_info.assume_N_nodes_arg > 0 && args_info.stringIDs_flag) {
		cerr << endl << "Usage error: --stringIDs and --assume_N_nodes are not allowed together. Exiting." << endl;
		exit(1);
	}
	if(args_info.latentspace_flag && args_info.K_arg == -1) {
		cerr << endl << "Usage error: Currently, the latent space model requires the number of clusters (-K) to be specified." << endl;
		exit(1);
	}
	PP(sbm :: ObjectiveFunction_Poisson :: s);
	PP(sbm :: ObjectiveFunction_Poisson :: theta);

	std :: auto_ptr<const sbm :: ObjectiveFunction> obj ( args_info.weighted_flag
		? up_cast<sbm :: ObjectiveFunction*>(new sbm :: ObjectiveFunction_Poisson(args_info.selfloop_flag, args_info.directed_flag, args_info.weighted_flag)    )
		: up_cast<sbm :: ObjectiveFunction*>(new sbm :: ObjectiveFunction_Bernoulli(args_info.selfloop_flag, args_info.directed_flag, args_info.weighted_flag)  )
		);
	std :: auto_ptr<const graph :: NetworkInterfaceConvertedToStringWithWeights > network;
	// Try the new graph loader
	if(args_info.stringIDs_flag) {
		network = graph :: loading :: make_Network_from_edge_list_string(edgeListFileName, args_info.directed_flag, args_info.weighted_flag, false,
				0 // "../../../survey_74_names_participants.txt"
				);
	} else {
		network = graph :: loading :: make_Network_from_edge_list_int64(edgeListFileName, args_info.directed_flag, args_info.weighted_flag, false, args_info.assume_N_nodes_arg);
	}
	{ // check for self loops
		if(network->get_plain_graph()->number_of_self_loops() > 0 && !args_info.selfloop_flag) {
			cerr << "Usage error: You have some self loops (" << network->get_plain_graph()->number_of_self_loops() << "). You must specify -s if to support the self-loop model. Exiting." << endl;
			exit(1);
		}
	}

	vector<int> groundTruth;
	if(args_info.GT_vector_given) { //populate groundTruth based on the --GT.vector option
		assert(groundTruth.size()==0);
		std :: ifstream GTvectorStream(args_info.GT_vector_arg);
		int n=0;
		while(1) {
			int z_i;
			GTvectorStream >> z_i;
			if(GTvectorStream.fail())
				break;
			groundTruth.push_back(z_i);
			++n;
		}
		if((int) groundTruth.size() != network->numNodes() || GTvectorStream.eof()==false ) {
			cerr << endl << "Error: the GT.vector file \"" << args_info.GT_vector_arg << "\" has " << groundTruth.size() << " valid lines (positive integers only). "
				<< "I was expecting " << network->numNodes() << " lines."
				<< endl;
			groundTruth.clear();
			exit(1);
		}
	}

	srand48(args_info.seed_arg);
	gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set(r, args_info.seed_arg);
	if(args_info.algo_sbm_cem_flag) {
		runCEM(network.get(), args_info.K_arg
				, obj.get()
				, groundTruth.empty() ? NULL : &groundTruth, args_info.iterations_arg, args_info, r);
	} else if(args_info.model_scf_flag) {
		runSCF(network.get(), args_info.K_arg, args_info.initGT_flag, groundTruth.empty() ? NULL : &groundTruth, args_info.iterations_arg, r);
	} else {
		runSBM(network.get(), args_info.K_arg, obj.get(), args_info.initGT_flag, groundTruth.empty() ? NULL : &groundTruth, args_info.iterations_arg, args_info.algo_gibbs_arg, args_info.algo_m3_arg, args_info, r);
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
		const int newCluster = static_cast<int>(drand48() * s._k);
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

format_flag_stack :: FormatFlagStack stack;

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
		cout << stack.push << fixed << setw(4) << setprecision(1) << 100.0L * static_cast<double>(this->a)/this->n << " %" << stack.pop; // cout << static_cast<double>(this->a)/this->n;
		cout << "\t" << this->a << " / " << this->n;
		cout << endl;
		forEach(const typeof(pair<const int,int>) &x, amd :: mk_range(this->mostRecent)) {
			if(this->n >= x.first) {
				cout << stack.push
					<< '\t' << setw(7) << x.first << " // " << setw(7) << x.second << "\t "
					<< fixed
					<< setw(4)
					<< setprecision(1)
					<< 100.0L * static_cast<double>(x.second)/x.first << " %" << endl
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

static long double delta_P_z_x__1RowOfBlocks(const sbm :: State &s, const sbm :: ObjectiveFunction *obj, const int pre_k, const int t, const int isolatedClusterId, const long double isolatedNodesSelfLoop);

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
	const int left = static_cast<int>(drand48() * s._k);
	int right_; do { right_ = static_cast<int>(drand48() * s._k); } while (left==right_);
	const int right = right_;
#ifdef M3_debugPrinting
	PP3(s._k, left,right);
#endif

	vector<int> randomizedNodeIDs; // the nodes in both clusters. Will be considered in a random order
	const sbm :: Cluster * const lCluster = s.labelling.clusters.at(left);
	const sbm :: Cluster * const rCluster = s.labelling.clusters.at(right);
	forEach(int n, amd :: mk_range(lCluster->get_members())) { randomizedNodeIDs.push_back(n); }
	forEach(int n, amd :: mk_range(rCluster->get_members())) { randomizedNodeIDs.push_back(n); }

	const int M = (int)randomizedNodeIDs.size();
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
			{ const int old_cl_id = s.moveNode(nToRemove, left); assert(old_cl_id == isolatedClusterId); }
			const long double post_z_K_left = s.P_z_orders();
			s.moveNode(nToRemove, isolatedClusterId);
			assert(VERYCLOSE(pre_z_K, s.P_z_orders()));
			const long double delta_z_K_left = post_z_K_left - pre_z_K;
			assert(isfinite(delta_z_K_left));
			s.moveNode(nToRemove, right);
			const long double post_z_K_right = s.P_z_orders();
			s.moveNode(nToRemove, isolatedClusterId);
			assert(VERYCLOSE(pre_z_K, s.P_z_orders()));
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
			assert(VERYCLOSE(pre_z_K , s.P_z_orders()));
			const long double delta_z_K_left = post_z_K_left - pre_z_K;
			assert(isfinite(delta_z_K_left));
			s.moveNode(n, right);
			const long double post_z_K_right = s.P_z_orders();
			s.moveNode(n, isolatedClusterId);
			assert(VERYCLOSE(pre_z_K , s.P_z_orders()));
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
	const int n = static_cast<int>(drand48() * s._N);
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
				assert(VERYCLOSE(pre_z_K, s.P_z_orders()));
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

	vector<long double> bits_latent_space_for_each_cluster;
	if(!s.cluster_to_points_map.empty()) {
		for(int k = 0; k < pre_k; ++k) {
			bits_latent_space_for_each_cluster.push_back(my_likelihood(n , s, obj, k));
		}
		assert(bits_latent_space_for_each_cluster.size() == s.cluster_to_points_map.size());
		assert((int)bits_latent_space_for_each_cluster.size() == pre_k);
	}

	// now to select cluster with probability proportional to exp(delta_P_x_z_IfIMoveIntoClusterT + delta_P_z_K_IfIMoveIntoClusterT)
	vector<long double> combinedDelta;
	for(int t=0; t<pre_k; t++) {
		combinedDelta.push_back(delta_P_x_z_IfIMoveIntoClusterT.at(t) + delta_P_z_K_IfIMoveIntoClusterT.at(t)
			+ ( s.cluster_to_points_map.empty() ? 0 : (bits_latent_space_for_each_cluster.at(t) )  )
		);
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

#if 0
	{
		cout << endl << "estimated changes:" << endl;
		PP2(delta_P_x_z_IfIMoveIntoClusterT.at(newCluster), delta_P_x_z_IfIMoveIntoClusterT.at(origClusterID));
		PP2(delta_P_z_K_IfIMoveIntoClusterT.at(newCluster), delta_P_z_K_IfIMoveIntoClusterT.at(origClusterID));
		PP(bits_latent_space_for_each_cluster.at(newCluster));
		PP(bits_latent_space_for_each_cluster.at(origClusterID));
		PP(bits_latent_space_for_each_cluster.at(newCluster)-bits_latent_space_for_each_cluster.at(origClusterID));
	}
#endif
	return + delta_P_x_z_IfIMoveIntoClusterT.at(newCluster) + delta_P_z_K_IfIMoveIntoClusterT.at(newCluster)
	       - delta_P_x_z_IfIMoveIntoClusterT.at(origClusterID) - delta_P_z_K_IfIMoveIntoClusterT.at(origClusterID)
	       + ( s.cluster_to_points_map.empty() ? 0 : (bits_latent_space_for_each_cluster.at(newCluster)-bits_latent_space_for_each_cluster.at(origClusterID) )  )
		;
}

#if 0
// for a given node, go through all the *other* nodes in that cluster, and check whether they are connected or not
vector<int> cluster_mates(sbm :: State &s, const int n) {
	const int z_n = s.labelling.cluster_id.at(n);
	const sbm :: Cluster *CL = s.labelling.clusters.at(z_n);
	const std :: list<int> * mem = &CL->get_members();
	std :: vector<int> in_this_cluster;
	For(m, (*mem)) {
		if(n!= *m)
			in_this_cluster.push_back(*m);
	}
	sort(in_this_cluster.begin(), in_this_cluster.end());
	assert((int)in_this_cluster.size() == CL->order() - 1);
	return in_this_cluster; // I hope it does RVO (return value optimization) with this!
}
#endif

struct ANormalDistribution {
	sbm :: State :: point_type mean;
	long double variance;
	sbm :: State :: point_type draw(gsl_rng *r) const {
		sbm :: State :: point_type proposed_new_location;
		for(int d = 0; d < sbm :: State :: point_type :: dimensionality; ++d) {
			proposed_new_location.at(d) = gsl_ran_gaussian(r, sqrt(variance));
			proposed_new_location.at(d) += mean.at(d);
		}
		return proposed_new_location;
	}
	long double pdf_prop( sbm :: State :: point_type x ) {
		const long double dist_2 = x.dist_2(this->mean);
		return M_LOG2E *  ( - dist_2 / ( 2 * variance ));
	}
   ANormalDistribution(const int n, const sbm :: State &s, const sbm :: ObjectiveFunction *obj) {
	// In theory, this can be pretty arbitrary. We just need something we can draw from, and be able to
	// ask it for the density at various points
	assert(!obj->weighted);
	assert(!obj->selfloops);
	this->mean.zero();
	long double sum_weight_in_this_cluster = 0.0L;
	const int k = s.labelling.cluster_id.at(n);
	const sbm :: Cluster *CL = s.labelling.clusters.at(k);
	for(auto m : CL->get_members()) {
		if(n!=m) {
			assert(s.labelling.cluster_id.at(m) == k);
			const long double sum_weight_on_this_rel = s.sum_weights_BOTH_directions(n,m);
			sum_weight_in_this_cluster += sum_weight_on_this_rel;
			sbm :: State :: point_type neighbours_position = s.cluster_to_points_map.at(k).at(m);
			assert(sbm :: is_integer(sum_weight_on_this_rel));
			switch((int)sum_weight_on_this_rel) {
				break; case 2: assert(obj->directed); this->mean += neighbours_position; this->mean += neighbours_position;
				break; case 1: this->mean += neighbours_position;
				break; case 0:
				break; default: assert(1==2); // shouldn't get here!  assert(sum_weight_on_this_rel == 0 || sum_weight_on_this_rel == 1 || sum_weight_on_this_rel == 2);
			}
		}
	}
	// now, we incorporate something like the prior into our proposal, by dividing by ( 1+ num_neighbours_in_same_cluster)
	for(int d = 0; d < sbm :: State :: point_type :: dimensionality; ++d) {
		this->mean.at(d) /= ((1.0/(2*sbm :: ls_prior_sigma_2))+sum_weight_in_this_cluster);
	}
	for(int d = 0; d < sbm :: State :: point_type :: dimensionality; ++d) {
		assert(isfinite(this->mean.at(d)));
	}

	this->variance = 1.0 /(  1.0/sbm :: ls_prior_sigma_2 + 2.0*sum_weight_in_this_cluster);
	assert(isfinite(this->variance));
    }
};


static long double my_likelihood(const int n, sbm :: State &s, const sbm :: ObjectiveFunction *obj, int k /* = -1 */) {
	// given this one node, what is the likelihood of the edges/non-edges in its cluster?
	if(k==-1)
		k = s.labelling.cluster_id.at(n);
	const sbm :: Cluster *CL = s.labelling.clusters.at(k);
	sbm :: State :: point_type current_position = s.cluster_to_points_map.at(k).at(n);
	double l2_bits = 0.0L;
	if(1){ // P( position | sigma) - the prior on the positions
		sbm :: State :: point_type pzero;
		pzero.zero();
		const double dist_2 = pzero.dist_2(current_position);
		const double ln_prior_position = - dist_2 / (2* sbm :: ls_prior_sigma_2);
		const double l2_prior_position = M_LOG2E * ln_prior_position;
		l2_bits += l2_prior_position;
	}
	for(auto m : CL->get_members()) {
		if(n==m)
			continue;
		sbm :: State :: point_type neighbours_position = s.cluster_to_points_map.at(k).at(m);
		const long double sum_weight_on_this_rel = s.sum_weights_BOTH_directions(n,m);
		assert(sbm :: is_integer(sum_weight_on_this_rel));
		switch(obj->directed) {
			break; case true:
				switch((int)sum_weight_on_this_rel) {
					break; case 0: l2_bits += 2.0L * l2_likelihood( current_position, neighbours_position, false);
					break; case 1: l2_bits += l2_likelihood( current_position, neighbours_position, false) + l2_likelihood( current_position, neighbours_position, true);
					break; case 2: l2_bits += 2.0L * l2_likelihood( current_position, neighbours_position, true);
					break; default: assert(1==2);
				}
			break; case false:
				switch((int)sum_weight_on_this_rel) {
					break; case 0: l2_bits += l2_likelihood( current_position, neighbours_position, false);
					break; case 1: l2_bits += l2_likelihood( current_position, neighbours_position, true);
					break; default: assert(1==2);
				}
		}
	}
	assert(isfinite(l2_bits));
	return l2_bits;
}
#if 0
long double update_one_nodes_position(const int n, sbm :: State &s, const sbm :: ObjectiveFunction *obj, AcceptanceRate *AR, gsl_rng * r) {
	const int k = s.labelling.cluster_id.at(n);

	sbm :: State :: point_type current_position = s.cluster_to_points_map.at(k).at(n);

	ANormalDistribution normpdf(n, s, obj); // density will be near the neighbours of n

	sbm :: State :: point_type proposed_new_location = normpdf.draw(r);
	// need to calculate the proposal probabilities at both locations
	// and the posterior densities
	//
	//
	// I should do the posterior densities first. This means checking this node against every other node in this cluster
	const long double current_likelihood = my_likelihood(n, s , obj);
	s.cluster_to_points_map.at(k).at(n) = proposed_new_location;
	const long double     new_likelihood = my_likelihood(n, s , obj);
	s.cluster_to_points_map.at(k).at(n) = current_position;

	const long double new_proposal_prob = normpdf.pdf(proposed_new_location);
	const long double old_proposal_prob = normpdf.pdf(current_position);

	const long double acceptance_prob = new_likelihood - current_likelihood - new_proposal_prob + old_proposal_prob;
	assert(isfinite(acceptance_prob));

	if( log2l(gsl_ran_flat(r,0,1)) < acceptance_prob ) {
		s.cluster_to_points_map.at(k).at(n) = proposed_new_location;
		AR->notify(true);
		return new_likelihood - current_likelihood;
	} else {
		AR->notify(false);
		return 0;
	}


}
#endif
long double gibbs_update_one_nodes_position(const int n, sbm :: State &s, const sbm :: ObjectiveFunction *obj, AcceptanceRate *, gsl_rng * r) {
	// cout << "gibbs_update_one_nodes_position" << endl;
	// PP(n);
	const int k = s.labelling.cluster_id.at(n);
	const long double current_likelihood = my_likelihood(n, s , obj, k);

	ANormalDistribution normpdf(n, s, obj); // density will be near the neighbours of n, and incorporates the prior
	// for(int d = 0; d<DIMENSIONALITY; ++d) { PP(normpdf.mean.at(d)); }
	int count_attempts = 0;
	// cout << "start the attempts to find a new position" << endl;
	while(1) {
		++count_attempts;
		if(count_attempts > 10)
			PP(count_attempts);
		const sbm :: State :: point_type proposed_new_location = normpdf.draw(r);

		long double acceptance_prob = 0.0L;
		assert(acceptance_prob <= 0.0);

		const sbm :: Cluster *CL = s.labelling.clusters.at(k);
		for(auto m : CL->get_members()) {
			if(n==m)
				continue;
			const sbm :: State :: point_type neighbours_position = s.cluster_to_points_map.at(k).at(m);
			const long double dist_2 = proposed_new_location.dist_2(neighbours_position);
			const long double this_is_what_is_actually_subtracted = -log2l(1+expl(sbm :: ls_alpha_k - dist_2));
			acceptance_prob += this_is_what_is_actually_subtracted;
			assert(this_is_what_is_actually_subtracted <= 0);
		}
		assert(acceptance_prob <= 0.0);

		if( log2l(gsl_ran_flat(r,0,1)) < acceptance_prob ) {
			s.cluster_to_points_map.at(k).at(n) = proposed_new_location;
			break;
		}
	}
	const long double new_likelihood = my_likelihood(n, s , obj, k);
	// PP3(new_likelihood, current_likelihood, count_attempts);
	return new_likelihood - current_likelihood;
}

long double update_ls_positions(sbm :: State &s, const sbm :: ObjectiveFunction *obj, AcceptanceRate *ar, gsl_rng * r) {
	// - take each cluster in turn.
	//   - take each node in that cluster in turn
	//     - propose a new position based on it cluster-mate-NEIGHBOURS
	//     - calculate proposal prob, again based on cluster-mate-NEIGHBOURS
	//     - calculate the posterior, based on ALL-cluster-mates

	// const long double pre = s.pmf(obj);
	assert(s._k == (int)s.cluster_to_points_map.size());
	assert(!obj->weighted);
	assert(!obj->selfloops); // I think I'll never like this!

	// update each node in turn, move it towards the cluster-mates it's connected to, and away from those it's not connected to
	long double l2_delta_bits = 0.0;
	for(int k = 0; k<s._k; ++k) {
		const sbm :: Cluster *CL = s.labelling.clusters.at(k);
		const std :: list<int> & mem = CL->get_members();
		for(auto n : mem) {
			//l2_delta_bits += update_one_nodes_position(n, s, obj, ar, r);
			l2_delta_bits += gibbs_update_one_nodes_position(n, s, obj, ar, r);
		}
	}

	// const long double post = s.pmf(obj);
	// PP2(post - pre , l2_delta_bits);
	// assert(VERYCLOSE(post - pre , l2_delta_bits));
	return l2_delta_bits;
}
pair<long double, long double> prepare_two_M3_ls_proposals(sbm :: State &s, const sbm :: ObjectiveFunction *obj, const int n, const int left, const int right) {
	assert(!obj->weighted && !obj->selfloops && !obj->directed);
	assert((int)s.cluster_to_points_map.size() == s._k);

	// given a node, and the two candidate communities:
	//    place it at the mean of its neighbours
	//    evaluate the overall pmf() for both options
	//    normalize the pmfs
	//
	const int temp_cluster_id = s._k-1;
	assert(s.labelling.cluster_id.at(n) == temp_cluster_id);
	s.moveNodeAndInformOfEdges(n, left);
	ANormalDistribution normpdf_left(n, s, obj); // density will be near the neighbours of n, and it incorporates the prior sigma
	s.cluster_to_points_map.at(left).at(n) = normpdf_left.mean;
	long double pmf_left = s.pmf(obj);

	s.moveNodeAndInformOfEdges(n, temp_cluster_id);

	assert(s.labelling.cluster_id.at(n) == temp_cluster_id);
	s.moveNodeAndInformOfEdges(n, right);
	ANormalDistribution normpdf_right(n, s, obj); // density will be near the neighbours of n, and it incorporates the prior sigma
	s.cluster_to_points_map.at(right).at(n) = normpdf_right.mean;
	long double pmf_right = s.pmf(obj);

	s.moveNodeAndInformOfEdges(n, temp_cluster_id);

	// PP2(pmf_left, pmf_right);
	const long double max_pmf = max(pmf_left, pmf_right);
	pmf_left -= max_pmf;
	pmf_right -= max_pmf;
	// PP2(pmf_left, pmf_right);
	assert(isfinite(pmf_left) && isfinite(pmf_right));
	assert(pmf_left <= 0);
	assert(pmf_right <= 0);

	const long double denominator = exp2l(pmf_left) + exp2l(pmf_right);
	return make_pair (exp2l(pmf_left) / denominator, exp2l(pmf_right) / denominator);
}
long double M3_LS(sbm :: State &s, const sbm :: ObjectiveFunction *obj, AcceptanceRate *AR, gsl_rng *r) {
	if(s._k < 2)
		return 0.0;
	assert(!obj->weighted);
	assert(!obj->selfloops);
	assert(!obj->directed);
	assert(!s.cluster_to_points_map.empty());
tryagain:
	const int left  = 0; //s._k * gsl_ran_flat(r,0,1);
	const int right = 1; //s._k * gsl_ran_flat(r,0,1);
	if(left==right) goto tryagain;

	const long double pre = s.pmf(obj);

	// store the old states
	const std :: vector< int > old_cluster_ids = s.labelling.cluster_id;
	const std :: vector< std :: vector<sbm :: State :: point_type> > old_cluster_to_points_map = s.cluster_to_points_map;

	// move into the temporary cluster
	vector<int> both_sets_of_nodes;
	const int temp_cluster_id = s._k;
	s.appendEmptyCluster();
	sbm :: State :: point_type pzero; pzero.zero();
	s.cluster_to_points_map.push_back( vector<sbm :: State :: point_type> (s._N, pzero) ) ; // give them a pretend location at the origin
	for(int i = 0; i<s._N; i++) {
		if(s.labelling.cluster_id.at(i) == left || s.labelling.cluster_id.at(i) == right) {
			s.moveNodeAndInformOfEdges(i, temp_cluster_id);
			both_sets_of_nodes.push_back(i);
		}
	}
	assert(s.labelling.clusters.at(left)->order()==0);
	assert(s.labelling.clusters.at(right)->order()==0);
	std :: random_shuffle(both_sets_of_nodes.begin(), both_sets_of_nodes.end());

	long double l2_reverse_prop = 0.0L;
	for(auto one_node : both_sets_of_nodes) {
		pair<long double, long double> two_options = prepare_two_M3_ls_proposals(s, obj, one_node, left, right);
		assert(VERYCLOSE(two_options.first + two_options.second, 1.0));
		const int orig_cluster_id = old_cluster_ids.at(one_node);
		if(orig_cluster_id == left) {
			l2_reverse_prop += log2l(two_options.first);
		} else {
			assert(orig_cluster_id == right);
			l2_reverse_prop += log2l(two_options.second);
		}
		s.moveNodeAndInformOfEdges(one_node, orig_cluster_id);
		ANormalDistribution normpdf_selected(one_node, s, obj);
		l2_reverse_prop += normpdf_selected.pdf_prop(s.cluster_to_points_map.at(orig_cluster_id).at(one_node));
	}
	// PP(l2_reverse_prop);

	{ // revert for verification, and then put things back
		s.deleteClusterFromTheEnd();
		s.cluster_to_points_map.pop_back();
		assert(s.cluster_to_points_map.size() == old_cluster_to_points_map.size());
		s.cluster_to_points_map = old_cluster_to_points_map;
		const long double verify_pre = s.pmf(obj);
		assert(VERYCLOSE(pre, verify_pre));
		s.cluster_to_points_map.push_back( vector<sbm :: State :: point_type> (s._N, pzero) ) ; // give them a pretend location at the origin
		s.appendEmptyCluster();
	}

	// we must move them out again, in order to do our random proposal
	for(auto one_node : both_sets_of_nodes) {
		s.moveNodeAndInformOfEdges(one_node, temp_cluster_id);
	}
	// at this stage, the positions are as before, except with the extra dummy positions for the temp cluster

	// finally, doing the random proposal
	long double l2_random_prop = 0.0L;
	for(auto one_node : both_sets_of_nodes) {
		pair<long double, long double> two_options = prepare_two_M3_ls_proposals(s, obj, one_node, left, right);
		assert(VERYCLOSE(two_options.first + two_options.second, 1.0));
		if(gsl_ran_flat(r,0,1) < two_options.first) {
			l2_random_prop += log2l(two_options.first);
			s.moveNodeAndInformOfEdges(one_node, left);
			assert(left == s.labelling.cluster_id.at(one_node));
		} else {
			l2_random_prop += log2l(two_options.second);
			s.moveNodeAndInformOfEdges(one_node, right);
			assert(right == s.labelling.cluster_id.at(one_node));
		}
		ANormalDistribution normpdf_selected(one_node, s, obj);
		sbm :: State :: point_type new_position = normpdf_selected.draw(r);
		l2_random_prop += normpdf_selected.pdf_prop(new_position);
		s.cluster_to_points_map.at(s.labelling.cluster_id.at(one_node)).at(one_node) = new_position;
		// PP2(__LINE__, l2_random_prop);
	}
	// PP2(__LINE__, l2_random_prop);
	// PP2( s.labelling.clusters.at(left)->order(), s.labelling.clusters.at(right)->order());

	// delete the temporary community, we don't need it any more
	s.deleteClusterFromTheEnd();
	s.cluster_to_points_map.pop_back();
	assert(s.cluster_to_points_map.size() == old_cluster_to_points_map.size());

	// now to decide whether to keep this (accept) or to revert (reject)
	const long double pmf_random = s.pmf(obj);
	// PP2(pre, pmf_random);
	// PP2(l2_reverse_prop, l2_random_prop);
	const long double l2_acceptance_probability = pmf_random - pre + l2_reverse_prop - l2_random_prop;
	// PP(l2_acceptance_probability);

	if( log2l(gsl_ran_flat(r,0,1)) < l2_acceptance_probability ) {
		// keep it
		// cout << "ACCEPTED M3_LS" << endl;
		AR->notify(true);
		return pmf_random - pre;
	} else {
		s.cluster_to_points_map = old_cluster_to_points_map;
		// restore everything
		for(auto i : both_sets_of_nodes) {
			const int orig_cluster_id = old_cluster_ids.at(i);
			if(s.labelling.cluster_id.at(i) != orig_cluster_id) {
				s.moveNodeAndInformOfEdges(i, orig_cluster_id);
			}
		}
		const long double post = s.pmf(obj);
		assert(VERYCLOSE(pre,post));
		AR->notify(false);
		return 0.0;
	}
}

static long double delta_P_z_x__1RowOfBlocks(const sbm :: State &s, const sbm :: ObjectiveFunction *obj, const int pre_k, const int t, const int isolatedClusterId, const long double isolatedNodesSelfLoop) {
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

		// NOTE: We ignore this on-diagonal data when using the latentspace model
		if(s.cluster_to_points_map.empty())
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

static
long double MoneNode(sbm :: State &s, const sbm :: ObjectiveFunction *obj, AcceptanceRate *AR) {
	if(s._k == 1)
	       return 0.0L;	// can't move a node unless there exist other clusters
	assert(s._k > 1); // can't move a node unless there exist other clusters
	const int n = static_cast<int>(drand48() * s._N);
	const int oldClusterID = s.labelling.cluster_id.at(n);
	int newClusterID;
	do {
		newClusterID = static_cast<int>(drand48() * s._k);
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
static long double MetropolisOnK(sbm :: State &s, const sbm :: ObjectiveFunction *obj __attribute__((unused)), AcceptanceRate *AR) {
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

			// + LOG2GAMMA(s._alpha)  // the change to SumOfLog2Facts
			// I belatedly noticed that this last line can be cancelled against some expressions in the next line.

			// log2(preK) - log2(s._N+preK)
			+(LOG2GAMMA(postK * s._alpha) - LOG2GAMMA(postK * s._alpha + s._N) /*- postK*LOG2GAMMA(s._alpha)*/ )
			-(LOG2GAMMA(preK * s._alpha) - LOG2GAMMA(preK  * s._alpha + s._N) /*- preK*LOG2GAMMA(s._alpha)*/ )
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
		const int clusterToProposeDelete = static_cast<int>(drand48() * s._k);
		if(s.labelling.clusters.at(clusterToProposeDelete)->order()>0) {
			// can't delete this. Time to get out
			assert(s._k==preK);
			AR->notify(false);
			return 0.0L;
		}
		if(clusterToProposeDelete != s._k-1) {
			// cout << "swapped during deletion "; PP2(clusterToProposeDelete , s._k-1);
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

static const long double a = 1.0; // this is 'a' as described in Nobile & Fearnside in section 3.2 ('Moves that change k'), except I'm just hardcoding at 1 for the moment, Uniform(0,1)
static long double EjectAbsorb_prop_prob_split(const int small_k, const int n_jBoth, const int n_j1, const int n_j2) {
	assert(small_k >= 1);
	const double proposalProbForMerge = small_k == 1 ? 0.0 : 0.5;
	return
			log2l(1.0-proposalProbForMerge) // probability of trying to split
			+ log2l(small_k)           // probability of selecting j1 to be the target of the split
			+ LOG2GAMMA(2 * a) - LOG2GAMMA(a) - LOG2GAMMA(a) // from section 3.2 of N+F 2001
			+ LOG2GAMMA(a + n_j1) + LOG2GAMMA(a + n_j2) - LOG2GAMMA(2*a + n_jBoth)
			;
}
static long double EjectAbsorb_prop_prob_merge(const int small_k) {
	assert(small_k >= 1);
	return
	-1 // log2l(0.5) , proposing a merge when k >= 2
	+ log2l(small_k)
	;
}
static long double EjectAbsorb(sbm :: State &s, const sbm :: ObjectiveFunction *obj, AcceptanceRate *AR, gsl_rng * r) {
	assert(obj);
	// cout << "EjectAbsorb" << endl;
	// we propose either to split, or to merge, two communities.
	const long double pre = s.pmf(obj);
	const double proposalProbForMerge = s._k == 1 ? 0.0 : 0.5;
	// PP2(s._k, proposalProbForMerge);
	if(drand48() < proposalProbForMerge) {
		assert(proposalProbForMerge == 0.5);
		assert(s._k >= 2);
		// cout << "Absorb" << endl;
		const int j2 = s._k - 1;
		const int j1 = drand48() * (s._k - 1);
		assert(j1 >= 0 && j1 < j2);
		const int n_j1 = s.labelling.clusters.at(j1)->order();
		const int n_j2 = s.labelling.clusters.at(j2)->order();
		const std :: vector<int> j2_members(s.labelling.clusters.at(j2)->get_members().begin(), s.labelling.clusters.at(j2)->get_members().end());
		// apply a merge
		For(i, j2_members) {
			const int shouldBej2 = s.moveNodeAndInformOfEdges(*i, j1);
			assert(j2 == shouldBej2);
		}
		const int n_jBoth = s.labelling.clusters.at(j1)->order();
		assert(n_jBoth == n_j1 + n_j2);
		s.deleteClusterFromTheEnd();
		const long double post = s.pmf(obj);
		// PP3(pre,post, post-pre);

		const long double full_proposal_probability = EjectAbsorb_prop_prob_merge(s._k);
		const long double reverse_proposal_probability = EjectAbsorb_prop_prob_split(s._k, n_jBoth, n_j1, n_j2);

		const long double acceptance_probability = post-pre - full_proposal_probability + reverse_proposal_probability;
		if(log2l(drand48()) < acceptance_probability) { // accept
			AR->notify(true);
			return post-pre;
		} else { // reject, split 'em again
			AR->notify(false);
			const int new_cluster_id = s.appendEmptyCluster();
			assert(new_cluster_id == j2);
			For(j2_member, j2_members) {
				const int shouldBej1 = s.moveNodeAndInformOfEdges(*j2_member, new_cluster_id);
				assert(j1 == shouldBej1);
			}
			return 0.0;
		}

		return post-pre;
	} else {
		// cout << "Eject" << endl;
		// apply a split
		// choose a cluster at random to split in two
		const int j1 = drand48() * s._k;
		const int j2 = s._k;
		assert(j1 >= 0 && j1 < j2);
		const std :: vector<int> j1_members(s.labelling.clusters.at(j1)->get_members().begin(), s.labelling.clusters.at(j1)->get_members().end());
		const int n_jBoth = j1_members.size();
		s.appendEmptyCluster();
		const long double p_E = gsl_ran_beta(r, a, a);
		int n_j2 = 0;
		int n_j1 = 0;
		For(j1_member, j1_members) {
			if(drand48() < p_E) {
				const int shouldBej1 = s.moveNodeAndInformOfEdges(*j1_member, j2);
				assert(j1 == shouldBej1);
				++ n_j2;
			} else
				++ n_j1;
		}
		assert(n_jBoth == n_j1 + n_j2);
		// gotta calculate the proposal probability properly
		const long double full_proposal_probability = EjectAbsorb_prop_prob_split(s._k - 1, n_jBoth, n_j1, n_j2);
		const long double reverse_proposal_probability = EjectAbsorb_prop_prob_merge(s._k - 1);
			;

		const long double post = s.pmf(obj);
		// PP3(pre,post, post-pre);

		const long double acceptance_probability = post-pre - full_proposal_probability + reverse_proposal_probability;
		if(log2l(drand48()) < acceptance_probability) { // accept
			AR->notify(true);
			return post-pre;
		} else { // reject
			AR->notify(false);
			const std :: vector<int> j2_members(s.labelling.clusters.at(j2)->get_members().begin(), s.labelling.clusters.at(j2)->get_members().end());
			For(j2_member, j2_members) {
				const int shouldBej2 = s.moveNodeAndInformOfEdges(*j2_member, j1);
				assert(j2 == shouldBej2);
			}
			s.deleteClusterFromTheEnd();
			return 0.0;
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
		assert(s._N == this->N);
		cout << "  CountSharedCluster(" << this->denominator << ")" << endl;
		for(int n=0; n<this->N; n++) {
			const std :: string node_name = s._g->node_name_as_string(n);
			for(int m=0; m<this->N; m++) {
				if(n==m)
					cout << "    ";
				else
					cout << stack.push << setw(4) << static_cast<int>(100.0 * static_cast<double>(shared.at(n*N+m)) / denominator + 0.5) << stack.pop;
			}
			const int my_cluster_id = s.labelling.cluster_id.at(n);
			cout
				<< "\t" << setw(1 + my_cluster_id) << my_cluster_id
				<< setw(s._k - my_cluster_id) << ""
				<< " \"" << node_name << "\"";
			cout << endl;
		}
	}
};

#define CHECK_PMF_TRACKER(track, actual) do { const long double _actual = (actual); long double & _track = (track); if(VERYCLOSE(_track,_actual)) { track = _actual; } else { PP(_actual - track); } assert(_track == _actual); } while(0)

static void runSBM(const graph :: NetworkInterfaceConvertedToStringWithWeights *g, const int commandLineK, const sbm :: ObjectiveFunction * const obj, const bool initializeToGT, const vector<int> * const groundTruth, const int iterations, const bool algo_gibbs, const bool algo_m3 , const  gengetopt_args_info &args_info, gsl_rng *r) {
	PP2(g->numNodes(), g->numRels());
	if(g->get_plain_graph()->number_of_self_loops() > 0 && !obj->selfloops ){
		cerr << endl << "Error: You must specify the -s flag to fully support self-loops. Your network has " << g->get_plain_graph()->number_of_self_loops() << " self-loops." << endl;
		exit(1);
	}
	sbm :: State s(g, args_info.mega_flag, args_info.alpha_arg);

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
	if(args_info.latentspace_flag) {
		assert(s._k == commandLineK);
		assert(s.cluster_to_points_map.empty());
		for(int k = 0; k < s._k; ++k) {
			vector<sbm :: State :: point_type> initial_random_locations(s._N);
			For(dbl, initial_random_locations) {
				for(int d = 0; d < sbm :: State :: point_type :: dimensionality; ++d) {
					dbl->at(d) = gsl_ran_gaussian(r, 1);
				}
			}
			s.cluster_to_points_map.push_back(initial_random_locations);
		}
		assert((int)s.cluster_to_points_map.size() == s._k);
		cout << "assigned initial random positions for the latent space model" << endl;
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
	AcceptanceRate AR_ea  ("EjectAbsorb");
	AcceptanceRate AR_lspos  ("LSSBM positions");
	AcceptanceRate AR_M3lspos  ("LSSBM M3");

	// some variables to check the PMP, i.e. the single most-visited state
	map< pair<int, vector<int> >, int> pmp_table;
	pair< pair<int, vector<int> >, int> best_pmp_so_far(make_pair(0, vector<int>()), 0);
	int64_t num_states_checked_for_pmp = 0;

	ofstream * save_z_fstream = NULL;
	if(args_info.save_z_arg[0]) {
		save_z_fstream = new ofstream(args_info.save_z_arg);
	}
	for(int i=1; i<=iterations; i++) {
		cout
			<< "Iteration:\t" << i
			<< "\tk0\t" << s._k
			<< "\tk1\t" << s.labelling.NonEmptyClusters
			<< "\t"
			;
		s.KandClusterSizes();
		if(0) { /// more swapping
			if(s._k > 1) // && drand48() < 0.01)
			{
				const int cl1 = static_cast<int>(s._k * drand48());
				const int cl2 = static_cast<int>(s._k * drand48());
				if(cl1 != cl2) {
					s.swapClusters(cl1,cl2);
					if(!s.cluster_to_points_map.empty()) {
						swap (s.cluster_to_points_map.at(cl1),s.cluster_to_points_map.at(cl2));
					}
				}
			}
		}

		int random_move = static_cast<int>(drand48() * 7);
		switch( random_move ) {
			break; case 0:
				if(commandLineK == -1) {
					if(args_info.algo_metroK_arg) {
						pmf_track += MetropolisOnK(s, obj, &AR_metroK);
					}
				} else
					assert(commandLineK == s._k);
			break; case 1:
				if(algo_gibbs) {
#if 0
					{
						CHECK_PMF_TRACKER(pmf_track, s.pmf(obj));
						cout << " == before Gibbs == " << endl;
						PP(s.P_edges_given_z_slow(obj)); // also includes *all* the latent space stuff
						PP(s.P_z_slow());
						PP(pmf_track);
						cout << " == about to Gibbs == " << endl;
					}
#endif

					pmf_track += gibbsOneNode(s, obj, &AR_gibbs);
#if 0

					{
						cout << " == done Gibbs == " << endl;
						PP(s.P_edges_given_z_slow(obj)); // also includes *all* the latent space stuff
						PP(s.P_z_slow());
						PP2(pmf_track, s.pmf(obj));
						CHECK_PMF_TRACKER(pmf_track, s.pmf(obj));
					}
#endif
				}
			break; case 2:
				if(s.cluster_to_points_map.empty()) {
				if(algo_m3)
					pmf_track += M3(s, obj, &AR_M3, &AR_M3little, &AR_M3very);
				}
			break; case 3:
				if(s.cluster_to_points_map.empty()) {
				if(commandLineK == -1) {
					if(args_info.algo_ejectabsorb_arg)
						pmf_track += EjectAbsorb(s, obj, &AR_ea, r);
				}
				}
			break; case 4:
				if(s.cluster_to_points_map.empty()) {
				if(args_info.algo_1node_arg)
					pmf_track += MoneNode(s, obj, &AR_metro1Node);
				}
			break; case 5:
				if(args_info.algo_lspos_arg) {
					assert(!s.cluster_to_points_map.empty());
					assert(commandLineK == s._k);
					pmf_track += update_ls_positions(s, obj, &AR_lspos, r);
				}
			break; case 6:
				if(args_info.algo_lsm3_arg) {
					assert(!s.cluster_to_points_map.empty());
					assert(commandLineK == s._k);
					pmf_track += M3_LS(s, obj, &AR_M3lspos, r);
				}
		}
		if(i > 30000) {
			if(count_shared_cluster)
				count_shared_cluster->consume(s.labelling.cluster_id);
		}
		if(i > 30000) {
			vector<int> clustering_copy( s.labelling.cluster_id );
			// let's put them in order before checking for PMP. this is a cheap form of label-switching
			int neg_id = -1;
			// find the first non-neg id, rename them all, and then try again
			for(int i=0; i<s._N; i++) {
				const int nonneg_id = clustering_copy.at(i);
				if(nonneg_id >= 0) {
					for(int j=0; j<s._N; j++) {
						if(clustering_copy.at(j)==nonneg_id)
							clustering_copy.at(j)=neg_id;
					}
					-- neg_id;
				}
			}
			assert(-neg_id == s.labelling.NonEmptyClusters+1);
			{
					pmp_table[ make_pair(s._k, clustering_copy) ]++;
					const pair< pair<int, vector<int> >, int> &ref_to_latest = *pmp_table.find( make_pair(s._k, clustering_copy) );
					if(ref_to_latest.second > best_pmp_so_far.second)
						best_pmp_so_far = ref_to_latest;
					++num_states_checked_for_pmp;
			}
		}

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
			PP4(num_states_checked_for_pmp, best_pmp_so_far.second, 100.0*best_pmp_so_far.second/num_states_checked_for_pmp , best_pmp_so_far.first.first);
			// s.summarizeEdgeCounts();
			AR_metroK.dump();
			AR_metro1Node.dump();
			AR_gibbs.dump();
			AR_ea.dump();
			AR_lspos.dump();
			AR_M3lspos.dump();
			AR_M3.dump();
			AR_M3little.dump();
			AR_M3very.dump();
			s.blockDetail(obj);
			if(count_shared_cluster)
				count_shared_cluster->dump(s);
			cout << " end of check at i==" << i << ". (" << double(clock()) / CLOCKS_PER_SEC << " seconds)" << endl;
			CHECK_PMF_TRACKER(pmf_track, s.pmf(obj));
			s.internalCheck();
		}
		if(i*2 >= iterations && args_info.save_z_arg[0]) {
			PP(i);
			*save_z_fstream << s._k << ':';
			for(int n=0; n<s._N; n++) {
				*save_z_fstream << s.labelling.cluster_id.at(n);
				if(n+1 < s._N)
					*save_z_fstream << ',';
			}
			*save_z_fstream << endl;
		}
	}
	if(save_z_fstream)
	       save_z_fstream->close();
	s.shortSummary(obj, groundTruth); /*s.summarizeEdgeCounts();*/ s.blockDetail(obj);
	s.internalCheck();
}

typedef vector<long double > theta_t;
typedef vector<vector<long double > > pi_t;
typedef vector<int> z_t;

static void CEM_update_theta(theta_t &theta, const z_t &z) {
	const int K = theta.size();
	const int N = z.size();
	for(int k=0; k<K; k++)
		theta.at(k) = 0;
	for(int n=0; n<N; n++) {
		const int z_n = z.at(n);
		theta.at(z_n) += 1.0L / N;
	}
	for(int k=0; k<K; k++)
		PP(theta.at(k));
}

static void CEM_update_pi(pi_t &pi, const z_t &z, const graph :: NetworkInterfaceConvertedToStringWithWeights *g
		, const sbm :: ObjectiveFunction * const obj
		) {
	assert(g->get_edge_weights()->is_weighted() == obj->weighted);
	assert(g->get_edge_weights()->is_directed() == obj->directed);
	assert(!g->get_edge_weights()->is_weighted());
	const graph :: VerySimpleGraphInterface *graph = g->get_plain_graph();
	const int K = pi.size();
	const int N = z.size();
	for(int k=0; k<K; k++) {
		vector<long double> &pi_k = pi.at(k);
		assert(int(pi_k.size())==K);
		for(int l=0; l<K; l++) {
			pi_k.at(l) = 0.0L;
		}
	}
	PP2(N,g->numNodes());
	// first, set up pi such that it stores the count of edges in that block.
	// then, we'll divide all the elements of pi appropriately
	int total_num_edges = 0;
	for(int rel = 0; rel < graph->numRels(); rel++) {
		const std :: pair<int32_t, int32_t> & eps = graph->EndPoints(rel);
		assert(eps.first <= eps.second);
		int z_1 = z.at(eps.first);
		int z_2 = z.at(eps.second);
		// if undirected, the lower triangle is pi must be empty
		if(g->get_edge_weights()->is_directed()) {
			pi.at(z_1).at(z_2) += g->get_edge_weights()->getl2h(rel);
			pi.at(z_2).at(z_1) += g->get_edge_weights()->geth2l(rel);
			total_num_edges += g->get_edge_weights()->getl2h(rel);
			total_num_edges += g->get_edge_weights()->geth2l(rel);
		} else {
			if(z_1 > z_2)
				swap(z_1, z_2);
			pi.at(z_1).at(z_2) += g->get_edge_weights()->getl2h(rel);
			assert(1 == g->get_edge_weights()->getl2h(rel));
			total_num_edges++;
		}
		// a self loop is only reported in l2h
	}
	PP("checking");
	if(g->get_edge_weights()->is_directed()) {
		PP("directed");
		assert(total_num_edges >= g->numRels());
		// for(int k=0; k<K; k++) { for(int j=0; j<K; j++) { PP3(k,j,pi.at(k).at(j)); } }
	} else {
		PP("undirected");
		assert(total_num_edges == g->numRels());
		for(int k=0; k<K; k++) {
			for(int j=0; j<k; j++) {
				assert(pi.at(k).at(j)==0);
			}
		}
	}

	vector<int> z_size(K, 0);
	for(int i=0; i<N; i++)
		z_size.at(z.at(i))++;

	// for each block, divide pi by the number of pairs of nodes
	int blocks_considered = 0;
	for(int k=0; k<K; k++) {
		for(int j=0; j<K; j++) {
			if(!g->get_edge_weights()->is_directed() && k>j) {
				assert(pi.at(k).at(j)==0);
				continue;
			}
			++blocks_considered;
			int pairs = z_size.at(k) * z_size.at(j);
			// BUT
			if(k==j) {
				if(obj->selfloops)
					pairs = z_size.at(k) * (z_size.at(k)+1) / 2;
				else
					pairs = z_size.at(k) * (z_size.at(k)-1) / 2;
			}
			assert(pairs >= 0);
			if(pairs == 0) {
				assert(pi.at(k).at(j) == 0);
				pi.at(k).at(j) = 0.5L;
			} else 
				pi.at(k).at(j) /= pairs;
		}
	}
	if(g->get_edge_weights()->is_directed()) {
		assert(blocks_considered == K * K);
	} else {
		assert(blocks_considered == K * (K+1) / 2);
	}
	
	for(int k=0; k<K; k++) {
		for(int j=0; j<K; j++) {
			PP3(k,j,pi.at(k).at(j));
		}
	}
}

static void recurse(z_t &z, vector<int> &z_size, const int n, const int N, const int K
		, const theta_t &theta, const pi_t &pi, const graph :: NetworkInterfaceConvertedToStringWithWeights *g
		, const sbm :: ObjectiveFunction * const obj
		, const long double score
		, long double & best_score_so_far
		, long double log2_max_theta
		, long double log2_max_pi
		) {
	assert(obj->directed);
	assert(obj->selfloops);
	assert(n >= 0 && n <= N);
	if(n==12)
		PP2(z.at(n-1), best_score_so_far);
	if(n==N) { // we've got a partition
		{ // double check the final score
			// PP(score);
#ifdef verify_recurse
			vector<int> verify_z_size(K,0);
			for(int i=0; i<N; i++) {
				verify_z_size.at(z.at(i)) ++;
			}
			for(int k=0; k<K; k++) {
				assert(z_size.at(k) == verify_z_size.at(k));
			}
			long double verify_score = 0.0;
			for(int i = 0; i<N; i++)
				verify_score += log2l(theta.at(z.at(i)));
			for(int i=0;i<N;i++) {
				for(int j=0;j<N;j++) {
					int z_1 = z.at(i);
					int z_2 = z.at(j);
					const long double pi_12 = pi.at(z_1).at(z_2);
					verify_score += log2l(1.0L-pi_12);
				}
			}
			for(int rel = 0; rel<g->numRels(); rel++) {
				const std :: pair <int32_t, int32_t> eps = g->get_plain_graph()->EndPoints(rel);
				int z_1 = z.at(eps.first);
				int z_2 = z.at(eps.second);
				if(g->get_edge_weights()->getl2h(rel)>0) { // edge from first to second
					const long double pi_12 = pi.at(z_1).at(z_2);
					assert(pi_12 > 0.0L);
					assert(pi_12 < 1.0L); // TODO: This'll definitely break down sooner or later
					verify_score += log2l(pi_12) - log2l(1.0L-pi_12);
				}
				if(g->get_edge_weights()->geth2l(rel)>0) { // edge from second to first
					const long double pi_21 = pi.at(z_2).at(z_1);
					assert(pi_21 > 0.0L);
					assert(pi_21 < 1.0L); // TODO: This'll definitely break down sooner or later
					verify_score += log2l(pi_21) - log2l(1.0L-pi_21);
				}

			}
			// PP2(verify_score, best_score_so_far);
			// PP(verify_score - score);
			assert(VERYCLOSE(verify_score , score));
#endif

		}
		if(score > best_score_so_far)
			best_score_so_far = score;
		// exit(0);
		return;
	}
	const int pairs_left = N * N - n * n;
	// PP3(N,n,pairs_left);
	if(score
			+ (N-n)*log2_max_theta + pairs_left * log2_max_pi
		< best_score_so_far) {
		// PP2(log2_max_theta, log2_max_pi);
		// PP3(n, score, best_score_so_far);
		// exit(2);
		return;
	}
	for(int k=0; k<K; k++) {
		z.at(n) = k;
		long double new_score = score;
		new_score += log2l(theta.at(k));
		assert(new_score < score);
		assert(new_score <= score + log2_max_theta);
		// consider just the neighbours (and pairs of nodes) at node n
		for(int l =0; l < K; l++) { // assume it's *disconnected* from every one of the clusters
			new_score += z_size.at(l) * log2l(1.0L-pi.at(k).at(l));
			new_score += z_size.at(l) * log2l(1.0L-pi.at(l).at(k));
		}
		new_score += log2l(1.0L-pi.at(k).at(k)); // this accounts for the self loop (assumed missing)
		assert(new_score < score);

		// now, to correct for the edges that actually are present
		const vector<int32_t> & neigh_rels = g->get_plain_graph()->neighbouring_rels_in_order(n);
		For(neigh_rel, neigh_rels) {
			const std :: pair<int32_t, int32_t> & eps = g->get_plain_graph()->EndPoints(*neigh_rel);
			if(eps.first <= n && eps.second <=n) {
				if(g->get_edge_weights()->getl2h(*neigh_rel) > 0) { // link from first to second
					const int l1 = z.at(eps.first);
					const int l2 = z.at(eps.second);
					new_score += log2l(pi.at(l1).at(l2)) - log2l(1.0L-pi.at(l1).at(l2));
				}
				if(eps.first != eps.second && g->get_edge_weights()->geth2l(*neigh_rel) > 0) { // link from second to first
					const int l2 = z.at(eps.first);
					const int l1 = z.at(eps.second);
					new_score += log2l(pi.at(l1).at(l2)) - log2l(1.0L-pi.at(l1).at(l2));
				}
			}
		}
		assert(new_score < score);

		z_size.at(k) ++;
		recurse(z, z_size, n+1, N, K, theta, pi, g, obj, new_score, best_score_so_far, log2_max_theta, log2_max_pi);
		z_size.at(k) --;
	}
}
static void CEM_update_z(z_t &z, const theta_t &theta, const pi_t &pi, const graph :: NetworkInterfaceConvertedToStringWithWeights *g
		, const sbm :: ObjectiveFunction * const obj
		) {
	const int N = z.size();
	const int K = pi.size();
	vector<int> z_size(K, 0);

	// estimate the 'best-case scenario for theta
	const long double max_theta = *max_element(theta.begin(), theta.end());
	vector<long double> pis;
	for(int k=0; k<K;k++) {
		for(int l=0; l<K;l++) {
			pis.push_back(pi.at(k).at(l));
			pis.push_back(1.0L-pi.at(k).at(l));
		}
	}
	assert(int(pis.size()) == 2 * K * K);
	const long double max_pi = *max_element(pis.begin(), pis.end());
	long double best_score_so_far = -DBL_MAX;
	recurse(z, z_size, 0, N, K, theta, pi, g, obj, 0, best_score_so_far, log2l(max_theta), log2l(max_pi));
}

static void runCEM(const graph :: NetworkInterfaceConvertedToStringWithWeights *g, const int commandLineK
		, const sbm :: ObjectiveFunction * const obj
		, const vector<int> * const groundTruth, const int iterations, const  gengetopt_args_info &args_info, gsl_rng *r) {
	// Given K and x, alternate
	// - maximize P(z,x|theta,pi,K) wrt theta and pi
	// - maximize P(z,x|theta,pi,K) wrt z
	cout << " == Classification EM == " << endl;
	PP2(g->numNodes(), g->numRels());

	const int N = g->numNodes();
	const int K = commandLineK;
	
	{
		assert(!args_info.weighted_flag); // weights aren't allowed with this yet.
	}

	// initialize z randomly
	z_t z(N);
	for(int i=0; i<N; i++)
		z.at(i) = gsl_ran_flat(r,0,1) * K;
	theta_t theta(K); // no need to initialize, this'll be overwritten first
	pi_t pi(K, vector<long double>(K) ); // no need to initialize, this'll be overwritten first

	for(int iter=0; iter<iterations; iter++) {
		if(groundTruth) {
			const double nmi = sbm :: State :: NMI(*groundTruth, z);
			PP(nmi);
		}
		PP(iter);
		CEM_update_theta(theta, z);
		CEM_update_pi(pi, z, g, obj);
		CEM_update_z(z, theta, pi, g, obj);
	}
}

