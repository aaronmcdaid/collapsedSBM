#include "graph/network.hpp"
#include "graph/loading.hpp"
using namespace std;
#include <algorithm>
#include <iomanip>
#include <vector>
#include <deque>
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
#include <climits>
#include <signal.h>

#include "Range.hpp"
#include "format_flag_stack/format_flag_stack.hpp"
#include "macros.hpp"
#include "sbm_state.hpp"
#include "scf.hpp"
#include "sbm.hpp"

#include "gitstatus.hpp"
#include "cmdline.h"
#include "args_info.hpp"


template<typename T,typename V> T up_cast(V x) { return x; }

using namespace std;

struct UsageMessage {
};

/*
static void memory_usage() {
	int ret = system("cat /proc/$PPID/status | egrep 'VmRSS|VmSize'");
	ret  = ret;
}
*/


static void runSBM(const graph :: NetworkInterfaceConvertedToStringWithWeights *g, const int commandLineK, const sbm :: ObjectiveFunction * const obj, const bool initializeToGT, const vector<int> * const groundTruth, const int iterations, const  gengetopt_args_info &args_info, gsl_rng *r) ;

static long double my_likelihood(const int n, sbm :: State &s, const sbm :: ObjectiveFunction *obj, int k = -1);
static bool drawPiAndTest(const sbm :: State &s, const sbm :: ObjectiveFunction *obj, gsl_rng *r);

gengetopt_args_info args_info;

int main(int argc, char **argv) {
	if (cmdline_parser (argc, argv, &args_info) != 0)
		exit(1) ;
	if(args_info.git_version_flag) {
		cout << "== comment.txt and gitstatus.txt ==" << endl;
		cout << gitstatus;
		cout << "====" << endl;
		cout << "args:";
		for (int i=0; i<argc; i++) {
			cout << ' ' << argv[i];
		}
		cout << endl;
		cout << "=====" << endl << endl;
	}
	if(args_info.inputs_num != 1) {
		cmdline_parser_print_help();
		exit(1);
	}

	const char * edgeListFileName   = args_info.inputs[0];
	// const char * directoryForOutput = args_info.inputs[1];
	if(args_info.verbose_flag) {
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
	PP(args_info.algo_sm_arg);
	PP(args_info.algo_ejectabsorb_arg);
	PP(args_info.initGT_flag);
	PP(args_info.stringIDs_flag);
	PP(args_info.mega_flag);
	PP(args_info.alpha_arg);
	PP2(args_info.beta1_arg, args_info.beta2_arg);;
	if(args_info.GT_vector_given)
		PP(args_info.GT_vector_arg);
	else
		cout << "args_info.GT_vector_arg:<undefined>" << endl;
	PP(args_info.model_scf_flag);
	PP(args_info.scf_flag);
	PP(args_info.assume_N_nodes_arg);
	PP(args_info.save_z_arg);
	PP(args_info.latentspace_flag);
	PP(args_info.lsalpha_arg);
	PP(args_info.algo_lspos_arg);
	PP(args_info.algo_lsm3_arg);
	PP(args_info.uniformK_flag);
	PP(args_info.geometricK_flag);
	PP(args_info.save_lsz_arg);
	PP(args_info.labels_arg);
	PP(args_info.keep_arg);
	}
	//PP(args_info.gamma_s_arg);
	//PP(args_info.gamma_phi_arg);
	sbm :: ObjectiveFunction_Poisson :: s     = args_info.gamma_s_arg;
	sbm :: ObjectiveFunction_Poisson :: theta = args_info.gamma_phi_arg;
	sbm :: ls_alpha_k = args_info.lsalpha_arg;
	sbm :: ObjectiveFunction_Bernoulli :: beta_1     = args_info.beta1_arg;
	sbm :: ObjectiveFunction_Bernoulli :: beta_2     = args_info.beta2_arg;


	if(args_info.model_scf_flag) {
		unless(args_info.K_arg == 2) {
			cerr << "Usage error: Currently, the stochastic community finding model (uncollapsed) requires -K 2 as an argument. Exiting." << endl;
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
	if(args_info.uniformK_flag && args_info.geometricK_flag) {
		cerr << endl << "Usage error: You can't have a Geometric prior (-g) and a Uniform prior (-u). Just one, or the other, or neither." << endl;
		exit(1);
	}
	if(args_info.verbose_flag) {
	PP(sbm :: ObjectiveFunction_Poisson :: s);
	PP(sbm :: ObjectiveFunction_Poisson :: theta);
	}

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
			if(GTvectorStream.peek()==',') {
				char c;
				GTvectorStream >> c;
			}
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
	if(args_info.model_scf_flag) {
		runSCF(network.get(), args_info.K_arg, args_info.initGT_flag, groundTruth.empty() ? NULL : &groundTruth, args_info.iterations_arg, r);
	} else {
		runSBM(network.get(), args_info.K_arg, obj.get(), args_info.initGT_flag, groundTruth.empty() ? NULL : &groundTruth, args_info.iterations_arg, args_info, r);
	}
}

void randomize(sbm :: State &s, const int K) { // randomize the partition and have K clusters in it
	assert(s._k <= K);
	assert(K >= 1);
	while(s._k < K) {
		s.appendEmptyCluster();
	}
	// cout << "Randomizing.. ";
	for(int n=0; n<s._N; n++) {
		const int newCluster = static_cast<int>(drand48() * s._k);
		if(newCluster != s.labelling.cluster_id.at(n))
			s.moveNodeAndInformOfEdges2(n, newCluster);
	}
	// PP2(s._k, s.labelling.NonEmptyClusters);
	assert(s._k == K);
	s.internalCheck();
}

static bool fiftyfifty() {
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
	if(log2l(drand48()) < delta)
		b = true;
	else
		b = false;
	if(AR)
		AR->notify(b);
	return b;
}

static long double delta_P_z_x__1RowOfBlocks(const sbm :: State &s, const sbm :: ObjectiveFunction *obj, const int pre_k, const int t, const int isolatedClusterId, const long double isolatedNodesSelfLoop);

enum POSSIBLE_MOVES {
		POS_MetroK
		,POS_Gibbs
		,POS_M3
		,POS_SM_Merge_M3
		,POS_SM_Split_M3
		,POS_SM_Merge_CF
		,POS_SM_Split_CF
		,POS_AE
		,POS_MoneNode
		,POS_update_ls_positions
		,POS_M3_LS
};

template<typename T>
static inline void Ignore(const T&) { }

enum Strategy { STRATEGY_M3, STRATEGY_CF };

long double SM_worker(sbm :: State &s, const sbm :: ObjectiveFunction *obj
		, gsl_rng * //r
		, const vector<int> &all_nodes
		, const int left
		, const int right
		, vector<int> &z
		, enum Strategy strategy
		, bool reject_some_simple_proposals_immediately = true// e.g. if num_nodes == 0, or if first two nodes are assigned to the same cluster
		) {
	// given a set of nodes, all currently unassigned,
	// assign them to one of two clusters randomly.
	// Return *this* part of the proposal probability.
	// The z vector is used to return, *or to input*,
	// the assignments of each node.
	assert(left>=0  &&  left<s._k);
	assert(right>=0 && right<s._k);
	assert(left != right);
	const size_t num = all_nodes.size();
	assert(num == z.size());
	if(reject_some_simple_proposals_immediately && num == 0) {
		return 0.0L;
	}
	long double this_prop_prob = 0.0L;
	pair<int,int> justTheseClusters(left,right);
	int size_of_left_so_far = 0;
	int size_of_right_so_far = 0;
	for(size_t ii = 0; ii < num; ++ii) {
		const int n = all_nodes.at(ii);
		assert(n>=0 && n<s._N);
		assert(-1 == s.labelling.cluster_id.at(n));
// #define Paranoid


		long double left_score = 0.0;
		long double right_score = 0.0;
switch(strategy) {
break; case STRATEGY_M3 :
	{
#ifdef Paranoid
		const long double miss_scoreF= s.P_all_fastish(obj);
#endif
		const long double miss_score = s.P_all_fastish(obj, justTheseClusters);

		s.insertNodeAndInformOfEdges(n, left);
#ifdef Paranoid
		long double left_scoreF= s.P_all_fastish(obj);
#endif
		left_score = s.P_all_fastish(obj, justTheseClusters);


		s.removeNodeAndInformOfEdges(n);
#ifdef Paranoid
		// assert(VERYCLOSE(miss_scoreF, s.P_all_fastish(obj)));
#endif

		s.insertNodeAndInformOfEdges(n, right);
#ifdef Paranoid
		long double right_scoreF= s.P_all_fastish(obj);
#endif
		right_score = s.P_all_fastish(obj, justTheseClusters);
#ifdef Paranoid
		assert(VERYCLOSE(left_scoreF  - left_score  , miss_scoreF - miss_score));
		assert(VERYCLOSE(right_scoreF - right_score , miss_scoreF - miss_score));
#endif
		s.removeNodeAndInformOfEdges(n);
#ifdef Paranoid
		// assert(VERYCLOSE(miss_scoreF, s.P_all_fastish(obj)));
#endif

		if(!args_info.weighted_flag) {
			assert(miss_score > left_score);
			assert(miss_score > right_score);
		}

		const long double max_score = (left_score > right_score) ? left_score : right_score;
		left_score -= max_score;
		right_score -= max_score;
		left_score = exp2l(left_score);
		right_score = exp2l(right_score);
		// PP2(left_score, right_score);
	}
		if(reject_some_simple_proposals_immediately)
		switch(ii) {
			break; case 0: left_score = 1; right_score = 0.000000;
			break; case 1:   left_score = 0.000000; right_score = 1;
		}
break; case STRATEGY_CF :
		// This is the 'community-finding' strategy: be attracted
		// to one of the two clusters simply depending on which one
		// you have the most connections to.
		const std :: vector<int32_t> & neighs_of_n = s.vsg->neighbouring_nodes_in_order(n);
		int already_a_neighbour_on_the_left = 0;
		int already_a_neighbour_on_the_right = 0;
		For(neigh, neighs_of_n) {
			const int neigh_z = s.labelling.cluster_id.at(*neigh);
			if(neigh_z == left)
				++already_a_neighbour_on_the_left;
			if(neigh_z == right)
				++already_a_neighbour_on_the_right;
		}
		left_score = (already_a_neighbour_on_the_left+0.5) / (size_of_left_so_far+1);
		right_score = (already_a_neighbour_on_the_right+0.5) / (size_of_right_so_far+1);
		// PP2(left_score, right_score);
}
		const long double total = left_score + right_score;
		left_score /= total;
		right_score /= total;
		// PP2(left_score, right_score);
		// assert(VERYCLOSE(left_score , 0.5L));
		// assert(VERYCLOSE(right_score , 0.5L));
		assert(isfinite(left_score));
		assert(isfinite(right_score));

		assert(VERYCLOSE(left_score + right_score , 1.0L));

		if(z.at(ii) == -1) {
			const long double unif = drand48();
			if(unif < left_score)
				z.at(ii) = left;
			else
				z.at(ii) = right;
		}
		s.insertNodeAndInformOfEdges(n, z.at(ii));

		if(z.at(ii) == left) {
			++ size_of_left_so_far;
			if(left_score == 0.0)
				return -LDBL_MAX;
			this_prop_prob += log2l(left_score);
		}
		else if(z.at(ii) == right) {
			++ size_of_right_so_far;
			if(right_score == 0.0)
				return -LDBL_MAX;
			this_prop_prob += log2l(right_score);
		} else
			assert(1==2);
		assert(isfinite(this_prop_prob));
	}
	assert(size_of_left_so_far == s.labelling.clusters.at(left)->order());
	assert(size_of_right_so_far == s.labelling.clusters.at(right)->order());
	return this_prop_prob;
}

long double SM_Split(sbm :: State &s, const sbm :: ObjectiveFunction *obj
		, gsl_rng *r
		, enum Strategy strategy
		) {
	/* Split
	 * 1. select a cluster at random and shuffle its nodes
	 * 2. remove all those nodes from that cluster
	 * 3. randomly reassign, recording the decisions made and remembering this part of the proposal probablity
	 * 4. calculate the two prob-ratios and execute
	 */
	/* Merge [implemented in SM_Merge(), not here]
	 * 1. select two clusters at random, and shuffle their nodes.
	 * 2. remove all those nodes from that cluster
	 * 3. NON-randomly reassign, FORCING the decisions made and remembering this part of the proposal probablity
	 * 4. calculate the two prob-ratios and execute
	 */
	if(args_info.maxK_arg > 0 && s._k >= args_info.maxK_arg) {
		return 0.0L;
	}
	const int pre_k = s._k;
	const int left = static_cast<int>(drand48() * pre_k);
	const int right = s._k;
	vector<int> all_nodes;
	{
		const sbm :: Cluster * const lCluster = s.labelling.clusters.at(left);
		For(node, lCluster->get_members()) { all_nodes.push_back(*node); }
		random_shuffle(all_nodes.begin(), all_nodes.end());
	}
	const int num = all_nodes.size();
	if(num==0) {
		return 0.0L; // change of plan again!  We will not split an empty cluster in
				// two, nor will we merge *two* *empty* clusters into one.
	}

	const long double pre_fast = s.P_all_fastish(obj);
	s.appendEmptyCluster();
	assert(s.labelling.clusters.at(s._k-1)->order() == 0);

	// remove all the nodes from their cluster
	// reassign randomly
	// calculate FULL proposal probability
	// Accept/Reject:
	// 	accept: return
	// 	reject: undo

	For(node, all_nodes) {
		const int n = *node;
		const int oldcl = s.labelling.removeNode(n);
		s.informNodeMove(n, oldcl, -1);
		assert(oldcl == left);
	}

	// all the nodes have now been removed, let's start using SM_worker
	vector<int> z(num, -1);
	const long double partial_prop_prob = SM_worker(s, obj, r, all_nodes, left, right, z, strategy);
	assert(partial_prop_prob != -LDBL_MAX);
	// Ignore(partial_prop_prob);
	const long double new_fast = s.P_all_fastish(obj);

	// const int small_k = pre_k;
	// const int big_k = s._k;
	assert(s._k == pre_k+1);
	const long double partial_acceptance_prob = new_fast - pre_fast
		// divide by this prop prob
		- partial_prop_prob
		// - (log2l(small_k) - log2l(big_k)) // select one of K clusters to be split, select one of K+1 slots into which to insert the new cluster
		// multiply by the reverse prop prob
		// - log2l(big_k) - log2l(small_k)
		; // select one of bigK clusters to be embiggened, select one of smallK to be the victim

	assert(!z.empty());
	const double unif = drand48();
	if(log2l(unif) < partial_acceptance_prob) {
		// Important: The new cluster shouldn't just be put on the end
		assert(s._k-1 == pre_k); // the new cluster's id is 'pre_k'
		const int empty_cluster_id = static_cast<int>(drand48() * s._k);
		if(empty_cluster_id != pre_k) {
			s.swapClusters(pre_k, empty_cluster_id);
		}
		return new_fast - pre_fast;
	}

	// Not accepted.  We should undo everything and return zero

	For(node, all_nodes) {
		const int n = *node;
		s.removeNodeAndInformOfEdges(n);
		s.insertNodeAndInformOfEdges(n,left);
	}

	s.deleteClusterFromTheEnd();
	const long double post2_fast = s.P_all_fastish(obj);
	assert(VERYCLOSE(pre_fast, post2_fast));

	return 0.0L;
}
long double SM_Merge(sbm :: State &s, const sbm :: ObjectiveFunction *obj
		, gsl_rng *r
		, enum Strategy strategy
		) {
	// cout << endl << "Attempt a merge" << endl;
	/* Merge
	 * 1. select two clusters at random, and shuffle their nodes.
	 * 2. record the pre_score_still_split
	 * 3. remove all those nodes from that cluster
	 * 4. NON-randomly reassign, FORCING the decisions made and remembering this part of the proposal probablity
	 * 5. verify that the score is back to the pre_score_still_split
	 * 6. force them to merge, and calculate post_score_now_merged
	 * 7. calculate the two prob-ratios and execute
	 */
	const int pre_k = s._k;
	if(pre_k < 2)
		return 0.0;

must_be_different_try_again:
	const int left = static_cast<int>(drand48() * pre_k);
	const int right = static_cast<int>(drand48() * pre_k);
	if(left == right)
		goto must_be_different_try_again;

	vector<int> all_nodes;
	{
		const sbm :: Cluster * const lCluster = s.labelling.clusters.at(left);
		For(node, lCluster->get_members()) { all_nodes.push_back(*node); }
		const sbm :: Cluster * const rCluster = s.labelling.clusters.at(right);
		For(node, rCluster->get_members()) { all_nodes.push_back(*node); }
		random_shuffle(all_nodes.begin(), all_nodes.end());
	}
	const int num = all_nodes.size();
	if(strategy == STRATEGY_M3) {
		bool will_fail = false;
		if(num>=2) {
			const int node1 = all_nodes.at(0);
			const int node2 = all_nodes.at(1);
			unless(s.labelling.cluster_id.at(node1) == left && s.labelling.cluster_id.at(node2) == right) {
				will_fail = true;
			}
		}
		if(num==1) {
			const int node1 = all_nodes.at(0);
			unless(s.labelling.cluster_id.at(node1) == left) {
				will_fail = true;
			}
		}
		if(will_fail) return 0.0L;
	}
	if(num==0) {
		return 0.0L; // change of plan again!  We will not split an empty cluster in
				// two, nor will we merge *two* *empty* clusters into one.
	}

	const long double pre_score_still_split_up = s.P_all_fastish(obj);

	// remove all the nodes from their cluster
	// reassign randomly
	// calculate FULL proposal probability
	// Accept/Reject:
	// 	accept: return
	// 	reject: undo

	vector<int> z(num, -1);
	for(int ii = 0; ii < num; ++ii) {
		const int n = all_nodes.at(ii);
		const int oldcl = s.labelling.removeNode(n);
		z.at(ii) = oldcl;
		s.informNodeMove(n, oldcl, -1);
		assert(oldcl == left || oldcl == right);
	}

	// all the nodes have now been removed, let's start using SM_worker
	const long double pre_worker = s.P_all_fastish(obj);
	const long double partial_prop_prob = SM_worker(s, obj, r, all_nodes, left, right, z, strategy);

	assert(partial_prop_prob != -LDBL_MAX);

	const long double post_worker = s.P_all_fastish(obj);
	assert(VERYCLOSE(pre_score_still_split_up, post_worker));
	if(!args_info.weighted_flag)
		assert(pre_worker >= post_worker);

	// Now, we need to force them to merge and calculate the score there
	for(int ii = 0; ii < num; ++ii) {
		const int n = all_nodes.at(ii);
		const int oldcl = s.labelling.cluster_id.at(n);
		assert(oldcl == z.at(ii));
		if(oldcl != left) {
			assert(oldcl == right);
			// const int oldcl = s.labelling.removeNode(n);
			// s.informNodeMove(n, oldcl, -1);
			s.moveNodeAndInformOfEdges(n, left);
		}
	}

	assert(s.labelling.clusters.at(right)->order() == 0);
	if(right != s._k-1) s.swapClusters(right, s._k-1);
	s.deleteClusterFromTheEnd();
	const long double post_score_now_merged = s.P_all_fastish(obj);

	// const int big_k = pre_k;
	// const int small_k = s._k;
	assert(s._k == pre_k-1);

	s.appendEmptyCluster();
	if(right != s._k-1) s.swapClusters(right, s._k-1);

	const long double partial_acceptance_prob = post_score_now_merged - pre_score_still_split_up
		// divide by this prop prob
		// + log2l(big_k) + log2l(small_k)
		// multiply by the reverse prop prob
		+ partial_prop_prob
		//+( - log2l(small_k)-log2l(big_k))
		;

	assert(!z.empty());
	const double unif = drand48();
	if(log2l(unif) < partial_acceptance_prob) {
		// cout << "     Merge "; PP3(partial_acceptance_prob ,+ log2l(pre_k) + log2l(pre_k-1) ,partial_prop_prob - log2l(pre_k-1));
		assert(s.labelling.clusters.at(right)->order() == 0);
		if(right != s._k-1) s.swapClusters(right, s._k-1);
		s.deleteClusterFromTheEnd();
		// PP2 (post_score_now_merged , pre_score_still_split_up);
		return post_score_now_merged - pre_score_still_split_up;
	}

	// Not accepted.  We should undo everything, splitting them, and return zero

	for(int ii = 0; ii < num; ++ii) {
		const int n = all_nodes.at(ii);
		const int oldcl = s.labelling.cluster_id.at(n);
		assert(oldcl == left);
		if(z.at(ii) != oldcl) {
			assert(z.at(ii) == right);
			s.moveNodeAndInformOfEdges(n, z.at(ii) );
		}
	}

	const long double post2_fast = s.P_all_fastish(obj);
	assert(VERYCLOSE(pre_score_still_split_up, post2_fast));

	return 0.0L;
}

static long double M3(sbm :: State &s, const sbm :: ObjectiveFunction *obj
		// , AcceptanceRate * const AR, AcceptanceRate * const AR_alittleConservative, AcceptanceRate * const AR_veryConservative
		, gsl_rng *r
		) {
	if(s._k < 2) {
		return 0.0L; // should this be recorded as a rejection for the purpose of the acceptance rate?
	}

	// Select two distinct clusters
	// Record current state.
	// Make random proposal, recording proposal probability along the way
	// Calculate new cpmf
	// Make forced proposal, recording proposal probability along the way

gotta_be_distinct:
	const int left  = static_cast<int>(drand48() * s._k);
	const int right = static_cast<int>(drand48() * s._k);
	if(left==right)
		goto gotta_be_distinct;
	const pair<int,int> justTheseClusters = make_pair(left,right);
#define dump_clustering2() do {\
PP2(__LINE__, s.P_all_fastish(obj)); for(int i=0; i<s._N; ++i) { cout << setw(3) << s.labelling.cluster_id.at(i); } cout << endl; \
} while(0)

#define dump_clustering() do{}while(0)

	dump_clustering();

	vector<int> all_nodes;
	{
		const sbm :: Cluster * const lCluster = s.labelling.clusters.at(left);
		For(node, lCluster->get_members()) { all_nodes.push_back(*node); }
		const sbm :: Cluster * const rCluster = s.labelling.clusters.at(right);
		For(node, rCluster->get_members()) { all_nodes.push_back(*node); }
		random_shuffle(all_nodes.begin(), all_nodes.end());
	}
	const int num = all_nodes.size();

	// The various states the nodes will be in
	vector<int> z_original(num, -1);
	vector<int> z_random(num, -1);

	// Record the original state, while removing all the nodes from both clusters
	const long double pre_score = s.P_all_fastish(obj);
	for(int ii = 0; ii < num; ++ii) {
		const int n = all_nodes.at(ii);
		const int oldcl = s.labelling.removeNode(n);
		z_original.at(ii) = oldcl;
		s.informNodeMove(n, oldcl, -1);
		assert(oldcl == left || oldcl == right);
	}
	const long double unassigned_score = s.P_all_fastish(obj);
	Ignore(unassigned_score);
	dump_clustering();

	// Make random proposal
	const long double partial_prop_prob_random = SM_worker(s, obj, r, all_nodes, left, right, z_random, STRATEGY_M3, false);
	Ignore(partial_prop_prob_random);
	const long double randomized_score = s.P_all_fastish(obj);
	Ignore(randomized_score);
	dump_clustering();

	// Empty the two clusters again
	for(int ii = 0; ii < num; ++ii) {
		const int n = all_nodes.at(ii);
		const int oldcl = s.labelling.removeNode(n);
		assert(z_random.at(ii) == oldcl);
		s.informNodeMove(n, oldcl, -1);
		assert(oldcl == left || oldcl == right);
	}
	dump_clustering();

	// Make forced proposal
	const long double partial_prop_prob_forced = SM_worker(s, obj, r, all_nodes, left, right, z_original, STRATEGY_M3, false);
	Ignore(partial_prop_prob_forced);

	const long double undone_score = s.P_all_fastish(obj);
	dump_clustering();

	//PP3( pre_score, unassigned_score , undone_score) ;
	assertVERYCLOSE( undone_score , pre_score) ;
	assertEQ( undone_score , pre_score) ;

	// Finally ready to decide
	const long double acceptance_prob = randomized_score - pre_score - partial_prop_prob_random + partial_prop_prob_forced;
	if( log2l(gsl_ran_flat(r,0,1)) < acceptance_prob) {
		// Accept
		for(int ii = 0; ii < num; ++ii) {
			const int n = all_nodes.at(ii);
			const int cl = z_random.at(ii);
			assert(cl == left || cl == right);
			s.removeNodeAndInformOfEdges(n);
			s.insertNodeAndInformOfEdges(n,cl);
		}
		return randomized_score - pre_score;
	} else {
		return undone_score - pre_score;
	}
}

static long double beta_draw(const int p_kl, const int y_kl, gsl_rng *r) {
	assert(args_info.scf_flag);
	const int successes = y_kl;
	assert(successes >= 0);
	const int failures = p_kl - y_kl;
	assert(failures >= 0);
	const long double draw = gsl_ran_beta(r, successes + sbm :: ObjectiveFunction_Bernoulli :: beta_1, failures + sbm :: ObjectiveFunction_Bernoulli :: beta_2) ;
	return draw;
}
static long double gamma_draw(const int p_kl, const int y_kl, gsl_rng *r) {
	const long double s_post = sbm :: ObjectiveFunction_Poisson :: s + y_kl;
	const long double phi_post = 1.0L / (1.0L/sbm :: ObjectiveFunction_Poisson :: theta + p_kl);
	const long double draw = gsl_ran_gamma(r, s_post, phi_post);
	return draw;
}

static
bool drawPiAndTest(const sbm :: State &s, const sbm :: ObjectiveFunction *obj, gsl_rng *r) {
	assert(args_info.scf_flag);
	if(s._k == 1)
		return true;

	// draw for each of the diagonal entries
	long double min_on_diagonal = 1000000000.0;
	for(int k=0; k<s._k; ++k) {
		const long int y_kk = obj->relevantWeight      (k, k, &s._edgeCounts);
		const long int p_kk = obj->numberOfPairsInBlock(k, k, &s.labelling);
		const long double pi_kk = obj->weighted ? gamma_draw(p_kk, y_kk, r) : beta_draw(p_kk, y_kk, r);
		//PP3(p_kk, y_kk, pi_kk);
		min_on_diagonal = min(min_on_diagonal, pi_kk);
	}
	long double max_off_diagonal = -1.0;
	for(int k=0; k<s._k; ++k) {
		for(int l=0; l<s._k; ++l) {
			if(k==l)
				continue;
			if(l>k && !obj->directed)
				break;
			const long int y_kl = obj->relevantWeight      (k, l, &s._edgeCounts);
			const long int p_kl = obj->numberOfPairsInBlock(k, l, &s.labelling);
			const long double pi_kl = obj->weighted ? gamma_draw(p_kl, y_kl, r) : beta_draw(p_kl, y_kl, r);
			//PP3(p_kl, y_kl, pi_kl);
			max_off_diagonal = max(max_off_diagonal, pi_kl);
		}
	}
	assert(max_off_diagonal > 0.0L);
	assert(min_on_diagonal < 1000000000.0);
	// PP2(min_on_diagonal, max_off_diagonal);
	// assert(min_on_diagonal <= max_off_diagonal);
	return min_on_diagonal > max_off_diagonal;
}

static
long double gibbsOneNode(sbm :: State &s, const sbm :: ObjectiveFunction *obj, AcceptanceRate *AR, gsl_rng *r, const int n) {
	if(args_info.scf_flag) {
		assert(!args_info.latentspace_flag);
		assert(s.cluster_to_points_map.empty());
	}

	if(s._k == 1) {
		AR->notify(false);
		return 0.0L;
	}
	const int pre_k = s._k;
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

	bool has_moved_to_new_cluster = newCluster != origClusterID;

// Now for the Stochastic Community Finding support.
// - draw from the posterior of \pi, and if it doesn't satisfy the constraint,
//   then we reject and revert the node to the origClusterID
	if(has_moved_to_new_cluster && args_info.scf_flag) {
		const bool SCF_posterior_test = drawPiAndTest(s, obj, r);
		if(!SCF_posterior_test) {
			has_moved_to_new_cluster = false;
			s.moveNodeAndInformOfEdges(n, origClusterID);
		}
	}

	AR->notify(has_moved_to_new_cluster);

	if (has_moved_to_new_cluster) {
		return + delta_P_x_z_IfIMoveIntoClusterT.at(newCluster) + delta_P_z_K_IfIMoveIntoClusterT.at(newCluster)
	       - delta_P_x_z_IfIMoveIntoClusterT.at(origClusterID) - delta_P_z_K_IfIMoveIntoClusterT.at(origClusterID)
	       + ( s.cluster_to_points_map.empty() ? 0 : (bits_latent_space_for_each_cluster.at(newCluster)-bits_latent_space_for_each_cluster.at(origClusterID) )  )
		;
	} else
		return 0;
}

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
   ANormalDistribution(): mean() {
	this->variance = sbm :: ls_prior_sigma_2;
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
	For(mp, CL->get_members()) {
		const int32_t m = *mp;
		if(n!=m) {
			assert(s.labelling.cluster_id.at(m) == k);
			const long double sum_weight_on_this_rel = s.sum_weights_BOTH_directions(n,m);
			sum_weight_in_this_cluster += sum_weight_on_this_rel;
			sbm :: State :: point_type neighbours_position = s.cluster_to_points_map.at(m);
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
	sbm :: State :: point_type current_position = s.cluster_to_points_map.at(n);
	double l2_bits = 0.0L;
	if(1){ // P( position | sigma) - the prior on the positions
		sbm :: State :: point_type pzero;
		pzero.zero();
		const double dist_2 = pzero.dist_2(current_position);
		const double ln_prior_position = - dist_2 / (2* sbm :: ls_prior_sigma_2);
		const double l2_prior_position = M_LOG2E * ln_prior_position;
		l2_bits += l2_prior_position;
		// PP2("prior ", l2_bits);
	}
	For(mp, CL->get_members()) {
		const int32_t m = *mp;
		if(n==m) {
			continue;
		}
		sbm :: State :: point_type neighbours_position = s.cluster_to_points_map.at(m);
		const long double sum_weight_on_this_rel = s.sum_weights_BOTH_directions(n,m);
		// PP3(n,m,sum_weight_on_this_rel);
		assert(sbm :: is_integer(sum_weight_on_this_rel));
		long double delta_l2_bits;
		switch(obj->directed) {
			break; case true:
				switch((int)sum_weight_on_this_rel) {
					break; case 0: delta_l2_bits = 2.0L * l2_likelihood( current_position, neighbours_position, false);
					break; case 1: delta_l2_bits = l2_likelihood( current_position, neighbours_position, false) + l2_likelihood( current_position, neighbours_position, true);
					break; case 2: delta_l2_bits = 2.0L * l2_likelihood( current_position, neighbours_position, true);
					break; default: assert(1==2);
				}
			break; case false:
				switch((int)sum_weight_on_this_rel) {
					break; case 0: delta_l2_bits = l2_likelihood( current_position, neighbours_position, false);
					break; case 1: delta_l2_bits = l2_likelihood( current_position, neighbours_position, true);
					break; default: assert(1==2);
				}
		}
		// PP(delta_l2_bits);
		l2_bits += delta_l2_bits;
	}
	assert(isfinite(l2_bits));
	// PP(l2_bits);
	return l2_bits;
}
static long double MH_update_one_nodes_position(const int n, sbm :: State &s, const sbm :: ObjectiveFunction *obj, AcceptanceRate *AR, gsl_rng * r) {
	// cout << "MH_update_one_nodes_position" << endl;
	// PP(n);
	const int k = s.labelling.cluster_id.at(n);

	// ANormalDistribution normpdf(n, s, obj); // density will be near the neighbours of n, and incorporates the prior
	ANormalDistribution prior_only;
	for(int i = 0; i < DIMENSIONALITY; ++i) {
		assert(prior_only.mean.at(i) == 0.0);
	}
	assert(prior_only.variance == 1.0);
	// for(int d = 0; d<DIMENSIONALITY; ++d) { PP(normpdf.mean.at(d)); }
	// cout << "start the attempts to find a new position" << endl;
	const sbm :: State :: point_type original_location = s.cluster_to_points_map.at(n);

	const sbm :: State :: point_type proposed_new_location = prior_only.draw(r);
	const long double new_prior_with_edges = prior_only.pdf_prop(proposed_new_location);
	const long double old_prior_with_edges = prior_only.pdf_prop(original_location);

	long double acceptance_prob = 0.0L;
	long double sort_of_reverse_acceptance_prob = 0.0L;
	// assert(acceptance_prob <= 0.0);

	const sbm :: Cluster *CL = s.labelling.clusters.at(k);
	For(mp, CL->get_members()) {
		const int32_t m = *mp;
		if(n==m) {
			assert(!obj->selfloops);
			continue;
		}
		const sbm :: State :: point_type neighbours_position = s.cluster_to_points_map.at(m);
		const long double old_dist_2 = original_location.dist_2(neighbours_position);
		const long double dist_2 = proposed_new_location.dist_2(neighbours_position);
		const long double sum_weight_on_this_rel = s.sum_weights_BOTH_directions(n,m);

		// I'll first make the calculations as if there are no edges, then correct for that if necessary
		const long double disconnected_prob = -log2l(1+exp2l(sbm :: ls_alpha_k - dist_2));
		acceptance_prob += disconnected_prob * (obj->directed ? 2 : 1);
		const long double old_disconnected_prob = -log2l(1+exp2l(sbm :: ls_alpha_k - old_dist_2));
		sort_of_reverse_acceptance_prob += old_disconnected_prob * (obj->directed ? 2 : 1);

		if(sum_weight_on_this_rel > 0) {
			const long double this_is_what_is_actually_subtracted = -log2l(1+exp2l(dist_2-sbm :: ls_alpha_k)) - disconnected_prob;
			acceptance_prob += this_is_what_is_actually_subtracted * sum_weight_on_this_rel;

			const long double old_this_is_what_is_actually_subtracted = -log2l(1+exp2l(old_dist_2-sbm :: ls_alpha_k)) - old_disconnected_prob;
			sort_of_reverse_acceptance_prob += old_this_is_what_is_actually_subtracted * sum_weight_on_this_rel;
		}
	}

	if(1)
	{ // sanity check
		assert(s.cluster_to_points_map.at(n) == original_location);
		const long double current_likelihood = my_likelihood(n, s , obj, k);
		s.cluster_to_points_map.at(n) = proposed_new_location;
		const long double new_likelihood = my_likelihood(n, s , obj, k);
		assert(VERYCLOSE(new_likelihood    , new_prior_with_edges + acceptance_prob));
		assert(VERYCLOSE(current_likelihood, old_prior_with_edges + sort_of_reverse_acceptance_prob));
		assert(VERYCLOSE(new_likelihood - current_likelihood, new_prior_with_edges - old_prior_with_edges + acceptance_prob - sort_of_reverse_acceptance_prob));
		s.cluster_to_points_map.at(n) = original_location;
	}
	// PP(acceptance_prob);

	// NOTE on Metropolis-Hastings here.
	// The proposal probability is the prior, and hence instead of
	// dividing by the prior, we just ignore it when calculating the posterior density
	assert(isfinite(acceptance_prob - sort_of_reverse_acceptance_prob));
	if( log2l(gsl_ran_flat(r,0,1)) < acceptance_prob - sort_of_reverse_acceptance_prob ) {
		s.cluster_to_points_map.at(n) = proposed_new_location;
		AR->notify(true);
		return new_prior_with_edges - old_prior_with_edges + acceptance_prob - sort_of_reverse_acceptance_prob;
	} else {
		AR->notify(false);
		return 0.0;
	}
}

static long double update_ls_positions(sbm :: State &s, const sbm :: ObjectiveFunction *obj, AcceptanceRate *ar, gsl_rng * r) {
	assert(args_info.scf_flag == 0);
	// - take each cluster in turn.
	//   - take each node in that cluster in turn
	//     - propose a new position based on it cluster-mate-NEIGHBOURS
	//     - calculate proposal prob, again based on cluster-mate-NEIGHBOURS
	//     - calculate the posterior, based on ALL-cluster-mates

	// const long double pre = s.pmf(obj);
	assert(s._N == (int)s.cluster_to_points_map.size());
	assert(!obj->weighted);
	assert(!obj->selfloops); // I think I'll never like this!

	// update each node in turn, move it towards the cluster-mates it's connected to, and away from those it's not connected to
	long double l2_delta_bits = 0.0;
	for(int k = 0; k<s._k; ++k) {
		const sbm :: Cluster *CL = s.labelling.clusters.at(k);
		const std :: list<int> & mem = CL->get_members();
		For(n, mem) {
			//l2_delta_bits += update_one_nodes_position(n, s, obj, ar, r);
			l2_delta_bits += MH_update_one_nodes_position(*n, s, obj, ar, r);
		}
	}

	// const long double post = s.pmf(obj);
	// PP2(post - pre , l2_delta_bits);
	// assert(VERYCLOSE(post - pre , l2_delta_bits));
	return l2_delta_bits;
}
static pair<long double, long double> prepare_two_M3_ls_proposals(sbm :: State &s, const sbm :: ObjectiveFunction *obj, const int n, const int left, const int right) {
	assert(!obj->weighted && !obj->selfloops && !obj->directed);
	assert((int)s.cluster_to_points_map.size() == s._N);

	// given a node, and the two candidate communities:
	//    place it at the mean of its neighbours
	//    evaluate the overall pmf() for both options
	//    normalize the pmfs
	//
	const int temp_cluster_id = s._k-1;
	assert(s.labelling.cluster_id.at(n) == temp_cluster_id);
	s.moveNodeAndInformOfEdges(n, left);
	ANormalDistribution normpdf_left(n, s, obj); // density will be near the neighbours of n, and it incorporates the prior sigma
	s.cluster_to_points_map.at(n) = normpdf_left.mean;
	long double pmf_left = s.pmf(obj);

	s.moveNodeAndInformOfEdges(n, temp_cluster_id);

	assert(s.labelling.cluster_id.at(n) == temp_cluster_id);
	s.moveNodeAndInformOfEdges(n, right);
	ANormalDistribution normpdf_right(n, s, obj); // density will be near the neighbours of n, and it incorporates the prior sigma
	s.cluster_to_points_map.at(n) = normpdf_right.mean;
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
static long double M3_LS(sbm :: State &s, const sbm :: ObjectiveFunction *obj, AcceptanceRate *AR, gsl_rng *r) {
	assert(args_info.scf_flag == 0);
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
	const std :: vector<sbm :: State :: point_type> old_cluster_to_points_map = s.cluster_to_points_map;

	// move into the temporary cluster
	vector<int> both_sets_of_nodes;
	const int temp_cluster_id = s._k;
	s.appendEmptyCluster();
	sbm :: State :: point_type pzero; pzero.zero();
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
	For(one_node, both_sets_of_nodes) {
		pair<long double, long double> two_options = prepare_two_M3_ls_proposals(s, obj, *one_node, left, right);
		assert(VERYCLOSE(two_options.first + two_options.second, 1.0));
		const int orig_cluster_id = old_cluster_ids.at(*one_node);
		if(orig_cluster_id == left) {
			l2_reverse_prop += log2l(two_options.first);
		} else {
			assert(orig_cluster_id == right);
			l2_reverse_prop += log2l(two_options.second);
		}
		s.moveNodeAndInformOfEdges(*one_node, orig_cluster_id);
		ANormalDistribution normpdf_selected(*one_node, s, obj);
		l2_reverse_prop += normpdf_selected.pdf_prop(s.cluster_to_points_map.at(*one_node));
	}
	// PP(l2_reverse_prop);

	{ // revert for verification, and then put things back
		s.deleteClusterFromTheEnd();
		assert(s.cluster_to_points_map.size() == old_cluster_to_points_map.size());
		s.cluster_to_points_map = old_cluster_to_points_map;
		const long double verify_pre = s.pmf(obj);
		assert(VERYCLOSE(pre, verify_pre));
		s.appendEmptyCluster();
	}

	// we must move them out again, in order to do our random proposal
	For(one_node, both_sets_of_nodes) {
		s.moveNodeAndInformOfEdges(*one_node, temp_cluster_id);
	}
	// at this stage, the positions are as before, except with the extra dummy positions for the temp cluster

	// finally, doing the random proposal
	long double l2_random_prop = 0.0L;
	For(one_node, both_sets_of_nodes) {
		pair<long double, long double> two_options = prepare_two_M3_ls_proposals(s, obj, *one_node, left, right);
		assert(VERYCLOSE(two_options.first + two_options.second, 1.0));
		if(gsl_ran_flat(r,0,1) < two_options.first) {
			l2_random_prop += log2l(two_options.first);
			s.moveNodeAndInformOfEdges(*one_node, left);
			assert(left == s.labelling.cluster_id.at(*one_node));
		} else {
			l2_random_prop += log2l(two_options.second);
			s.moveNodeAndInformOfEdges(*one_node, right);
			assert(right == s.labelling.cluster_id.at(*one_node));
		}
		ANormalDistribution normpdf_selected(*one_node, s, obj);
		sbm :: State :: point_type new_position = normpdf_selected.draw(r);
		l2_random_prop += normpdf_selected.pdf_prop(new_position);
		s.cluster_to_points_map.at(*one_node) = new_position;
		// PP2(__LINE__, l2_random_prop);
	}
	// PP2(__LINE__, l2_random_prop);
	// PP2( s.labelling.clusters.at(left)->order(), s.labelling.clusters.at(right)->order());

	// delete the temporary community, we don't need it any more
	s.deleteClusterFromTheEnd();
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
		For(i, both_sets_of_nodes) {
			const int orig_cluster_id = old_cluster_ids.at(*i);
			if(s.labelling.cluster_id.at(*i) != orig_cluster_id) {
				s.moveNodeAndInformOfEdges(*i, orig_cluster_id);
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
	assert(args_info.scf_flag == 0);
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
static long double MetropolisOnK(sbm :: State &s, const sbm :: ObjectiveFunction *obj __attribute__((unused)), gsl_rng* r, AcceptanceRate *AR) {
	assert(s.cluster_to_points_map.empty());
	/// const long double prePMF = s.pmf(obj);
	/// const long double prePMF12 = s.P_z_K();
	const int preK = s._k;
	if(fiftyfifty()) { // propose increase in K
		if(args_info.maxK_arg > 0 && s._k >= args_info.maxK_arg) {
			return 0.0;
		}
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
			( args_info.uniformK_flag ? 0 : args_info.geometricK_flag ? -1 : - log2(preK+1) )

			// + LOG2GAMMA(s._alpha)  // the change to SumOfLog2Facts
			// I belatedly noticed that this last line can be cancelled against some expressions in the next line.

			// log2(preK) - log2(s._N+preK)
			+(LOG2GAMMA(postK * s._alpha) - LOG2GAMMA(postK * s._alpha + s._N) /*- postK*LOG2GAMMA(s._alpha)*/ )
			-(LOG2GAMMA(preK * s._alpha) - LOG2GAMMA(preK  * s._alpha + s._N) /*- preK*LOG2GAMMA(s._alpha)*/ )
			;
		bool provisional_accept = acceptTest(presumed_delta);
		if(provisional_accept && args_info.scf_flag) {
			provisional_accept = drawPiAndTest(s, obj, r);
		}
		if(provisional_accept) {
			// cout << "k: acc inc" << endl;
			assert(s._k>preK);
			// PP(s.pmf(obj) - prePMF);
			// PP2(preK, s._k);
			// PP(s._N);
			// PP( presumed_delta );
			// PP(postPMF12 - prePMF12);
			// assert(VERYCLOSE(presumed_delta, s.pmf(obj) - prePMF));

			AR->notify(true);

			// 2012-12-02: We shouldn't just add the new one at the end.
			const int new_cluster_id = static_cast<int>(drand48() * s._k);
			if(new_cluster_id != s._k-1) {
				s.swapClusters(new_cluster_id, s._k-1);
			}

			return presumed_delta;
		} else {
			// cout << "k: rej inc" << endl;
			s.deleteClusterFromTheEnd();
			// assert(s.pmf(obj)==prePMF);
			// assert(s.P_z_K()==prePMF12);
			assert(s._k==preK);
			AR->notify(false);
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
				( args_info.uniformK_flag ? 0 : args_info.geometricK_flag ?  1 : + log2(preK) )

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
			bool provisional_accept = acceptTest(presumed_delta);
			if(provisional_accept && args_info.scf_flag) {
				provisional_accept = drawPiAndTest(s, obj, r);
			}
			if(provisional_accept) {
				assert(s._k<preK);
				AR->notify(true);
				return presumed_delta;
			} else {
				// actually, I think we'll never get here, except with weird priors on K and alpha maybe? !
				s.appendEmptyCluster();
				assert(s._k==preK);
				AR->notify(false);
				if(clusterToProposeDelete != s._k-1) {
					s.swapClusters(clusterToProposeDelete , s._k-1);
				}
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
		const int j2 = static_cast<int>(drand48() * s._k);
must_be_different:
		const int j1 = static_cast<int>(drand48() * s._k);
		if(j1==j2)
			goto must_be_different;
		const int n_j1 = s.labelling.clusters.at(j1)->order();
		const int n_j2 = s.labelling.clusters.at(j2)->order();
		const std :: vector<int> j2_members(s.labelling.clusters.at(j2)->get_members().begin(), s.labelling.clusters.at(j2)->get_members().end());
		// apply a merge
		For(i, j2_members) {
			const int shouldBej2 = s.moveNodeAndInformOfEdges(*i, j1);
			assert(j2 == shouldBej2);
		}
		assert(0==s.labelling.clusters.at(j2)->order());
		const int n_jBoth = s.labelling.clusters.at(j1)->order();
		assert(n_jBoth == n_j1 + n_j2);
		// 2012-12-02: We shouldn't just add the new one at the end.
		if(j2 != s._k-1) { s.swapClusters(j2, s._k-1); }
		assert(0==s.labelling.clusters.at(s._k-1)->order());
		s.deleteClusterFromTheEnd();
		const long double post = s.pmf(obj);
		// PP3(pre,post, post-pre);

		const long double full_proposal_probability = EjectAbsorb_prop_prob_merge(s._k);
		const long double reverse_proposal_probability = EjectAbsorb_prop_prob_split(s._k, n_jBoth, n_j1, n_j2);

		const long double acceptance_probability = post-pre - full_proposal_probability + reverse_proposal_probability;

		bool provisional_accept = acceptTest(acceptance_probability);
		if(provisional_accept && args_info.scf_flag) {
			provisional_accept = drawPiAndTest(s, obj, r);
		}
		if(provisional_accept) { // accept
			AR->notify(true);
			return post-pre;
		} else { // reject, split 'em again
			AR->notify(false);
			s.appendEmptyCluster();
			if(j2 != s._k-1) { s.swapClusters(j2, s._k-1); }
			For(j2_member, j2_members) {
				const int shouldBej1 = s.moveNodeAndInformOfEdges(*j2_member, j2);
				assert(j1 == shouldBej1);
			}
			return 0.0;
		}
		assert(1==2); // should never get here
	} else {
		if(args_info.maxK_arg > 0 && s._k >= args_info.maxK_arg) {
			return 0.0L;
		}
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
		bool provisional_accept = acceptTest(acceptance_probability);
		if(provisional_accept && args_info.scf_flag) {
			provisional_accept = drawPiAndTest(s, obj, r);
		}
		if(provisional_accept) { // accept
			AR->notify(true);
			{
				// 2012-12-02: We shouldn't just add the new one at the end.
				const int new_cluster_id = static_cast<int>(drand48() * s._k);
				if(new_cluster_id != s._k-1) {
					s.swapClusters(new_cluster_id, s._k-1);
				}
			}
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
		assert(1==2); // should never get here
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

static void recursive(const int deciding, const int K
		, vector<int> &new_names, vector<bool> &already_taken
		, const vector<bool> &is_currently_empty
		, const vector< vector<int> > & kbyk
		, int & best_score_so_far
		, vector<int> & best_relabelling_so_far
		, const int partial_score
		) {
	if(partial_score >= best_score_so_far) {
		// never gonna improve, may as well quit now
		return;
	}
	assert((int)new_names.size() == K);
	assert((int)kbyk.size() == K);
	if(deciding == K) {
		int score = 0;
		// cout << "leaf: ";
		if(best_score_so_far > score) {
			best_score_so_far = score;
			best_relabelling_so_far = new_names;
		}
		return;
	}
	if(is_currently_empty.at(deciding)) {
		// we don't care about finding a label for the empty clusters
		new_names.at(deciding) = -1;
		recursive(deciding+1, K, new_names, already_taken, is_currently_empty, kbyk, best_score_so_far, best_relabelling_so_far, partial_score);
		return;
	}
	for(int new_k = 0; new_k < K; ++new_k) {
		if(already_taken.at(new_k)) // could this be made more efficient? keeping a linked list of available ids?
			continue;
		new_names.at(deciding) = new_k;
		already_taken.at(new_k) = true;
		recursive(deciding + 1, K, new_names, already_taken, is_currently_empty, kbyk
				, best_score_so_far, best_relabelling_so_far
				, partial_score + kbyk.at(deciding).at(new_k)
				);
		already_taken.at(new_k) = false;
	}
}
static vector<int> calculate_best_relabelling(const vector<int> & z, const vector< vector<int> > & relab_freq) {
	const size_t N = relab_freq.size();
	assert(N>0);
	const size_t K = relab_freq.front().size();
	vector<bool> is_currently_empty(K,true);

	// Consider the clustering in z, and construct a K-by-K matrix of similarity
	// THEN, convert it to a *distance* matrix
	vector< vector<int> > kbyk(K, vector<int>(K,0));
	for(size_t n=0; n<N; ++n) {
		const int current_z_n = z.at(n);
		is_currently_empty.at(current_z_n) = false;
		const vector<int> & relab_freq_n = relab_freq.at(n);
		for(size_t new_k=0; new_k<K; ++new_k) {
			kbyk.at(current_z_n).at(new_k) += relab_freq_n.at(new_k);
		}
	}
	/*
	For(row, kbyk) {
		For(cell, *row) {
			cout << ' '
				<< stack.push << fixed << setw(6)
				<< *cell
				<< stack.pop;
		}
		cout << endl;
	}
	*/

	// First, we'll try the optimistic approach,
	// setting each cluster to its best match,
	// where we hope that they all map to distinct clusters

	{
		vector<int> new_names(K, -1);
		vector<bool> already_taken(K);
		bool optimism_failed = false;
		for(size_t old_k = 0; old_k<K; ++old_k) {
			if(!is_currently_empty.at(old_k)) {
				const vector<int> & row = kbyk.at(old_k);
				const size_t best_fit = max_element(row.begin()        , row.end()        ) - row.begin();
				new_names.at(old_k) = best_fit;
				if( already_taken.at(best_fit) ) {
					optimism_failed = true;
					break; // gotta give up on the optimistic route
				}
				already_taken.at(best_fit) = true;
			}
		}
		if(!optimism_failed) {
			// the optimistic route succeeded, applying the relabelling and return.
			return new_names;
		}
	}

	{ // we're here because the optimistic route didn't work,
	  // so we have to try the hard way
	  // instead of maximizing similarity, let's minimize dissimilarity
	  //
	  // the minimum dissimilarity is zero, so let's subtract from the rows
		for(size_t k=0; k<K; ++k) {
			vector<int> * row = &kbyk.at(k);
			const int largest_element = *max_element(row->begin(), row->end());
			if(is_currently_empty.at(k)) {
				assert(largest_element==0);
				For(cell, *row) {
					assert(*cell == 0);
				}
			} else {
				assert(largest_element>0);
			For(cell, *row) {
				*cell = largest_element - *cell;
			}
			}
		}
		vector<int> new_names(K);
		vector<bool> already_taken(K);
		int best_score_so_far = INT_MAX;
		vector<int> best_relabelling_so_far;
		recursive(0, K, new_names, already_taken, is_currently_empty, kbyk, best_score_so_far, best_relabelling_so_far, 0);
		assert(best_relabelling_so_far.size() == K);
		return best_relabelling_so_far;
	}
}
static void assert_valid_relabelling(const int K, const vector<int> &relabelling, int lineno) {
	assert((int)relabelling.size() == K);
	set<int> used_ids;
	For(r, relabelling) {
		const int id = *r;
		if(id!=-1) {
			assert(id>=0);
			assert(id< K);
			bool was_inserted =
				used_ids.insert(id)
				.second
				;
			if(!was_inserted) {
				PP2(lineno, id);
				assert(was_inserted);
			}
		}
	}
}

static void label_switch(
		const graph :: NetworkInterfaceConvertedToStringWithWeights *g
		, const sbm :: State &
		, const size_t N, const size_t K
		, deque< pair< pair<int,int> , vector<int> > > & all_burned_in_z
		, const vector<int> * const groundTruth
		) {
// TODO
//  - subtract column-wise, instead of row-wise?

	// Note: K includes the empty clusters.  This is necessary as the z_n may be as high as K-1
	assert(N>=1);
	assert(K>=1);
	assert(!all_burned_in_z.empty());
	// already sorted by the number of non-empty clusters
	vector< vector<int> > relab_freq(N, vector<int>(K, 0) );
	int num_states_processed = 0;
	For(one_state, all_burned_in_z) {
		if(num_states_processed % 1000 == 0)
			cout
				<< " .. label-switching: "
				<< num_states_processed
				<< "/" << all_burned_in_z.size()
				<< "\t(after " << ELAPSED() << " seconds)"
				<< endl;
		++num_states_processed;
		{ // verify the number of non-empty clusters
			set<int> non_empty_clusters;
			for(size_t n=0; n<N; ++n) {
				assert(one_state->second.size() == N);
				const size_t z_n = one_state->second.at(n);
				assert(z_n < K);
				non_empty_clusters.insert(z_n);
			}
			assert( (int) non_empty_clusters.size() == one_state -> first . second);
		}
		if(one_state == all_burned_in_z.begin()) {
			// leave the first one as-is
		} else {
			const vector<int> relabelling = calculate_best_relabelling(one_state->second, relab_freq);
			assert_valid_relabelling(K, relabelling, __LINE__);
			vector<int> & z = one_state -> second;
			for(size_t n=0; n<N; ++n) {
				const int current_z_n = z.at(n);
				const int new_z_n = relabelling.at(current_z_n);
				assert(new_z_n >= 0);
				assert(new_z_n <  (int)K);
				z.at(n) = new_z_n;
			}
		}
		// Now, the 'current' kz has been relabelled.  We must add its details to the counts
		for(size_t n=0; n<N; ++n) {
			assert(one_state->second.size() == N);
			const size_t z_n = one_state->second.at(n);
			assert(z_n < K);
			++ relab_freq.at(n).at(   z_n  );
		}
	}

	// now, we'll relabel either to the ground truth, or,
	//   if there's no ground truth, put the largest clusters first
	if(groundTruth) {
		assert(groundTruth->size() == N);
		vector<int> relabelling = calculate_best_relabelling(*groundTruth, relab_freq);
		{ // if there are any -1s (for empty clusters) we should fill
		  // them in with unused ids
			assert_valid_relabelling(K, relabelling, __LINE__);
			set<int> unused_ids;
			for(size_t k=0; k<K; ++k)
				unused_ids.insert(k);
			assert(unused_ids.size() == K);
			for(size_t k=0; k<K; ++k) {
				if(relabelling.at(k) != -1) {
					const bool was_erased = unused_ids.erase(relabelling.at(k));
					assert(was_erased);
				}
			}
			for(size_t k=0; k<K; ++k) {
				if(relabelling.at(k) == -1) {
					assert(!unused_ids.empty());
					const int use_this_id_for_empty = *unused_ids.begin();
					unused_ids.erase(unused_ids.begin());
					relabelling.at(k) = use_this_id_for_empty;
				}
			}
			assert(unused_ids.empty());
		}
		// GT cluster k is closest to found cluster relabelling.at(k)
		for(size_t n=0; n<N; ++n) {
			vector<int> & one_node = relab_freq.at(n);
			vector<int> new_row_from_relab_freq(K);
			for(size_t k=0; k<K; ++k) {
				assert(relabelling.at(k) >= 0);
				assert(relabelling.at(k) <  (int)K);
				new_row_from_relab_freq.at( k ) = one_node.at( relabelling.at(k) );
			}
			new_row_from_relab_freq.swap( one_node );
		}
	} else {
		vector< pair<int64_t,size_t> > size_of_clusters(K);
		for(size_t k=0; k < K; ++k) {
			size_of_clusters.at(k).second = k;
		}
		For(row, relab_freq) {
			for(size_t k=0; k < K; ++k) {
				size_of_clusters.at(k).first += row->at(k);
			}
		}
		sort(size_of_clusters.begin(), size_of_clusters.end());
		reverse(size_of_clusters.begin(), size_of_clusters.end());
		for(size_t n=0; n<N; ++n) {
			vector<int> & one_node = relab_freq.at(n);
			vector<int> new_row_from_relab_freq(K);
			for(size_t k=0; k<K; ++k) {
				new_row_from_relab_freq.at( k ) = one_node.at( size_of_clusters.at(k).second );
			}
			new_row_from_relab_freq.swap( one_node );
		}
	}

	cout << " label-switching complete. (after " << ELAPSED() << " seconds)" << endl;

	{ // print out the summary output from label-switching
#define FIELD_WIDTH 9
			cout << endl;
			cout << "Each node, and its favourite cluster(s):" << endl;
			cout
					<< stack.push << fixed << setw(FIELD_WIDTH)
					<< "nodeid"
					<< stack.pop;
			cout
					<< stack.push << fixed << setw(FIELD_WIDTH)
					<< "cluster"
					<< stack.pop;
			cout
					<< ' ' << stack.push << fixed << setw(FIELD_WIDTH)
					<< "nodename"
					<< stack.pop;
			cout
					<< "  favourite clusters, i.e. P(z_ik > 0.01)";
			cout << endl;

		// the last few lines printed the column names. Next, we print the data

		for(size_t n=0; n<N; ++n) {
			// first, create a nicely-formatted version of the
			// node*name*, the string/int that appeared in the edge-list file
			ostringstream node_name;
			if(g->was_string_data_not_int() ) {
				node_name << '"' << g->node_name_as_string(n) << '"';
			} else {
				node_name << g->node_name_as_string(n);
			}
			cout
					<< stack.push << fixed << setw(FIELD_WIDTH)
					<< n
					<< stack.pop;
			cout
					<< stack.push << fixed << setw(FIELD_WIDTH)
					<< max_element(relab_freq.at(n).begin(), relab_freq.at(n).end()) - relab_freq.at(n).begin()
					<< stack.pop;
			cout
					<< ' ' << stack.push << fixed << setw(FIELD_WIDTH)
					<< node_name.str()
					<< stack.pop;
			cout << "  ";
			bool first = true;
			vector< pair<double, int> > ziks;
			for(size_t k=0; k<K; ++k) {
				const double zik = double(relab_freq.at(n).at(k)) / all_burned_in_z.size();
				ziks.push_back( make_pair( zik, k) );
			}
			sort(ziks.begin(), ziks.end());
			reverse(ziks.begin(), ziks.end());
			For(zikk, ziks) {
				const double zik = zikk->first; // double(relab_freq.at(n).at(k)) / all_burned_in_z.size();
				const int k = zikk->second;
				if(zik > 0.01) {
					ostringstream oss;
					oss << k << "(" << fixed << setprecision(1) << 100.0 * zik << "%)";
					// cout << stack.push << setw(10) << left << oss.str() << stack.pop;
					if(!first)
						cout << ',';
					first = false;
					cout << oss.str();
				}
			}
			cout << endl;
		}
	}
	{ // print out the summary output from label-switching
		cout << endl;
		cout << "Best single clustering:" << endl;
		for(size_t n=0; n<N; ++n) {
			const vector<int> & relab_freq_n = relab_freq.at(n);
			if(n>0)
				cout << ',';
			cout << max_element(relab_freq_n.begin(), relab_freq_n.end()) - relab_freq_n.begin();
		}
		cout << endl;
		cout << endl;
		cout << "Summary of cluster sizes across all label-switched states:" << endl;
		vector<double> cluster_sizes(K,0.0);
		for(size_t k=0; k<K; ++k) {
			size_t sum_z_k = 0;
			for(size_t i=0; i<N; ++i) {
				sum_z_k += relab_freq.at(i).at(k);
			}
			cout << k << "\t"
				<< stack.push << fixed << setw(7) << setprecision(3)
				<< double(sum_z_k) / double(all_burned_in_z.size())
				<< stack.pop
				<< endl;
			cluster_sizes.at(k) = double(sum_z_k) / double(all_burned_in_z.size());
		}
		// For each block, we want to know how many edges are in it
		vector< vector<double> > rate_of_edges_z(K, vector<double>(K,0) );
#define MAX_CLUSTERS_FOR_BLOCK_SUMMARY 15
		for(int rel = 0; rel<g->numRels(); rel++) {
			const std :: pair <int32_t, int32_t> eps = g->get_plain_graph()->EndPoints(rel);
			const int32_t left = eps.first;
			const int32_t right = eps.second;
			assert(left<=right);
			const double l2r = g->get_edge_weights()->getl2h(rel);
			const double r2l = g->get_edge_weights()->geth2l(rel);
			if(left==right)
				assert(r2l==0);
			assert(l2r + r2l > 0.0);
			for(size_t k=0; k<K; ++k) {
				if(relab_freq.at(left).at(k) == 0 && relab_freq.at(right).at(k) == 0) continue;
			for(size_t l=0; l<K; ++l) {
				if(relab_freq.at(left).at(l) == 0 && relab_freq.at(right).at(l) == 0) continue;
				const double l2r_k_to_l = double(relab_freq.at(left).at(k))
					* double(relab_freq.at(right).at(l))
					/ double(all_burned_in_z.size())
					/ double(all_burned_in_z.size());
				double r2l_k_to_l = double(relab_freq.at(right).at(k))
					* double(relab_freq.at(left).at(l))
					/ double(all_burned_in_z.size())
					/ double(all_burned_in_z.size());
				if(args_info.directed_flag) {
					rate_of_edges_z.at(k).at(l) += l2r * l2r_k_to_l;
					rate_of_edges_z.at(k).at(l) += r2l * r2l_k_to_l;
				} else {
					if(k != l) {
						rate_of_edges_z.at(k).at(l) += l2r * l2r_k_to_l;
						rate_of_edges_z.at(k).at(l) += r2l * r2l_k_to_l;
						rate_of_edges_z.at(l).at(k) += l2r * l2r_k_to_l;
						rate_of_edges_z.at(l).at(k) += r2l * r2l_k_to_l;
					}
					if(k==l) {
						rate_of_edges_z.at(k).at(l) += l2r * l2r_k_to_l;
						rate_of_edges_z.at(k).at(l) += r2l * r2l_k_to_l;
					}
				}
			}
			}
		}
		cout << endl;
		cout << "Summary of block densities:" << endl;
		vector< vector<double> > num_pairs(MAX_CLUSTERS_FOR_BLOCK_SUMMARY, vector<double>(MAX_CLUSTERS_FOR_BLOCK_SUMMARY,0) );
		{
			// First, we'll do calculations where k!=l - these are easy.
			for(size_t k=0; k<K; ++k) {
				if(k >= MAX_CLUSTERS_FOR_BLOCK_SUMMARY) break;
				for(size_t l=0; l<K; ++l) {
					if(l >= MAX_CLUSTERS_FOR_BLOCK_SUMMARY) break;
					if(k!=l)
						num_pairs.at(k).at(l) = cluster_sizes.at(k) * cluster_sizes.at(l);
				}
			}
			// Now, we'll do the diagonal blocks - where we need to be careful and self loops and direction.
			for(size_t k=0; k<K; ++k) {
				if(k >= MAX_CLUSTERS_FOR_BLOCK_SUMMARY) break;
				double self_pairs_in_this_cluster = 0.0;
				for(size_t n=0; n<N; ++n) {
					self_pairs_in_this_cluster +=
						double(relab_freq.at(n).at(k))
						/ double(all_burned_in_z.size())
						* double(relab_freq.at(n).at(k))
						/ double(all_burned_in_z.size())
					;
				}
				const long double non_self_pairs = (cluster_sizes.at(k) * cluster_sizes.at(k) - self_pairs_in_this_cluster) / (args_info.directed_flag ? 1.0 : 2.0);
				num_pairs.at(k).at(k) = non_self_pairs + (args_info.selfloop_flag ? self_pairs_in_this_cluster : 0.0);
			}
		}
		double verify_edges = 0;
		for(size_t k=0; k<K; ++k) {
			if(k>=MAX_CLUSTERS_FOR_BLOCK_SUMMARY) {
					cout << " ... too many blocks to print all of them. I won't print any more. ..." << endl;
					break;
			}
			cout << k << "\t";
			for(size_t l=0; l<K; ++l) {
				if(l>=MAX_CLUSTERS_FOR_BLOCK_SUMMARY) {
					cout << " ..." << endl;
					break;
				}
				if(args_info.directed_flag || l<=k)
					verify_edges += rate_of_edges_z.at(k).at(l);
				cout << " |";
				cout
					<< stack.push << fixed << setw(5) << setprecision(1)
					<< rate_of_edges_z.at(k).at(l) << " /" <<  setw(5) << num_pairs.at(k).at(l);
				if(isfinite(100.0 * rate_of_edges_z.at(k).at(l) /  num_pairs.at(k).at(l))) {
					if(!args_info.weighted_flag && rate_of_edges_z.at(k).at(l) >  num_pairs.at(k).at(l)) {
						assert(VERYCLOSE(rate_of_edges_z.at(k).at(l) , num_pairs.at(k).at(l)));
						rate_of_edges_z.at(k).at(l) = num_pairs.at(k).at(l);
					}
					cout << " =";
					if(args_info.weighted_flag) {
						cout << setw(6);
						cout <<  rate_of_edges_z.at(k).at(l) /  num_pairs.at(k).at(l);
					} else {
						cout << setw(5);
						cout <<  100.0 * rate_of_edges_z.at(k).at(l) /  num_pairs.at(k).at(l);
						cout << "%";
					}
				} else
					cout << "        ";
				cout << stack.pop;
			}
			cout << endl;
		}
	}
}

sighandler_t original_ctrl_C_handler = NULL;
bool was_ctrl_C_caught_in_MCMC = false;
void sig_ctrl_C_caught_in_MCMC(int) {
	cerr << "Ctrl-C received. Exiting out of MCMC early." << endl;
	was_ctrl_C_caught_in_MCMC = true;

}

long double entropy_of_a_clustering(const sbm :: State &s) {
	// summing over the non-empty clusters
	long double entropy = 0.0L;
	for(int k=0; k<s._k; ++k) {
		const int o = s.labelling.clusters.at(k)->order();
		if(o>0) {
			const long double proportion = (long double)o / (long double)s._N;
			entropy -= proportion * logl(proportion);
		}
	}
	return entropy;
}

#define CHECK_PMF_TRACKER(track, actual) do { const long double _actual = (actual); long double & _track = (track); if(VERYCLOSE(_track,_actual)) { track = _actual; } else { PP(_actual - track); } assert(_track == _actual); } while(0)

static void runSBM(const graph :: NetworkInterfaceConvertedToStringWithWeights *g, const int commandLineK, const sbm :: ObjectiveFunction * const obj, const bool initializeToGT, const vector<int> * const groundTruth, const int iterations, const  gengetopt_args_info &args_info, gsl_rng *r) {
	if(g->get_plain_graph()->number_of_self_loops() > 0 && !obj->selfloops ){
		cerr << endl << "Error: You must specify the -s flag to fully support self-loops. Your network has " << g->get_plain_graph()->number_of_self_loops() << " self-loops." << endl;
		exit(1);
	}
	sbm :: State s(g, args_info.mega_flag, args_info.alpha_arg, args_info.latentspace_flag);

	s.internalCheck();
	cout
		<< "nodes: " << g->numNodes()
		<< "\t\tedges: " << s.total_edge_weight
		<< endl;

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
		while(commandLineK > s._k) {
			s.appendEmptyCluster();
		}
		if(commandLineK != -1 && commandLineK != s._k) {
			cerr << endl
				<< "Error: Problem using the --initGT and -K args. The GT.vector file \""
				<< args_info.GT_vector_arg
				<< "\" has " << s._k << " clusters, "
				<< "but you requested -K " << commandLineK << " clusters."
				<< endl;
			exit(1);

		}
	} else {
		if(commandLineK != -1)
			randomize(s, commandLineK);
	}
	if(args_info.verbose_flag) {
		s.shortSummary(obj, groundTruth); /*s.summarizeEdgeCounts();*/ s.blockDetail(obj);
	}
	s.internalCheck();
	if(args_info.latentspace_flag) {
		assert(s._k == commandLineK);
		assert(s.cluster_to_points_map.empty());
		s.cluster_to_points_map.resize(s._N);

		For(dbl, s.cluster_to_points_map) {
			for(int d = 0; d < sbm :: State :: point_type :: dimensionality; ++d) {
				dbl->at(d) = gsl_ran_gaussian(r, 1);
			}
		}
		cout << "assigned initial random positions for the latent space model" << endl;
	}
	if(args_info.verbose_flag) {
		s.shortSummary(obj, groundTruth); /*s.summarizeEdgeCounts();*/ s.blockDetail(obj);
	}
	s.internalCheck();

	long double pmf_track = s.pmf(obj);

	AcceptanceRate AR_metroK("metroK");
	AcceptanceRate AR_metro1Node("metro1Node");
	AcceptanceRate AR_gibbs("gibbs");
	AcceptanceRate AR_ea  ("EjectAbsorb");
	AcceptanceRate AR_lspos  ("LSSBM positions");
	AcceptanceRate AR_M3lspos  ("LSSBM M3");

	// some variables to check the PMP, i.e. the single most-visited state
/*
	map< pair<int, vector<int> >, int> pmp_table;
	pair< pair<int, vector<int> >, int> best_pmp_so_far(make_pair(0, vector<int>()), 0);
	int64_t num_states_checked_for_pmp = 0;
*/

	ofstream * save_z_fstream = NULL;
	ofstream * save_lsz_fstream = NULL;
	if(args_info.save_z_arg[0]) {
		save_z_fstream = new ofstream(args_info.save_z_arg);
	}
	if(args_info.save_lsz_arg[0]) {
		save_lsz_fstream = new ofstream(args_info.save_lsz_arg);
		*save_lsz_fstream << "N:\t" << s._N << "\t\t" << endl;
	}
	deque< pair< pair<int,int>, vector<int> > > all_burned_in_z; // store *half* of the states, for label-switching after
	cout << endl << " = Starting MCMC =  (after " << ELAPSED() << " seconds)" << endl << endl;
	long double lagging_time = ELAPSED();
	int iteration;
	original_ctrl_C_handler = signal(SIGINT, sig_ctrl_C_caught_in_MCMC);
	vector<int> all_nodes_random_order;
	int current_gibbs_position = 0;
	{
		for(int n=0; n<s._N; ++n) {
			all_nodes_random_order.push_back(n);
		}
		assert(s._N == (int)all_nodes_random_order.size());
		random_shuffle(all_nodes_random_order.begin(), all_nodes_random_order.end());
	}
	size_t pushed_back_z=0;
	for(iteration = 0; /*see the 'break's below */; iteration++) {
		if(was_ctrl_C_caught_in_MCMC)
			break;
		if(commandLineK != -1)
			assert(commandLineK == s._k);
		if(args_info.maxK_arg != -1) {
			assert(args_info.maxK_arg > 0);
			assert(s._k <= args_info.maxK_arg);
		}
		{ // should we print the really-short summary?
			// It's every 10 iterations and every minute
			bool has_a_minute_elapsed = false;
			if( int(lagging_time / 60) < int(ELAPSED() / 60) ) {
				lagging_time = ELAPSED();
				// another minute has elapsed
				has_a_minute_elapsed = true;
			}
			if(has_a_minute_elapsed || iteration%10 == 0 || iteration == iterations) {
				cout
					<< " .. iteration: " << iteration
					<< "/" << iterations
					<< "\tnonEmpty: " << s.labelling.NonEmptyClusters
					<< "\tK: " << s._k;
				if(groundTruth) {
					vector<int> z_vector(s._N);
					for (int n=0; n < s._N; n++) {
						const int id_of_cluster = s.labelling.cluster_id.at(n);
						z_vector.at(n) = id_of_cluster;
					}
					assert(z_vector.size() == groundTruth->size());
					cout
						<< "\tnmi: "
						<< sbm :: State :: NMI(z_vector, *groundTruth) * 100;
				}
				cout << "\tentropy: " << entropy_of_a_clustering(s);
				cout << "\t(after " << ELAPSED() << " seconds)";
				cout << endl;
			}
		}
		if(iteration>=iterations) {
			break;
		}
		/*
		cout
			<< "iteration\t" << i
			<< "\telapsed\t" << ELAPSED()
			<< endl;
		cout
			<< "Iteration:\t" << i
			<< "\tk0\t" << s._k
			<< "\tk1\t" << s.labelling.NonEmptyClusters
			<< "\t"
			;
		s.KandClusterSizes();
			*/

		vector<POSSIBLE_MOVES> possible_moves;
				if(commandLineK == -1 && args_info.algo_metroK_arg) {
					possible_moves.push_back(POS_MetroK);
				}
				if(args_info.algo_gibbs_arg) {
					possible_moves.push_back(POS_Gibbs);
				}
				if(s.cluster_to_points_map.empty() && args_info.algo_m3_arg) {
					possible_moves.push_back(POS_M3);
				}
				if(s.cluster_to_points_map.empty() && commandLineK == -1  && args_info.algo_sm_arg          ) {
					possible_moves.push_back(POS_SM_Merge_M3);
					possible_moves.push_back(POS_SM_Split_M3);
				}
				if(s.cluster_to_points_map.empty() && commandLineK == -1  && args_info.algo_cf_arg          ) {
					possible_moves.push_back(POS_SM_Merge_CF);
					possible_moves.push_back(POS_SM_Split_CF);
				}
				if(s.cluster_to_points_map.empty() && commandLineK == -1 && args_info.algo_ejectabsorb_arg) {
					possible_moves.push_back(POS_AE);
				}
				if(s.cluster_to_points_map.empty() && args_info.algo_1node_arg) {
					possible_moves.push_back(POS_MoneNode);
				}
				if(args_info.algo_lspos_arg) {
					possible_moves.push_back(POS_update_ls_positions);
				}
				if(args_info.algo_lsm3_arg) {
					possible_moves.push_back(POS_M3_LS);
				}
		assert(!possible_moves.empty());
		int random_move = static_cast<int>(drand48() * possible_moves.size());
		assert(s.labelling.missing_nodes == 0);
		switch( possible_moves.at(random_move) ) {
			break; case POS_MetroK: // can NOT handle LSSBM
				if(commandLineK == -1 && args_info.algo_metroK_arg) {
					for(int rep = 0; rep<args_info.algo_metroK_arg; ++rep) {
						pmf_track += MetropolisOnK(s, obj, r, &AR_metroK);
					}
				} else {
					assert(1==2);
				}
			break; case POS_Gibbs: // CAN handle LSSBM
				if(args_info.algo_gibbs_arg) {
					for(int rep = 0; rep<args_info.algo_gibbs_arg; ++rep) {
						++ current_gibbs_position;
						while(current_gibbs_position >= s._N)
							current_gibbs_position -= s._N;
						const int n = all_nodes_random_order.at(current_gibbs_position);
						pmf_track += gibbsOneNode(s, obj, &AR_gibbs, r, n);
					}
				}
			break; case POS_M3: // can NOT handle LSSBM
				if(s.cluster_to_points_map.empty() && args_info.algo_m3_arg) {
					for(int rep = 0; rep<args_info.algo_m3_arg; ++rep) {
						pmf_track += M3(s, obj, r);
					}
				} else
					assert(1==2);
			break; case POS_SM_Merge_M3 : // can NOT handle LSSBM
				if(s.cluster_to_points_map.empty() && args_info.algo_sm_arg) {
					for(int rep = 0; rep<args_info.algo_sm_arg; ++rep)
						pmf_track += SM_Merge(s, obj, r, STRATEGY_M3);
				} else
					assert(1==2);
			break; case POS_SM_Split_M3 : // can NOT handle LSSBM
				if(s.cluster_to_points_map.empty() && commandLineK == -1  && args_info.algo_sm_arg) {
					for(int rep = 0; rep<args_info.algo_sm_arg; ++rep)
						pmf_track += SM_Split(s, obj, r, STRATEGY_M3);
				} else
					assert(1==2);
			break; case POS_SM_Merge_CF : // can NOT handle LSSBM
				if(s.cluster_to_points_map.empty() && args_info.algo_cf_arg) {
					for(int rep = 0; rep<args_info.algo_cf_arg; ++rep)
						pmf_track += SM_Merge(s, obj, r, STRATEGY_CF);
				} else
					assert(1==2);
			break; case POS_SM_Split_CF : // can NOT handle LSSBM
				if(s.cluster_to_points_map.empty() && commandLineK == -1  && args_info.algo_cf_arg) {
					for(int rep = 0; rep<args_info.algo_cf_arg; ++rep)
						pmf_track += SM_Split(s, obj, r, STRATEGY_CF);
				} else
					assert(1==2);
			break; case POS_AE: // can NOT handle LSSBM
				if(s.cluster_to_points_map.empty() && commandLineK == -1 && args_info.algo_ejectabsorb_arg) {
					for(int rep = 0; rep<args_info.algo_ejectabsorb_arg; ++rep)
						pmf_track += EjectAbsorb(s, obj, &AR_ea, r);
				} else
					assert(1==2);
			break; case POS_MoneNode: // can NOT handle LSSBM
				if(s.cluster_to_points_map.empty() && args_info.algo_1node_arg) {
					for(int rep = 0; rep<args_info.algo_1node_arg; ++rep)
						pmf_track += MoneNode(s, obj, &AR_metro1Node);
				} else
					assert(1==2);
			break; case POS_update_ls_positions: // CAN handle LSSBM
				for(int rep = 0; rep<args_info.algo_lspos_arg; ++rep) {
					assert(!s.cluster_to_points_map.empty());
					assert(commandLineK == s._k);
					pmf_track += update_ls_positions(s, obj, &AR_lspos, r);
				}
			break; case POS_M3_LS: // CAN handle LSSBM
				for(int rep = 0; rep<args_info.algo_lsm3_arg; ++rep) {
					assert(!s.cluster_to_points_map.empty());
					assert(commandLineK == s._k);
					pmf_track += M3_LS(s, obj, &AR_M3lspos, r);
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
#if 0
		if(iteration%100 == 0 && iteration*2>=iterations) {
			cout << "paranoid";
			cout << "\t" << s._k;
			cout << "\t" << s.labelling.NonEmptyClusters;
			cout << "\t:";
			for(int k=0; k<s._k; ++k) {
				cout << "\t" << s.labelling.clusters.at(k)->order();
			}
			cout << endl;
		}
#endif

		if(iteration % 10 == 0 && args_info.save_current_state_arg) {
			ofstream file(args_info.save_current_state_arg);
			for(int i=0; i<s._N; ++i) {
				if(i>0)
					file << ',';
				file << s.labelling.cluster_id.at(i);
			}
			file << endl;
			file.close();
		}

		if(iteration%args_info.keep_arg == 0)
		{ // store the states, even if we don't do label-switching, as it is needed to calculate highest_K_sampled and so on
			all_burned_in_z.push_back( make_pair( make_pair(s._k, s.labelling.NonEmptyClusters), s.labelling.cluster_id ) );
			++ pushed_back_z;
			// Every *second* iteration, we'll drop the first state.
			// This is to ensure that, at the end of the algorithm we only have the last half of the iterations
			// This might seem like a strange design, but it does enable us to drop out early if the user types Ctrl-C.
			if(pushed_back_z % 2 == 0)
				all_burned_in_z.pop_front();
		}
		if(args_info.verbose_flag && iteration % args_info.printEveryNIters_arg == 0) {
			cout << endl;
			PP(iteration);
			s.shortSummary(obj, groundTruth);
			// PP4(num_states_checked_for_pmp, best_pmp_so_far.second, 100.0*best_pmp_so_far.second/num_states_checked_for_pmp , best_pmp_so_far.first.first);
			// s.summarizeEdgeCounts();
			AR_metroK.dump();
			AR_metro1Node.dump();
			AR_gibbs.dump();
			AR_ea.dump();
			AR_lspos.dump();
			AR_M3lspos.dump();
			s.blockDetail(obj);
			cout << " end of check at iteration==" << iteration << ". (" << double(clock()) / CLOCKS_PER_SEC << " seconds)" << endl;
			CHECK_PMF_TRACKER(pmf_track, s.pmf(obj));
			s.internalCheck();
		}
		if(iteration*2 >= iterations && args_info.save_z_arg[0]) {
			PP(iteration);
			*save_z_fstream << s._k << ':';
			for(int n=0; n<s._N; n++) {
				*save_z_fstream << s.labelling.cluster_id.at(n);
				if(n+1 < s._N)
					*save_z_fstream << ',';
			}
			*save_z_fstream << endl;
		}
		if(iteration*2 >= iterations
				&& save_lsz_fstream
				&& iteration % args_info.printEveryNIters_arg == 0
			) {
			// print the iteration number and also, for each node,
			//   its colour and its position.
			// Note: The first line of this file will record the number of nodes.
			for(int n=0; n < s._N; n++) {
				const std :: string node_name = s._g->node_name_as_string(n);
				*save_lsz_fstream
					<< '\"' << node_name << '\"'
					<< '\t' << s.labelling.cluster_id.at(n);
				sbm :: State :: point_type position = s.cluster_to_points_map.at(n);
				for(int d = 0; d<DIMENSIONALITY; d++)
					*save_lsz_fstream << '\t' << position.at(d);
				*save_lsz_fstream
					<< endl;
			}
		}
	}
	{ // no need for our custom Ctrl-C handler any more, because the MCMC is done.
		if(original_ctrl_C_handler != SIG_ERR) {
			// restore the original
			signal(SIGINT, original_ctrl_C_handler);
		}
	}
	cout << endl << " = MCMC complete =  (after " << ELAPSED() << " seconds)" << endl << endl;

	vector<int> K_freq(10000); // we're not expecting more than 10000 clusters :-)
	vector<int> KnonEmpty_freq(10000);
	int highest_K_sampled = 0;
	int highest_KnonEmpty_sampled = 0;
	{ // check the burned-in iterations, and calculate some summaries
		For(burned_in, all_burned_in_z) {
			const int k_then = burned_in->first.first;
			const int non_empty_k_then = burned_in->first.second;
			K_freq.at( k_then                        ) ++;
			KnonEmpty_freq.at( non_empty_k_then) ++;
			if(highest_K_sampled < k_then)
			   highest_K_sampled = k_then;
			if(highest_KnonEmpty_sampled < non_empty_k_then)
			   highest_KnonEmpty_sampled = non_empty_k_then;
		}
	}
	const size_t burned_in_iters = all_burned_in_z.size();
	assert(highest_K_sampled > 0);
	assert(highest_K_sampled >= highest_KnonEmpty_sampled);
	if(save_z_fstream)
	       save_z_fstream->close();
	if(save_lsz_fstream)
	       save_lsz_fstream->close();
	if(args_info.verbose_flag) {
		s.shortSummary(obj, groundTruth); /*s.summarizeEdgeCounts();*/ s.blockDetail(obj);
	}
	s.internalCheck();
	{
		PP2(highest_K_sampled, highest_KnonEmpty_sampled);
		for(int k=1; k<=highest_K_sampled; ++k) {
			const double PKk = double(K_freq.at(k)) / burned_in_iters;
			if(PKk > 0.001)
				cout << "P(K=" << k << ")=\t"
					<< stack.push << fixed << setw(6) << setprecision(4)
					<< PKk
					<< stack.pop
					<< endl;
		}
		for(int k=1; k<=highest_K_sampled; ++k) {
			const double PKnonemptyk = double(K_freq.at(k)) / burned_in_iters;
			if(PKnonemptyk > 0.001)
				cout << "P(KnonEmpty=" << k << ")=\t"
					<< stack.push << fixed << setw(6) << setprecision(4)
					<< PKnonemptyk
					<< stack.pop
					<< endl;
		}
		cout << "modalK=\t"         << max_element(K_freq.begin()        , K_freq.end()        ) - K_freq.begin() << endl;
		cout << "modalNonEmptyK=\t" << max_element(KnonEmpty_freq.begin(), KnonEmpty_freq.end()) - KnonEmpty_freq.begin() << endl;
	}
	if(args_info.labels_arg) {
		cout << endl << "Completed " << iteration << " iterations.  Using the last half (" << all_burned_in_z.size() << " iterations) for label switching." << endl << endl;
		assert(!all_burned_in_z.empty());
		cout << " label-switching ...";
		//cout.flush(); cout << " sorted ..."
		cout << endl;
		sort( all_burned_in_z.begin(), all_burned_in_z.end() );
		size_t max_K_to_consider_in_label_switching = highest_K_sampled;
		if(groundTruth) {
			const size_t max_K_in_ground_truth = 1 + *max_element(groundTruth->begin(), groundTruth->end());
			if(max_K_to_consider_in_label_switching < max_K_in_ground_truth) {
				max_K_to_consider_in_label_switching = max_K_in_ground_truth;
			}
		}
		label_switch(s._g, s, s._N, max_K_to_consider_in_label_switching, all_burned_in_z, groundTruth);
	}
	cout << "SBM complete. (after " << ELAPSED() << " seconds)" << endl;
}

typedef vector<long double > theta_t;
typedef vector<vector<long double > > pi_t;
typedef vector<int> z_t;

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
