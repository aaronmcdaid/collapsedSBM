using namespace std;
#include <algorithm>
#include <vector>

#include <getopt.h>
#include <unistd.h>
#include <libgen.h>
#include <float.h>
#include <gsl/gsl_sf.h>

#include "aaron_utils.hpp"
#include "shmGraphRaw.hpp"
#include "sbm_state.hpp"
#include "mmsb_state.hpp"

const char gitstatus[] = 
#include "comment.txt"
#include "gitstatus.txt"
;
#include "cmdline.h"



struct UsageMessage {
};

template<bool selfloops, bool directed, bool weighted>
void runSBM(const sbm::GraphType *g, const int commandLineK);
void runMMSB(const sbm::GraphType *g, const int commandLineK);

static
void dumpGraph(shmGraphRaw::ReadableShmGraphTemplate<shmGraphRaw::PlainMem> *g, const shmGraphRaw:: EdgeDetails< shmGraphRaw:: NoDetails > & edge_details) {
	PP(g->numNodes());
	PP(g->numRels());
	for(int rel=0; rel<g->numRels(); rel++) {
		std::pair<int,int> eps = g->EndPoints(rel);
		std::pair<const char*, const char*> epsNames = g->EndPointsAsStrings(rel);
		cout << rel
			<< '\t' << eps.first << '"' << epsNames.first << '"'
			<< '\t' << eps.second << '"' << epsNames.second << '"'
			<< endl;
	}
}
static
void dumpGraph(shmGraphRaw::ReadableShmGraphTemplate<shmGraphRaw::PlainMem> *g, const shmGraphRaw:: EdgeDetails< shmGraphRaw:: DirectedLDoubleWeights > & edge_details) {
	PP(g->numNodes());
	PP(g->numRels());
	for(int rel=0; rel<g->numRels(); rel++) {
		std::pair<int,int> eps = g->EndPoints(rel);
		std::pair<const char*, const char*> epsNames = g->EndPointsAsStrings(rel);
		cout << rel
			<< '\t' << eps.first << '"' << epsNames.first << '"'
			<< '\t' << eps.second << '"' << epsNames.second << '"'
			<< '\t' << edge_details.dw.at(rel).first << ',' << edge_details.dw.at(rel).second
			<< endl;
	}
}

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
	PP(args_info.mmsb_flag);
	PP(args_info.K_arg);
	PP(args_info.directed_flag);
	PP(args_info.weighted_flag);
	PP(args_info.selfloop_flag);

	if(!args_info.directed_flag && !args_info.weighted_flag) { // UNdir UNwei
		shmGraphRaw:: EdgeDetails< shmGraphRaw:: NoDetails > edge_details;
		auto_ptr<shmGraphRaw::ReadableShmGraphTemplate<shmGraphRaw::PlainMem> > g (shmGraphRaw::loadEdgeList<shmGraphRaw::PlainMem>(edgeListFileName, edge_details));
		dumpGraph(g.get(), edge_details);
		if(args_info.selfloop_flag)
			runSBM<true ,false,false>(g.get(), args_info.K_arg);
		else
			runSBM<false,false,false>(g.get(), args_info.K_arg);
	}
	if( args_info.directed_flag &&  args_info.weighted_flag) { // UNdir UNwei
		shmGraphRaw:: EdgeDetails< shmGraphRaw:: DirectedLDoubleWeights > edge_details;
		auto_ptr<shmGraphRaw::ReadableShmGraphTemplate<shmGraphRaw::PlainMem> > g (shmGraphRaw::loadEdgeList<shmGraphRaw::PlainMem>(edgeListFileName, edge_details));
		dumpGraph(g.get(), edge_details);
		if(args_info.selfloop_flag)
			runSBM<true ,true,true>(g.get(), args_info.K_arg);
		else
			runSBM<false,true,true>(g.get(), args_info.K_arg);
	}
	exit(0);
	// if(args_info.mmsb_flag)
		// runMMSB(g.get(), commandLineK);
	// else
		// runSBM(g.get(), commandLineK);
}

void randomize(sbm::State &s, const int K) { // randomize the partition and have K clusters in it
	assert(s._k <= K);
	assert(K >= 1);
	for(int i=0; i<1000 || s._k < K; i++) {
		int n;
		do {
			n = drand48() * s._N;
			//cout << endl << "Moving node: " << n << " move# " << i << endl;
			s.isolateNode(n);
		} while (s._k <= K);
		const int newClusterID = drand48() * (s._k-1); // -1 because we must move the node from the "temporary" community
		// s.shortSummary(); s.summarizeEdgeCounts();

		s.unIsolateTempNode(n, newClusterID);
		assert(s._k <= K);
		// s.shortSummary(); s.summarizeEdgeCounts();
		s.internalCheck();
	}
	assert(s._k == K);
	s.internalCheck();
}

bool fiftyfifty() {
	if(drand48() < 0.5)
		return true;
	else
		return false;
}
bool acceptTest(const long double delta) {
	if(log2(drand48()) < delta)
		return true;
	else
		return false;
}

long double MoneNode(sbm::State &s) {
	if(s._k == 1)
	       return 0.0L;	// can't move a node unless there exist other clusters
	assert(s._k > 1); // can't move a node unless there exist other clusters
	const long double pre = s.pmf();
	const int n = drand48() * s._N;
	const int oldClusterID = s.labelling.cluster_id.at(n);
	int newClusterID;
	do {
		newClusterID = drand48() * s._k;
	} while (newClusterID == oldClusterID);
	assert(newClusterID != oldClusterID);
	// PP(oldClusterID);
	// PP(newClusterID);
	s.moveNode(n, newClusterID);
	s.informNodeMove(n, oldClusterID, newClusterID);
	const long double post = s.pmf();
	const long double delta = post - pre;
	// PP(pre);
	// PP(post);
	// PP(delta);
	if(acceptTest(delta)) {
		// cout << " + ";
		return delta;
	} else {
		// cout << "   ";
		s.moveNode(n, oldClusterID);
		s.informNodeMove(n, newClusterID, oldClusterID);
		// assert(VERYCLOSE(s.pmf(), pre)); // make sure it has undone it properly
		return 0.0L;
	}
}
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
		const long double LeftOverRight = exp2(log2LeftOverRight);
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
static OneChoice M3_oneNode(sbm::State &s, const int n, const int candCluster) {
	// const long double pre = s.pmf(); // TODO remove this, and the corresponding check at the end
	// given an isolated node, and a candidate cluster to join, what's the delta in the fitness?
	
	const int isolatedClusterID = s.labelling.cluster_id.at(n);
	const sbm:: Cluster * isoCl = s.labelling.clusters.at(isolatedClusterID);
	assert(isoCl->order() == 1);
	const sbm:: Cluster * candCl = s.labelling.clusters.at(candCluster);
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
void M3(sbm::State &s) {
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
	const sbm:: Cluster * CL1 = s.labelling.clusters.at(cl1);
	const sbm:: Cluster * CL2 = s.labelling.clusters.at(cl2);
	vector<int> allNodes;
	boost::unordered_map<int, int> statusQuoClustering;
	allNodes.insert(allNodes.end(), CL1->members.begin(), CL1->members.end());
	allNodes.insert(allNodes.end(), CL2->members.begin(), CL2->members.end());
	// PP(CL1->members.size());
	// PP(CL2->members.size());
	// PP(allNodes.size());
	assert(CL1->members.size() + CL2->members.size() == allNodes.size());
	forEach(int x, amd::mk_range(CL1->members)) { statusQuoClustering[x] = cl1; }
	forEach(int y, amd::mk_range(CL2->members)) { statusQuoClustering[y] = cl2; }
	// forEach(int x, amd::mk_range(CL1->members)) { PPt(x); } cout << endl;
	// forEach(int y, amd::mk_range(CL2->members)) { PPt(y); } cout << endl;
	random_shuffle(allNodes.begin(), allNodes.end());
	// forEach(int z2, amd::mk_range(allNodes    )) { PPt(z2); } cout << endl;
	long double deltaSumOfTheStatusQuo = 0.0L;
	long double log2ProductOfProposalProbabilitiesForStatusQuo = 0.0;
	for(vector<int>::const_reverse_iterator remover = allNodes.rbegin(); remover != allNodes.rend(); ++remover) {
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
		const sbm:: Cluster * old_cluster = s.labelling.clusters.at(old_clusterID);
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
		for(vector<int>::const_iterator adder = allNodes.begin(); adder != allNodes.end(); ++adder) {
			// cout << endl << "  random proposal for M3" << endl << endl;
			// const long double preM3OneRandom = s.pmf();
			const int node_to_Add = *adder;
			// PP(node_to_Add);
			const int clID = s.labelling.cluster_id.at(node_to_Add);
			assert(clID + 1 == s._k);
			const sbm:: Cluster * clIsolated = s.labelling.clusters.at(clID);
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

	// Now we either accept it with probability exp2(acceptanceLog2), or reject it

	if(log2(drand48()) < acceptanceLog2 /*|| IsRandomProposalIdenticalToStatusQuo*/) {
		// accepting.
		// if(verbose)
			cout << "M3 ACCEPT: " << acceptanceLog2 << endl;
	} else {
		// reject. let's put them all back
		for(vector<int>::const_iterator reAdder = allNodes.begin(); reAdder != allNodes.end(); ++reAdder) {
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
void MetropolisOnK(sbm::State &s) {
	// const long double prePMF = s.pmf();
	const long double prePMF12 = s.P_z_K();
	const int preK = s._k;
	if(fiftyfifty()) { // propose increase in K
		s.appendEmptyCluster();
		const long double postPMF12 = s.P_z_K();
		// assert(VERYCLOSE(s.pmf(), prePMF - prePMF12 + postPMF12));
		// const long double postPMF = prePMF - prePMF12 + postPMF12; 
		// assert(VERYCLOSE(s.pmf(), postPMF));
		// assert(postPMF < prePMF);
		assert(postPMF12 < prePMF12);
		// assert(VERYCLOSE(postPMF - prePMF, postPMF12 - prePMF12));
		if(acceptTest(postPMF12 - prePMF12)) {
			// cout << "k: acc inc" << endl;
			assert(s._k>preK);
		} else {
			// cout << "k: rej inc" << endl;
			s.deleteClusterFromTheEnd();
			// assert(s.pmf()==prePMF);
			assert(s.P_z_K()==prePMF12);
			assert(s._k==preK);
		}
	} else { // propose decrease
		if(s._k >= 1 && s.labelling.clusters.back()->order()==0) {
			s.deleteClusterFromTheEnd();
			const long double postPMF12 = s.P_z_K();
			// assert(VERYCLOSE(s.pmf(), prePMF - prePMF12 + postPMF12));
			// const long double postPMF = prePMF - prePMF12 + postPMF12; 
			// assert(VERYCLOSE(s.pmf(), postPMF));
			// assert(postPMF > prePMF);
			if(acceptTest(postPMF12 - prePMF12)) {
				// cout << "k: acc dec" << endl;
				assert(s._k<preK);
			} else {
				assert(1==2); // it'll always like decreases, except maybe if the prior on K is an increasing function.
				// cout << "k: rej dec" << endl;
				s.appendEmptyCluster();
				// assert(s.pmf()==postPMF);
				assert(s.P_z_K()==prePMF12);
				assert(s._k==preK);
			}
		} else {// otherwise, not possible to remove it
			// cout << "k: rej-dec" << endl;
			assert(s._k==preK);
		}
	}
}

template<bool selfloops, bool directed, bool weighted>
void runSBM(const sbm::GraphType *g, const int commandLineK) {
	sbm::State s(g);

	s.shortSummary(); s.summarizeEdgeCounts(); s.blockDetail();
	s.internalCheck();

	/*
	s.isolateNode(0);
	s.isolateNode(1); // to bring us up to three clusters

	s.shortSummary(); s.summarizeEdgeCounts(); s.blockDetail();
	s.internalCheck();
	PP(s.pmf());

	*/
	if(commandLineK != -1)
		randomize(s, commandLineK);
	// s.shortSummary(); s.summarizeEdgeCounts(); s.blockDetail();

	PP(s.pmf());

	for(int i=1; i<=400000000; i++) {
		if(commandLineK == -1)
			MetropolisOnK(s);
		// PP(i);
		MoneNode(s);
		// if(i%50 == 0)
			// M3(s);
	
		// PP(s.pmf());
		// cout << endl;
		s.internalCheck();
		if(i%100 == 0) {
			cout << endl;
			PP(i);
			s.shortSummary(); s.summarizeEdgeCounts(); s.blockDetail();
			cout << " end of check at i==" << i << endl;
		}
	}
	s.shortSummary(); s.summarizeEdgeCounts(); s.blockDetail();
	s.internalCheck();
}

void runMMSB(const sbm::GraphType *g, const int commandLineK) {
	assert(commandLineK > 1);
	PP2(g->numNodes(), g->numRels());
	sbm:: MMSBstate s(g);
	PP(s._N);
	// s.P_zs_given_K();
	// cout << endl;

	while(s._k < commandLineK)
		s.appendEmptyCluster();
	PP(s._k);

	s.P_zs_given_K();
	cout << endl;

	if(0) {
		/* 18 3 4 5 7
		 * 10 12 17 2 8
		 * 1 11 13 14 6
		 * 15 16 19 20 9
		 */
		const int comm1[] = {18,3,4,5,7,-1};
		const int comm2[] = {10,12,17,2,8,-1};
		const int comm3[] = {1,11,13,14,6,-1};
		const int comm4[] = {15,16,19,20,9,-1};
		const int * comms[] = {comm1,comm2,comm3,comm4,NULL};
		for(int commid=1; commid<4;commid++) {
			const int * comm = comms[commid];
			for(const int *commW = comm; *commW != -1; commW++) {
				const int W = *commW;
				PP(W);
				const int w = g->StringToNodeId(printfstring("%d", W).c_str());
				for(int v = 0; v<s._N;v++) {
					if(w!=v)
						s.performMoveAndUpdateEdges(w,v,commid);
				}
			}
		}
	}
	// with SIMPLER_Z==0.5, pmf(GT) == 2271.4
	// with SIMPLER_Z==2.0, pmf(GT) ==  536.319
	if(1) { // random initialization
		for(int w=0; w<s._N; w++) {
			for(int v=0; v<s._N; v++) {
				if(w!=v) {
					const int randomID = drand48()*s._k;
					if(randomID>0)
						s.performMoveAndUpdateEdges(w,v,randomID);
				}
			}
		}
	}

	s.P_zs_given_K();
	cout << endl;

	long double deltas = 0.0L;

	for(int i=0; i<=100000000; i++) {
		const int w=g->numNodes() * drand48();
		const int v=g->numNodes() * drand48();
		if(w!=v) {
			const long double pre = s.pmf_slow();
			const long double delta = s.MetropolisMoveOnePair(w,v,4*drand48());
			const long double post = s.pmf_slow();
			PP3(pre,delta,post);
			PP(pre+delta - post);
			assert(VERYCLOSE(pre+delta , post));
			deltas += delta;
		}
		if(i%1==0) {
			PP(i);
			s.P_zs_given_K();
			PP(deltas);
			cout << endl;
		}
	}

	s.P_zs_given_K();
	PP(deltas);
	cout << endl;
}
