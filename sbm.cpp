using namespace std;
#include <algorithm>
#include <vector>

#include <getopt.h>
#include <unistd.h>
#include <libgen.h>
#include <gsl/gsl_sf.h>

#include "aaron_utils.hpp"
#include "shmGraphRaw.hpp"
#include "sbm_state.hpp"

const char gitstatus[] = 
#include "comment.txt"
#include "gitstatus.txt"
;



struct UsageMessage {
};

void runSBM(const sbm::GraphType *g);

int main(int argc, char **argv) {
	PP(gitstatus);
	for (int i=0; i<argc; i++) {
		PP(argv[i]);
	}
	{ int c, option_index; while (1)
		{
      static const struct option long_options[] = {
        // {"seed", required_argument,       0, 22},
        {0, 0, 0, 0}
      };
      /* getopt_long stores the option index here. */
      c = getopt_long (argc, argv, "t:T:k:", long_options, &option_index);
      if (c == -1) break; /* Detect the end of the options. */
     
      switch (c) {
        case '?': /* getopt_long already printed an error message. */ break;
        default: abort (); break;
        case 0:
          /* If this option set a flag, do nothing else now. */
          if (long_options[option_index].flag != 0)
            break;
          printf ("option %s", long_options[option_index].name);
          if (optarg) printf (" with arg %s", optarg);
          printf ("\n");
          break;
        // case 'T':
					// option_inclusiveThreshold = true;
      }
    }
	}

	unless (argc - optind == 2) {
		throw UsageMessage();
		// cout << "Usage: edge_list directory_for_output" << endl;
		// exit(1);
	}
	const char * edgeListFileName = argv[optind];
	const char * directoryForOutput = argv[optind+1];
	PP(edgeListFileName);
	PP(directoryForOutput);

	auto_ptr<shmGraphRaw::ReadableShmGraphTemplate<shmGraphRaw::PlainMem> > g (shmGraphRaw::loadEdgeList<shmGraphRaw::PlainMem>(edgeListFileName));
	runSBM(g.get());
}

void randomize(sbm::State &s, const int K) { // randomize the partition and have K clusters in it
	assert(s._k <= K);
	for(int i=0; i<10000 || s._k < K; i++) {
		int n;
		do {
			n = drand48() * s._N;
			//cout << endl << "Moving node: " << n << " move# " << i << endl;
			s.isolateNode(n);
		} while (s._k <= K);
		const int newClusterID = drand48() * (s._k-1); // -1 because we must move the node from the "temporary" community
		// s.shortSummary(); s.summarizeEdgeCounts();
		s.internalCheck();

		s.unIsolateTempNode(n, newClusterID);
		assert(s._k <= K);
		// s.shortSummary(); s.summarizeEdgeCounts();
		s.internalCheck();
	}
	assert(s._k == K);
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
	const long double pre = s.pmf();
	const int n = drand48() * s._N;
	const int oldClusterID = s.cluster_id.at(n);
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
		cout << " + ";
		return delta;
	} else {
		cout << "   ";
		s.moveNode(n, oldClusterID);
		s.informNodeMove(n, newClusterID, oldClusterID);
		assert(VERYCLOSE(s.pmf(), pre)); // make sure it has undone it properly
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
		PP2(left_deltaSum, right_deltaSum);
		PP2(Pleft, Pright);
		assert(Pleft >= 0.0L);
		assert(Pleft <= 1.0L);
		assert(Pright >= 0.0L);
		assert(Pright <= 1.0L);
		assert(isfinite(Pleft));
		assert(isfinite(Pright));
	}
};
static OneChoice M3_oneNode(sbm::State &s, const int n, const int candCluster) {
	const long double pre = s.pmf(); // TODO remove this, and the corresponding check at the end
	// given an isolated node, and a candidate cluster to join, what's the delta in the fitness?
	
	const int isolatedClusterID = s.cluster_id.at(n);
	const sbm::State:: Cluster * isoCl = s.clusters.at(isolatedClusterID);
	assert(isoCl->order() == 1);
	const sbm::State:: Cluster * candCl = s.clusters.at(candCluster);
	assert(isolatedClusterID != candCluster);
	const int old_order = candCl->order();
	const int new_order = old_order+1;
	const int old_NonEmpty = s.NonEmptyClusters;
	const int new_NonEmpty = old_order==0 ? s.NonEmptyClusters : (s.NonEmptyClusters-1);
	const long double preSumOfLog2l = s.SumOfLog2LOrders;
	const long double postReMergeSumOfLog2l = preSumOfLog2l + log2l(new_order) - ( old_order>1 ? log2l(old_order) : 0.0L);
	PP(postReMergeSumOfLog2l);


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

	PP(delta2);
	PP(delta3);
	PP(delta4);
	
	const long double postTempMove = s.pmf();
	assert(VERYCLOSE(postTempMove , pre + delta2 + delta3 + delta4));

	s.moveNodeAndInformOfEdges(n, isolatedClusterID); // move the node back again

	const long double undone = s.pmf();
	PP(undone);
	assert(VERYCLOSE(pre, undone)); // to ensure that we've undone the change

	return OneChoice(delta2 , delta3 , delta4);
}
void M3(sbm::State &s) {
	cout << endl << "     ========== M3 =========" << endl;
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
	const sbm::State::Cluster * CL1 = s.clusters.at(cl1);
	const sbm::State::Cluster * CL2 = s.clusters.at(cl2);
	vector<int> allNodes;
	boost::unordered_map<int, int> statusQuoClustering;
	allNodes.insert(allNodes.end(), CL1->members.begin(), CL1->members.end());
	allNodes.insert(allNodes.end(), CL2->members.begin(), CL2->members.end());
	PP(CL1->members.size());
	PP(CL2->members.size());
	PP(allNodes.size());
	assert(CL1->members.size() + CL2->members.size() == allNodes.size());
	forEach(int x, amd::mk_range(CL1->members)) { statusQuoClustering[x] = cl1; }
	forEach(int y, amd::mk_range(CL2->members)) { statusQuoClustering[y] = cl2; }
	forEach(int x, amd::mk_range(CL1->members)) { PPt(x); } cout << endl;
	forEach(int y, amd::mk_range(CL2->members)) { PPt(y); } cout << endl;
	random_shuffle(allNodes.begin(), allNodes.end());
	forEach(int z2, amd::mk_range(allNodes    )) { PPt(z2); } cout << endl;
	long double deltaSumOfTheStatusQuo = 0.0L;
	long double log2ProductOfProposalProbabilitiesForStatusQuo = 0.0;
	for(vector<int>::const_reverse_iterator remover = allNodes.rbegin(); remover != allNodes.rend(); ++remover) {
		const int node_to_remove = *remover;
		PP(node_to_remove);
		const long double pre = s.pmf();
		const long double pre1 = s.P_z_K();
		const long double pre2 = s.P_z_orders();
		const long double pre3 = s.P_edges_given_z_baseline();
		const long double pre4 = s.P_edges_given_z_correction();
		assert(pre == pre1+pre2+pre3+pre4);
		const long double preSumOfLog2l = s.SumOfLog2LOrders;
		PP(preSumOfLog2l);
		const long double preNonEmpty = s.NonEmptyClusters;
		assert(pre == pre1 + pre2 + pre3 + pre4);

		const int old_clusterID = s.cluster_id.at(node_to_remove);
		const sbm::State:: Cluster * old_cluster = s.clusters.at(old_clusterID);
		const int old_order = old_cluster->order();
		assert(old_order>=1);
		long double delta2 = -LOG2FACT(old_order);
		long double delta4 = - s.P_edges_given_z_correction_JustOneCluster(old_clusterID); // this must be remembered before the move

		const int tempClusterID = s.isolateNode(node_to_remove);

		const long double postSumOfLog2l = s.SumOfLog2LOrders;
		const long double postNonEmpty = s.NonEmptyClusters;
		const int new_order = old_order-1;
		if(new_order>1)
			delta2 += LOG2FACT(new_order);
		assert(new_order >= 0);
		const long double post2 = s.P_z_orders();
		assert(VERYCLOSE(post2, pre2 + delta2));

		const long double post3 = s.P_edges_given_z_baseline();
		PP(post3 - pre3);
		long double delta3 = 0.0;
		PP2(old_order, new_order);
		delta3 += postSumOfLog2l * (postNonEmpty-1);
		delta3 -= preSumOfLog2l  * (preNonEmpty-1);
		if(new_order >= 2)
			delta3 += log2l((new_order * new_order - new_order)/2);
		if(old_order >= 2)
			delta3 -= log2l((old_order * old_order - old_order)/2);
		delta3 = -delta3;
		PP(delta3);
		assert(VERYCLOSE(delta3 , post3 - pre3));

		const long double post4 = s.P_edges_given_z_correction();
		delta4 += s.P_edges_given_z_correction_JustOneCluster(old_clusterID);
		delta4 += s.P_edges_given_z_correction_JustOneCluster(tempClusterID);
		// Any edges between tempClusterID and old_clusterID (i.e. old internal edges) will have been double counted.
		const int doubleCounted_edges = s._edgeCounts.get(tempClusterID, old_clusterID);
		delta4 += M_LOG2E * gsl_sf_lnchoose(new_order/*doubleCounted_pairs*/, doubleCounted_edges);

		PP2(post4,pre4);
		PP(post4 - pre4);
		PP(delta4);
		assert(post4 - pre4 == delta4);

		const long double post1 = s.P_z_K();
		long double delta1 = post1 - pre1;

		PP(delta1);
		PP(delta2);
		PP(delta3);
		PP(delta4);

		cout << " ==  M3_oneNode ==" << endl;
		long double left  = M3_oneNode(s, node_to_remove, cl1).deltaSum();
		long double right = M3_oneNode(s, node_to_remove, cl2).deltaSum();
		TwoChoices two_choices(M3_oneNode(s, node_to_remove, cl1),M3_oneNode(s, node_to_remove, cl2));
		assert(cl1 == old_clusterID || cl2 == old_clusterID);
		const long double statusQuo = cl1 == old_clusterID ? left : right;
		PP2(left, right);
		PP2(two_choices.Pleft, two_choices.Pright);
		cout << " == ~M3_oneNode ==" << endl;
		assert(VERYCLOSE(-statusQuo , delta2 + delta3 + delta4));
		const long double post = s.pmf();
		assert(VERYCLOSE(post, post1+post2+post3+post4));
		assert(VERYCLOSE(post, pre + delta1 + delta2 + delta3 + delta4));

		const long double prpsl = cl1 == old_clusterID ? two_choices.Pleft : two_choices.Pright;
		PP(prpsl);
		log2ProductOfProposalProbabilitiesForStatusQuo += log2(prpsl);

		deltaSumOfTheStatusQuo += cl1 == old_clusterID ? two_choices.left_deltaSum : two_choices.right_deltaSum;
	}

	const long double midM3 = s.pmf();
	const long double midM3_1 = s.P_z_K();
	// const long double midM3_2 = s.P_z_orders();
	// const long double midM3_3 = s.P_edges_given_z_baseline();
	// const long double midM3_4 = s.P_edges_given_z_correction();
	PP2(-deltaSumOfTheStatusQuo, midM3 - midM3_1 - preM3 + preM3_1);
	assert(VERYCLOSE(-deltaSumOfTheStatusQuo , midM3 - midM3_1 - preM3 + preM3_1));

	long double log2ProductOfProposalProbabilitiesForNewProposal = 0.0L;
	{ // random proposal
		cout << endl << "  random proposals for M3" << endl << endl;
		for(vector<int>::const_iterator adder = allNodes.begin(); adder != allNodes.end(); ++adder) {
			cout << endl << "  random proposal for M3" << endl << endl;
			const long double preM3OneRandom = s.pmf();
			const int node_to_Add = *adder;
			PP(node_to_Add);
			const int clID = s.cluster_id.at(node_to_Add);
			assert(clID + 1 == s._k);
			const sbm::State:: Cluster * clIsolated = s.clusters.at(clID);
			assert(clIsolated->order()==1);
			assert(clIsolated->members.front()==node_to_Add);
			// which of the two to add to?
			long double left  = M3_oneNode(s, node_to_Add, cl1).deltaSum();
			long double right = M3_oneNode(s, node_to_Add, cl2).deltaSum();
			PP2(left,right);
			assert(cl1 != clID && cl2 != clID);
			assert(cl1 != cl2);
			TwoChoices two_choices(M3_oneNode(s, node_to_Add, cl1),M3_oneNode(s, node_to_Add, cl2));
			PP2(left , two_choices.left.deltaSum());
			assert(VERYCLOSE(left , two_choices.left.deltaSum()));
			assert(VERYCLOSE(right, two_choices.right.deltaSum()));
			PP2(two_choices.left.deltaSum(),two_choices.right.deltaSum());
			PP2(two_choices.Pleft,two_choices.Pright);
			long double prpsl;
			if(drand48() < two_choices.Pright) {
				cout << " go right" << endl;
				s.moveNodeAndInformOfEdges(node_to_Add, cl2);
				assert(VERYCLOSE(s.pmf() , preM3OneRandom + two_choices.right.deltaSum()));
				// the above assert is (correctly) using an 'out-of-date' value of _k. Hence we don't delete this temporary (now empty) cluster until the next line
				s.deleteClusterFromTheEnd();
				prpsl = two_choices.Pright;
			} else {
				cout << " go left" << endl;
				s.moveNodeAndInformOfEdges(node_to_Add, cl1);
				assert(VERYCLOSE(s.pmf() , preM3OneRandom + two_choices.left.deltaSum()));
				// the above assert is (correctly) using an 'out-of-date' value of _k. Hence we don't delete this temporary (now empty) cluster until the next line
				s.deleteClusterFromTheEnd();
				prpsl = two_choices.Pleft;
			}

			log2ProductOfProposalProbabilitiesForNewProposal += log2(prpsl);

			cout << endl << " ~random proposal for M3" << endl << endl;
		}
	}

	cout << "Accept or Reject? Not yet implemented." << endl;
	PP2(log2ProductOfProposalProbabilitiesForStatusQuo, log2ProductOfProposalProbabilitiesForNewProposal);

	// let's put them all back
	assert(preM3_k == s._k);
	for(vector<int>::const_iterator reAdder = allNodes.begin(); reAdder != allNodes.end(); ++reAdder) {
		// these should all be isolated nodes, and at the end of the list of clusters
		const int node_to_reAdd = *reAdder;
		PP(node_to_reAdd);
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
	PP2(preM3, postM3);
	cout << endl << "     ========= ~M3 =========" << endl;
	assert(preM3_1 == postM3_1);
	assert(preM3_2 == postM3_2);
	assert(VERYCLOSE(preM3_3, postM3_3));
	assert(preM3_4 == postM3_4);
	assert(VERYCLOSE(preM3 , postM3));
}
void MetropolisOnK(sbm::State &s) {
	const long double pre = s.pmf();
	const int preK = s._k;
	if(fiftyfifty()) { // propose increase in K
		s.appendEmptyCluster();
		const long double post = s.pmf();
		assert(post < pre);
		if(acceptTest(post - pre)) {
			// cout << "k: acc inc" << endl;
			assert(s._k>preK);
		} else {
			// cout << "k: rej inc" << endl;
			s.deleteClusterFromTheEnd();
			assert(s.pmf()==pre);
			assert(s._k==preK);
		}
	} else { // propose decrease
		if(s._k >= 1 && s.clusters.back()->order()==0) {
			s.deleteClusterFromTheEnd();
			const long double post = s.pmf();
			assert(post > pre);
			if(acceptTest(post - pre)) {
				// cout << "k: acc dec" << endl;
				assert(s._k<preK);
			} else {
				assert(1==2); // it'll always like decreases, except maybe if the prior on K is an increasing function.
				// cout << "k: rej dec" << endl;
				s.appendEmptyCluster();
				assert(s.pmf()==post);
				assert(s._k==preK);
			}
		} else {// otherwise, not possible to remove it
			// cout << "k: rej-dec" << endl;
			assert(s._k==preK);
		}
	}
}

void runSBM(const sbm::GraphType *g) {
	sbm::State s(g);

	s.shortSummary(); s.summarizeEdgeCounts(); s.blockDetail();
	s.internalCheck();

	s.isolateNode(0);
	s.isolateNode(1); // to bring us up to three clusters

	s.shortSummary(); s.summarizeEdgeCounts(); s.blockDetail();
	s.internalCheck();
	PP(s.pmf());

	randomize(s, 3);
	s.shortSummary(); s.summarizeEdgeCounts(); s.blockDetail();
	PP(s.pmf());

	for(int i=1; i<=20000000; i++) {
		// PP(i);
		MoneNode(s);
		// PP(s.pmf());
		// cout << endl;
		s.internalCheck();
		if(i%100 == 0) {
			cout << endl;
			PP(i);
			s.shortSummary(); s.summarizeEdgeCounts(); s.blockDetail();
			M3(s);
		}
		MetropolisOnK(s);
	}
	s.shortSummary(); s.summarizeEdgeCounts(); s.blockDetail();
	s.internalCheck();
}
