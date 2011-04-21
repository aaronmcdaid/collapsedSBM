#include <iostream>
#include "aaron_utils.hpp"
#include "scf.hpp"
using namespace std;
struct SCFreals { // the real valued parameters, the community densities
	long double pi_0; // background
	long double pi_1; // community 1
	long double pi_2; // community 2
	SCFreals() : pi_0(drand48()), pi_1(drand48()), pi_2(drand48()) {
	}
};

static void SCFiteration(sbm :: State &s, SCFreals &reals) {
	PP3(reals.pi_0, reals.pi_1, reals.pi_2);
}

void runSCF(const sbm::GraphType *g, const int commandLineK, const shmGraphRaw:: EdgeDetailsInterface * const edge_details, const bool initializeToGT, const vector<int> * const groundTruth, const int iterations) {
	cout << endl << "Stochastic Community Finding" << endl << endl;
	assert(commandLineK == 2);

	sbm::State s(g, edge_details);

	sbm:: ObjectiveFunction *obj = new sbm:: ObjectiveFunction_Bernoulli(false, false, false);
	s.shortSummary(obj, groundTruth); s.summarizeEdgeCounts(); s.blockDetail(obj);
	s.internalCheck();

	if(groundTruth && initializeToGT) {
		cout << "Loading GT" << endl;
		for(int v = 0; v < g->numNodes(); v++ ) {
			const int z_i = groundTruth->at(v);
			while(z_i >= s._k) {
				s.appendEmptyCluster();
			}
			s.moveNodeAndInformOfEdges2(v, z_i);
		}
	}

	s.shortSummary(obj, groundTruth); s.summarizeEdgeCounts(); s.blockDetail(obj); s.internalCheck();

	SCFreals reals;
	for(int iters=0; iters<1; iters++) {
		SCFiteration(s, reals);
	}
}
