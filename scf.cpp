#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "aaron_utils.hpp"
#include "scf.hpp"
#include "sbm.hpp"
using namespace std;
struct SCFreals { // the real valued parameters, the community densities
	long double pi_0; // background
	long double pi_1; // community 1
	long double pi_2; // community 2
	SCFreals() : pi_0(drand48()), pi_1(drand48()), pi_2(drand48()) {
	}
};

static void SCFiteration(gsl_rng * r, sbm :: State &s, const sbm:: ObjectiveFunction *obj, SCFreals &reals, AcceptanceRate *AR_metro);
static void newSCFreals(const gsl_rng * r, const sbm :: State &s, const sbm:: ObjectiveFunction *obj, SCFreals &reals);
static bool metroNode(const gsl_rng * r, sbm :: State &s, const sbm:: ObjectiveFunction *obj, const SCFreals &reals, AcceptanceRate *AR_metro);
static long double pmf_scf_x_given_z(const sbm :: State &s, const sbm:: ObjectiveFunction *obj, const SCFreals &reals);

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
	} else {
		if(commandLineK != -1)
			randomize(s, commandLineK);
	}

	assert(s._k == 2);

	s.shortSummary(obj, groundTruth); s.summarizeEdgeCounts(); s.blockDetail(obj); s.internalCheck();

	gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);

	SCFreals reals;
	AcceptanceRate AR_metro("metro");
	for(int iter=0; iter<10000; iter++) {
		cout << endl;
		PP(iter);
		PP(pmf_scf_x_given_z(s, obj, reals) + s.P_z_slow());
		s.shortSummary(obj, groundTruth); s.summarizeEdgeCounts(); s.blockDetail(obj); s.internalCheck();
		PP3(reals.pi_0, reals.pi_1, reals.pi_2);
		SCFiteration(r, s, obj, reals, &AR_metro);
		AR_metro.dump();
	}
}

static void SCFiteration(gsl_rng * r, sbm :: State &s, const sbm:: ObjectiveFunction *obj, SCFreals &reals, AcceptanceRate *AR_metro) {
	newSCFreals(r, s, obj, reals);
	bool accepted = metroNode(r, s, obj, reals, AR_metro);
	PP(accepted);
}

static void newSCFreals(const gsl_rng * r, const sbm :: State &s, const sbm:: ObjectiveFunction *obj, SCFreals &reals) {
	// PP3(reals.pi_0, reals.pi_1, reals.pi_2);
	const int bg_pairs = obj->numberOfPairsInBlock(0,1, &s.labelling);
	const int bg_edges = s._edgeCounts.read(0,1) + s._edgeCounts.read(1,0);
	const int c1_pairs = obj->numberOfPairsInBlock(0,0, &s.labelling);
	const int c1_edges = s._edgeCounts.read(0,0);
	const int c2_pairs = obj->numberOfPairsInBlock(1,1, &s.labelling);
	const int c2_edges = s._edgeCounts.read(1,1);
	// PP2(bg_edges, bg_pairs);
	// PP2(c1_edges, c1_pairs);
	// PP2(c2_edges, c2_pairs);
	// draw new values from the posteriors, *independently* of each other.
	// BUT then reject them all if the constraint isn't satisfied.

	while(1) {
		const double bBG = gsl_ran_beta(r, 1+bg_edges, 1+bg_pairs-bg_edges);
		// const double b1  = gsl_ran_beta(r, 1+c1_edges, 1+c1_pairs-c1_edges);
		// const double b2  = gsl_ran_beta(r, 1+c2_edges, 1+c2_pairs-c2_edges);
		const double b12  = gsl_ran_beta(r, 1+c1_edges+c2_edges, 1+c1_pairs+c2_pairs-c1_edges-c2_edges);
		assert(isfinite(bBG));
		// assert(isfinite(b1));
		// assert(isfinite(b2));
		assert(isfinite(b12));
		// PP2(bBG, b12);
		// if(bBG < b1 && bBG < b2)
		if(bBG < b12)
		{
			// PP3(bBG, b1, b2);
			reals.pi_0 = bBG;
			reals.pi_1 = b12;
			reals.pi_2 = b12;
			return;
		}
		else {
			// cout << "reject ";
			// PP3(bBG, b1, b2);
			// assert(1==2);
			continue;
		}
	}
}

static long double pmf_scf_x_given_z(const sbm :: State &s, const sbm:: ObjectiveFunction *obj, const SCFreals &reals) {
	const bool verbose = false;
	assert(s._k==2);
	const int bg_pairs = obj->numberOfPairsInBlock(0,1, &s.labelling);
	const int bg_edges = s._edgeCounts.read(0,1) + s._edgeCounts.read(1,0);
	const int c1_pairs = obj->numberOfPairsInBlock(0,0, &s.labelling);
	const int c1_edges = s._edgeCounts.read(0,0);
	const int c2_pairs = obj->numberOfPairsInBlock(1,1, &s.labelling);
	const int c2_edges = s._edgeCounts.read(1,1);
	if(verbose) PP2(bg_edges, bg_pairs);
	if(verbose) PP2(c1_edges, c1_pairs);
	if(verbose) PP2(c2_edges, c2_pairs);
	const long double x0_z = bg_edges * log2l(reals.pi_0) + (bg_pairs-bg_edges) * log2l(1.0L-reals.pi_0);
	const long double x1_z = c1_edges * log2l(reals.pi_1) + (c1_pairs-c1_edges) * log2l(1.0L-reals.pi_1);
	const long double x2_z = c2_edges * log2l(reals.pi_2) + (c2_pairs-c2_edges) * log2l(1.0L-reals.pi_2);
	if(verbose) PP(x0_z);
	if(verbose) PP(x1_z);
	if(verbose) PP(x2_z);
	return x0_z + x1_z + x2_z;
}

static bool metroNode(const gsl_rng * r, sbm :: State &s, const sbm:: ObjectiveFunction *obj, const SCFreals &reals, AcceptanceRate *AR_metro) {
	assert(s._k==2);
	const long double pre = pmf_scf_x_given_z(s, obj, reals) + s.P_z_slow();
	const int randomNode = drand48() * s._N;
	const int oldCluster = s.labelling.cluster_id.at(randomNode);
	const int newCluster = 1 - oldCluster;
	// PP(pre);
	// PP(randomNode);
	// PP2(oldCluster, newCluster);
	s.moveNodeAndInformOfEdges(randomNode, newCluster);
	const long double post = pmf_scf_x_given_z(s, obj, reals) + s.P_z_slow();
	// PP(post);

	// const long double u = drand48();
	if(log2l(drand48()) < post - pre) {
		cout << "Accept metroNode" << endl;
		AR_metro->notify(true);
		return true;
	} else {
		s.moveNodeAndInformOfEdges(randomNode, oldCluster);
		assert(pre == pmf_scf_x_given_z(s, obj, reals) + s.P_z_slow()); // TODO VERYCLOSE
		AR_metro->notify(false);
		return false;
	}
}
