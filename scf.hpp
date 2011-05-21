#include "sbm_state.hpp"

void runSCF(const sbm :: GraphType *g, const int commandLineK, const graph :: weights :: EdgeDetailsInterface * const edge_details, const bool initializeToGT, const std :: vector<int> * const groundTruth, const int iterations, const gsl_rng *r);
