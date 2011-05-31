#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include "macros.hpp"
#include "scf.hpp"
#include "sbm.hpp"
using namespace std;
struct SCFreals { // the real valued parameters, the community densities
	long double pi_0; // background
	long double pi_1; // community 1
	long double pi_2; // community 2
	SCFreals() : pi_0(drand48()), pi_1(drand48()), pi_2(pi_1) {
	}
};

static void SCFiteration(const gsl_rng * r, sbm :: State &s, const sbm :: ObjectiveFunction *obj, SCFreals &reals, AcceptanceRate *AR_metro);
static void newSCFreals(const gsl_rng * r, const sbm :: State &s, const sbm :: ObjectiveFunction *obj, SCFreals &reals);
static bool metroNode(const gsl_rng * r, sbm :: State &s, const sbm :: ObjectiveFunction *obj, const SCFreals &reals, AcceptanceRate *AR_metro);
static long double pmf_scf_x_given_z(const sbm :: State &s, const sbm :: ObjectiveFunction *obj, const SCFreals &reals);
static long double my_gsl_cdf_beta_Pinv (const long double P, const long double a, const long double b);
static long double my_gsl_cdf_beta_Qinv (const long double P, const long double a, const long double b);

void runSCF(const sbm :: GraphType *g, const int commandLineK, const bool initializeToGT, const std :: vector<int> * const groundTruth, const int iterations, const gsl_rng *r) {
	cout << endl << "Stochastic Community Finding" << endl << endl;
	PP2(g->numNodes(), g->numRels());
	assert(commandLineK == 2);

	sbm :: State s(g, false); // why true here?

	sbm :: ObjectiveFunction *obj = new sbm :: ObjectiveFunction_Bernoulli(false, false, false);
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

	SCFreals reals;
	AcceptanceRate AR_metro("metro");
	for(int iter=1; iter<=iterations; iter++) {
		SCFiteration(r, s, obj, reals, &AR_metro);
		if(iter % 10 == 0) {
			cout << endl;
			PP(iter);
			s.shortSummary(obj, groundTruth); s.summarizeEdgeCounts(); s.blockDetail(obj); s.internalCheck();
			PP3(reals.pi_0, reals.pi_1, reals.pi_2);
			PP(pmf_scf_x_given_z(s, obj, reals) + s.P_z_slow());
			AR_metro.dump();
		}
	}
}

static void SCFiteration(const gsl_rng * r, sbm :: State &s, const sbm :: ObjectiveFunction *obj, SCFreals &reals, AcceptanceRate *AR_metro) {
	newSCFreals(r, s, obj, reals);
	metroNode(r, s, obj, reals, AR_metro);
}

template <typename T>
static inline int32_t to_integer(T d) {
	assert(isfinite(d));
	assert(d == floor(d));
	return floor(d);
}

static void newSCFreals(const gsl_rng * r, const sbm :: State &s, const sbm :: ObjectiveFunction *obj, SCFreals &reals) {
	// PP3(reals.pi_0, reals.pi_1, reals.pi_2);
	const int64_t bg_pairs = obj->numberOfPairsInBlock(0,1, &s.labelling);
	const int bg_edges = to_integer(s._edgeCounts.read(0,1) + s._edgeCounts.read(1,0));
	const int64_t c1_pairs = obj->numberOfPairsInBlock(0,0, &s.labelling);
	const int c1_edges = to_integer(s._edgeCounts.read(0,0));
	const int64_t c2_pairs = obj->numberOfPairsInBlock(1,1, &s.labelling);
	const int c2_edges = to_integer(s._edgeCounts.read(1,1));
	// PP2(bg_edges, bg_pairs);
	// PP2(c1_edges, c1_pairs);
	// PP2(c2_edges, c2_pairs);
	// draw new values from the posteriors, *independently* of each other.
	// BUT then reject them all if the constraint isn't satisfied.

	if(1) { // Gibbs on reals.pi_0
		// we want a beta from gsl_ran_beta(r, 1+bg_edges, 1+bg_pairs-bg_edges);
		// BUT restricted to 0 < bBG < reals.pi_1
		// PP(__LINE__);
		// PP(reals.pi_0);
		assert(reals.pi_1 == reals.pi_2);
		const double unif_lbound = gsl_cdf_beta_P(0         , 1+bg_edges, 1+bg_pairs-bg_edges);
		const double unif_ubound = gsl_cdf_beta_P(reals.pi_1, 1+bg_edges, 1+bg_pairs-bg_edges);
		// PP2(unif_lbound, unif_ubound);
		assert(unif_lbound == 0.0);
		const double unif_cdf = gsl_ran_flat(r, unif_lbound, unif_ubound);
		// PP(unif_cdf);
		// PP2(1+bg_edges, 1+bg_pairs-bg_edges);
		const long double newbBG_ = my_gsl_cdf_beta_Pinv(unif_cdf, 1+bg_edges, 1+bg_pairs-bg_edges);
		// PP(newbBG_);
		// const long double newbBG = gsl_cdf_beta_Pinv(unif_cdf, 1+bg_edges, 1+bg_pairs-bg_edges);
		// PP(newbBG);
		// assert(newbBG_ == newbBG);
		reals.pi_0 = newbBG_;
	}
	if(1) { // Gibbs on reals.pi_1
		// we want a beta from gsl_ran_beta(r, 1+c1_edges+c2_edges, 1+c1_pairs+c2_pairs-c1_edges-c2_edges);
		// BUT restricted to        reals.pi_0  < newpi_1 < 1
		// PP(__LINE__);
		// PP(reals.pi_1);
		assert(reals.pi_1 == reals.pi_2);
		const double unif_lbound = gsl_cdf_beta_P(reals.pi_0, 1+c1_edges+c2_edges, 1+c1_pairs+c2_pairs-c1_edges-c2_edges);
		const double unif_ubound = gsl_cdf_beta_P(1         , 1+c1_edges+c2_edges, 1+c1_pairs+c2_pairs-c1_edges-c2_edges);
		// PP2(unif_lbound, unif_ubound);
		assert(unif_ubound == 1.0);
		const double unif_cdf = gsl_ran_flat(r, unif_lbound, unif_ubound);
		// PP(unif_cdf);
		const long double new_pi_1 = my_gsl_cdf_beta_Pinv(unif_cdf, 1+c1_edges+c2_edges, 1+c1_pairs+c2_pairs-c1_edges-c2_edges);
		// PP(new_pi_1);
		reals.pi_2 = reals.pi_1 = new_pi_1;
	}
}

static long double pmf_scf_x_given_z(const sbm :: State &s, const sbm :: ObjectiveFunction *obj, const SCFreals &reals) {
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

static bool metroNode(const gsl_rng *  , sbm :: State &s, const sbm :: ObjectiveFunction *obj, const SCFreals &reals, AcceptanceRate *AR_metro) {
	assert(s._k==2);
	const long double pre = pmf_scf_x_given_z(s, obj, reals) + s.P_z();
	const int randomNode = drand48() * s._N; // should use gsl maybe, to be consistent in the source of randomness?
	const int oldCluster = s.labelling.cluster_id.at(randomNode);
	const int newCluster = 1 - oldCluster;
	// PP(pre);
	// PP(randomNode);
	// PP2(oldCluster, newCluster);
	s.moveNodeAndInformOfEdges(randomNode, newCluster);
	const long double post = pmf_scf_x_given_z(s, obj, reals) + s.P_z();
	// PP(post);

	// const long double u = drand48();
	if(log2l(drand48()) < post - pre) {
		AR_metro->notify(true);
		return true;
	} else {
		s.moveNodeAndInformOfEdges(randomNode, oldCluster);
		// assert(pre == pmf_scf_x_given_z(s, obj, reals) + s.P_z()); // TODO VERYCLOSE
		AR_metro->notify(false);
		return false;
	}
}

static long double 
bisect (long double x, long double P, long double a, long double b, long double xtol, long double Ptol)
{
  long double x0 = 0, x1 = 1, Px;

  while (fabsl(x1 - x0) > xtol) {
    Px = gsl_cdf_beta_P (x, a, b);
    if (fabsl(Px - P) < Ptol) {
      /* return as soon as approximation is good enough, including on
         the first iteration */
      return x;  
    } else if (Px < P) {
      x0 = x;
    } else if (Px > P) {
      x1 = x;
    }
    x = 0.5L * (x0 + x1);
  }
  return x;
}  
#define CDF_ERROR(s, y) do { cerr << s << endl; exit(1); } while (0)

#include <gsl/gsl_math.h>
static long double my_gsl_cdf_beta_Pinv (const long double P, const long double a, const long double b)
{
	// PP3(P,a,b);
  long double x, mean;

  if (P < 0.0L || P > 1.0L)
    {
      CDF_ERROR ("P must be in range 0 < P < 1", GSL_EDOM);
    }

  if (a < 0.0L)
    {
      CDF_ERROR ("a < 0", GSL_EDOM);
    }

  if (b < 0.0L)
    {
      CDF_ERROR ("b < 0", GSL_EDOM);
    }

  if (P == 0.0L)
    {
      return 0.0L;
    }

  if (P == 1.0L)
    {
      return 1.0L;
    }

  if (P > 0.5L)
    {
      return my_gsl_cdf_beta_Qinv (1.0L - P, a, b);
    }

  mean = a / (a + b);

  if (P < 0.1)
    {
      /* small x */

      long double lg_ab = gsl_sf_lngamma (a + b);
      long double lg_a = gsl_sf_lngamma (a);
      long double lg_b = gsl_sf_lngamma (b);

      long double lx = (logl (a) + lg_a + lg_b - lg_ab + logl (P)) / a;
      if (lx <= 0) {
        x = expl (lx);             /* first approximation */
        x *= pow (1 - x, -(b - 1) / a);   /* second approximation */
      } else {
        x = mean;
      }

      if (x > mean)
        x = mean;
    }
  else
    {
      /* Use expected value as first guess */
      x = mean;
    }

  /* Do bisection to get closer */
  x = bisect (x, P, a, b, 0.01, 0.01);

  {
    long double lambda, dP, phi;
    unsigned int n = 0;

  start:
    dP = P - gsl_cdf_beta_P (x, a, b);
    phi = gsl_ran_beta_pdf (x, a, b);

    if (dP == 0.0L || n++ > 64) {
      goto end;
    }

    lambda = dP / GSL_MAX (2 * fabsl (dP / x), phi);

    {
      long double step0 = lambda;
      long double step1 = -((a - 1) / x - (b - 1) / (1 - x)) * lambda * lambda / 2;

      long double step = step0;

      if (fabsl (step1) < fabsl (step0))
        {
          step += step1;
        }
      else
        {
          /* scale back step to a reasonable size when too large */
          step *= 2 * fabsl (step0 / step1);
        };

      if (x + step > 0 && x + step < 1)
        {
          x += step;
        }
      else
        {
          x = sqrtl (x) * sqrtl (mean);   /* try a new starting point */
        }

      if (fabsl (step0) > 1e-10L * x)
        goto start;
    }

  end:

    if (0 && fabsl(dP) > GSL_SQRT_DBL_EPSILON * P)
      {
	      // PP3(dP, n, x);
	      // PP(gsl_cdf_beta_P(x, a, b));
        GSL_ERROR_VAL("inverse failed to converge", GSL_EFAILED, GSL_NAN);
      }

    return x;
  }
}

static long double
my_gsl_cdf_beta_Qinv (const long double Q, const long double a, const long double b)
{
	// PP3(Q,a,b);

  if (Q < 0.0L || Q > 1.0L)
    {
      CDF_ERROR ("Q must be inside range 0 < Q < 1", GSL_EDOM);
    }

  if (a < 0.0L)
    {
      CDF_ERROR ("a < 0", GSL_EDOM);
    }

  if (b < 0.0L)
    {
      CDF_ERROR ("b < 0", GSL_EDOM);
    }

  if (Q == 0.0L)
    {
      return 1.0L;
    }

  if (Q == 1.0L)
    {
      return 0.0L;
    }

  if (Q > 0.5L)
    {
      return my_gsl_cdf_beta_Pinv (1 - Q, a, b);
    }
  else
    {
      return 1 - my_gsl_cdf_beta_Pinv (Q, b, a);
    };
}
