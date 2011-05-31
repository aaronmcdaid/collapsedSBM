#ifndef _SBM_STATE_HPP_
#define _SBM_STATE_HPP_

#include <list>
#include <vector>
#include <set>
#include <gsl/gsl_sf.h>
#include "graph/network.hpp"
#include <cmath>


static inline double LOG2GAMMA(const double x) {
	assert(x>0);
	return M_LOG2E * gsl_sf_lngamma(x);
}
static inline double LOG2GAMMA(const long double x) {
	return LOG2GAMMA(static_cast<double>(x));
}
static inline double LOG2FACT(const int x) {
	assert(x>0);
	return M_LOG2E * gsl_sf_lnfact(x);
}

namespace sbm {
typedef graph :: NetworkInterfaceConvertedToStringWithWeights GraphType;

struct Cluster {
private:
		int _order;
		std :: list<int> members; // nodes ids of the members
public:
		int order() const;
		std :: list<int>:: iterator newMember(const int n);
		void eraseMember(const std :: list<int>:: iterator it);
		const std :: list<int> & get_members() const;
		Cluster();
}; // each cluster to know which nodes are in it

struct Labelling {
	const int _N;
	int _k;
	const long double _alpha;
	std :: vector< Cluster* > clusters; // numbered 0 to k-1
	std :: vector< int > cluster_id; // the cluster that each node is in
	std :: vector< std :: list<int>:: iterator > its; // an iterator into the relevant part of Cluster :: members
	int NonEmptyClusters;
	mutable long double SumOfLog2Facts;
	std :: vector<double> log2GammaAlphaPlus;
	Labelling(const int _N, const long double _alpha);


	int appendEmptyCluster();
	void deleteClusterFromTheEnd();
	void deltaSumOfLog2Facts(const int oldClOrder, const int clOrder);
	void fixUpIterators(const int n, Cluster *cl, Cluster *oldcl);
	int moveNode(const int n, const int newClusterID);
	void swapClusters(const int cl1, const int cl2);
};

struct ObjectiveFunction;

struct State {
	const GraphType * const _g; // the graph
	const graph :: VerySimpleGraphInterface * const vsg;
	const graph :: weights :: EdgeDetailsInterface * const _edge_details;
	long double total_edge_weight;
	const int _N; // the number of nodes in the graph
	const long double _alpha; // the parameter to the Dirichlet prior on z
	const bool _mega; // if this is true, the print less
	explicit State(const GraphType * const g, const bool mega = false);

	Labelling	labelling;
	// the clustering
	int _k; // the number of clusters (including empty ones)

	int appendEmptyCluster();
	void deleteClusterFromTheEnd() ;
	int moveNode(const int n, const int newClusterID);
	int isolateNode(const int n); // create a new (probably temporary) cluster to hold this one node
	void unIsolateTempNode(const int n, const int newClusterID); // move a node from its 'temporary' cluster to an existing cluster
	void swapClusters(const int cl1, const int cl2);

	// internalcheck
	void internalCheck() const;

	// summaries
	void shortSummary(const ObjectiveFunction *obj, const std :: vector<int> *groundTruth) const;
	void summarizeEdgeCounts() const ;
	void blockDetail(const ObjectiveFunction *obj) const ;

	// counting the edges between each pair of clusters
	struct EdgeCounts {
		const graph :: weights :: EdgeDetailsInterface * const _edge_details;
		mutable std :: vector< std :: vector<long double> > counts;
		long double externalEdgeWeight;
		EdgeCounts(const graph :: weights :: EdgeDetailsInterface *edge_details);
		void   inform(const int cl1, const int cl2, int relId) ; // inform us of an edge between cl1 and cl2
		void uninform(const int cl1, const int cl2, int relId) ; // UNinform us of an edge between cl1 and cl2
		long double & at(int i,int j) ;
		long double read(int i, int j) const ;
		private:
		friend void State :: summarizeEdgeCounts() const;
		friend void State :: internalCheck() const;
		friend void State :: swapClusters(int,int);
	};
	void informNodeMove(const int n, const int oldcl, const int newcl); // a node has just moved from one cluster to another. We must consider it's neighbours for _edgeCounts
	int moveNodeAndInformOfEdges(const int n, const int newcl);
	int moveNodeAndInformOfEdges2(const int n, const int newcl);
	EdgeCounts _edgeCounts;

	// pmf - probability mass function to guide sampling
	// The overall objective can be divided into four parts:
	// 1. The K-dependant part, including its prior
	// 2. The bits for each cluster's size
	// 3. The edges "baseline", where there are assumed to be no edges
	// 4. The edges correction over the baseline.
	long double P_z_K() const; // 1.
	long double P_z_orders() const; // 2.
	long double P_z_slow() const;
	long double P_z() const;

	long double P_edges_given_z_slow(const ObjectiveFunction *obj) const;
	long double P_edges_given_z(const ObjectiveFunction *obj) const;

	long double pmf_slow(const ObjectiveFunction *obj) const;
	long double pmf(const ObjectiveFunction *obj) const;
};
struct ObjectiveFunction {
	const bool selfloops;
	const bool directed;
	const bool weighted;
	// these first three bools decide which blocks to loop over, and how to calculate the valid pairs.
	// But give the total edge weight of a block, which form of function will actually be used?
	ObjectiveFunction(const bool s, const bool d, const bool w);
	virtual long double log2OneBlock(const long double edge_total, const long int pairs, bool isDiagonal) const = 0;
	bool isValidBlock(const int i, const int j) const;
	bool isTwoSided(const int i, const int j) const; // should the other direction be included when counting the edges?
	long double relevantWeight(const int i, const int j, const State :: EdgeCounts *edge_counts) const; // this might include the other direction, if it's undirected
	long int numberOfPairsInBlock(const int i, const int j, const Labelling *edge_counts) const;
	virtual ~ ObjectiveFunction() {}
};
struct ObjectiveFunction_Bernoulli : public ObjectiveFunction {
	ObjectiveFunction_Bernoulli(const bool s, const bool d, const bool w);
	virtual long double log2OneBlock(const long double edge_total, const long int pairs, bool isDiagonal) const;
};
struct ObjectiveFunction_Poisson : public ObjectiveFunction {
	static const long double s = 1.0L;
	static const long double theta = 1000.0L;
	ObjectiveFunction_Poisson(const bool s, const bool d, const bool w);
	virtual long double log2OneBlock(const long double edge_total, const long int pairs, bool isDiagonal) const;
};
	long double assertNonPositiveFinite_line(const long double x, const int lineno);

} // namespace sbm

#endif
