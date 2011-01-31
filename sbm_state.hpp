#include <list>
#include <vector>
#include "shmGraphRaw.hpp"

static inline double LOG2GAMMA(double x) {
	return M_LOG2E * gsl_sf_lngamma(x);
}
static inline double LOG2FACT(double x) {
	return M_LOG2E * gsl_sf_lnfact(x);
}

namespace sbm {
typedef shmGraphRaw::ReadableShmGraphTemplate<shmGraphRaw::PlainMem> GraphType;
struct State {
	const GraphType * const _g; // the graph
	const int _N; // the number of nodes in the graph
	explicit State(const GraphType * const g);

	// the clustering
	int _k; // the number of clusters (including empty ones)
	struct Cluster {
		std::list<int> members; // nodes ids of the members
		const int order() const;
		std::list<int>::iterator newMember(const int n);
	}; // each cluster to know which nodes are in it
	std::vector< Cluster* > clusters; // numbered 0 to k-1
	std::vector< int > cluster_id; // the cluster that each node is in
	std::vector< std::list<int>::iterator > its; // an iterator into the relevant part of Cluster::members
	int NonEmptyClusters;
	mutable long double SumOfLog2LOrders;

	int appendEmptyCluster();
	void deleteClusterFromTheEnd() ;
	void moveNode(const int n, const int newClusterID);
	int isolateNode(const int n); // create a new (probably temporary) cluster to hold this one node
	void unIsolateTempNode(const int n, const int newClusterID); // move a node from its 'temporary' cluster to an existing cluster

	// internalcheck
	void internalCheck() const;

	// summaries
	void shortSummary() const;
	void summarizeEdgeCounts() const ;
	void blockDetail() const ;

	// counting the edges between each pair of clusters
	struct EdgeCounts {
		typedef boost::unordered_map< int , boost::unordered_map<int,int> >::value_type outer_value_type;
		typedef                             boost::unordered_map<int,int>  ::value_type inner_value_type;
		boost::unordered_map< int , boost::unordered_map<int,int> > counts;
		void   inform(const int cl1, const int cl2) ; // inform us of an edge between cl1 and cl2
		void uninform(const int cl1, const int cl2) ; // UNinform us of an edge between cl1 and cl2
		int       get(const int cl1, const int cl2) const throw() ;
		private:
		void partialUnInform(const int cl1, const int cl2);
	};
	void informNodeMove(const int n, const int oldcl, const int newcl); // a node has just moved from one cluster to another. We must consider it's neighbours for _edgeCounts
	void moveNodeAndInformOfEdges(const int n, const int newcl);
	EdgeCounts _edgeCounts;

	// pmf - probability mass function to guide sampling
	// The overall objective can be divided into four parts:
	// 1. The K-dependant part, including its prior
	// 2. The bits for each cluster's size
	// 3. The edges "baseline", where there are assumed to be no edges
	// 4. The edges correction over the baseline.
	long double P_z_K() const; // 1.
	long double P_z_orders() const; // 2.
	long double P_z() const;
	long double P_edges_given_z_slow() const;
	long double P_edges_given_z() const;
	long double P_edges_given_z_baseline() const; // 3.
	long double P_edges_given_z_correction() const; // 4.
	long double P_edges_given_z_correction_JustOneCluster(const int clusterID) const;
	long double pmf_slow() const;
	long double pmf() const;
};

} // namespace sbm

