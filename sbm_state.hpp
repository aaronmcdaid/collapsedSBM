#include <list>
#include <vector>
#include "shmGraphRaw.hpp"
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

	int appendEmptyCluster();
	void moveNode(const int n, const int newClusterID);
	int isolateNode(const int n); // create a new (probably temporary) cluster to hold this one node
	void unIsolateTempNode(const int n, const int newClusterID); // move a node from its 'temporary' cluster to an existing cluster

	// internalcheck
	void internalCheck() const;

	// summaries
	void shortSummary() const;

	struct EdgeCounts {
		typedef boost::unordered_map< int , boost::unordered_map<int,int> >::value_type outer_value_type;
		typedef                             boost::unordered_map<int,int>  ::value_type inner_value_type;
		boost::unordered_map< int , boost::unordered_map<int,int> > counts;
		void   inform(const int cl1, const int cl2) ; // inform us of an edge between cl1 and cl2
		void uninform(const int cl1, const int cl2) ; // UNinform us of an edge between cl1 and cl2
		private:
		void partialUnInform(const int cl1, const int cl2);
	};
	void informNodeMove(const int n, const int oldcl, const int newcl); // a node has just moved from one cluster to another. We must consider it's neighbours for _edgeCounts
	EdgeCounts _edgeCounts;
};

} // namespace sbm

