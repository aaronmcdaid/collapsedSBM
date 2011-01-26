#include <list>
#include <vector>
#include "shmGraphRaw.hpp"
namespace sbm {
struct State {
	const shmGraphRaw::ReadableShmGraphBase * const _g; // the graph
	const int _N; // the number of nodes in the graph
	explicit State(const shmGraphRaw::ReadableShmGraphBase * const g);

	// the clustering
	int _k; // the number of clusters (including empty ones)
	struct Cluster {
		std::list<int> members; // nodes ids of the members
	}; // each cluster to know which nodes are in it
	std::vector< Cluster > clusters; // numbered 0 to k-1
	std::vector< int > cluster_id; // the cluster that each node is in
	std::vector< std::list<int>::iterator > its; // an iterator into the relevant part of Cluster::members

	// internalcheck
	void internalCheck() const;
};

} // namespace sbm

