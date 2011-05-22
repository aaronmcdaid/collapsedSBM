#ifndef _GRAPH_NETWORK_
#define _GRAPH_NETWORK_

#include "graph.hpp"
#include "weights.hpp"
#include <string>
#include <memory>
namespace graph {

struct Network {
	/* A Network is a VerySimpleGraph with some extra attributes.
	 * The nodes will have string names, and the edges might have directionality and weights
	 */
	std :: auto_ptr<graph :: VerySimpleGraph> vsg;
	const std :: string & node_name(int node_id);
	std :: auto_ptr<graph :: weights :: EdgeDetailsInterface> edge_weights;
	Network(const bool directed, const bool weighted);
	virtual ~ Network() throw(); // this forces derivations to declare a destructor.
};

} // namespace graph

#endif
