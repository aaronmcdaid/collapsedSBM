#include "strings.hpp"
#include "graph.hpp"
#include "weights.hpp"
#include <memory>
namespace graph {

struct Network {
	/* A Network is a VerySimpleGraph with some extra attributes.
	 * The nodes will have string names, and the edges might have directionality and weights
	 */
	std :: auto_ptr<strings :: StringArray> string_array;
	std :: auto_ptr<graph :: VerySimpleGraph> vsg;
	std :: auto_ptr<graph :: weights :: EdgeDetailsInterface> edge_weights;
};

} // namespace graph
