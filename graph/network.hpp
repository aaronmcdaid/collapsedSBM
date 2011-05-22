#ifndef _GRAPH_NETWORK_
#define _GRAPH_NETWORK_

#include "graph.hpp"
#include "weights.hpp"
#include <string>
#include <memory>
#include <cstdlib>
namespace graph {

template <class NodeNameT>
struct Network;

struct NodeNameIsInt32; // the default node_name type. It's nice to have the names sorted by this, so that "10" comes after "2"
struct NodeNameIsString; // .. but if the user wants to specify arbitrary strings in their edge list, they can do so explicitly.
struct BadlyFormattedNodeName : public std :: exception { // if the text doesn't define an int properly when trying to use NodeNameIsInt32
};

typedef Network<NodeNameIsInt32> NetworkInt32;
typedef Network<NodeNameIsString> NetworkString;

struct NodeNameIsInt32 {
	typedef int32_t value_type;
	inline static value_type fromString(const std :: string &s) {
		assert(!s.empty());
		value_type i;
		char *end_ptr;
		i = strtol(s.c_str(), &end_ptr, 10);
		if(*end_ptr == '\0') // success, according to http://linux.die.net/man/3/strtol
			return i;
		else
			throw BadlyFormattedNodeName();
	}
};
struct NodeNameIsString {
	typedef std :: string value_type;
	inline static value_type fromString(const std :: string &s) {
		assert(!s.empty());
		return s;
	}
};
template <class NodeNameT>
struct Network {
	/* A Network is a VerySimpleGraph with some extra attributes.
	 * The nodes will have string names, and the edges might have directionality and weights
	 */
	std :: auto_ptr<const graph :: VerySimpleGraph> plain_graph;
	std :: auto_ptr<const graph :: weights :: EdgeDetailsInterface> edge_weights;
public: // I should make the above private some time!
	virtual ~ Network() throw(); // this forces derivations to declare a destructor.
	Network(const bool directed, const bool weighted);
	virtual int32_t numNodes() const { assert(plain_graph.get()); return plain_graph->numNodes(); }
	virtual int32_t numRels()  const { assert(plain_graph.get()); return plain_graph->numRels(); }
};

} // namespace graph

#endif
