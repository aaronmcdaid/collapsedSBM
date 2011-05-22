#ifndef _GRAPH_HPP_
#define _GRAPH_HPP_
namespace graph {

class VerySimpleGraph { // consecutive ints. No attributes, or directions, or anything
public:
	virtual int numNodes() const = 0;
	virtual int numRels() const = 0;
	virtual ~ VerySimpleGraph() {}
};

} // namespace graph
#endif
