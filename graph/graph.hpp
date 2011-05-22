#ifndef _GRAPH_HPP_
#define _GRAPH_HPP_

#include <utility>

namespace graph {

class VerySimpleGraph { // consecutive ints. No attributes, or directions, or anything
public:
	virtual int numNodes() const = 0;
	virtual int numRels() const = 0;
	virtual const std :: pair<int, int> & EndPoints(int relId) const = 0;
	virtual ~ VerySimpleGraph() {}
};

} // namespace graph
#endif
