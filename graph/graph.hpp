#ifndef _GRAPH_HPP_
#define _GRAPH_HPP_

#include <utility>
#include <vector>
#include <stdint.h>

namespace graph {

class VerySimpleGraphInterface { // consecutive ints. No attributes, or directions, or anything
public:
	virtual int numNodes() const = 0;
	virtual int numRels() const = 0;
	virtual const std :: pair<int, int> & EndPoints(int relId) const = 0;
	virtual ~ VerySimpleGraphInterface() {}
	virtual const std :: vector<int32_t> & neighbouring_rels_in_order(const int32_t node_id) = 0;
};

} // namespace graph
#endif
