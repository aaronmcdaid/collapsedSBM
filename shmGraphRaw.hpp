#ifndef SHMGRAPHRAW_H
#define SHMGRAPHRAW_H

#include <vector> 
#include <sstream> 

#include <boost/unordered_set.hpp>


#include "Range.hpp"

#include "graph/graph.hpp"
#include "graph/weights.hpp"
namespace shmGraphRaw {

class ReadableShmGraphBase /*: public graph :: VerySimpleGraphInterface*/ {
public:
	bool hasASelfLoop;
	virtual ~ReadableShmGraphBase();
	virtual int32_t numNodes() const = 0;
	virtual int32_t numRels() const = 0;
	virtual const std :: pair<int32_t, int32_t> & EndPoints(int32_t relId) const = 0;
	virtual bool are_connected(int32_t v1, int32_t v2) const = 0;
	virtual int32_t oppositeEndPoint(int32_t relId, int32_t oneEnd) const; // impure function.

	virtual std :: pair<const char*, const char*> EndPointsAsStrings(int32_t relId) const = 0;
	virtual const char * NodeAsString(int32_t v) const = 0;
	virtual int32_t StringToNodeId(const char *s) const = 0;
	virtual std :: string WhichNode(int32_t v) const; // impure function

	virtual int32_t degree(int32_t v) const = 0; // implemented in ReadableShmGraphTemplate<T>
};

class ReadableShmGraphTemplate : public ReadableShmGraphBase { // this is mostly just an interface, but note that oppositeEndPoint is defined in this class
public:
	virtual int degree(int v) const { return this->myRels(v).size(); }
	virtual const boost :: unordered_set<int> & myRels(int n) const = 0;
};

ReadableShmGraphTemplate * loadEdgeList(const char *graphTextFileName, const bool selfloops_allowed, graph :: weights :: EdgeDetailsInterface *edge_details);

struct SelfLoopsNotSupported : public std :: exception {
};
} // namespace shmGraphRaw 

#endif 
