#ifndef SHMGRAPHRAW_H
#define SHMGRAPHRAW_H

#include "strings.hpp"

#include <set> 
#include <map> 
#include <vector> 
#include <sstream> 

#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/hashed_index.hpp>
namespace bmi = boost :: multi_index;

#include "Range.hpp"

#include "graph/weights.hpp"
namespace shmGraphRaw {


class VerySimpleGraph { // consecutive ints. No attributes, or directions, or anything
public:
	virtual int numNodes() const = 0;
	virtual int numRels() const = 0;
	virtual ~ VerySimpleGraph() {}
};




class ReadableShmGraphBase : public VerySimpleGraph {
public:
	bool hasASelfLoop;
	virtual ~ReadableShmGraphBase();
	virtual int numNodes() const = 0;
	virtual int numRels() const = 0;
	virtual const std :: pair<int, int> & EndPoints(int relId) const = 0;
	virtual bool are_connected(int v1, int v2) const = 0;
	virtual int oppositeEndPoint(int relId, int oneEnd) const; // impure function.

	virtual std :: pair<const char*, const char*> EndPointsAsStrings(int relId) const = 0;
	virtual const char * NodeAsString(int v) const = 0;
	virtual int StringToNodeId(const char *s) const = 0;
	virtual std :: string WhichNode(int v) const; // impure function

	virtual int degree(int v) const = 0; // implemented in ReadableShmGraphTemplate<T>
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
