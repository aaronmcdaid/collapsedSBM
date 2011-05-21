#ifndef SHMGRAPHRAW_H
#define SHMGRAPHRAW_H

#include <set> 
#include <map> 
#include <vector> 
#include <sstream> 


#include <boost/interprocess/managed_mapped_file.hpp>
#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/allocators/allocator.hpp>
namespace bip = boost :: interprocess;

typedef boost :: interprocess :: managed_mapped_file MMapType; // typedef boost :: interprocess :: managed_shared_memory MMapType;

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

struct idT{}; // dummy tag type for use in boost :: interprocess index names
struct nameT{}; // dummy tag type for use in boost :: interprocess index names
struct nodeIdsT{}; // dummy tag type for use in boost :: interprocess index names

struct relationship
{
	int         relId;
	typedef std :: pair<int,int> relPairType;
	std :: pair<int,int>		nodeIds;
	relationship( int id_ , const std :: pair<int,int> &nodes_);
};
typedef bmi :: multi_index_container<
 		relationship,
  		bmi :: indexed_by<
	 		bmi :: hashed_unique  <bmi :: tag<idT>,  BOOST_MULTI_INDEX_MEMBER(relationship,int,relId)>,
	 		bmi :: hashed_unique  <bmi :: tag<nodeIdsT>,BOOST_MULTI_INDEX_MEMBER(relationship,relationship :: relPairType,nodeIds)>
		>
> relationship_set;

typedef boost :: unordered_map < int, boost :: unordered_set<int> > neighbours_to_relationships_map;


class StrH { // string handle. It just wraps an int that refers to the memory mapped file
	int i;
public:
	explicit StrH(int _i);
	bool operator == (StrH r) const;
	struct hasher;
	int get_underlying_id() const; // Shouldn't really call this too much.
};

class StringArray {
public:
	virtual const char * operator[] (StrH s) const = 0;
	virtual StrH StringToStringId(const char *s) const = 0;
};

class ReadableShmGraphBase {
public:
	bool hasASelfLoop;
	virtual ~ReadableShmGraphBase();
	virtual int numNodes() const = 0;
	virtual int numRels() const = 0;
	virtual int numNodesWithAtLeastOneRel() const = 0;
	virtual std :: pair<const char*, const char*> EndPointsAsStrings(int relId) const = 0;
	virtual const char * NodeAsString(int v) const = 0;
	virtual int StringToNodeId(const char *s) const = 0;
	virtual const std :: pair<int, int> & EndPoints(int relId) const = 0;
	virtual bool are_connected(int v1, int v2) const = 0;
	virtual int oppositeEndPoint(int relId, int oneEnd) const; // impure function.
	virtual std :: string WhichNode(int v) const; // impure function

	virtual int degree(int v) const = 0; // implemented in ReadableShmGraphTemplate<T>
	virtual const std :: set<int> & neighbours(int v) const = 0; // implemented in ReadableShmGraphTemplate<T>
};

class ReadableShmGraphTemplate : public ReadableShmGraphBase { // this is mostly just an interface, but note that oppositeEndPoint is defined in this class
private:
	mutable std :: map<int, std :: set<int> > neighbours_cache;
public:
	virtual int degree(int v) const { return this->myRels(v).size(); }
	virtual const std :: set<int> & neighbours(int v) const { // sorted list of neighbours. Sorted by internal int id, not by the original string name
		std :: set<int> &  neighs = neighbours_cache[v]; // Will create an empty one, if it hasn't been requested before
		if(neighs.size() == 0 && this->degree(v) != 0) {
			forEach(int rel, amd :: mk_range(this->myRels(v))) {
				int otherEnd = this->oppositeEndPoint(rel, v);
				neighs.insert(otherEnd);
			}
		}
		return neighs;
	}
	virtual const boost :: unordered_set<int> & myRels(int n) const = 0;
};

ReadableShmGraphTemplate * loadEdgeList(const char *graphTextFileName, const bool selfloops_allowed, graph :: weights :: EdgeDetailsInterface *edge_details);

struct SelfLoopsNotSupported : public std :: exception {
};
} // namespace shmGraphRaw 

#endif 
