#ifndef SHMGRAPHRAW_H
#define SHMGRAPHRAW_H

#include <set> 
#include <map> 


#include <boost/interprocess/managed_mapped_file.hpp>
#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/allocators/allocator.hpp>

typedef boost::interprocess::managed_mapped_file MMapType;
// typedef boost::interprocess::managed_shared_memory MMapType;

#include <boost/unordered_set.hpp>

#include "Range.hpp"

namespace shmGraphRaw {

typedef boost::unordered_set<int, boost::hash<int>,  std::equal_to<int>, boost::interprocess::allocator< int, MMapType::segment_manager> > mmap_uset_of_ints;

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

class ReadableShmGraph { // this is mostly just an interface, but not that oppositeEndPoint is defined in this class
public:
	virtual ~ReadableShmGraph() {};
	virtual int numNodes() const = 0;
	virtual int numRels() const = 0;
	virtual int numNodesWithAtLeastOneRel() const = 0;
	virtual const mmap_uset_of_ints & myRels(int n) const = 0;
	virtual std::pair<const char*, const char*> EndPointsAsStrings(int relId) const = 0;
	virtual const char * NodeAsString(int v) const = 0;
	virtual int StringToNodeId(const char *s) const = 0;
	virtual const std::pair<int, int> & EndPoints(int relId) const = 0;
	virtual bool are_connected(int v1, int v2) const = 0;
	virtual int oppositeEndPoint(int relId, int oneEnd) const; // impure function.
	virtual std::string WhichNode(int v) const; // impure function
	virtual int degree(int v) const { return this->myRels(v).size(); }
	mutable std::map<int, std::set<int> > neighbours_cache;
	virtual const std::set<int> & neighbours(int v) const { // sorted list of neighbours. Sorted by internal int id, not by the original string name
		std::set<int> &  neighs = neighbours_cache[v]; // Will create an empty one, if it hasn't been requested before
		if(neighs.size() == 0 && this->degree(v) != 0) {
			forEach(int rel, amd::mk_range(this->myRels(v))) {
				int otherEnd = this->oppositeEndPoint(rel, v);
				neighs.insert(otherEnd);
			}
		}
		return neighs;
	}
};


ReadableShmGraph * loadMmapFile(const char *directory, const char * graphTextFileName);

} // namespace shmGraphRaw 

#endif 
