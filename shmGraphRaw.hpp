#ifndef SHMGRAPHRAW_H
#define SHMGRAPHRAW_H

#include <boost/interprocess/managed_mapped_file.hpp>
#include <boost/interprocess/allocators/allocator.hpp>

#include <boost/unordered_set.hpp>

namespace shmGraphRaw {

typedef boost::unordered_set<int, boost::hash<int>,  std::equal_to<int>, boost::interprocess::allocator< int, boost::interprocess::managed_mapped_file::segment_manager> > mmap_uset_of_ints;

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
	virtual int oppositeEndPoint(int relId, int oneEnd) const; // impure function.
	virtual std::string WhichNode(int v) const; // impure function
};


ReadableShmGraph * loadMmapFile(const char *directory, const char * graphTextFileName);

} // namespace shmGraphRaw 

#endif 
