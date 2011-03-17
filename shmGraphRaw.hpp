#ifndef SHMGRAPHRAW_H
#define SHMGRAPHRAW_H

#include <set> 
#include <map> 
#include <vector> 
#include <sstream> 


#include <boost/interprocess/managed_mapped_file.hpp>
#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/allocators/allocator.hpp>
namespace bip = boost::interprocess;

typedef boost::interprocess::managed_mapped_file MMapType; // typedef boost::interprocess::managed_shared_memory MMapType;

#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/hashed_index.hpp>
namespace bmi = boost::multi_index;

#include "Range.hpp"

namespace shmGraphRaw {

struct idT{}; // dummy tag type for use in boost::interprocess index names
struct nameT{}; // dummy tag type for use in boost::interprocess index names
struct nodeIdsT{}; // dummy tag type for use in boost::interprocess index names

struct relationship
{
	int         relId;
	typedef std::pair<int,int> relPairType;
	std::pair<int,int>		nodeIds;
	relationship( int id_ , const std::pair<int,int> &nodes_);
};
typedef bmi::multi_index_container<
 		relationship,
  		bmi::indexed_by<
	 		bmi::hashed_unique  <bmi::tag<idT>,  BOOST_MULTI_INDEX_MEMBER(relationship,int,relId)>,
	 		bmi::hashed_unique  <bmi::tag<nodeIdsT>,BOOST_MULTI_INDEX_MEMBER(relationship,relationship::relPairType,nodeIds)>
		>
		, MMapType  ::allocator<relationship>::type
> relationship_set_MapMem;
typedef bmi::multi_index_container<
 		relationship,
  		bmi::indexed_by<
	 		bmi::hashed_unique  <bmi::tag<idT>,  BOOST_MULTI_INDEX_MEMBER(relationship,int,relId)>,
	 		bmi::hashed_unique  <bmi::tag<nodeIdsT>,BOOST_MULTI_INDEX_MEMBER(relationship,relationship::relPairType,nodeIds)>
		>
		// , MMapType  ::allocator<relationship>::type
> relationship_set_PlainMem;

struct MapMem {
	typedef boost::unordered_set<int, boost::hash<int>,  std::equal_to<int>, boost::interprocess::allocator< int, MMapType::segment_manager> > mmap_uset_of_ints;
	typedef mmap_uset_of_ints neighbouring_relationship_set;
	typedef std::pair<const int, neighbouring_relationship_set> valtype;
	typedef bip::allocator< valtype, MMapType::segment_manager> ShmemAllocator;
	typedef boost::unordered_map
    		< int               , neighbouring_relationship_set
    		, boost::hash<int>  ,std::equal_to<int>
    		, ShmemAllocator>
			neighbours_to_relationships_map;
	typedef relationship_set_MapMem relationship_set;
	static const int isMapMem = 1;
	typedef MMapType segment_type;
};
struct PlainMem {
	typedef boost::unordered_set<int, boost::hash<int>,  std::equal_to<int> > mmap_uset_of_ints;
	typedef mmap_uset_of_ints neighbouring_relationship_set;
	typedef std::pair<const int, neighbouring_relationship_set> valtype;
	// typedef bip::allocator< valtype, MMapType::segment_manager> ShmemAllocator;
	typedef boost::unordered_map
    		< int               , neighbouring_relationship_set
    		, boost::hash<int>  ,std::equal_to<int>
    		>
			neighbours_to_relationships_map;
	typedef relationship_set_PlainMem relationship_set;
	static const int isMapMem = 0;
	enum nil{};
	typedef nil segment_type;
};


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
	virtual std::pair<const char*, const char*> EndPointsAsStrings(int relId) const = 0;
	virtual const char * NodeAsString(int v) const = 0;
	virtual int StringToNodeId(const char *s) const = 0;
	virtual const std::pair<int, int> & EndPoints(int relId) const = 0;
	virtual bool are_connected(int v1, int v2) const = 0;
	virtual int oppositeEndPoint(int relId, int oneEnd) const; // impure function.
	virtual std::string WhichNode(int v) const; // impure function

	virtual int degree(int v) const = 0; // implemented in ReadableShmGraphTemplate<T>
	virtual const std::set<int> & neighbours(int v) const = 0; // implemented in ReadableShmGraphTemplate<T>
};

template <class T>
class ReadableShmGraphTemplate : public ReadableShmGraphBase { // this is mostly just an interface, but note that oppositeEndPoint is defined in this class
private:
	mutable std::map<int, std::set<int> > neighbours_cache;
public:
	virtual int degree(int v) const { return this->myRels(v).size(); }
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
	virtual const typename T::mmap_uset_of_ints & myRels(int n) const = 0;
};

struct EdgeDetailsInterface {
	virtual long double getl2h(const int relId) = 0;
	virtual long double geth2l(const int relId) = 0;
};
template <class W>
struct EdgeDetails : public EdgeDetailsInterface {
	std:: vector < typename W::datumT > dw; // the directions and weights of all the edges
	int size() const {
		return this->dw.size();
	}
	void new_rel(int relId, std::pair<int,int> nodeIds, std::string &weight) {
		if(relId == (int)this->dw.size())
			this->dw.push_back( typename W::datumT() );
		assert(relId+1 == (int)this->dw.size());
		this->dw.at(relId).inform(nodeIds.first > nodeIds.second, weight);
		// this->dw.back().inform(nodeIds.first > nodeIds.second, weight);
	}
	virtual long double getl2h(const int relId) { // this returns the value in the undirected case, and it handles self loops
		return this->dw.at(relId).getl2h();
	}
	virtual long double geth2l(const int relId) { // this is only relevant in directed graphs.
		return this->dw.at(relId).geth2l();
	}
};
struct NoDetails { // unweighted, undirected
	typedef struct nil {
		void inform(const bool, const std::string) const {
		}
		long double getl2h() const {
			return 1.0L;
		}
		long double geth2l() const {
			return 0.0L;
		}
	} datumT;
};
struct DirectedLDoubleWeights {
	typedef struct LdblPair : public std:: pair<long double,long double> {
		LdblPair() {
			this->first = this->second = 0.0L;
		}
		struct DuplicateWeightedEdge : public std::exception {
		};
		void inform(const bool highToLow, const std::string weight) { // the weight string might be empty. I suppose we let that default to 0
			if( highToLow ? this->second : this->first != 0) {
				throw DuplicateWeightedEdge();
			}
			std:: istringstream oss(weight);
			long double w = 0;
			assert(oss.peek() != EOF);
			oss >> w;
			assert(oss.peek() == EOF);
			if(highToLow) {
				this->second = w;
			} else {
				this->first = w;
			}
		}
		long double getl2h() const {
			return this->first;
		}
		long double geth2l() const {
			return this->second;
		}
	} datumT; // the type needed to store the weights in each direction.
};
struct DirectedNoWeights {
	typedef struct IntPair : public std:: pair<int,int> {
		IntPair() {
			this->first = this->second = 0;
		}
		void inform(const bool highToLow, const std::string) { // the weight string might be empty. I suppose we let that default to 0
			if(highToLow) {
				this->second = 1;
			} else {
				this->first = 1;
			}
		}
		long double getl2h() const {
			return this->first;
		}
		long double geth2l() const {
			return this->second;
		}
	} datumT; // the type needed to store the weights in each direction.
};
struct WeightNoDir {
	typedef struct LdblPair {
		long double weight;
		LdblPair() {
			this->weight = 0.0L;
		}
		struct DuplicateWeightedEdge : public std::exception {
		};
		void inform(const bool , const std::string weight) { // the weight string might be empty. I suppose we let that default to 0
			if( this->weight != 0) {
				throw DuplicateWeightedEdge();
			}
			std:: istringstream oss(weight);
			long double w = 0;
			assert(oss.peek() != EOF);
			oss >> w;
			assert(oss.peek() == EOF);
			this->weight = w;
		}
		long double getl2h() const {
			return this->weight;
		}
		long double geth2l() const {
			return 0.0L;
		}
	} datumT; // the type needed to store the weights in each direction.
};

template<class T, class W>
ReadableShmGraphTemplate<T> * loadEdgeList(const char * graphTextFileName, const bool selfloops_allowed, EdgeDetails<W> &edge_details, const char * directory = NULL);

//typedef ReadableShmGraphTemplate<MapMem> ReadableShmGraph; // TODO: Deprecate this

struct SelfLoopsNotSupported : public std::exception {
};
} // namespace shmGraphRaw 

#endif 
