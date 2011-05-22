#include "shmGraphRaw.hpp"

#include <vector>
#include <memory>
#include <fstream>
#include <string.h>

#include "aaron_utils.hpp"

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

ReadableShmGraphBase ::~ReadableShmGraphBase() {
}

int ReadableShmGraphBase :: oppositeEndPoint(int relId, int oneEnd) const {
		const std :: pair<int, int> & eps = this->EndPoints(relId);
		if(eps.first == oneEnd) {
			return eps.second;
		} else {
			assert(eps.second == oneEnd);
			return eps.first;
		}
}
std :: string ReadableShmGraphBase :: WhichNode(int v) const {
		ostringstream s;
		s << "\"";
		s << this->NodeAsString(v); // (*strings_wrapRO)[nodesRO->get<idT>().find(v )->string_h];
		s << "\"";
		return s.str();
}

struct fileDoesNotExistException {};
struct invalidGraphFileFormatException {};

/*
 * strings. Every string, including but not limited to node names, to be stored in a mmaped file
 */


class ModifiableStringArray : public strings :: StringArray {
public:
	virtual size_t size() const = 0;
	virtual strings :: StrH insert(const char * s) = 0;
        virtual strings :: StrH insert(std :: string &s) {
			return insert(s.c_str());
	}
	virtual const char * operator[] (strings :: StrH sh) const = 0;
	virtual strings :: StrH StringToStringId(const char *s) const = 0;
};
class ModifiableStringArrayInPlainMemory : public ModifiableStringArray {
private:
	std :: vector<std :: string> the_strings;
	boost :: unordered_map<std :: string, int> string2id; // TODO: map by const char *
public:
	explicit ModifiableStringArrayInPlainMemory() {
			this->insert("");
	}
	virtual size_t size() const {
		return the_strings.size();
	}
	virtual strings :: StrH insert(const char * s) {
		std :: string str(s);
		if(string2id.count(str)==0) {
			int nextID = the_strings.size();
			string2id[str] = nextID;
			assert((int)string2id.size() == nextID+1);
			the_strings.push_back(str);
			return strings :: StrH(nextID);
		} else {
			return strings :: StrH(string2id.at(str));
		}
	}
	virtual const char * operator[] (strings :: StrH sh) const {
		return the_strings.at(sh.get_underlying_id()).c_str();
	}
	virtual strings :: StrH StringToStringId(const char *s) const {
		std :: string str(s);
		return strings :: StrH(string2id.at(str));
	}
};

/*
 * nodes and relationships
 * neighbours. For each node, a set of the ids of its relationships.
 */
struct nodeWithName {
	int         id;
	strings :: StrH  string_h;
	nodeWithName( int _id
	        , strings :: StrH _string_h
	        )
	   : id(_id), string_h(_string_h)
	{}
};
typedef bmi :: multi_index_container<
		nodeWithName,
		bmi :: indexed_by<
			bmi :: hashed_unique  <bmi :: tag<idT>,  BOOST_MULTI_INDEX_MEMBER(nodeWithName,int,id)>,
			bmi :: hashed_unique  <bmi :: tag<nameT>,BOOST_MULTI_INDEX_MEMBER(nodeWithName,strings :: StrH,string_h) ,strings :: StrH :: hasher>
		>
	> nodeWithName_Set;


relationship :: relationship( int id_ , const std :: pair<int,int> &nodes_) : relId(id_), nodeIds(nodes_) {
		assert(nodes_ . first <= nodes_ . second); // make sure each relationship, including self-loops, appears only once
}



class DumbGraphReadableTemplate : public ReadableShmGraphTemplate {
protected:
	const shmGraphRaw :: nodeWithName_Set *nodesRO;
	const shmGraphRaw :: relationship_set *relationshipsRO;
	const shmGraphRaw :: neighbours_to_relationships_map *neighbouring_relationshipsRO;
	std :: auto_ptr<const strings :: StringArray> strings_wrapRO;
	const boost :: unordered_set<int> *empty_set_for_neighboursRO;
public:
	virtual int numNodes() const;
	virtual int numRels()  const;
	virtual const boost :: unordered_set<int> & myRels(int n) const;
	virtual pair<const char*, const char*> EndPointsAsStrings(int relId) const;
	virtual const char * NodeAsString(int v) const;
	virtual int StringToNodeId(const char *s) const;
	virtual const std :: pair<int, int> & EndPoints(int relId) const;
	virtual bool are_connected(int v1, int v2) const;
};
/* virtual */ int DumbGraphReadableTemplate :: numNodes() const { return nodesRO->size(); }
/* virtual */ int DumbGraphReadableTemplate :: numRels()  const { return relationshipsRO->size(); }
/* virtual */ const boost :: unordered_set<int> & DumbGraphReadableTemplate :: myRels(int n) const {
		shmGraphRaw :: neighbours_to_relationships_map :: const_iterator i = neighbouring_relationshipsRO->find(n);
		if(i == neighbouring_relationshipsRO->end()) { // it's possible to add a node without any corresponding edges. Must check for this.
			return *empty_set_for_neighboursRO;
		}
		else
			return i->second;
	}
/* virtual */ pair<const char*, const char*> DumbGraphReadableTemplate :: EndPointsAsStrings(int relId) const {
		assert(relId >=0 && relId < this->numRels() );
		const pair<int,int> &endpoints = relationshipsRO->get<idT>().find(relId)->nodeIds;
		const char *l = this->NodeAsString(endpoints.first); // (*strings_wrap)[nodes->get<idT>().find(endpoints.first )->string_h];
		const char *r = this->NodeAsString(endpoints.second);// (*strings_wrap)[nodes->get<idT>().find(endpoints.second)->string_h];
		return make_pair(l,r);
	}
/* virtual */ const char * DumbGraphReadableTemplate :: NodeAsString(int v) const {
		return (*strings_wrapRO)[nodesRO->get<idT>().find(v )->string_h];
	}
/* virtual */ int DumbGraphReadableTemplate :: StringToNodeId(const char *s) const {
		strings :: StrH string_handle = strings_wrapRO->StringToStringId(s);
		// assert (0==strcmp(s , (*strings_wrapRO)[string_handle]));
		shmGraphRaw :: nodeWithName_Set :: index_iterator<nameT>:: type i = nodesRO->get<nameT>().find(string_handle );
		assert(i != nodesRO->get<nameT>().end());
		// assert (0==strcmp(s , (*strings_wrapRO)[i->string_h]));
		return i->id;
	}
/* virtual */ const std :: pair<int, int> & DumbGraphReadableTemplate :: EndPoints(int relId) const {
		assert(relId >=0 && relId < this->numRels() );
		return relationshipsRO->get<idT>().find(relId)->nodeIds;
	}
/* virtual */ bool DumbGraphReadableTemplate :: are_connected(int v1, int v2) const {
		if(v1 > v2)
			swap(v1,v2);
		assert (v1 <= v2);
		return (relationshipsRO->get<nodeIdsT>().end() != relationshipsRO->get<nodeIdsT>().find(make_pair(v1,v2)) );
	}

class DumbGraphRaw : public DumbGraphReadableTemplate {

	const boost :: unordered_set<int> empty_set_for_neighbours;

	shmGraphRaw :: nodeWithName_Set *nodes;
	shmGraphRaw :: relationship_set *relationships;
	shmGraphRaw :: neighbours_to_relationships_map *neighbouring_relationships;
public:
	ModifiableStringArray *strings_wrap;
public:
	virtual ~DumbGraphRaw() {
		assert(this->strings_wrap == this->strings_wrapRO.get()); // delete strings_wrap; // DO NOT delete here, the auto_ptr already has it
		// delete nodes; delete relationships; delete neighbouring_relationships; // seems like you can't/shouldn't delete objects like this
	}
	explicit DumbGraphRaw();
	int insertNode(strings :: StrH node_name) {
		shmGraphRaw :: nodeWithName_Set :: index_iterator<nameT>:: type i = nodes->get<nameT>().find(node_name);
		if(i == nodes->get<nameT>().end()) {
			int proposedNewId = nodes->size();
			std :: pair<shmGraphRaw :: nodeWithName_Set :: iterator, bool> insertionResult = nodes->insert(nodeWithName(proposedNewId, node_name));
			assert(proposedNewId == insertionResult.first->id        );
			assert(node_name     == insertionResult.first->string_h  );
			assert(insertionResult.second);
			return insertionResult.first->id ;
		} else {
			return i->id ;
		}
	}
	int insertRel(pair<int,int> p) {
		if(p.first > p.second)
			swap(p.first, p.second);
		assert(p.first <= p.second);
		shmGraphRaw :: relationship_set :: index_iterator<nodeIdsT>:: type i = relationships->get<nodeIdsT>().find(p);
		if(i == relationships->get<nodeIdsT>().end()) {
			int relId = relationships->size();
			std :: pair<shmGraphRaw :: relationship_set :: iterator, bool> insertionResult = relationships->insert(relationship(relId, p));
			assert(relId == insertionResult.first->relId        );
			assert(p.first       == insertionResult.first->nodeIds.first  );
			assert(p.second      == insertionResult.first->nodeIds.second  );
			assert(insertionResult.second);
			{ // add this to the neighbouring relationships object too
				shmGraphRaw :: neighbours_to_relationships_map :: iterator i;
				{
					i = neighbouring_relationships->find(p.first);
					if(i == neighbouring_relationships->end())
						i = neighbouring_relationships->insert( make_pair(p.first, empty_set_for_neighbours) ) .first ;
					i->second.insert(relId);
				}
				{
					i = neighbouring_relationships->find(p.second);
					if(i == neighbouring_relationships->end())
						i = neighbouring_relationships->insert( make_pair(p.second, empty_set_for_neighbours) ) .first ;
					i->second.insert(relId);
				}
			}
			return relId;
		} else {
			// no need to update the neighbouring_relationships object here. The relationship was already entered
			return i->relId ;
		}
	}
};
DumbGraphRaw :: DumbGraphRaw()
{
		nodes         = new shmGraphRaw :: nodeWithName_Set();
		relationships = new shmGraphRaw :: relationship_set();
	 	neighbouring_relationships
		              = new shmGraphRaw :: neighbours_to_relationships_map();
		strings_wrap = new ModifiableStringArrayInPlainMemory();
		this->nodesRO = nodes;
		this->relationshipsRO = relationships;
		this->neighbouring_relationshipsRO = neighbouring_relationships;
		this->strings_wrapRO.reset(strings_wrap);
		this->empty_set_for_neighboursRO = &empty_set_for_neighbours;
}



ReadableShmGraphTemplate * loadEdgeList(const char *graphTextFileName, const bool selfloops_allowed, graph :: weights :: EdgeDetailsInterface *edge_details) {
	assert(graphTextFileName);
	DumbGraphRaw *nodes_and_rels_wrap = new DumbGraphRaw();

	assert(nodes_and_rels_wrap);
	nodes_and_rels_wrap->hasASelfLoop = false; // this will be changed if/when a self loop is found

	ifstream graphTextStream(graphTextFileName); 
	if(!graphTextStream.is_open())
		throw fileDoesNotExistException();
	forEach(const std :: string &s, amd :: rangeOverStream(graphTextStream)) {
			istringstream oneLine(s);
			amd :: RangeOverStream fields(oneLine, " \t");
			!fields.empty() || ({ throw invalidGraphFileFormatException(); 0; });
			std :: string l = fields.front();
			fields.popFront();
			!fields.empty() || ({ throw invalidGraphFileFormatException(); 0; });
			std :: string r = fields.front();
			
			fields.popFront();
			std :: string weight = fields.front();
			if(weight.length()>0 && *weight.rbegin() == '\r') {
				weight.resize(weight.length()-1);
			}

			pair<int,int> edgeAsIds = make_pair(
					nodes_and_rels_wrap->insertNode(nodes_and_rels_wrap->strings_wrap->insert(l))
					,nodes_and_rels_wrap->insertNode(nodes_and_rels_wrap->strings_wrap->insert(r))
				);
			if(edgeAsIds.first == edgeAsIds.second) {
				nodes_and_rels_wrap->hasASelfLoop = true;
				if(! selfloops_allowed)
					throw shmGraphRaw :: SelfLoopsNotSupported();
			}
			int relId = nodes_and_rels_wrap->insertRel(edgeAsIds);

			edge_details->new_rel(relId, edgeAsIds, weight);

	}
	return nodes_and_rels_wrap;
}

} // namespace shmGraphRaw {
