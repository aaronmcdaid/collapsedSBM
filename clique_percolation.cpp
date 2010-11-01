/*
 * Copyright Aaron McDaid (2010). aaronmcdaid@gmail.com
 * Available under GNU Affero General Public License, version 3
 * (I should put a proper disclaimer and notice on this)
 */
#include <iostream>
#include <fstream>
#include <iterator>
#include <list>
#include <map>
#include <set>
#include <algorithm>
#include "Range.hpp"
#include "clique_percolation.hpp"

#include <sys/stat.h>
#include <cerrno>

using namespace std;
using namespace cliques;
using namespace amd;

/* TODO
 * Foreach for Containers
 * Something to speed up the scan across all communities for a parent. Presumably should keep the lists ordered.
 * delete the empty communities after the merging
 * Put the call to lower level built into merge
 * Non-maximal overlaps?
 */

/* Storage for the cliques */

typedef	vector<V> Clique;
typedef	unsigned int uint;
struct Community;
typedef map<unsigned int, set<Community*> > MapType;
vector< pair<double,bool> > option_thresholds; // this'll be set by an external module, like acp.cpp

struct CliqueStore : public CliqueSink {
	vector<Clique > all_cliques;
	CliqueStore(const SimpleIntGraph &_g, const string& fileName) : CliqueSink(_g, fileName) {}
	void operator () (Clique Compsub) { // We MUST sort this, the clique percolation relies on sorted cliques to help with checking of overlaps.
		sort(Compsub.begin(), Compsub.end());
		if(Compsub.size() >= 3) {
			CliqueSink::operator() (Compsub); // Write them out to the cliques file.
			all_cliques.push_back(Compsub);
		}
	}
};

string toString(const SimpleIntGraph &g, const Clique &c) {
	ContainerRange< const Clique > r (c);
	ostringstream s;
	bool firstItem = true;
	Foreach (V v, r) {
		if(!firstItem) s << ",";
		s << g->NodeAsString(v);
		firstItem=false;
	}
	return s.str();
}

int64 hash_false_positives = 0;
int64 hash_true_positives = 0;
int64 hash_true_negative = 0;
int64 hash_1_true_negatives = 0;
int64 hash_1_false_positives = 0;
int64 hash_1_true_positives = 0;
map < pair<uint, uint> , pair<uint, uint> > overlapSizesAndResults;

/* Storing Communities */
typedef pair<uint64,uint64> HashMembersType;
HashMembersType & operator |= (HashMembersType &l, HashMembersType &r) {
	l.first  |= r.first;
	l.second |= r.second;
	return l;
}

// http://www.concentric.net/~Ttwang/tech/inthash.htm
uint64_t hash64shift(uint64_t key)
{
  key = (~key) + (key << 21); // key = (key << 21) - key - 1;
  key = key ^ (key >> 24);
  key = (key + (key << 3)) + (key << 8); // key * 265
  key = key ^ (key >> 14);
  key = (key + (key << 2)) + (key << 4); // key * 21
  key = key ^ (key >> 28);
  key = key + (key << 31);
  return key;
}

void hashClique (const Clique *clique, HashMembersType &hash) {
		ForeachContainer(V v , *clique) {
			HashMembersType x;

			unsigned int h1 = ( (v * 2654435761U)      ) % (sizeof(HashMembersType::first_type)*8);
			x.first = 1;
			x.first  <<= h1;

			// x.first = bits(v);
			// assert(bitCount(x.first)==N);

			unsigned int h2 = (unsigned int) (   ( hash64shift(v) ) % (sizeof(HashMembersType::first_type)*8)  );
			x.second = 1;
			x.second <<= h2;

			//x.second = bits(hash64shift(v));
			//assert(bitCount(x.second)==N);

			hash.first |= x.first;
			hash.second |= x.second;
		}
}

class Community {
public:
	typedef ::HashMembersType HashMembersType;
private:
	unsigned int k; // cliques will be (at least) this size, and will have (at least) K-1 overlap with all others in this community, or its parent.
	HashMembersType hash_members; // Each community with have a hash of cliquesAtThisLevel
public:
	//const vector<Clique> & CliquesAtThisLevel() const { return cliquesAtThisLevel; } // allow read-only access only
	const vector<const Clique *> & CliquesAtThisLevel;
	const HashMembersType & Hash_members; // Each community with have a hash of cliquesAtThisLevel
private:
	vector<const Clique *> cliquesAtThisLevel; // .. so these must have EXACTLY overlap k-1 with cliques here in this comm. // TODO: Base optimization on this? No overlap with first clique, means no overlap is possible with the one it overlaps with
public:
#define hashes1
#ifdef hashes1
	vector<HashMembersType> hashes;
#endif
	vector<Community*> parents;   // an overlap >=k with anything in this comm implies it should be in one of the parents
	Community *child; // Will be non-NULL, except if k==3 of course.

	Community(unsigned int K) : k(K), hash_members(make_pair(0,0)), CliquesAtThisLevel(cliquesAtThisLevel), Hash_members(hash_members), child(NULL) {
		if(k>3) {
			child = new Community(k-1);
			child->parents.push_back(this);
		}
	}
	unsigned int K() const { return k; }
	void newCliqueAtThisLevel(const Clique *clique, const HashMembersType &hash) {
		cliquesAtThisLevel.push_back(clique);
#ifdef hashes1
		hashes.push_back(hash);
#endif
		hash_members.first  |= hash.first;
		hash_members.second |= hash.second;
	}
	void newCliqueAtThisLevel(const Clique *clique) {
		HashMembersType hash(0,0); // The total hash for this one clique
		hashClique(clique, hash);

		newCliqueAtThisLevel(clique, hash);
	}
	void merge(Community *pull_this_in, MapType &commsForEachK) {
		if(this == pull_this_in)
			return;
		assert(K() == pull_this_in->K());
		ForeachContainer(Community *p, pull_this_in->parents) { // The parents of pull_this_in are to have a new child
			assert(p->child == pull_this_in);
			p->child = this; // It was equal to pull_this_in
		}
		parents.insert(parents.end(), pull_this_in->parents.begin(), pull_this_in->parents.end());

		cliquesAtThisLevel.insert(cliquesAtThisLevel.end(), pull_this_in->cliquesAtThisLevel.begin(), pull_this_in->cliquesAtThisLevel.end());
#ifdef hashes1
		hashes.insert(hashes.end(), pull_this_in->hashes.begin(), pull_this_in->hashes.end());
#endif

		pull_this_in->parents.clear();
		pull_this_in->cliquesAtThisLevel.clear();
#ifdef hashes1
		pull_this_in->hashes.clear();
#endif
		hash_members |= pull_this_in->hash_members;
		pull_this_in->hash_members.first = 0;
		pull_this_in->hash_members.second = 0;
		// This'll leave some empty communities in memory. We'll just delete them later when they're removed from the map.

		if(k>3) { // TODO: Double check this section. I'm not sure I like it!
			Community * child_of_pull_this_in = pull_this_in->child;
			pull_this_in->child = NULL;
			for(unsigned int i = 0; i < child_of_pull_this_in->parents.size(); i++) {
				if(child_of_pull_this_in->parents[i] == pull_this_in) {
					child_of_pull_this_in->parents.erase(child_of_pull_this_in->parents.begin()+i);
					// cout << "REMOVED" << endl;
					break; // Gotta break out now regardless, the iterators in the Range are invalid.
				}
			}
			if(child != child_of_pull_this_in)
				child->merge(child_of_pull_this_in, commsForEachK);
		}

		assert(pull_this_in->CliquesAtThisLevel.empty() && pull_this_in->parents.empty());
		assert(pull_this_in->child==NULL);
		commsForEachK[k].erase(pull_this_in); // we could delete it earlier, in the middle of the merge operation, but I think it's a little tidier here; It's more obvious that it's safe to delete.
		delete(pull_this_in);
	}
};




bool emptyComm(const Community *c) { // TODO: make this a method of Community
	if(!c->CliquesAtThisLevel.empty()) return false;
	ContainerRange< const vector<Community*> > parentsRange(c->parents);   // an overlap >=k with anything in this comm implies it should be in one of the parents
	Foreach(const Community *parentComm, parentsRange) {
		if(!emptyComm(parentComm)) return false;
	}
	return true;
}

void dumpComm(const SimpleIntGraph &g, const Community * c) {
		for(int i=c->K(); i; --i) { cout << ' '; } cout << c->K();
		if(emptyComm(c)) cout << "~";
		cout << endl;

		//ContainerRange< const vector < const Clique *> > cliquesInThisComm(c->CliquesAtThisLevel); //Foreach(const Clique *cl, cliquesInThisComm)
		ForeachContainer(const Clique *cl, c->CliquesAtThisLevel)
		{
			cout << '\t' << toString(g, *cl) << endl;
		}
		ContainerRange< const vector<Community*> > parentsRange(c->parents);   // an overlap >=k with anything in this comm implies it should be in one of the parents
		Foreach(const Community *parentComm, parentsRange) {
			dumpComm(g, parentComm);
		}
		// for(int i=c->K; i; --i) { cout << '.'; } cout << c->K << endl;
}

unsigned int sizeOfBiggestOverlap(const Clique &clique, const Community &comm, const Community::HashMembersType &hash_of_this_clique, bool recursive_call = false);
void extractAllNodes(set<V> &thisComm, const Community *c);

bool moreBySize(const vector<V> & l, const vector<V> & r) { return l.size() > r.size(); }
bool checkHashForOverlap(const Community::HashMembersType &hash_of_this_clique, const Community *c) { // merge this into sizeOfBiggestOverlap
	if( (hash_of_this_clique.first & c->Hash_members.first)  && (hash_of_this_clique.second & c->Hash_members.second) )
		return true; // There might well be at least one overlap
	ForeachContainer(const Community * parent, c->parents)
		if(checkHashForOverlap(hash_of_this_clique, parent))
			return true;
	return false;
}

void cliquePercolation(const SimpleIntGraph &g_, const string &outputDirectory, unsigned int minimumSize) {
	create_directory(outputDirectory);

	string cliquesFileName(outputDirectory + "/cliques");
	CliqueStore cliquesOut(g_, cliquesFileName);

	findCliques<CliqueStore>(g_, cliquesOut, minimumSize);

	cout << cliquesOut.all_cliques.size() << " cliques found" << endl;
	if(cliquesOut.all_cliques.size() == 0) {
		cout << "There were no cliques of the minimum size found" << endl;
		return;
	}

	sort(cliquesOut.all_cliques.begin(), cliquesOut.all_cliques.end(), moreBySize );
	
	const vector< vector<V> > &all_cliques = cliquesOut.all_cliques;


	int sizeOfBiggestClique = all_cliques.front().size();
	PP(sizeOfBiggestClique);

	// Create the map from K to the-communities-for-K
	typedef map<unsigned int, set<Community*> > MapType;
	MapType commsForEachK;

	
	ContainerRange< const vector< vector<V> > > all_cliquesR(all_cliques);
	uint current_clique_size = 0;
	Foreach (const Clique &clique, all_cliquesR) {
		// PP(toString(g, clique));
		{
			uint cliqueNo = &clique - &(all_cliques.front());
			if(cliqueNo % 10000 == 0) {
				// extern StopWatch beginLoading;
				ostringstream s;
				s << setw(6) << right << fixed << setprecision(2) << (double) cliqueNo / all_cliques.size() * 100;
				// beginLoading.laptime(s.str() + "% of cliques");
			}
		}
		const unsigned int k = clique.size();
		if(current_clique_size!=k) {
			current_clique_size = k;
			PP(current_clique_size);
		}
		assert(k>=3);
		Community * newComm = new Community(k);
		newComm->newCliqueAtThisLevel(&clique);
		const Community::HashMembersType hash_of_this_clique = newComm->Hash_members;
		{ // insert this community (and its children) into our map
			Community *c = newComm;
			for(unsigned int j=k; j>=3; j--) {
				assert(c);
				assert(c->K() == j);
				commsForEachK[j].insert(c); // Now every Community will have a pointer from commsForEachK as well as the parent/child links in the main structure.
				c = c->child;
				if(j==3) assert(c==NULL);
			}
		}
		{ // now, the difficult bit. Checking for overlaps against all existing communities and merging if necessary.
			// invariant: newComm might move, but will always point to the k-community containing the currently-being-inserted clique.
			// invariant: we know that we have no cliques smaller than k, but there may be interesting communities for less-than-k.
			// We just need to check all the existing k+ communities for overlaps against the current clique. For each one, the size of the largest overlap is sufficient.

			// cout << "new clique is: " << toString(g, clique) << endl;
/*
			bool followThisClique = g.name_of_one_node(clique.front()) == 22945;
			if(followThisClique) {
				cout << "missing clique: " << toString(g, clique) << endl;
				const vector<Community*> &k4 = commsForEachK[4];
				dumpComm(g, k4.back());
			}
*/
			MapType::mapped_type &commsForEachKk = commsForEachK[k]; // efficiency, to avoid the map extracting [k] repeatedly
			for(MapType::mapped_type::iterator it = commsForEachKk.begin(); it!=commsForEachKk.end(); ++it) {
				// we sometimes erase+delete newComm, but the iterator it is always valid.
				Community *c = *it;
				if(c == newComm)
					continue;
				assert(!emptyComm(c));
				//cout << "checking for overlap in" << endl;
				unsigned int overlap;
				//if(!checkHashForOverlap(hash_of_this_clique, c)) {
					//overlap = 0;
				//} else {
					overlap = sizeOfBiggestOverlap(clique, *c, hash_of_this_clique);
				//}
				if(overlap>=2) {
					// PP((hash_of_this_clique & c->Hash_members));
				}
				//dumpComm(c);
				//PP(overlap);
				// overlap must be < k. If overlap+1 == k, then merge, otherwise move down a level.
/*
				if(overlap>=2 && followThisClique) {
					cout << "our special clique overlaps with an existing community" << endl;
					dumpComm(g, c);
				}
*/
				if(overlap >=2 && overlap+1 == k) {
					c->merge(newComm, commsForEachK);
					// The below is the only line that changes newComm
					newComm = c; // newComm is empty, ignore it for now and clean it up later.
				}
				if(overlap >=2 && overlap+1 <  k) { // must be 2<>4 for the authors400 data
					Community *c2 = c;
					Community *new2 = newComm;
					int skipDown = k - overlap - 1;
					while(skipDown-- > 0) {
						c2 = c2->child;
						new2 = new2->child;
					}
					c2->merge(new2, commsForEachK);
				}
			}
/*
			if(followThisClique) {
				const vector<Community*> &k4 = commsForEachK[4]; dumpComm(g, k4.back());
				const vector<Community*> &k3 = commsForEachK[3]; dumpComm(g, k3.back());
			}
*/
		}
	}
	
	if(1) { // dump the final communities
		PP(commsForEachK.size());
		IteratorRange<MapType::iterator> allCommsRange(commsForEachK.begin(), commsForEachK.end());
		Foreach(const MapType::value_type &v, allCommsRange) {
			ostringstream cliquesFileName;
			cliquesFileName << outputDirectory << "/comms" << v.first;
			ofstream commsKFile(cliquesFileName.str().c_str());
			unsigned int countComms=0;
			ForeachContainer(const Community *c, v.second) {
				if(!emptyComm(c)) {
					countComms++;
					set<V> thisComm;
					extractAllNodes(thisComm, c);
					ForeachContainer(V v, thisComm) {
						commsKFile << g_->NodeAsString(v) << " ";
					}
					commsKFile << endl;
				}
			}
			cout << "#comms for k=" << v.first << ": " << countComms << endl;
		}
/*
		PP(commsForEachK.size());
		IteratorRange<MapType::iterator> allCommsRange2(commsForEachK.begin(), commsForEachK.end());
		Foreach(const MapType::value_type &v, allCommsRange2) {
			cout << "#comms for k=" << v.first << ": " << v.second.size() << endl;
		}
*/
	}


	// map < pair<uint, uint> , pair<uint, uint> > overlapSizesAndResults;
	typedef map < pair<uint, uint> , pair<uint, uint> > ::value_type value_type;
	cout << endl;
	ForeachContainer(const value_type &v, overlapSizesAndResults) {
		cout << "" << setw(4) << v.first.first << ',' << setw(4) << v.first.second << ": " << setw(10) << v.second.first << ',' << setw(10) << v.second.second << endl;
	}
	PPdec(hash_false_positives);
	PPdec(hash_true_positives);
	PPdec(hash_true_negative);
#ifdef hashes1
	PPdec(hash_1_false_positives);
	PPdec(hash_1_true_positives);
	PPdec(hash_1_true_negatives);
#endif
	PP(sizeof(Community::HashMembersType));
}

size_t max(size_t l, size_t r) { return l<r ? r : l; }

/*
struct count_inserter {
This isn't working right. There should be a nice way to count the size of an intersection, ignoring the contents.
	unsigned int count;
	V v;
	count_inserter() : count(0) , v(0) {}
	V & operator * () { return v; }
	void operator ++ () { ++count; }
};
*/ 

unsigned int sizeOfBiggestOverlap(const Clique &clique, const Community &comm, const Community::HashMembersType &hash_of_this_clique, bool recursive_call /*= false*/) { // TODO: This is a clear candidate for optimization, with the hash function thing for example.
	recursive_call = false; // Cheat, look at each level as a separate attempt to match based on hash
	const unsigned int k = clique.size();
	unsigned int overlap = 0;
//#define gatherOverlapPerfStats
#ifdef gatherOverlapPerfStats
	pair<uint, uint> & counterForTheseSizes = overlapSizesAndResults[make_pair(comm.CliquesAtThisLevel.size(), comm.CliquesAtThisLevel.size())];
#endif
	if(!recursive_call) {
#ifdef gatherOverlapPerfStats
		counterForTheseSizes.first++;
#endif
		//overlapSizesAndResults[make_pair(comm.CliquesAtThisLevel.size(), comm.hashes.size())].first ++;
	}
	if( (hash_of_this_clique.second & comm.Hash_members.second) && (hash_of_this_clique.first & comm.Hash_members.first) ) {
		int i=0;
		ForeachContainer(const Clique *cl, comm.CliquesAtThisLevel) {
#ifdef hashes1
			if( (hash_of_this_clique.second & comm.hashes[i].second) && (hash_of_this_clique.first & comm.hashes[i].first) )
				{} // There might be an overlap
			else
			{
				hash_1_true_negatives++;
				i++;
				continue; // definitely no overlap
			}
#endif
			i++;
			uint overlap_ = 0;
			{
				Clique::const_iterator b1 = clique.begin();
				Clique::const_iterator b2 = cl    ->begin();
				Clique::const_iterator e1 = clique.end();
				Clique::const_iterator e2 = cl    ->end();
				while (b1!=e1 && b2!=e2)
				{
					if (*b1<*b2) ++b1;
					else if (*b2<*b1) ++b2;
					else { ++overlap_; b1++; b2++; }
					if(overlap_ + 1 == k)
						break;
				}
			}
			// Clique o;
			// set_intersection(clique.begin(), clique.end(), cl.begin(), cl.end(), back_inserter(o));
			// assert(o.size()==overlap_);
			overlap = max(overlap, overlap_);
			if(overlap_==0) {
				hash_1_false_positives++;
				//PP(bitCount(hash_of_this_clique.first & comm.hashes[i].first));
				//PP(bitCount(comm.hashes[i].first));
				//PP(cl.size());
				//PP(bitCount(hash_of_this_clique.first));
				//PP(clique.size());
			} else {
				hash_1_true_positives++;
			}
			assert(overlap!=k);
			if(overlap+1 == k) {
				++hash_true_positives;
				return overlap;
			}
		}
		if(overlap==0) {
			++hash_false_positives;
			if(!recursive_call) {
#ifdef gatherOverlapPerfStats
				counterForTheseSizes.second++;
#endif
			}
			// PPhex(hash_of_this_clique);
			// printf("  "); PPhex(comm.Hash_members);
			// PP(comm.CliquesAtThisLevel.size());
		} else {
			++hash_true_positives;
		}
	} else {
		if(comm.CliquesAtThisLevel.empty())
			assert(comm.Hash_members.first==0 && comm.Hash_members.second==0);
		else
			++hash_true_negative;
	}
	ForeachContainer(const Community *c, comm.parents) {
		overlap = max(overlap, sizeOfBiggestOverlap(clique, *c, hash_of_this_clique, true));
		assert(overlap!=k);
		if(overlap+1 == k) return overlap;
	}
	return overlap;
}

void extractAllNodes(set<V> &thisComm, const Community *c) {
	ForeachContainer(const Clique *cl, c->CliquesAtThisLevel)
		thisComm.insert(cl->begin(), cl->end());
	ForeachContainer(const Community *c2, c->parents)
		extractAllNodes(thisComm, c2);
}

/*
authors150
k=7 correct                                   <StopWatch> 13.032 (+13.032) Clique Percolation
authors100
k=7 correct                                   <StopWatch> 109.931 (+109.931) Clique Percolation

after basic hash usage...
           7.9
						82


too many false negatives
authors 150 8.4s
hash_false_positives:23,926,559
hash_true_positives:     33,963
hash_true_negative: 134,355,939

// 16005532  >>5
// 15939520 >> 8
//  4459857 >> 8 (with less messing with v)
//  4472956  // hash64shift

using 128 bits of hashing...
	authors100   73.7s (versus 43 for CFinder) false_pos  27m / 846m
	authors100   84.7s (versus 43 for CFinder) false_pos 136m / 737m // Knuth-only

*/

struct CliquesVector /*: public CliqueSink*/ {
	vector<Clique > all_cliques;
	CliquesVector() /*: CliqueSink(_g, fileName)*/ {}
	void operator () (Clique Compsub) { // We MUST sort this, the clique percolation relies on sorted cliques to help with checking of overlaps.
		sort(Compsub.begin(), Compsub.end());
		if(Compsub.size() >= 3) {
			all_cliques.push_back(Compsub);
		}
	}
};

struct ConnectedComponents {
	int C; // the number of cliques
	vector<int> component;
	vector<int> next;
	vector<int> prev;
	vector<int> sizes;
	ConnectedComponents() : C(-1) {}
	void setNumCliques(int _C) {
		assert(this->C==-1);
		this->C = _C;
		this->component.resize(C);
		this->next.resize(C);
		this->prev.resize(C);
		this->sizes.resize(C,1);
		assert(this->C>0);
		for(int i=0; i<C; i++) {
			component.at(i) = i;
			next     .at(i) = i;
			prev     .at(i) = i;
		}
	}
	void joinNodesIntoSameComponent(int cl1, int cl2) {
		assert(this->C>0);
		const int comp1 = this->component.at(cl1);
		const int comp2 = this->component.at(cl2);
		{ // this'd be faster if comp2 is smaller
			if(this->sizes.at(comp1) < this->sizes.at(comp2)) {
				this->joinNodesIntoSameComponent(cl2,cl1);
				return;
			}
		}
		assert(comp1 != comp2); // TODO: 
#ifdef checkCompSizes
		int sizeA = 0;
		int sizeB = 0;
		{
			int cl = comp1;
			do {
				assert(this->component.at(cl) == comp1);
				sizeA++;
				// cout << "clique " << cl << " is in component " << this->component.at(cl) << endl;
				assert(cl == this->prev.at(this->next.at(cl)));


				cl = this->next.at(cl);
			} while (cl != comp1);
		}
		{
			int cl = comp2;
			do {
				assert(this->component.at(cl) == comp2);
				sizeB++;
				// cout << "clique " << cl << " is in component " << this->component.at(cl) << endl;
				assert(cl == this->prev.at(this->next.at(cl)));


				cl = this->next.at(cl);
			} while (cl != comp2);
		}
			// PP(sizeA);
			// PP(sizeB);
#endif
		assert(comp1 == this->component.at(comp1));
		assert(comp2 == this->component.at(comp2));
		// abolish comp2, renaming all of its to comp1
		int size2 = 0;
		for(int cl = comp2; this->component.at(cl)==comp2; cl = this->next.at(cl) ) {
			this->component.at(cl) = comp1;
			size2++;
		}
		assert(size2 == this->sizes.at(comp2));
		this->sizes.at(comp1) += size2;
		this->sizes.at(comp2) = 0;
		assert(comp1 == this->component.at(comp2));
		const int comp1SndLast = this->prev.at(comp1);
		const int comp2SndLast = this->prev.at(comp2);
		this->next.at(this->prev.at(comp1)) = comp2;
		this->next.at(this->prev.at(comp2)) = comp1;
		// we modified the nexts in terms of the prevs
		// hence the nexts are correct, but the prevs aren't
		this->prev.at(comp1) = comp2SndLast;
		this->prev.at(comp2) = comp1SndLast;

#ifdef checkCompSizes
		{
			int cl = comp1;
			int size = 0;
			do {
				assert(this->component.at(cl) == comp1);
				size++;
				// cout << "clique " << cl << " is in component " << this->component.at(cl) << endl;
				assert(cl == this->prev.at(this->next.at(cl)));


				cl = this->next.at(cl);
			} while (cl != comp1);
			assert(size == sizeA + sizeB);
		}
#endif
	}
};

void myAdjacentCliques(vector<int> &cliquesIShareANodeWith, const int cliqueID, /*vector<ConnectedComponents> &cpms, */const vector< vector<int> > &nodeToCliquesMap, const vector<Clique> &all_cliques, const SimpleIntGraph &g) {
	const Clique &clique = all_cliques.at(cliqueID);
	forEach(const int v, mk_range(clique)) {
		const size_t split = cliquesIShareANodeWith.size();
		vector<int>::const_iterator i     = nodeToCliquesMap.at(v).begin();
		vector<int>::const_iterator i_end = nodeToCliquesMap.at(v).end();
		i_end = lower_bound(i, i_end, cliqueID);
		for(;i!=i_end;i++){
			const int adjClique = *i;
			assert(adjClique < cliqueID);
			if(adjClique < cliqueID) { // clique[cliqueID] should only be compared against larger cliques
				// if(cpm4.component.at(cliqueID) != cpm4.component.at(adjClique) ) {
					// PP(adjClique);
					cliquesIShareANodeWith.push_back(adjClique);
				// }
			}
		}
		// The cliques before the split, will be sorted, as will those after. But we need to merge them.
		if(split > 0 && split < cliquesIShareANodeWith.size()) {
			vector<int>::iterator splitIter = cliquesIShareANodeWith.begin() + split;
			// printf("%p\n", &*cliquesIShareANodeWith.begin());
			// printf("%p\n", &*splitIter);
			// printf("%p\n", &*cliquesIShareANodeWith.end());
			inplace_merge(cliquesIShareANodeWith.begin(), splitIter, cliquesIShareANodeWith.end());
		}
	}
}
void percolateThis(const int cliqueID, vector<ConnectedComponents> &cpms, vector<ConnectedComponents> &byRelative, const vector< vector<int> > &nodeToCliquesMap, const vector<Clique> &all_cliques, const SimpleIntGraph &g) {
/*
	cout << "                                                                  "; PP(cliqueID);
	cout << "                                                                  "; PP(cpm4.next.at(cliqueID));
	cout << "                                                                  "; PP(cpm4.prev.at(cliqueID));
	cout << "                                                                  "; PP(cpm4.component.at(cliqueID));
	cout << "                                                                  "; PP(cpm4.component.at(142));
*/
	vector<int> cliquesIShareANodeWith;
	// cliquesIShareANodeWith.reserve(1000000);
	// cliquesIShareANodeWith.clear();
	// assert(cliquesIShareANodeWith.size()==0);

	myAdjacentCliques(cliquesIShareANodeWith, cliqueID, /*cpms, */nodeToCliquesMap, all_cliques, g);

	if(cliquesIShareANodeWith.size()==0)
		return;

	int currentClique = -1;
	vector<int>::iterator next = cliquesIShareANodeWith.begin();
	assert( next != cliquesIShareANodeWith.end());
	while(currentClique < (int)all_cliques.size()) {
		int consecutiveLikeThis = 0;
		while(next != cliquesIShareANodeWith.end() && currentClique == *next) {
			++consecutiveLikeThis;
			++next;
		}
		if(currentClique > -1) {
			assert(consecutiveLikeThis < (int)all_cliques.at(cliqueID).size());
			assert(consecutiveLikeThis < (int)all_cliques.at(currentClique).size());
			// cout << consecutiveLikeThis << " instances of clique #" << currentClique << endl;
			if(consecutiveLikeThis >= 2) {
				// PP(consecutiveLikeThis);
				for(int k = 3; k <= consecutiveLikeThis+1; k++) {
					if((int)all_cliques.at(cliqueID).size() >= k) {
						ConnectedComponents &cpmk = cpms.at(k);
						if(cpmk.component.at(cliqueID) != cpmk.component.at(currentClique)) {
							// PP(currentClique);
							// PP(cpmk.component.at(currentClique));
							assert((int)all_cliques.at(cliqueID).size() >= k);
							assert((int)all_cliques.at(currentClique).size() >= k);
							assert((int)all_cliques.at(cliqueID).size() <= (int)all_cliques.at(currentClique).size());
							assert(consecutiveLikeThis >= k-1);
							cpmk.joinNodesIntoSameComponent(cliqueID, currentClique);
							// PP(cpmk.component.at(currentClique));
						}
					}
				}
			}
			{ // no restriction on consecutiveLikeThis. Even a single node in the overlap is enough.
				for(int k = 0; k < (int)option_thresholds.size(); k++) {
					const int s_small = (int)all_cliques.at(cliqueID)     .size();
					const int s_big   = (int)all_cliques.at(currentClique).size();
					assert(consecutiveLikeThis > 0);
					assert(consecutiveLikeThis < s_small);
					assert(consecutiveLikeThis < s_big);
					assert(s_small <= s_big);
					if(
							   ( option_thresholds.at(k).second && consecutiveLikeThis >= s_small * 0.01*option_thresholds.at(k).first)
							|| (!option_thresholds.at(k).second && consecutiveLikeThis >  s_small * 0.01*option_thresholds.at(k).first)
						) { // always comparing smaller to bigger, so we can just apply the relative threshold to the size of clique #cliqueID
						ConnectedComponents &cpmk = byRelative.at(k);
						if(cpmk.component.at(cliqueID) != cpmk.component.at(currentClique)) {
							// PP(currentClique);
							// PP(cpmk.component.at(currentClique));
							// assert((int)all_cliques.at(cliqueID).size() >= k);
							// assert((int)all_cliques.at(currentClique).size() >= k);
							assert((int)all_cliques.at(cliqueID).size() <= (int)all_cliques.at(currentClique).size());
							// assert(consecutiveLikeThis >= k-1);
							cpmk.joinNodesIntoSameComponent(cliqueID, currentClique);
							// PP(cpmk.component.at(currentClique));
						}
					}
				}
			}
		}
		if(next == cliquesIShareANodeWith.end()) {
			// all done.
			currentClique = all_cliques.size();
		} else {
			assert(currentClique < *next);
			currentClique = *next;
			assert(currentClique < (int) all_cliques.size());
		}
	}
}

int printCommsToFile(ofstream &cpm4Results, const ConnectedComponents &one_set_of_comms, const int numCliques, const CliquesVector &cliques, int k, const SimpleIntGraph &g_) {
		int numComps = 0;
		for(int comp=0; comp<numCliques; comp++) {
			if(comp == one_set_of_comms.component.at(comp) && one_set_of_comms.sizes.at(comp)==1 && (int)cliques.all_cliques.at(comp).size() < k) {
				continue; // This clique is too small to be relevant at this level of k. It will not have been merged into any communities.
			}
			if(comp == one_set_of_comms.component.at(comp)) {
				//PP(comp);
				numComps++;
				set<int> nodesInThisCPMComm;
				{
					int cliquesInThisCPMComm = 0;
					int cl = comp;
					do {
						assert(one_set_of_comms.component.at(cl) == comp);
						assert((int)cliques.all_cliques.at(cl).size() >= k);
						cliquesInThisCPMComm++;
						// cout << "clique " << cl << " is in component " << one_set_of_comms.component.at(cl) << endl;
						assert(cl == one_set_of_comms.prev.at(one_set_of_comms.next.at(cl)));
						forEach(const int v, mk_range(cliques.all_cliques.at(cl))) {
							nodesInThisCPMComm.insert(v);
						}
						cl = one_set_of_comms.next.at(cl);
					} while (cl != comp);
					assert(cliquesInThisCPMComm == one_set_of_comms.sizes.at(comp)); // TODO: Shouldn't allow that clique into the community anyway.
				}
				{ // output to the file
					forEach(const int v, mk_range(nodesInThisCPMComm)) {
						cpm4Results << ' ' << atoi(g_->NodeAsString(v));
					}
					cpm4Results << endl;
				}
			} else
				assert(0 == one_set_of_comms.sizes.at(comp));
		}
		cpm4Results.close();
		return numComps;
}

int printCommsAndCliquesToFile(ofstream &cpm4Results, const ConnectedComponents &one_set_of_comms, const int numCliques, const CliquesVector &cliques, int k, const SimpleIntGraph &g_) {
		int numComps = 0;
		for(int comp=0; comp<numCliques; comp++) {
			if(comp == one_set_of_comms.component.at(comp) && one_set_of_comms.sizes.at(comp)==1 && (int)cliques.all_cliques.at(comp).size() < k) {
				continue; // This clique is too small to be relevant at this level of k. It will not have been merged into any communities.
			}
			if(comp == one_set_of_comms.component.at(comp)) {
				//PP(comp);
				numComps++;
				{
					int cliquesInThisCPMComm = 0;
					set<int> cliqueIDsInThisComm;
					{
						int cl = comp;
						do {
							cliqueIDsInThisComm.insert(cl);
							cl = one_set_of_comms.next.at(cl);
						} while (cl != comp);
					}
					forEach(const int cl, mk_range(cliqueIDsInThisComm)) {
						assert(one_set_of_comms.component.at(cl) == comp);
						assert((int)cliques.all_cliques.at(cl).size() >= k);
						cliquesInThisCPMComm++;
						// cout << "clique " << cl << " is in component " << one_set_of_comms.component.at(cl) << endl;
						assert(cl == one_set_of_comms.prev.at(one_set_of_comms.next.at(cl)));
						cpm4Results << " {" << cl << "}";
						set<string> nodeNamesInThisClique;
						forEach(const int v, mk_range(cliques.all_cliques.at(cl))) {
							nodeNamesInThisClique.insert(g_->NodeAsString(v));
						}
						forEach(const string &s, mk_range(nodeNamesInThisClique)) {
							cpm4Results << ' ' << s;
						}
					}
					cpm4Results << endl;
					assert(cliquesInThisCPMComm == one_set_of_comms.sizes.at(comp));
					assert(cliquesInThisCPMComm == (int)cliqueIDsInThisComm.size());
				}
			} else
				assert(0 == one_set_of_comms.sizes.at(comp));
		}
		cpm4Results.close();
		return numComps;
}

void cliquePercolation3(const SimpleIntGraph &g_, const string &outputDirectory, unsigned int minimumSize) {
	assert(minimumSize >= 3);
	create_directory(outputDirectory);

	CliquesVector cliques;

	{ Timer timer (printfstring("find cliques of at least size %d", minimumSize));
		findCliques<CliquesVector>(g_, cliques, minimumSize);
	}
	const int numCliques = cliques.all_cliques.size();
	PP(numCliques);
	if(numCliques == 0) {
		Die("No cliques of the required size!");
	}
	{ Timer timer (printfstring("put the biggest cliques to the front of the list"));
		sort(cliques.all_cliques.begin(), cliques.all_cliques.end(), moreBySize );
	}
	const int maxCliqueSize = cliques.all_cliques.at(0).size();
	PP(maxCliqueSize);
	{
		map<int,int> cliqueSizeFrequencies;
		forEach(const Clique &cl, mk_range(cliques.all_cliques)) {
			cliqueSizeFrequencies[cl.size()]++;
		}
		cout << "clique-size frequencies" << endl;
		forEach(const typeof(pair<int,int>) &freq, mk_range(cliqueSizeFrequencies)) {
			cout << "#" << freq.first << "\t" << freq.second << endl;
		}
	}

	vector< vector<int> > nodeToCliquesMap(g_->numNodes());

	{ Timer timer("building nodeToCliquesMap");
		for(int cliqueID = 0; cliqueID < numCliques; cliqueID++) {
			const Clique &clique = cliques.all_cliques.at(cliqueID);
			forEach(const int v, mk_range(clique)) {
				assert(v < g_->numNodes());
				nodeToCliquesMap.at(v).push_back(cliqueID);
			}
		}
	}
	vector<ConnectedComponents> cpms(1+maxCliqueSize);
	vector<ConnectedComponents> byRelative(option_thresholds.size());
	for(int i=3; i<=maxCliqueSize; i++)
		cpms.at(i).setNumCliques(numCliques);
	for(int i=0; i<(int)option_thresholds.size(); i++)
		byRelative.at(i).setNumCliques(numCliques);
	{ Timer timer("do the clique percolation");
		for(int cliqueID = 0; cliqueID < numCliques; cliqueID++) {
			percolateThis(cliqueID, cpms, byRelative, nodeToCliquesMap, cliques.all_cliques, g_);
		}
	}
	{	Timer timer("print the results");
		for(int k=3; k<=maxCliqueSize; k++) {
			ofstream cpm4Results_((outputDirectory + printfstring("/comm%d_cliques", k)).c_str());
			ofstream cpm4Results((outputDirectory + printfstring("/comm%d", k)).c_str());
			ConnectedComponents &one_set_of_comms = cpms.at(k);
                printCommsAndCliquesToFile(cpm4Results_, one_set_of_comms, numCliques, cliques, k, g_);
			const int numComps = printCommsToFile(cpm4Results, one_set_of_comms, numCliques, cliques, k, g_);
			cout << "k=" << k
				<< "\t#communities " << numComps
				<< endl;
		}
		for(int t=0; t< (int)option_thresholds.size(); t++) {
			assert(t < (int)option_thresholds.size());
			assert(t < (int)byRelative       .size());
			const char * suffix = option_thresholds.at(t).second ? "inc" : "" ;
			ofstream fileForOutput_((outputDirectory + printfstring("/thresh_%g%s_cliques", option_thresholds.at(t).first, suffix)).c_str());
			ofstream fileForOutput((outputDirectory + printfstring("/thresh_%g%s", option_thresholds.at(t).first, suffix)).c_str());
			           printCommsAndCliquesToFile(fileForOutput_,byRelative.at(t), numCliques, cliques, 3, g_);
			const int numComps = printCommsToFile(fileForOutput, byRelative.at(t), numCliques, cliques, 3, g_);
			cout << "threshold=" << option_thresholds.at(t).first << '%' << suffix
				<< "\t#communities " << numComps
				<< endl;
		}
	}
}
