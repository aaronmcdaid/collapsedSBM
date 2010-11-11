#include "clique_percolation3.hpp"
#include "graph_utils.hpp"
#include <queue>
typedef pair< vector<int>::const_iterator, vector<int>::const_iterator > ListItem;

struct OptimisiticHeap {
	static bool cmp (const ListItem &l, const ListItem &r) {
		if(l.first == l.second) {
			if(r.first == r.second)
				return false;
			else
				return true;
		} else
			if(r.first == r.second) {
				return false;
			}
		assert(l.first < l.second);
		assert(r.first < r.second);
		return *l.first > *r.first; // the smaller cliqueIDs (i.e. the larger cliques) are brought to the top of the queue.
	}
	// assume that pushed items should probably be near the top of the heap.
	vector<ListItem> h;
	void push(const ListItem &l) {
		h.push_back(l);
		this->bubbleUp(h.size()-1);
	}
	void bubbleUp(const size_t offset) {
		// the node at offset may be better than it's parent.
		if(offset==0)
			return;
		const size_t parentoffset = (offset-1)/2;
		assert(parentoffset < offset);
		if(cmp(this->h.at(parentoffset), this->h.at(offset))) {
			swap(this->h.at(parentoffset), this->h.at(offset));
			this->bubbleUp(parentoffset);
		}
	}
	void bubbleDown(size_t offset) { // tail-recursion trick going on here, as it spends a lot of time here (recursively)
	while(1) {
		// the node at offset may be worse than it's children.
		const size_t child1 = offset*2+1;
		const size_t child2 = child1+1;
		// assert(offset < child1);
		if(child1 >= h.size())
			return; // there are no children
		if(child2 == h.size()) { // just one child location, check it and then return
			if(cmp(this->h.at(offset), this->h.at(child1))) {
				swap(this->h.at(offset), this->h.at(child1));
			}
			return;
		}
		// assert(child2 < h.size()); // both children are in scope.
		const size_t bestChild = cmp(this->h.at(child1), this->h.at(child2)) ? child2 : child1;
		// assert(*this->h.at(bestChild).first <= *this->h.at(child1).first);
		// assert(*this->h.at(bestChild).first <= *this->h.at(child2).first);
		if(cmp(this->h.at(offset), this->h.at(bestChild))) {
			swap(this->h.at(offset), this->h.at(bestChild));
			offset = bestChild;
			continue;
		}
		break;
	}
	}
	int size() const {
		return this->h.size();
	}
	ListItem & top () {
		return h.front();
	}
	bool empty() const {
		return h.empty();
	}
};

static int64 pushes = 0;

void advanceAsFarAsPossible(ListItem &l, const int newTop, const int cliqueID, const vector<cliques::Clique> &all_cliques, vector<amd::ConnectedComponents> &cpms
		) {
	// const int currentComponent = cpmk.component.at(cliqueID);
			int fastforward = 0;
			do {
				++ fastforward;
				++l.first;
				if(l.first == l.second)
					break;
				const int ahead_clique = *l.first;
				if(ahead_clique < newTop)
					continue;
				const int current_size = all_cliques.at(ahead_clique).size();
				const vector<int> &cpmk = cpms.at(current_size).component;
				if(cpmk.at(ahead_clique) == cpmk.at(cliqueID))
					continue;
				break;
			} while(1);
			// if(fastforward>200) PP(fastforward);
}
// ca-AstroPh
//   35,752,927 pushes
// cit-HepPh
//   1,329,152,018 pushes
//   612,387,598 down to 101 seconds
//   search from mid to smallest.
//     1,351,901,215         178s(133s in cp)
//     632,423,636    126s(81s in cp)
//   my Heap
//     1,685,956,731 139s(95s in cp)
//     1,685,956,731 111s(67s in cp)

static void myAdjacentCliques(const int cliqueID, const vector< vector<int> > &nodeToCliquesMap, const vector<cliques::Clique> &all_cliques, vector<amd::ConnectedComponents> &cpms) {
	// this function is called repeatedly, once for each clique. Starting at the largest clique (smallest cliqueID)
	// this function then compares its clique against all the smaller cliques (i.e. larger cliqueID).
	// TODO: reorder the cliqueID such that the last two lines are not so confusing!
	
	vector<int> overlaps(all_cliques.size());
	OptimisiticHeap qo;

	const cliques::Clique &clique = all_cliques.at(cliqueID);
	// const int clique_size = clique.size();
	// PP2(cliqueID, clique_size);
	forEach(const int v, amd::mk_range(clique)) {
		vector<int>::const_iterator i_mid = upper_bound(nodeToCliquesMap.at(v).begin(),nodeToCliquesMap.at(v).end(),cliqueID); // this seems to take no time. Happy days
		vector<int>::const_iterator i_end = nodeToCliquesMap.at(v).end();
		if(i_mid != i_end) {
			assert(*i_mid > cliqueID);
			//q.push(ListItem(i_mid, i_end));
			qo.push(ListItem(i_mid, i_end));
		}
		while(i_mid != i_end) {
			++ overlaps.at(*i_mid);
			++ i_mid;
		}
	}
	// PP(q.size());
	/*
	PP(qo.size());
	forEach(const ListItem &l, amd::mk_range(qo.h)) {
		PP2(l.second - l.first, *l.first);
	}
	*/
	if(!qo.empty())
	while(1) {
		assert(!qo.empty());
		ListItem &l_top = qo.top();
		if(l_top.first == l_top.second)
			break; // the best list left is empty, hence they're all empty
		assert(l_top.first < l_top.second);
		const int adjClique = *l_top.first;
		int overlap=0;
		do{
			++l_top.first;
			++overlap;
			qo.bubbleDown(0);
			// l_top might no longer point at the same thing, after the bubbleDown.
			++pushes;
		} while(l_top.first != l_top.second && *l_top.first == adjClique);
		assert(overlap == overlaps.at(adjClique));
		int k = overlap + 1;
		while (k >= 3) {
			amd::ConnectedComponents &cpmk = cpms.at(k);
			bool wasAnActualMerging = cpmk.joinNodesIntoSameComponent(cliqueID, adjClique);
			if(!wasAnActualMerging) break; // we know that we always mark all the way down to three. Hence, if these two cliques have already been merged, then we needn't do that again.
			k--;
		}
	}

#if 0
	if(0) while(!q.empty()) {
		const int adjClique = *q.top().first;
		int overlap=0;
		assert(adjClique > cliqueID);
		const int adjClique_size = all_cliques.at(adjClique).size();
		assert(adjClique_size <= clique_size);
		while(!q.empty() && *q.top().first == adjClique) {
			overlap++;
			ListItem l = q.top();
			q.pop();
			assert(l.first != l.second);
			// we can see what's on the top of the heap now. We should keep incrementing until we're at it
			const int newTop = *q.top().first;
			assert(newTop >= adjClique);
			// amd::ConnectedComponents &cpmk = cpms.at(adjClique_size);
			advanceAsFarAsPossible(l, newTop, cliqueID, all_cliques, cpms);
			if(l.first != l.second) {
				q.push(l);
				++pushes;
			}
		}
		assert(overlap > 0);
		// PP2(adjClique, overlap);
		assert(clique.size() >= all_cliques.at(adjClique).size());
		assert(overlap < adjClique_size);
		int k = overlap + 1;
		while (k >= 3) {
			amd::ConnectedComponents &cpmk = cpms.at(k);
			bool wasAnActualMerging = cpmk.joinNodesIntoSameComponent(cliqueID, adjClique);
			if(!wasAnActualMerging) break; // we know that we always mark all the way down to three. Hence, if these two cliques have already been merged, then we needn't do that again.
			k--;
		}
	}
	// assert(q.empty());
#endif
}
#if 0
static void percolateThis(const int cliqueID, vector<amd::ConnectedComponents> &cpms, vector<amd::ConnectedComponents> &byRelative, const vector< vector<int> > &nodeToCliquesMap, const vector<cliques::Clique> &all_cliques, const SimpleIntGraph &g, const vector<pair<double,bool> > &option_thresholds) {
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

	vector<int>::const_iterator cliquesIShareANodeWith_begin = cliquesIShareANodeWith.begin();
	vector<int>::const_iterator cliquesIShareANodeWith_end = cliquesIShareANodeWith.end();

	const bool do_asserts = false;

	int currentClique = -1;
	const int cliqueID_size = (int)all_cliques.at(cliqueID).size();
	vector<int>::const_iterator next = cliquesIShareANodeWith_begin;
	assert( next != cliquesIShareANodeWith_end);
	while(currentClique < (int)all_cliques.size()) {
		int consecutiveLikeThis = 0;
		while(next != cliquesIShareANodeWith_end && currentClique == *next) {
			++consecutiveLikeThis;
			++next;
		}
		if(currentClique > -1) {
			if(do_asserts) assert(consecutiveLikeThis < cliqueID_size);
			if(do_asserts) assert(consecutiveLikeThis < (int)all_cliques.at(currentClique).size());
			// cout << consecutiveLikeThis << " instances of clique #" << currentClique << endl;
			if(consecutiveLikeThis >= 2) {
				// PP(consecutiveLikeThis);
				for(int k = 3; k <= consecutiveLikeThis+1 && k <= cliqueID_size; k++) {
					// if(cliqueID_size >= k)
					{
						// assert(k <= cliqueID_size);
						amd::ConnectedComponents &cpmk = cpms.at(k);
						if(cpmk.component.at(cliqueID) != cpmk.component.at(currentClique)) {
							// PP(currentClique);
							// PP(cpmk.component.at(currentClique));
							if(do_asserts) assert(cliqueID_size >= k);
							if(do_asserts) assert((int)all_cliques.at(currentClique).size() >= k);
							if(do_asserts) assert(cliqueID_size <= (int)all_cliques.at(currentClique).size());
							if(do_asserts) assert(consecutiveLikeThis >= k-1);
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
					if(do_asserts) assert(consecutiveLikeThis > 0);
					if(do_asserts) assert(consecutiveLikeThis < s_small);
					if(do_asserts) assert(consecutiveLikeThis < s_big);
					if(do_asserts) assert(s_small <= s_big);
					if(
							   ( option_thresholds.at(k).second && consecutiveLikeThis >= s_small * 0.01*option_thresholds.at(k).first)
							|| (!option_thresholds.at(k).second && consecutiveLikeThis >  s_small * 0.01*option_thresholds.at(k).first)
						) { // always comparing smaller to bigger, so we can just apply the relative threshold to the size of clique #cliqueID
						amd::ConnectedComponents &cpmk = byRelative.at(k);
						if(cpmk.component.at(cliqueID) != cpmk.component.at(currentClique)) {
							// PP(currentClique);
							// PP(cpmk.component.at(currentClique));
							// assert((cliqueID_size>= k);
							// assert((int)all_cliques.at(currentClique).size() >= k);
							if(do_asserts) assert(cliqueID_size <= (int)all_cliques.at(currentClique).size());
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
#endif

static int printCommsToFile(ofstream &cpm4Results, const amd::ConnectedComponents &one_set_of_comms, const int numCliques, const cliques::CliquesVector &cliques, int k, const SimpleIntGraph &g_) {
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
						forEach(const int v, amd::mk_range(cliques.all_cliques.at(cl))) {
							nodesInThisCPMComm.insert(v);
						}
						cl = one_set_of_comms.next.at(cl);
					} while (cl != comp);
					assert(cliquesInThisCPMComm == one_set_of_comms.sizes.at(comp)); // TODO: Shouldn't allow that clique into the community anyway.
				}
				{ // output to the file
					forEach(const int v, amd::mk_range(nodesInThisCPMComm)) {
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

static int printCommsAndCliquesToFile(ofstream &cpm4Results, const amd::ConnectedComponents &one_set_of_comms, const int numCliques, const cliques::CliquesVector &cliques, int k, const SimpleIntGraph &g_) {
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
					forEach(const int cl, amd::mk_range(cliqueIDsInThisComm)) {
						assert(one_set_of_comms.component.at(cl) == comp);
						assert((int)cliques.all_cliques.at(cl).size() >= k);
						cliquesInThisCPMComm++;
						// cout << "clique " << cl << " is in component " << one_set_of_comms.component.at(cl) << endl;
						assert(cl == one_set_of_comms.prev.at(one_set_of_comms.next.at(cl)));
						cpm4Results << " {" << cl << "}";
						set<string> nodeNamesInThisClique;
						forEach(const int v, amd::mk_range(cliques.all_cliques.at(cl))) {
							nodeNamesInThisClique.insert(g_->NodeAsString(v));
						}
						forEach(const string &s, amd::mk_range(nodeNamesInThisClique)) {
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

void cliquePercolation3(const SimpleIntGraph &g_, const string &outputDirectory, unsigned int minimumSize, const vector< pair<double,bool> > &option_thresholds) {
	assert(minimumSize >= 3);
	amd::create_directory(outputDirectory);

	cliques::CliquesVector cliques;

	{ Timer timer (printfstring("find cliques of at least size %d", minimumSize));
		cliques::findCliques<cliques::CliquesVector>(g_, cliques, minimumSize);
	}
	const int numCliques = cliques.all_cliques.size();
	PP(numCliques);
	if(numCliques == 0) {
		Die("No cliques of the required size!");
	}
	{ Timer timer (printfstring("put the biggest cliques to the front of the list"));
		sort(cliques.all_cliques.begin(), cliques.all_cliques.end(), cliques::moreBySize );
	}
	const int maxCliqueSize = cliques.all_cliques.at(0).size();
	PP(maxCliqueSize);
	{
		map<int,int> cliqueSizeFrequencies;
		forEach(const cliques::Clique &cl, amd::mk_range(cliques.all_cliques)) {
			cliqueSizeFrequencies[cl.size()]++;
		}
		cout << "clique-size frequencies" << endl;
		forEach(const typeof(pair<int,int>) &freq, amd::mk_range(cliqueSizeFrequencies)) {
			cout << "#" << freq.first << "\t" << freq.second << endl;
		}
	}

	vector< vector<int> > nodeToCliquesMap(g_->numNodes());

	{ Timer timer("building nodeToCliquesMap");
		for(int cliqueID = 0; cliqueID < numCliques; cliqueID++) {
			const cliques::Clique &clique = cliques.all_cliques.at(cliqueID);
			forEach(const int v, amd::mk_range(clique)) {
				assert(v < g_->numNodes());
				nodeToCliquesMap.at(v).push_back(cliqueID);
			}
		}
	}
	vector<amd::ConnectedComponents> cpms(1+maxCliqueSize);
	vector<amd::ConnectedComponents> byRelative(option_thresholds.size());
	assert(option_thresholds.size()==0); // option_thresholds won't be supported for a while. Breaking it until the heap is fixed.
	for(int i=3; i<=maxCliqueSize; i++)
		cpms.at(i).setNumCliques(numCliques);
	for(int i=0; i<(int)option_thresholds.size(); i++)
		byRelative.at(i).setNumCliques(numCliques);
	{ Timer timer("do the clique percolation");
		int current_size = -1;
		std::auto_ptr<Timer> timer2(NULL);
		for(int cliqueID = 0; cliqueID < numCliques; cliqueID++) {
			const int cliqueID_size = cliques.all_cliques.at(cliqueID).size();
			if (cliqueID_size != current_size) {
				PP(pushes);
				timer2.reset(new Timer(printfstring("do the %d-cliques ", cliqueID_size)));
				current_size = cliqueID_size;
				PP(cliqueID_size);
			}
			// percolateThis(cliqueID, cpms, byRelative, nodeToCliquesMap, cliques.all_cliques, g_, option_thresholds);
			myAdjacentCliques(cliqueID, nodeToCliquesMap, cliques.all_cliques, cpms);
		}
		PP(pushes);
		timer2.reset(NULL);
	}
	{	Timer timer("print the results");
		for(int k=3; k<=maxCliqueSize; k++) {
			ofstream cpm4Results_((outputDirectory + printfstring("/comm%d_cliques", k)).c_str());
			ofstream cpm4Results((outputDirectory + printfstring("/comms%d", k)).c_str());
			amd::ConnectedComponents &one_set_of_comms = cpms.at(k);
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
