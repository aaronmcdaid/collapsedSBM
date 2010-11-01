#include "clique_percolation3.hpp"
#include "graph_utils.hpp"
void myAdjacentCliques(vector<int> &cliquesIShareANodeWith, const int cliqueID, /*vector<amd::ConnectedComponents> &cpms, */const vector< vector<int> > &nodeToCliquesMap, const vector<cliques::Clique> &all_cliques, const SimpleIntGraph &g) {
	const cliques::Clique &clique = all_cliques.at(cliqueID);
	forEach(const int v, amd::mk_range(clique)) {
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
void percolateThis(const int cliqueID, vector<amd::ConnectedComponents> &cpms, vector<amd::ConnectedComponents> &byRelative, const vector< vector<int> > &nodeToCliquesMap, const vector<cliques::Clique> &all_cliques, const SimpleIntGraph &g) {
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
						amd::ConnectedComponents &cpmk = cpms.at(k);
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
						amd::ConnectedComponents &cpmk = byRelative.at(k);
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

int printCommsToFile(ofstream &cpm4Results, const amd::ConnectedComponents &one_set_of_comms, const int numCliques, const cliques::CliquesVector &cliques, int k, const SimpleIntGraph &g_) {
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

int printCommsAndCliquesToFile(ofstream &cpm4Results, const amd::ConnectedComponents &one_set_of_comms, const int numCliques, const cliques::CliquesVector &cliques, int k, const SimpleIntGraph &g_) {
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

void cliquePercolation3(const SimpleIntGraph &g_, const string &outputDirectory, unsigned int minimumSize) {
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
