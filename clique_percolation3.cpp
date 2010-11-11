#include "clique_percolation3.hpp"
#include "graph_utils.hpp"
#include <queue>

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

	const cliques::Clique &clique = all_cliques.at(cliqueID);
	forEach(const int v, amd::mk_range(clique)) {
		vector<int>::const_iterator i_mid = upper_bound(nodeToCliquesMap.at(v).begin(),nodeToCliquesMap.at(v).end(),cliqueID); // this seems to take no time. Happy days
		vector<int>::const_iterator i_end = nodeToCliquesMap.at(v).end();
		while(i_mid != i_end) {
			++ overlaps.at(*i_mid);
			++ i_mid;
		}
	}
	for( size_t adjClique = cliqueID+1; adjClique < all_cliques.size(); adjClique++) {
		int overlap = overlaps.at(adjClique);
		int k = overlap + 1;
		while (k >= 3) {
			amd::ConnectedComponents &cpmk = cpms.at(k);
			bool wasAnActualMerging = cpmk.joinNodesIntoSameComponent(cliqueID, adjClique);
			if(!wasAnActualMerging) break; // we know that we always mark all the way down to three. Hence, if these two cliques have already been merged, then we needn't do that again.
			k--;
		}
	}

}

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
				timer2.reset(new Timer(printfstring("do the %d-cliques ", cliqueID_size)));
				current_size = cliqueID_size;
				PP(cliqueID_size);
			}
			// percolateThis(cliqueID, cpms, byRelative, nodeToCliquesMap, cliques.all_cliques, g_, option_thresholds);
			myAdjacentCliques(cliqueID, nodeToCliquesMap, cliques.all_cliques, cpms);
		}
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
