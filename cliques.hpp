#ifndef _CLIQUES_HPP_
#define _CLIQUES_HPP_

#include <list>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sys/stat.h>
#include <cerrno>
#include <stdexcept>

#include "Range.hpp"
#include "aaron_utils.hpp"
#include "shmGraphRaw.hpp"

using namespace std;

typedef int V;
typedef int VertexIDType;
typedef const shmGraphRaw::ReadableShmGraphBase* SimpleIntGraph;

namespace cliques {

struct CliqueFunctionAdaptor {
	virtual void operator() (const std::vector<V> &clique) = 0;
	virtual ~CliqueFunctionAdaptor();
};

template <class T> void findCliques(const SimpleIntGraph &, T &cliquesOut, unsigned int minimumSize);
void cliquesForOneNode(const SimpleIntGraph &g, CliqueFunctionAdaptor &cliquesOut, int minimumSize, V v);
void cliquesToDirectory          (SimpleIntGraph, const string &outputDirectory, unsigned int minimumSize); // You're not allowed to ask for the 2-cliques

template<typename T>
struct CliqueFunctionAdaptor_ : CliqueFunctionAdaptor {
	T & callThis;
	const SimpleIntGraph &g;
	CliqueFunctionAdaptor_(T &c, const SimpleIntGraph &g_) : callThis(c), g(g_) {}
	void operator() (const vector<V> &clique) {
		vector<V> copy = clique;
		sort(copy.begin(), copy.end());
		callThis(copy);
	}
};

template <class T> void findCliques(const SimpleIntGraph &g, T & cliquesOut, unsigned int minimumSize) {
	unless(minimumSize >= 3) throw std::invalid_argument("the minimumSize for findCliques() must be at least 3");

	CliqueFunctionAdaptor_<T> cfa(cliquesOut, g);
	
	for(V v = 0; v < (V) g->numNodes(); v++) {
		// if(v && v % 100 ==0) PP(v);
		cliquesForOneNode(g, cfa, minimumSize, v);
	}
}


struct CliqueSink { // Dump the cliques to a file in the CFinder format
	const SimpleIntGraph &g;
	int n;
	ofstream out;
	CliqueSink(const SimpleIntGraph &_g, const string& fileName) : g(_g), n(0), out(fileName.c_str()) { }
	void operator () (const vector<V> & Compsub) {
		if(Compsub.size() >= 3) {
			out << n << ": ";
			ForeachContainer(V v, Compsub) {
				out << g->NodeAsString(v) << ' ';
			}
			out << endl;
			n++;
		}
	}
};
typedef	vector<V> Clique;
struct CliquesVector /*: public CliqueSink*/ {
	vector<Clique > all_cliques;
	CliquesVector();
	void operator () (Clique Compsub);
};
bool moreBySize(const vector<V> & l, const vector<V> & r);

} // namespace cliques


#endif
