#include "cliques.hpp"
extern vector< pair<double,bool> > option_thresholds; // this'll be set by an external module, like acp.cpp
void cliquePercolation(const SimpleIntGraph &, const string &outputDirectory, unsigned int minimumSize); // You're not allowed to ask for the 2-cliques
