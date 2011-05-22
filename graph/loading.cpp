#include"loading.hpp"
#include"../pp.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

using namespace std;

namespace graph {
namespace loading {

struct ModifiableNetwork;
static void read_edge_list_from_file(ModifiableNetwork *network, const string file_name);
typedef pair< pair<string, string> , string> ThreeStrings;
static ThreeStrings parseLine(const string &lineOrig);

struct ModifiableNetwork : public Network { // Network is the read-only interface we want to expose, but this is the derived class that will do the heavy lifting
	std :: auto_ptr< map<string, int> > ordered_node_names;
	ModifiableNetwork(const bool directed, const bool weighted) : Network(directed, weighted) {
	}
	virtual ~ ModifiableNetwork() throw() { }
};

static ThreeStrings parseLine(const string &lineOrig) {
	PP(lineOrig);
	string line(lineOrig); // copy the line. We want to keep the original in order to print error messages
	// line will not have any newlines in it
	// So I'll convert all the delimeters (space, tab, pipe) to newlines
	for(string :: iterator i = line.begin(); i != line.end(); i++) {
		assert(*i != '\n');
		if(*i == ',' || *i == ' ' || *i == '\t' || *i == '|')
			*i = '\n';
	}
	istringstream fields(line);
	string sourceNodeName;
	string targetNodeName;
	string weightAsString;

	getline(fields, sourceNodeName); 
	if( fields.fail() ) throw BadlyFormattedLine( 5, line);
	getline(fields, targetNodeName); 
	if( fields.fail() ) throw BadlyFormattedLine( 5, line);
	getline(fields, weightAsString); 
	// it's OK if this last one fails, the user doesn't have to specify a weight.

	ThreeStrings t;
	t.first.first = sourceNodeName;
	t.first.second = targetNodeName;
	t.second = weightAsString;
	return t;
}

static void read_edge_list_from_file(ModifiableNetwork *network, const string file_name) {
	PP(file_name);
	ifstream f(file_name.c_str());
	string line;
	while( getline(f, line) ) {
		// There might be a '\r' at the end of this line (dammit!)
		if(!line.empty() && *line.rbegin() == '\r') { line.erase( line.length()-1, 1); }
		ThreeStrings t = parseLine(line);
	}
}

std :: auto_ptr<Network> make_Network_from_edge_list(const std :: string file_name, const bool directed, const bool weighted) {
	ModifiableNetwork *network = new ModifiableNetwork(directed, weighted);
	read_edge_list_from_file(network, file_name);
	return auto_ptr<Network>(network);
}

BadlyFormattedLine :: BadlyFormattedLine(int _line_number, std :: string _bad_line) : line_number(_line_number), bad_line(_bad_line) {
}
const char* BadlyFormattedLine :: what() const throw() {
	ostringstream s;
	s << "Error: badly formatted line in edge list at line " << this->line_number << ". ExitinG. \"" << this->bad_line << "\"" << endl;
	return s.str().c_str();
}
BadlyFormattedLine :: ~ BadlyFormattedLine() throw() {
}

} // namespace loading
} // namespace graph
