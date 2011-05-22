#include"loading.hpp"
#include"../pp.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <algorithm>

using namespace std;

namespace graph {
namespace loading {

struct ModifiableNetwork;
static void read_edge_list_from_file(ModifiableNetwork *network, const string file_name);
typedef pair< pair<string, string> , string> ThreeStrings;
static ThreeStrings parseLine(const string &lineOrig);

struct ModifiableNetwork : public Network { // Network is the read-only interface we want to expose, but this is the derived class that will do the heavy lifting
	vector<string> ordered_node_names;
	ModifiableNetwork(const bool directed, const bool weighted) : Network(directed, weighted) {
	}
	virtual ~ ModifiableNetwork() throw() {
	}
	int find_ordered_node_names_offset(const string &node_name) {
		// node_name must be in this->ordered_node_names. Now to find *where* it is.
		const int offset = lower_bound(this->ordered_node_names.begin(), this->ordered_node_names.end(), node_name) - this->ordered_node_names.begin();
		assert( this->ordered_node_names.at(offset) == node_name);
		return offset;
	}
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

static void read_edge_list_from_file(ModifiableNetwork *modifiable_network, const string file_name) {
	assert(modifiable_network->ordered_node_names.empty());
	/*
	 * This will make *three* passes:
	 * - One pass to identify all the node names (strings in the text file) and then to sort them so that consecutive-integer IDs can be assigned to them (respecting the order of the strings)
	 * - A second pass to identify all the unique relationships that exist, then they will be sorted
	 * - Finally, now that the node_names, node_ids and relationship_ids are sorted, read the network into the prepared datastructures, including the weights
	 */
	PP(file_name);
	{ // first pass: just store the node names
		ifstream f(file_name.c_str());
		string line;
		set<string> set_of_node_names; // This will store all the node names. Re
		while( getline(f, line) ) {
			// There might be a '\r' at the end of this line (dammit!)
			if(!line.empty() && *line.rbegin() == '\r') { line.erase( line.length()-1, 1); }
			ThreeStrings t = parseLine(line);
			set_of_node_names.insert(t.first.first);
			set_of_node_names.insert(t.first.second);
		}
		PP(set_of_node_names.size());
		for( set<string> :: const_iterator i = set_of_node_names.begin(); i != set_of_node_names.end(); i++) {
			modifiable_network->ordered_node_names.push_back(*i);
		}
		PP(modifiable_network->ordered_node_names.size());
		assert(modifiable_network->ordered_node_names.size() == set_of_node_names.size());
	}
	{ // second pass. Find all the distinct relationships (node_id_1, node_id_2; where node_id_1 <= node_id_2)
		ifstream f(file_name.c_str());
		string line;
		while( getline(f, line) ) {
			if(!line.empty() && *line.rbegin() == '\r') { line.erase( line.length()-1, 1); }
			ThreeStrings t = parseLine(line);
			const int source_node_id = modifiable_network->find_ordered_node_names_offset(t.first.first);
			const int target_node_id = modifiable_network->find_ordered_node_names_offset(t.first.second);
			PP2(source_node_id, target_node_id);
		}
	}
}

std :: auto_ptr<Network> make_Network_from_edge_list(const std :: string file_name, const bool directed, const bool weighted) throw(BadlyFormattedLine) {
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
