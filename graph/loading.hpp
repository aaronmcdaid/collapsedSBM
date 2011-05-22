#include "network.hpp"

#include <string>
#include <exception>
#include <memory>

namespace graph {
namespace loading {

struct BadlyFormattedLine : public std :: exception {
	const int line_number;
	const std :: string bad_line;
	BadlyFormattedLine(int _line_number, std :: string _bad_line);
	virtual const char* what() const throw();
	virtual ~ BadlyFormattedLine() throw();
};
	
std :: auto_ptr<Network> make_Network_from_edge_list(const std :: string file_name, const bool directed, const bool weighted);

} // namespace loading
} // namespace graph
