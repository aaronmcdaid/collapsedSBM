#include "saving.hpp"
#include "../pp.hpp"

namespace graph {
namespace saving {

void print_Network_to_screen( const graph :: NetworkInterfaceConvertedToString * net ) {
		for(int n=0; n< net->numNodes(); n++) {
			PP2(n, net->get_plain_graph()->degree(n));
			graph :: neighbouring_rel_id_iterator it(net->get_plain_graph(), n);
			while (!it.at_end()) {
				const int rel = *it++;
				PP(rel);
			}
		}
}

} // namespace graph
} // namespace saving

