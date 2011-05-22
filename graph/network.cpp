#include "network.hpp"
namespace graph {

template <class NodeNameT>
Network<NodeNameT> :: Network(const bool directed, const bool weighted) {
	if(directed) {
		if(weighted) {
			edge_weights.reset(new graph :: weights :: EdgeDetails< graph :: weights :: DirectedLDoubleWeights >());
		} else {
			edge_weights.reset(new graph :: weights :: EdgeDetails< graph :: weights :: DirectedNoWeights >());
		}
	} else {
		if(weighted) {
			edge_weights.reset(new graph :: weights :: EdgeDetails< graph :: weights :: WeightNoDir >());
		} else {
			edge_weights.reset(new graph :: weights :: EdgeDetails< graph :: weights :: NoDetails >());
		}
	}
}
template <class NodeNameT>
Network<NodeNameT> :: ~ Network() throw() {
}

template Network<NodeNameIsInt32>  :: Network(const bool directed, const bool weighted);
template Network<NodeNameIsString> :: Network(const bool directed, const bool weighted);
template Network<NodeNameIsInt32>  :: ~Network();
template Network<NodeNameIsString> :: ~Network();
} // namespace graph
