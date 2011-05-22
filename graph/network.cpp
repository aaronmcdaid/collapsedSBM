#include "network.hpp"
namespace graph {
	Network :: Network(const bool directed, const bool weighted) {
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
	Network :: ~ Network() throw() {
	}
} // namespace graph
