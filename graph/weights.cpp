#include "weights.hpp"
namespace graph {
namespace weights {

	template <class W>
	long double EdgeDetails<W> :: geth2l(const int relId) const { // this is only relevant in directed graphs.
		return this->dw.at(relId).geth2l();
	}

// template  // <> here would forward-declare template (or something like that). We need to leave out the <> in order to ensure it's actually created
template long double EdgeDetails<NoDetails> :: geth2l(const int relId) const;
template long double EdgeDetails<DirectedLDoubleWeights> :: geth2l(const int relId) const;
template long double EdgeDetails<DirectedNoWeights> :: geth2l(const int relId) const;
template long double EdgeDetails<WeightNoDir> :: geth2l(const int relId) const;
} // namespace weights
} // namespace graph

