#include "weights.hpp"
namespace graph {
namespace weights {

	template <class W>
	int EdgeDetails<W> :: size() const {
		return this->dw.size();
	}
	template <class W>
	void EdgeDetails<W> :: new_rel(int relId, std :: pair<int,int> nodeIds, std :: string &weight) {
		if(relId == (int)this->dw.size())
			this->dw.push_back( typename W :: datumT() );
		assert(relId <  (int)this->dw.size());
		this->dw.at(relId).inform(nodeIds.first > nodeIds.second, weight);
		// this->dw.back().inform(nodeIds.first > nodeIds.second, weight);
	}
	template <class W>
	long double EdgeDetails<W> :: getl2h(const int relId) const { // this returns the value in the undirected case, and it handles self loops
		return this->dw.at(relId).getl2h();
	}
	template <class W>
	long double EdgeDetails<W> :: geth2l(const int relId) const { // this is only relevant in directed graphs.
		return this->dw.at(relId).geth2l();
	}

// template  // <> here would forward-declare template (or something like that). We need to leave out the <> in order to ensure it's actually created
template long double EdgeDetails<NoDetails>              :: geth2l(const int relId) const;
template long double EdgeDetails<DirectedLDoubleWeights> :: geth2l(const int relId) const;
template long double EdgeDetails<DirectedNoWeights>      :: geth2l(const int relId) const;
template long double EdgeDetails<WeightNoDir>            :: geth2l(const int relId) const;
template long double EdgeDetails<NoDetails>              :: getl2h(const int relId) const;
template long double EdgeDetails<DirectedLDoubleWeights> :: getl2h(const int relId) const;
template long double EdgeDetails<DirectedNoWeights>      :: getl2h(const int relId) const;
template long double EdgeDetails<WeightNoDir>            :: getl2h(const int relId) const;
template void        EdgeDetails<NoDetails>              :: new_rel(int relId, std :: pair<int,int> nodeIds, std :: string &weight);
template void        EdgeDetails<DirectedLDoubleWeights> :: new_rel(int relId, std :: pair<int,int> nodeIds, std :: string &weight);
template void        EdgeDetails<DirectedNoWeights>      :: new_rel(int relId, std :: pair<int,int> nodeIds, std :: string &weight);
template void        EdgeDetails<WeightNoDir>            :: new_rel(int relId, std :: pair<int,int> nodeIds, std :: string &weight);
template int         EdgeDetails<NoDetails>              :: size()  const;
template int         EdgeDetails<DirectedLDoubleWeights> :: size()  const;
template int         EdgeDetails<DirectedNoWeights>      :: size()  const;
template int         EdgeDetails<WeightNoDir>            :: size()  const;
} // namespace weights
} // namespace graph

