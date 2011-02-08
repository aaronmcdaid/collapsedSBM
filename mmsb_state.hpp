#include "sbm_state.hpp"
namespace sbm {
	struct PairCounts : public boost::unordered_map< std::pair<int,int>, int > {
		int & set(const int i, const int j);
		int   get(const int i, const int j) const;
	};
	struct MMSBstate {
		const int _N;
		int _k;
		std::vector<Labelling*> ls; // one labelling for each node
		PairCounts numPairs;
		PairCounts numEdges;

		explicit MMSBstate (const GraphType * const g);
		void appendEmptyCluster();
		long double P_z_K() const;
		void P_zs_given_K() const;
	};
} // namespace sbm
