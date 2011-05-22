namespace graph {
namespace weights {

struct EdgeDetailsInterface {
	virtual long double getl2h(const int relId) const = 0;
	virtual long double geth2l(const int relId) const = 0;
};

template <class W>
struct EdgeDetails : public EdgeDetailsInterface {
	std :: vector < typename W :: datumT > dw; // the directions and weights of all the edges
	int size() const {
		return this->dw.size();
	}
	void new_rel(int relId, std :: pair<int,int> nodeIds, std :: string &weight) {
		if(relId == (int)this->dw.size())
			this->dw.push_back( typename W :: datumT() );
		assert(relId <  (int)this->dw.size());
		this->dw.at(relId).inform(nodeIds.first > nodeIds.second, weight);
		// this->dw.back().inform(nodeIds.first > nodeIds.second, weight);
	}
	virtual long double getl2h(const int relId) const { // this returns the value in the undirected case, and it handles self loops
		return this->dw.at(relId).getl2h();
	}
	virtual long double geth2l(const int relId) const { // this is only relevant in directed graphs.
		return this->dw.at(relId).geth2l();
	}
};

struct NoDetails { // unweighted, undirected
	typedef struct {
		void inform(const bool, const std :: string) const {
		}
		long double getl2h() const {
			return 1.0L;
		}
		long double geth2l() const {
			return 0.0L;
		}
	} datumT;
};
struct DirectedLDoubleWeights {
	typedef struct LdblPair : public std :: pair<long double,long double> {
		LdblPair() {
			this->first = this->second = 0.0L;
		}
		struct DuplicateWeightedEdge : public std :: exception {
		};
		void inform(const bool highToLow, const std :: string weight) { // the weight string might be empty. I suppose we let that default to 0
			if( highToLow ? this->second : this->first != 0) {
				throw DuplicateWeightedEdge();
			}
			std :: istringstream oss(weight);
			long double w = 0;
			assert(oss.peek() != EOF);
			oss >> w;
			assert(oss.peek() == EOF);
			if(highToLow) {
				this->second = w;
			} else {
				this->first = w;
			}
		}
		long double getl2h() const {
			return this->first;
		}
		long double geth2l() const {
			return this->second;
		}
	} datumT; // the type needed to store the weights in each direction.
};
struct DirectedNoWeights {
	typedef struct IntPair : public std :: pair<int,int> {
		IntPair() {
			this->first = this->second = 0;
		}
		void inform(const bool highToLow, const std :: string) { // the weight string might be empty. I suppose we let that default to 0
			if(highToLow) {
				this->second = 1;
			} else {
				this->first = 1;
			}
		}
		long double getl2h() const {
			return this->first;
		}
		long double geth2l() const {
			return this->second;
		}
	} datumT; // the type needed to store the weights in each direction.
};
struct WeightNoDir {
	typedef struct LdblPair {
		long double weight;
		LdblPair() {
			this->weight = 0.0L;
		}
		struct DuplicateWeightedEdge : public std :: exception {
		};
		void inform(const bool , const std :: string weight) { // the weight string might be empty. I suppose we let that default to 0
			if( this->weight != 0) {
				throw DuplicateWeightedEdge();
			}
			std :: istringstream oss(weight);
			long double w = 0;
			assert(oss.peek() != EOF);
			oss >> w;
			assert(oss.peek() == EOF);
			this->weight = w;
		}
		long double getl2h() const {
			return this->weight;
		}
		long double geth2l() const {
			return 0.0L;
		}
	} datumT; // the type needed to store the weights in each direction.
};

} // namespace weights
} // namespace graph