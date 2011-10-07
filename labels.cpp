#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <limits>
#include <tr1/unordered_map>
#include "macros.hpp"
using namespace std;

typedef pair<int, vector<int> > kz_t;

struct sort_on_pair_first {
	bool operator() (const kz_t &l, const kz_t &r) const {
		return l.first < r.first;
	}
};
struct hash_pair {
	size_t operator() (const pair<int,int> &v ) const {
		return std :: tr1 :: hash<int>()(-1-v.first) ^ std :: tr1 :: hash<int>()(v.second) ;
	}
};

void recurse(const int K, const vector< vector<int> > &k_by_k, const int k, vector<bool> &already_taken, const int score_so_far, int & best_score_so_far) {
	if(score_so_far < best_score_so_far) {
		// it can only get worse! return
		return;
	}
	if(k==K) { // all assigned, print and return
		PP(score_so_far);
		if(score_so_far > best_score_so_far) {
			best_score_so_far = score_so_far;
		}
		return;
	}
	assert(k<K);
	// what should cluster k be relabelled to?
	for(int new_k=0; new_k<K; ++new_k) {
		if(!already_taken.at(new_k)) {
			already_taken.at(new_k) = true;
			assert(k_by_k.at(k).at(new_k) <= 0);
			const int new_score_so_far = score_so_far + k_by_k.at(k).at(new_k);
			PP3(k,new_k, new_score_so_far);
			recurse(K, k_by_k, k+1, already_taken, new_score_so_far, best_score_so_far);
			already_taken.at(new_k) = false;
		}
	}
}

void find_best_labellings(const int K, const vector< vector<int> > &k_by_k) {
	int best_score_so_far = numeric_limits<int>::min();
	assert(K == int(k_by_k.size()));
	vector<bool> already_taken(K, false);
	recurse(K, k_by_k, 0, already_taken, 0, best_score_so_far);
}

int main() {
	string line;
	int N = -1;
	vector<kz_t> kzs;
	while(getline(cin, line)) {
		istringstream fields(line);
		int K;
		char sink;
		fields >> K;
		fields >> sink;
		vector<int> z;
		int z_i;
		while(fields >> z_i) {
			fields >> sink;
			z.push_back(z_i);
		}
		assert(fields.eof());
		if(N == -1)
			N = z.size();
		assert(N == int(z.size()));
		kzs.push_back(make_pair(K, z));
	}
	PP2(N, kzs.size());

	stable_sort(kzs.begin(), kzs.end(), sort_on_pair_first());

	typedef tr1 :: unordered_map< pair<int,int> , int, hash_pair> node_k_counts_t; // each node,cluster pair and how often it was assigned

	node_k_counts_t node_k_counts;
	For(kz, kzs) {
		const int K = kz->first;
		PP2(K, kz->second.at(0));
		// decide which relabelling to use:
		// - try an optimisitc relabelling first
		// - else fail
		// store it
		//
		vector<int> relabelling;
		if(kz == kzs.begin()) { // just leave the first one as is
			for(int k=K-1; k>=0; k--) {
				relabelling.push_back(k);
			}
		} else {

			// but first, build the K*K table of relabelling scores.
			vector< vector<int> > k_by_k(K, vector<int>(K,0));
			// k_by_k is like a summary of node_k_counts, but where the nodes are summed into their clusters.
			for(int n=0; n<74; n++) {
				const int z_i = kz->second.at(n);
				for(int k=0; k<K; k++) {
					k_by_k.at(z_i).at(k) += node_k_counts[make_pair(n,k)];
					PP3(z_i, k, k_by_k.at(z_i).at(k));
				}
			}
			// the aim now is to find a nice relabelling that maximizes the sums in k_by_k

			// to enable optimization, subtract from each row of k its maximum.
			for(int k=0; k<K; k++) {
				vector<int> & this_row = k_by_k.at(k);
				const int max_in_this_row = *max_element(this_row.begin(),this_row.end());
				For(y, this_row) {
					*y -= max_in_this_row;
				}
			}
			find_best_labellings(K, k_by_k);
		}

		assert(int(relabelling.size()) == K);
		vector<int> relabelled_z;
		for(int n=0; n<N; n++) {
			relabelled_z.push_back(relabelling.at(kz->second.at(n)));
		}
		for(int n=0; n<N; n++) {
			const int relabelled_z_i = relabelled_z.at(n);
			node_k_counts[make_pair(n,relabelled_z_i)] ++;
		}
		if(kz != kzs.begin())
			exit(1);
	}
}
