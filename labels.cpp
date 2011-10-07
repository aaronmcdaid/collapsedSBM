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

void recurse(const int K, const vector< vector<int64_t> > &k_by_k, const int k, vector<int> &already_taken, const int64_t score_so_far, int64_t & best_score_so_far, vector< vector<int> > & best_relabellings_so_far) {
	if(score_so_far < best_score_so_far) {
		// it can only get worse! return
		return;
	}
	if(k==K) { // all assigned, print and return
		// PP(score_so_far);
		if(score_so_far > best_score_so_far) {
			best_score_so_far = score_so_far;
			best_relabellings_so_far.clear();
		}
		if(score_so_far == best_score_so_far) {
			vector<int> inverse(K, -1);
			for (int j=0; j<K; j++) {
				const int x = already_taken.at(j);
				// PP(x);
				inverse.at(x) = j;
			}
			// For(rel, inverse) { PP(*rel); }
			// best_relabellings_so_far.push_back(already_taken);
			best_relabellings_so_far.push_back(inverse);
		}
		return;
	}
	assert(k<K);
	// what should cluster k be relabelled to?
	for(int new_k=0; new_k<K; ++new_k) {
		if(already_taken.at(new_k) == -1) {
			already_taken.at(new_k) = k;
			assert(k_by_k.at(k).at(new_k) <= 0);
			const int64_t new_score_so_far = score_so_far + k_by_k.at(k).at(new_k);
			// PP3(k,new_k, new_score_so_far);
			recurse(K, k_by_k, k+1, already_taken, new_score_so_far, best_score_so_far, best_relabellings_so_far);
			already_taken.at(new_k) = -1;
		}
	}
}

const pair<int64_t, vector<int> > find_best_labellings(const int K, const vector< vector<int64_t> > &k_by_k) {
	int64_t best_score_so_far = numeric_limits<int>::min();
	assert(K == int(k_by_k.size()));
	vector<int> already_taken(K, -1);
	vector< vector<int> > best_relabellings_so_far;
	recurse(K, k_by_k, 0, already_taken, 0, best_score_so_far, best_relabellings_so_far);
	PP(best_relabellings_so_far.size());
	PP(best_score_so_far);
	assert(1<=best_relabellings_so_far.size());
	const int offset = drand48() * best_relabellings_so_far.size();
	return make_pair(best_score_so_far,  best_relabellings_so_far.at(offset));
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
	int num_kz_so_far = 0;
	For(kz, kzs) {
		const int K = kz->first;
		PP3(num_kz_so_far, K, kz->second.at(0));
		// decide which relabelling to use:
		// - try an optimisitc relabelling first
		// - else fail
		// store it
		//
		vector<int> relabelling;
		if(kz == kzs.begin()) { // just leave the first one as is
			// for(int k=K-1; k>=0; k--) {
			for(int k=0; k<K; k++) {
				relabelling.push_back(k);
			}
		} else {

			// but first, build the K*K table of relabelling scores.
			vector< vector<int64_t> > k_by_k(K, vector<int64_t>(K,0));
			// k_by_k is like a summary of node_k_counts, but where the nodes are summed into their clusters.
			for(int n=0; n<74; n++) {
				const int z_i = kz->second.at(n);
				for(int k=0; k<K; k++) {
					k_by_k.at(z_i).at(k) += node_k_counts[make_pair(n,k)];
					// PP3(z_i, k, k_by_k.at(z_i).at(k));
				}
			}
			// the aim now is to find a nice relabelling that maximizes the sums in k_by_k

			// to enable optimization, subtract from each row of k its maximum.
			for(int k=0; k<K; k++) {
				vector<int64_t> & this_row = k_by_k.at(k);
				const int64_t max_in_this_row = *max_element(this_row.begin(),this_row.end());
				For(y, this_row) {
					*y -= max_in_this_row;
				}
			}
			const pair<int64_t, vector<int> > one_of_the_best_labellings = find_best_labellings(K, k_by_k);
			const int64_t best_score_so_far = one_of_the_best_labellings.first;
			{
				vector<int> verify_relabelling(one_of_the_best_labellings.second);
				assert(int(verify_relabelling.size()) == K);
				sort(verify_relabelling.begin(), verify_relabelling.end());
				for(int k=0; k<K; k++) {
					assert(verify_relabelling.at(k) == k);
				}
			}
			PP(best_score_so_far);
			{ // print the relabelled form
				for(int n=0; n<N; n++) {
					const int z_i = kz->second.at(n);
					const int rel_z_i = one_of_the_best_labellings.second.at(z_i);
					cout << ' ' << rel_z_i;
				}
				cout << endl;
			}
			assert(best_score_so_far == 0);
			For(x, one_of_the_best_labellings.second)
				relabelling.push_back(*x);
		}

		assert(int(relabelling.size()) == K);
		vector<int> relabelled_z;
		int verification_score = 0;
		for(int n=0; n<N; n++) {
			const int z_i = kz->second.at(n);
			const int rel_z_i = relabelling.at(z_i);
			verification_score += node_k_counts[make_pair(n, rel_z_i)];
		}
		// PP3(K,verification_score, num_kz_so_far);
		// assert(verification_score >= 74 * num_kz_so_far);
		for(int n=0; n<N; n++) {
			relabelled_z.push_back(relabelling.at(kz->second.at(n)));
		}
		for(int n=0; n<N; n++) {
			const int relabelled_z_i = relabelled_z.at(n);
			node_k_counts[make_pair(n,relabelled_z_i)] ++;
		}
		++ num_kz_so_far;
		// if(kz != kzs.begin()) exit(1);
	}
}
