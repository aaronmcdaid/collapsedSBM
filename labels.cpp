#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cassert>
#include "macros.hpp"
using namespace std;

typedef pair<int, vector<int> > kz_t;

struct sort_on_pair_first {
	bool operator() (const kz_t &l, const kz_t &r) const {
		return l.first < r.first;
	}
};

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

	For(kz, kzs) {
		PP2(kz->first, kz->second.at(0));
	}
}
