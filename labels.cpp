#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cassert>
#include "macros.hpp"
using namespace std;

int main() {
	string line;
	int N = -1;
	while(getline(cin, line)) {
		PP(line);
		For(i, line) {
			if(*i==':' || *i == ',')
				*i='\n';
		}
		istringstream fields(line);
		int K;
		fields >> K;
		PP(K);
		vector<int> z;
		int z_i;
		while(fields >> z_i) {
			z.push_back(z_i);
		}
		PP(z.size());
		if(N == -1)
			N = z.size();
		assert(N == int(z.size()));
	}
}
