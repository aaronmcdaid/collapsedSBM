#include <map>
#include <vector>
void randomize(sbm :: State &s, const int K);
struct AcceptanceRate {
	std :: vector<bool> acc;
	std :: map<int,int> mostRecent; // in order to do a moving average
	int n;
	int a; // a/n is the acceptance rate
	const std :: string _name;
	AcceptanceRate(const char * name);
	void notify(bool accepted);
	void dump() const;
};
