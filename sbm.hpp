void randomize(sbm::State &s, const int K);
struct AcceptanceRate {
	vector<bool> acc;
	map<int,int> mostRecent; // in order to do a moving average
	int n;
	int a; // a/n is the acceptance rate
	const string _name;
	AcceptanceRate(const char * name);
	void notify(bool accepted);
	void dump() const;
};
