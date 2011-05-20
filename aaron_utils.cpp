#include "aaron_utils.hpp"

#include <ext/stdio_filebuf.h>
#include <sys/stat.h> // for stat(2)

int testForError(FILE * x) { return x==NULL; }
int testForError(int x) { return x==-1; }

string show(int64 x) { ostringstream s; s << x; return s.str(); }
string show(const char * x) { ostringstream s; s << x; return s.str(); }
string show(const string &s) { return s; }

runningAverage :: runningAverage() : total(0) , n(0) {}
void runningAverage :: operator() (int64 i) {
	total+=i;
	if(total<0) Die("total overflowed %Ld", i);
	n++;
}
string runningAverage :: operator() (void) const {
	ostringstream o;
	o << "Average Time (" << n << " points)" << ": " << (double)total/(double)n;
	return o.str();
}

DummyOutputStream& DummyOutputStream :: operator << (int) { return *this; }
DummyOutputStream& DummyOutputStream :: operator << (const char*) { return *this; }
DummyOutputStream dummyOutputStream;

string thousandsSeparated(uint64 x) {
	ostringstream s;
	if(x>=1000) {
		s << thousandsSeparated(x/1000) << ',';
		s << setfill('0');
	}
	s << setw(3) << x % 1000;
	x/=1000;
	return s.str();
}

istream *amd :: zcatThis(const char *gzippedFile) {
	ostringstream commandLine;
	const char * ZCATTHIS_ARG = "ZCAT_ARG_FOR_FORK";
	setenv(ZCATTHIS_ARG, gzippedFile, 1);
	commandLine << "zcat \"$" << ZCATTHIS_ARG << "\" "; //| head -n20";
	FILE * unzipped = popen(commandLine.str().c_str(), "r");
	__gnu_cxx :: stdio_filebuf<char> *myFileBuf = new __gnu_cxx :: stdio_filebuf<char>(unzipped,std :: ios_base :: in);
	return new istream(myFileBuf);
}
bool amd :: fileExists(const char *fileName) {
	struct stat buf;
	if(stat(fileName, &buf)==0) // The file already exists, we can grow it
		return true;
	else
		return false;
}
amd :: DebugCounter :: DebugCounter(const std :: string &_s, int _freq /* = 1000 */) : count(0), freq(_freq), s(_s) {}
void amd :: DebugCounter :: operator++ () {
	++count;
	if(count%freq == 0)
		cout << s << ":" << count << endl;
}

ostream & amd :: operator<< (ostream &str, const amd :: FormatFlagStack :: PushT & pusher) {
	pusher.parent_stack->do_push(str);
	return str;
}
ostream & amd :: operator<< (ostream &str, const amd :: FormatFlagStack :: PopT  & pusher) {
	pusher.parent_stack->do_pop (str);
	return str;
}
amd :: FormatFlagStack :: PushT :: PushT(FormatFlagStack *parent_stack_) : parent_stack(parent_stack_) {}
amd :: FormatFlagStack :: PopT  :: PopT (FormatFlagStack *parent_stack_) : parent_stack(parent_stack_) {}
amd :: FormatFlagStack :: FormatFlagStack() : push(this),pop(this) {}
void amd :: FormatFlagStack :: do_push(ostream &str) {
		this->the_stack.push_back( str.flags() );
		this->the_stack_of_precision.push_back( str.precision() );
}
void amd :: FormatFlagStack :: do_pop(ostream &str) {
		assert(!this->the_stack.empty());
		assert(this->the_stack.size() == this->the_stack_of_precision.size());
		str.flags     ( this->the_stack.back() );
		str.precision ( this->the_stack_of_precision.back() );
		this->the_stack.pop_back();
		this->the_stack_of_precision.pop_back();
}


