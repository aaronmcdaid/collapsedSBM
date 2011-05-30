#ifndef _FORMAT_FLAG_STACK_
#include <vector>
#include <ios>
#include <ostream>
namespace format_flag_stack {
struct FormatFlagStack {
	std :: vector< std :: ios_base :: fmtflags> the_stack;
	std :: vector< std      :: streamsize > the_stack_of_precision; // setprecision(...)
	struct PushT {
		FormatFlagStack * const parent_stack;
		PushT (FormatFlagStack *parent_stack_);
	} push;
	struct PopT {
		FormatFlagStack * const parent_stack;
		PopT  (FormatFlagStack *parent_stack_);
	} pop;
	FormatFlagStack();
	void do_push(std :: ostream &str);
	void do_pop(std :: ostream &str);
};
std :: ostream & operator<< (std :: ostream &str, const FormatFlagStack :: PushT & pusher);
std :: ostream & operator<< (std :: ostream &str, const FormatFlagStack :: PopT  & pusher);

} // namespace format_flag_stack
#endif
