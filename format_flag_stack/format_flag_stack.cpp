#include "format_flag_stack.hpp"
#include <cassert>

using namespace std;

namespace format_flag_stack {
ostream & operator<< (ostream &str, const FormatFlagStack :: PushT & pusher) {
	pusher.parent_stack->do_push(str);
	return str;
}
ostream & operator<< (ostream &str, const FormatFlagStack :: PopT  & pusher) {
	pusher.parent_stack->do_pop (str);
	return str;
}
FormatFlagStack :: PushT :: PushT(FormatFlagStack *parent_stack_) : parent_stack(parent_stack_) {}
FormatFlagStack :: PopT  :: PopT (FormatFlagStack *parent_stack_) : parent_stack(parent_stack_) {}
FormatFlagStack :: FormatFlagStack() : push(this),pop(this) {}
void FormatFlagStack :: do_push(ostream &str) {
		this->the_stack.push_back( str.flags() );
		this->the_stack_of_precision.push_back( str.precision() );
}
void  FormatFlagStack :: do_pop(ostream &str) {
		assert(!this->the_stack.empty());
		assert(this->the_stack.size() == this->the_stack_of_precision.size());
		str.flags     ( this->the_stack.back() );
		str.precision ( this->the_stack_of_precision.back() );
		this->the_stack.pop_back();
		this->the_stack_of_precision.pop_back();
}

} // namespace format_flag_stack
