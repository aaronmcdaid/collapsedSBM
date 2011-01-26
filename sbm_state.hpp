#include "shmGraphRaw.hpp"
namespace sbm {
struct State {
	const shmGraphRaw::ReadableShmGraphBase * const _g; // the graph
	explicit State(const shmGraphRaw::ReadableShmGraphBase * const g);

};

} // namespace sbm

