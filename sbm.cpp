using namespace std;
#include <getopt.h>
#include <unistd.h>
#include <libgen.h>

#include "aaron_utils.hpp"
#include "shmGraphRaw.hpp"
#include "sbm_state.hpp"

const char gitstatus[] = 
#include "comment.txt"
#include "gitstatus.txt"
;



struct UsageMessage {
};

void runSBM(const sbm::GraphType *g);

int main(int argc, char **argv) {
	PP(gitstatus);
	for (int i=0; i<argc; i++) {
		PP(argv[i]);
	}
	{ int c, option_index; while (1)
		{
      static const struct option long_options[] = {
        // {"seed", required_argument,       0, 22},
        {0, 0, 0, 0}
      };
      /* getopt_long stores the option index here. */
      c = getopt_long (argc, argv, "t:T:k:", long_options, &option_index);
      if (c == -1) break; /* Detect the end of the options. */
     
      switch (c) {
        case '?': /* getopt_long already printed an error message. */ break;
        default: abort (); break;
        case 0:
          /* If this option set a flag, do nothing else now. */
          if (long_options[option_index].flag != 0)
            break;
          printf ("option %s", long_options[option_index].name);
          if (optarg) printf (" with arg %s", optarg);
          printf ("\n");
          break;
        // case 'T':
					// option_inclusiveThreshold = true;
      }
    }
	}

	unless (argc - optind == 2) {
		throw UsageMessage();
		// cout << "Usage: edge_list directory_for_output" << endl;
		// exit(1);
	}
	const char * edgeListFileName = argv[optind];
	const char * directoryForOutput = argv[optind+1];
	PP(edgeListFileName);
	PP(directoryForOutput);

	auto_ptr<shmGraphRaw::ReadableShmGraphTemplate<shmGraphRaw::PlainMem> > g (shmGraphRaw::loadEdgeList<shmGraphRaw::PlainMem>(edgeListFileName));
	runSBM(g.get());
}

void randomize(sbm::State &s, const int K) { // randomize the partition and have K clusters in it
	assert(s._k <= K);
	for(int i=0; i<100000 || s._k < K; i++) {
		int n;
		do {
			n = drand48() * s._N;
			//cout << endl << "Moving node: " << n << " move# " << i << endl;
			s.isolateNode(n);
		} while (s._k <= K);
		const int newClusterID = drand48() * (s._k-1); // -1 because we must move the node from the "temporary" community
		// s.shortSummary(); s.summarizeEdgeCounts();
		s.internalCheck();

		s.unIsolateTempNode(n, newClusterID);
		assert(s._k <= K);
		// s.shortSummary(); s.summarizeEdgeCounts();
		s.internalCheck();
	}
	assert(s._k == K);
}

void runSBM(const sbm::GraphType *g) {
	sbm::State s(g);

	s.shortSummary(); s.summarizeEdgeCounts(); s.internalCheck();

	s.isolateNode(0);
	s.isolateNode(1); // to bring us up to three clusters
	s.shortSummary(); s.summarizeEdgeCounts(); s.internalCheck();
	PP(s.pmf());

	randomize(s, 3);
	s.shortSummary(); s.summarizeEdgeCounts();
	PP(s.pmf());
}
