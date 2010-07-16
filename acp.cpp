using namespace std;
#include <getopt.h>

#include "aaron_utils.hpp"
#include "shmGraphRaw.hpp"
#include "cliques.hpp"
#include "clique_percolation.hpp"

int main(int argc, char **argv) {
	// for (int i=0; i<argc; i++) {
		// PP(argv[i]);
	// }
	{ int c, option_index; while (1)
		{
      static const struct option long_options[] = {
        // {"seed", required_argument,       0, 22},
        {0, 0, 0, 0}
      };
      /* getopt_long stores the option index here. */
      c = getopt_long (argc, argv, "", long_options, &option_index);
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

/*
        case 21: // --saveMOSESscores= {"saveMOSESscores", required_argument,         0, 21},
					strcpy(option_saveMOSESscores, optarg);
					break;
        case 22: // --seed= {"seed", required_argument,         0, 22},
					option_seed = atol(optarg);
					break;
*/

      }
    }
	}

	if (argc - optind != 2) {
		cout << "Usage: edge_list directory_for_output" << endl;
		exit(1);
	}
	const char * edgeListFileName = argv[optind];
	const char * directoryForOutput = argv[optind+1];
	PP(edgeListFileName);
	PP(directoryForOutput);

	auto_ptr<shmGraphRaw::ReadableShmGraph> g (shmGraphRaw::loadMmapFile("./binaryBlob", edgeListFileName));
	PP(g->numNodes());
	PP(g->numRels());
	// cliques::cliquesToDirectory(g.get(), "acp_results", 3);
	cliquePercolation(g.get(), directoryForOutput, 3); // You're not allowed to ask for the 2-cliques
}
