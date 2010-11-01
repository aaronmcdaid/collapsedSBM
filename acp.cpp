using namespace std;
#include <getopt.h>
#include <unistd.h>

#include "aaron_utils.hpp"
#include "shmGraphRaw.hpp"
#include "cliques.hpp"
#include "clique_percolation.hpp"
#include "clique_percolation3.hpp"

int option_minCliqueSize = 3;

vector< pair<double,bool> > option_thresholds;

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
      c = getopt_long (argc, argv, "t:T:k:", long_options, &option_index);
      if (c == -1) break; /* Detect the end of the options. */
     
			bool option_inclusiveThreshold = false; // this is not a global fixed option. Remember that multiple -t and -T may be specified by the user
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
        case 'T':
					option_inclusiveThreshold = true;
        case 't':
					{
					char * endptr;
					double threshold = strtod(optarg, &endptr);
					unless(threshold >= 0.0 && threshold <= 100.0)
						Die("Percentage relative percolation threshold should be between 0 and 100: %g. Exiting", threshold);
					unless(*endptr == '\0')
						Die("Badly formed -t or -T relative threshold. Must be a real number between 0 and 100: '%s'", optarg);
					if(option_inclusiveThreshold)
						unless(threshold > 0.0)
							Die("No point having a relative threshold of zero if it is inclusive. Exiting");

					option_thresholds.push_back(make_pair(threshold, option_inclusiveThreshold));
					}
					break;
        case 'k':
					option_minCliqueSize = atoi(optarg);
					unless(option_minCliqueSize >= 3)
						Die("-k option must be at least three: %d. Exiting", option_minCliqueSize);
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

	const char * defaultPrefix = "/tmp";
	if(getenv("TMP"))
		defaultPrefix = getenv("TMP");
	string directoryForBinaryBlobString = defaultPrefix;

	if(directoryForBinaryBlobString.length()==0) directoryForBinaryBlobString = ".";
	if(directoryForBinaryBlobString.at(directoryForBinaryBlobString.length()-1)!='/') directoryForBinaryBlobString += "/";
	directoryForBinaryBlobString += "acp-graph.XXXXXX";
	char directoryForBinaryBlob[1000];
	strcpy(directoryForBinaryBlob, directoryForBinaryBlobString.c_str());
	if(NULL == mkdtemp(directoryForBinaryBlob)) {
		cerr << "Couldn't create temp file: " << directoryForBinaryBlob << endl;
		cerr << "You could specify an alternative folder with the TMP environment variable" << endl;
		exit(1);
	}

	auto_ptr<shmGraphRaw::ReadableShmGraph> g (shmGraphRaw::loadMmapFile(directoryForBinaryBlob, edgeListFileName));
	PP(g->numNodes());
	PP(g->numRels());
	// cliques::cliquesToDirectory(g.get(), "acp_results", 3);
	if(0)
		cliquePercolation(g.get(), directoryForOutput, option_minCliqueSize); // You're not allowed to ask for the 2-cliques
	else
		cliquePercolation3(g.get(), directoryForOutput, option_minCliqueSize, option_thresholds); // You're not allowed to ask for the 2-cliques

	UNUSED int ignore = system( (string("rm -r ") + directoryForBinaryBlob) .c_str() );
}
