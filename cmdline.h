/** @file cmdline.h
 *  @brief The header file for the command line option parser
 *  generated by GNU Gengetopt version 2.22.2
 *  http://www.gnu.org/software/gengetopt.
 *  DO NOT modify this file, since it can be overwritten
 *  @author GNU Gengetopt by Lorenzo Bettini */

#ifndef CMDLINE_H
#define CMDLINE_H

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h> /* for FILE */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef CMDLINE_PARSER_PACKAGE
/** @brief the program name (used for printing errors) */
#define CMDLINE_PARSER_PACKAGE "Stochastic Block Models"
#endif

#ifndef CMDLINE_PARSER_PACKAGE_NAME
/** @brief the complete program name (used for help and version) */
#define CMDLINE_PARSER_PACKAGE_NAME "Stochastic Block Models"
#endif

#ifndef CMDLINE_PARSER_VERSION
/** @brief the program version */
#define CMDLINE_PARSER_VERSION "0.1"
#endif

/** @brief Where the command line options are stored */
struct gengetopt_args_info
{
  const char *help_help; /**< @brief Print help and exit help description.  */
  const char *version_help; /**< @brief Print version and exit help description.  */
  int git_version_flag;	/**< @brief detailed version description (default=off).  */
  const char *git_version_help; /**< @brief detailed version description help description.  */
  int verbose_flag;	/**< @brief detailed debugging (default=off).  */
  const char *verbose_help; /**< @brief detailed debugging help description.  */
  int K_arg;	/**< @brief Number of clusters, K (default='-1').  */
  char * K_orig;	/**< @brief Number of clusters, K original value given at command line.  */
  const char *K_help; /**< @brief Number of clusters, K help description.  */
  int directed_flag;	/**< @brief directed (default=off).  */
  const char *directed_help; /**< @brief directed help description.  */
  int weighted_flag;	/**< @brief weighted (default=off).  */
  const char *weighted_help; /**< @brief weighted help description.  */
  int selfloop_flag;	/**< @brief selfloops allowed (default=off).  */
  const char *selfloop_help; /**< @brief selfloops allowed help description.  */
  int seed_arg;	/**< @brief seed to drand48() (default='0').  */
  char * seed_orig;	/**< @brief seed to drand48() original value given at command line.  */
  const char *seed_help; /**< @brief seed to drand48() help description.  */
  char * GT_vector_arg;	/**< @brief The ground truth. a file with N lines. Starts from ZERO..  */
  char * GT_vector_orig;	/**< @brief The ground truth. a file with N lines. Starts from ZERO. original value given at command line.  */
  const char *GT_vector_help; /**< @brief The ground truth. a file with N lines. Starts from ZERO. help description.  */
  int algo_metroK_arg;	/**< @brief Use the simple Metropolis move on K (default='1').  */
  char * algo_metroK_orig;	/**< @brief Use the simple Metropolis move on K original value given at command line.  */
  const char *algo_metroK_help; /**< @brief Use the simple Metropolis move on K help description.  */
  int algo_1node_arg;	/**< @brief Use the simple Metropolis on one node (default='1').  */
  char * algo_1node_orig;	/**< @brief Use the simple Metropolis on one node original value given at command line.  */
  const char *algo_1node_help; /**< @brief Use the simple Metropolis on one node help description.  */
  int algo_gibbs_arg;	/**< @brief Use the simple Gibbs in the algorithm (default='1').  */
  char * algo_gibbs_orig;	/**< @brief Use the simple Gibbs in the algorithm original value given at command line.  */
  const char *algo_gibbs_help; /**< @brief Use the simple Gibbs in the algorithm help description.  */
  int algo_m3_arg;	/**< @brief Use M3 in the algorithm (default='1').  */
  char * algo_m3_orig;	/**< @brief Use M3 in the algorithm original value given at command line.  */
  const char *algo_m3_help; /**< @brief Use M3 in the algorithm help description.  */
  int algo_ejectabsorb_arg;	/**< @brief Use N+F's eject-absorb move (default='1').  */
  char * algo_ejectabsorb_orig;	/**< @brief Use N+F's eject-absorb move original value given at command line.  */
  const char *algo_ejectabsorb_help; /**< @brief Use N+F's eject-absorb move help description.  */
  int iterations_arg;	/**< @brief How many iterations (default='120000').  */
  char * iterations_orig;	/**< @brief How many iterations original value given at command line.  */
  const char *iterations_help; /**< @brief How many iterations help description.  */
  int initGT_flag;	/**< @brief Initialize to the ground truth (default=off).  */
  const char *initGT_help; /**< @brief Initialize to the ground truth help description.  */
  int model_scf_flag;	/**< @brief Stochastic community finding (default=off).  */
  const char *model_scf_help; /**< @brief Stochastic community finding help description.  */
  int stringIDs_flag;	/**< @brief string IDs in the input (default=off).  */
  const char *stringIDs_help; /**< @brief string IDs in the input help description.  */
  int mega_flag;	/**< @brief dumb down the algorithm for *big* networks (default=off).  */
  const char *mega_help; /**< @brief dumb down the algorithm for *big* networks help description.  */
  int printEveryNIters_arg;	/**< @brief How often to print an update (default='10').  */
  char * printEveryNIters_orig;	/**< @brief How often to print an update original value given at command line.  */
  const char *printEveryNIters_help; /**< @brief How often to print an update help description.  */
  int assume_N_nodes_arg;	/**< @brief Pre-create N nodes (0 to N-1), which may be left with zero degree (default='0').  */
  char * assume_N_nodes_orig;	/**< @brief Pre-create N nodes (0 to N-1), which may be left with zero degree original value given at command line.  */
  const char *assume_N_nodes_help; /**< @brief Pre-create N nodes (0 to N-1), which may be left with zero degree help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int git_version_given ;	/**< @brief Whether git-version was given.  */
  unsigned int verbose_given ;	/**< @brief Whether verbose was given.  */
  unsigned int K_given ;	/**< @brief Whether K was given.  */
  unsigned int directed_given ;	/**< @brief Whether directed was given.  */
  unsigned int weighted_given ;	/**< @brief Whether weighted was given.  */
  unsigned int selfloop_given ;	/**< @brief Whether selfloop was given.  */
  unsigned int seed_given ;	/**< @brief Whether seed was given.  */
  unsigned int GT_vector_given ;	/**< @brief Whether GT.vector was given.  */
  unsigned int algo_metroK_given ;	/**< @brief Whether algo.metroK was given.  */
  unsigned int algo_1node_given ;	/**< @brief Whether algo.1node was given.  */
  unsigned int algo_gibbs_given ;	/**< @brief Whether algo.gibbs was given.  */
  unsigned int algo_m3_given ;	/**< @brief Whether algo.m3 was given.  */
  unsigned int algo_ejectabsorb_given ;	/**< @brief Whether algo.ejectabsorb was given.  */
  unsigned int iterations_given ;	/**< @brief Whether iterations was given.  */
  unsigned int initGT_given ;	/**< @brief Whether initGT was given.  */
  unsigned int model_scf_given ;	/**< @brief Whether model.scf was given.  */
  unsigned int stringIDs_given ;	/**< @brief Whether stringIDs was given.  */
  unsigned int mega_given ;	/**< @brief Whether mega was given.  */
  unsigned int printEveryNIters_given ;	/**< @brief Whether printEveryNIters was given.  */
  unsigned int assume_N_nodes_given ;	/**< @brief Whether assume_N_nodes was given.  */

  char **inputs ; /**< @brief unamed options (options without names) */
  unsigned inputs_num ; /**< @brief unamed options number */
} ;

/** @brief The additional parameters to pass to parser functions */
struct cmdline_parser_params
{
  int override; /**< @brief whether to override possibly already present options (default 0) */
  int initialize; /**< @brief whether to initialize the option structure gengetopt_args_info (default 1) */
  int check_required; /**< @brief whether to check that all required options were provided (default 1) */
  int check_ambiguity; /**< @brief whether to check for options already specified in the option structure gengetopt_args_info (default 0) */
  int print_errors; /**< @brief whether getopt_long should print an error message for a bad option (default 1) */
} ;

/** @brief the purpose string of the program */
extern const char *gengetopt_args_info_purpose;
/** @brief the usage string of the program */
extern const char *gengetopt_args_info_usage;
/** @brief all the lines making the help output */
extern const char *gengetopt_args_info_help[];

/**
 * The command line parser
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser (int argc, char * const *argv,
  struct gengetopt_args_info *args_info);

/**
 * The command line parser (version with additional parameters - deprecated)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use cmdline_parser_ext() instead
 */
int cmdline_parser2 (int argc, char * const *argv,
  struct gengetopt_args_info *args_info,
  int override, int initialize, int check_required);

/**
 * The command line parser (version with additional parameters)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param params additional parameters for the parser
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_ext (int argc, char * const *argv,
  struct gengetopt_args_info *args_info,
  struct cmdline_parser_params *params);

/**
 * Save the contents of the option struct into an already open FILE stream.
 * @param outfile the stream where to dump options
 * @param args_info the option struct to dump
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_dump(FILE *outfile,
  struct gengetopt_args_info *args_info);

/**
 * Save the contents of the option struct into a (text) file.
 * This file can be read by the config file parser (if generated by gengetopt)
 * @param filename the file where to save
 * @param args_info the option struct to save
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_file_save(const char *filename,
  struct gengetopt_args_info *args_info);

/**
 * Print the help
 */
void cmdline_parser_print_help(void);
/**
 * Print the version
 */
void cmdline_parser_print_version(void);

/**
 * Initializes all the fields a cmdline_parser_params structure 
 * to their default values
 * @param params the structure to initialize
 */
void cmdline_parser_params_init(struct cmdline_parser_params *params);

/**
 * Allocates dynamically a cmdline_parser_params structure and initializes
 * all its fields to their default values
 * @return the created and initialized cmdline_parser_params structure
 */
struct cmdline_parser_params *cmdline_parser_params_create(void);

/**
 * Initializes the passed gengetopt_args_info structure's fields
 * (also set default values for options that have a default)
 * @param args_info the structure to initialize
 */
void cmdline_parser_init (struct gengetopt_args_info *args_info);
/**
 * Deallocates the string fields of the gengetopt_args_info structure
 * (but does not deallocate the structure itself)
 * @param args_info the structure to deallocate
 */
void cmdline_parser_free (struct gengetopt_args_info *args_info);

/**
 * Checks that all the required options were specified
 * @param args_info the structure to check
 * @param prog_name the name of the program that will be used to print
 *   possible errors
 * @return
 */
int cmdline_parser_required (struct gengetopt_args_info *args_info,
  const char *prog_name);


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* CMDLINE_H */
