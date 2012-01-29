/*
  File autogenerated by gengetopt version 2.22.2
  generated with the following command:
  gengetopt --unamed-opts 

  The developers of gengetopt consider the fixed text that goes in all
  gengetopt output files to be in the public domain:
  we make no copyright claims on it.
*/

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef FIX_UNUSED
#define FIX_UNUSED(X) (void) (X) /* avoid warnings for unused params */
#endif

#include "getopt.h"

#include "cmdline.h"

const char *gengetopt_args_info_purpose = "Fit models to network data";

const char *gengetopt_args_info_usage = "Usage: Stochastic Block Models [OPTIONS]... [FILES]...";

const char *gengetopt_args_info_description = "";

const char *gengetopt_args_info_help[] = {
  "  -h, --help                  Print help and exit",
  "  -V, --version               Print version and exit",
  "      --git-version           detailed version description  (default=off)",
  "  -v, --verbose               detailed debugging  (default=off)",
  "  -K, --K=INT                 Number of clusters, K  (default=`-1')",
  "  -d, --directed              directed  (default=off)",
  "  -w, --weighted              weighted  (default=off)",
  "  -s, --selfloop              selfloops allowed  (default=off)",
  "      --seed=INT              seed to drand48()  (default=`0')",
  "      --GT.vector=STRING      The ground truth. a file with N lines. Starts \n                                from ZERO.",
  "      --algo.metroK=INT       Use the simple Metropolis move on K  \n                                (default=`1')",
  "      --algo.1node=INT        Use the simple Metropolis on one node  \n                                (default=`1')",
  "      --algo.gibbs=INT        Use the simple Gibbs in the algorithm  \n                                (default=`1')",
  "      --algo.m3=INT           Use M3 in the algorithm  (default=`1')",
  "      --algo.ejectabsorb=INT  Use N+F's eject-absorb move  (default=`1')",
  "  -i, --iterations=INT        How many iterations  (default=`120000')",
  "      --initGT                Initialize to the ground truth  (default=off)",
  "      --model.scf             Stochastic community finding  (default=off)",
  "      --algo.sbm.cem          Classification EM (CEM) for the SBM  \n                                (default=off)",
  "      --stringIDs             string IDs in the input  (default=off)",
  "      --mega                  dumb down the algorithm for *big* networks  \n                                (default=off)",
  "      --printEveryNIters=INT  How often to print an update  (default=`10')",
  "      --assume_N_nodes=INT    Pre-create N nodes (0 to N-1), which may be left \n                                with zero degree  (default=`0')",
  "  -a, --alpha=FLOAT           alpha. How uniform the cluster sizes  \n                                (default=`1')",
  "  -z, --save.z=STRING         save burnt-in z to this file  (default=`')",
  "      --gamma.s=FLOAT         (for weighted only). Shape of Gamma prior  \n                                (default=`1')",
  "      --gamma.phi=FLOAT       (for weighted only). Scale of Gamma prior  \n                                (default=`1')",
    0
};

typedef enum {ARG_NO
  , ARG_FLAG
  , ARG_STRING
  , ARG_INT
  , ARG_FLOAT
} cmdline_parser_arg_type;

static
void clear_given (struct gengetopt_args_info *args_info);
static
void clear_args (struct gengetopt_args_info *args_info);

static int
cmdline_parser_internal (int argc, char * const *argv, struct gengetopt_args_info *args_info,
                        struct cmdline_parser_params *params, const char *additional_error);


static char *
gengetopt_strdup (const char *s);

static
void clear_given (struct gengetopt_args_info *args_info)
{
  args_info->help_given = 0 ;
  args_info->version_given = 0 ;
  args_info->git_version_given = 0 ;
  args_info->verbose_given = 0 ;
  args_info->K_given = 0 ;
  args_info->directed_given = 0 ;
  args_info->weighted_given = 0 ;
  args_info->selfloop_given = 0 ;
  args_info->seed_given = 0 ;
  args_info->GT_vector_given = 0 ;
  args_info->algo_metroK_given = 0 ;
  args_info->algo_1node_given = 0 ;
  args_info->algo_gibbs_given = 0 ;
  args_info->algo_m3_given = 0 ;
  args_info->algo_ejectabsorb_given = 0 ;
  args_info->iterations_given = 0 ;
  args_info->initGT_given = 0 ;
  args_info->model_scf_given = 0 ;
  args_info->algo_sbm_cem_given = 0 ;
  args_info->stringIDs_given = 0 ;
  args_info->mega_given = 0 ;
  args_info->printEveryNIters_given = 0 ;
  args_info->assume_N_nodes_given = 0 ;
  args_info->alpha_given = 0 ;
  args_info->save_z_given = 0 ;
  args_info->gamma_s_given = 0 ;
  args_info->gamma_phi_given = 0 ;
}

static
void clear_args (struct gengetopt_args_info *args_info)
{
  FIX_UNUSED (args_info);
  args_info->git_version_flag = 0;
  args_info->verbose_flag = 0;
  args_info->K_arg = -1;
  args_info->K_orig = NULL;
  args_info->directed_flag = 0;
  args_info->weighted_flag = 0;
  args_info->selfloop_flag = 0;
  args_info->seed_arg = 0;
  args_info->seed_orig = NULL;
  args_info->GT_vector_arg = NULL;
  args_info->GT_vector_orig = NULL;
  args_info->algo_metroK_arg = 1;
  args_info->algo_metroK_orig = NULL;
  args_info->algo_1node_arg = 1;
  args_info->algo_1node_orig = NULL;
  args_info->algo_gibbs_arg = 1;
  args_info->algo_gibbs_orig = NULL;
  args_info->algo_m3_arg = 1;
  args_info->algo_m3_orig = NULL;
  args_info->algo_ejectabsorb_arg = 1;
  args_info->algo_ejectabsorb_orig = NULL;
  args_info->iterations_arg = 120000;
  args_info->iterations_orig = NULL;
  args_info->initGT_flag = 0;
  args_info->model_scf_flag = 0;
  args_info->algo_sbm_cem_flag = 0;
  args_info->stringIDs_flag = 0;
  args_info->mega_flag = 0;
  args_info->printEveryNIters_arg = 10;
  args_info->printEveryNIters_orig = NULL;
  args_info->assume_N_nodes_arg = 0;
  args_info->assume_N_nodes_orig = NULL;
  args_info->alpha_arg = 1;
  args_info->alpha_orig = NULL;
  args_info->save_z_arg = gengetopt_strdup ("");
  args_info->save_z_orig = NULL;
  args_info->gamma_s_arg = 1;
  args_info->gamma_s_orig = NULL;
  args_info->gamma_phi_arg = 1;
  args_info->gamma_phi_orig = NULL;
  
}

static
void init_args_info(struct gengetopt_args_info *args_info)
{


  args_info->help_help = gengetopt_args_info_help[0] ;
  args_info->version_help = gengetopt_args_info_help[1] ;
  args_info->git_version_help = gengetopt_args_info_help[2] ;
  args_info->verbose_help = gengetopt_args_info_help[3] ;
  args_info->K_help = gengetopt_args_info_help[4] ;
  args_info->directed_help = gengetopt_args_info_help[5] ;
  args_info->weighted_help = gengetopt_args_info_help[6] ;
  args_info->selfloop_help = gengetopt_args_info_help[7] ;
  args_info->seed_help = gengetopt_args_info_help[8] ;
  args_info->GT_vector_help = gengetopt_args_info_help[9] ;
  args_info->algo_metroK_help = gengetopt_args_info_help[10] ;
  args_info->algo_1node_help = gengetopt_args_info_help[11] ;
  args_info->algo_gibbs_help = gengetopt_args_info_help[12] ;
  args_info->algo_m3_help = gengetopt_args_info_help[13] ;
  args_info->algo_ejectabsorb_help = gengetopt_args_info_help[14] ;
  args_info->iterations_help = gengetopt_args_info_help[15] ;
  args_info->initGT_help = gengetopt_args_info_help[16] ;
  args_info->model_scf_help = gengetopt_args_info_help[17] ;
  args_info->algo_sbm_cem_help = gengetopt_args_info_help[18] ;
  args_info->stringIDs_help = gengetopt_args_info_help[19] ;
  args_info->mega_help = gengetopt_args_info_help[20] ;
  args_info->printEveryNIters_help = gengetopt_args_info_help[21] ;
  args_info->assume_N_nodes_help = gengetopt_args_info_help[22] ;
  args_info->alpha_help = gengetopt_args_info_help[23] ;
  args_info->save_z_help = gengetopt_args_info_help[24] ;
  args_info->gamma_s_help = gengetopt_args_info_help[25] ;
  args_info->gamma_phi_help = gengetopt_args_info_help[26] ;
  
}

void
cmdline_parser_print_version (void)
{
  printf ("%s %s\n",
     (strlen(CMDLINE_PARSER_PACKAGE_NAME) ? CMDLINE_PARSER_PACKAGE_NAME : CMDLINE_PARSER_PACKAGE),
     CMDLINE_PARSER_VERSION);
}

static void print_help_common(void) {
  cmdline_parser_print_version ();

  if (strlen(gengetopt_args_info_purpose) > 0)
    printf("\n%s\n", gengetopt_args_info_purpose);

  if (strlen(gengetopt_args_info_usage) > 0)
    printf("\n%s\n", gengetopt_args_info_usage);

  printf("\n");

  if (strlen(gengetopt_args_info_description) > 0)
    printf("%s\n\n", gengetopt_args_info_description);
}

void
cmdline_parser_print_help (void)
{
  int i = 0;
  print_help_common();
  while (gengetopt_args_info_help[i])
    printf("%s\n", gengetopt_args_info_help[i++]);
}

void
cmdline_parser_init (struct gengetopt_args_info *args_info)
{
  clear_given (args_info);
  clear_args (args_info);
  init_args_info (args_info);

  args_info->inputs = 0;
  args_info->inputs_num = 0;
}

void
cmdline_parser_params_init(struct cmdline_parser_params *params)
{
  if (params)
    { 
      params->override = 0;
      params->initialize = 1;
      params->check_required = 1;
      params->check_ambiguity = 0;
      params->print_errors = 1;
    }
}

struct cmdline_parser_params *
cmdline_parser_params_create(void)
{
  struct cmdline_parser_params *params = 
    (struct cmdline_parser_params *)malloc(sizeof(struct cmdline_parser_params));
  cmdline_parser_params_init(params);  
  return params;
}

static void
free_string_field (char **s)
{
  if (*s)
    {
      free (*s);
      *s = 0;
    }
}


static void
cmdline_parser_release (struct gengetopt_args_info *args_info)
{
  unsigned int i;
  free_string_field (&(args_info->K_orig));
  free_string_field (&(args_info->seed_orig));
  free_string_field (&(args_info->GT_vector_arg));
  free_string_field (&(args_info->GT_vector_orig));
  free_string_field (&(args_info->algo_metroK_orig));
  free_string_field (&(args_info->algo_1node_orig));
  free_string_field (&(args_info->algo_gibbs_orig));
  free_string_field (&(args_info->algo_m3_orig));
  free_string_field (&(args_info->algo_ejectabsorb_orig));
  free_string_field (&(args_info->iterations_orig));
  free_string_field (&(args_info->printEveryNIters_orig));
  free_string_field (&(args_info->assume_N_nodes_orig));
  free_string_field (&(args_info->alpha_orig));
  free_string_field (&(args_info->save_z_arg));
  free_string_field (&(args_info->save_z_orig));
  free_string_field (&(args_info->gamma_s_orig));
  free_string_field (&(args_info->gamma_phi_orig));
  
  
  for (i = 0; i < args_info->inputs_num; ++i)
    free (args_info->inputs [i]);

  if (args_info->inputs_num)
    free (args_info->inputs);

  clear_given (args_info);
}


static void
write_into_file(FILE *outfile, const char *opt, const char *arg, const char *values[])
{
  FIX_UNUSED (values);
  if (arg) {
    fprintf(outfile, "%s=\"%s\"\n", opt, arg);
  } else {
    fprintf(outfile, "%s\n", opt);
  }
}


int
cmdline_parser_dump(FILE *outfile, struct gengetopt_args_info *args_info)
{
  int i = 0;

  if (!outfile)
    {
      fprintf (stderr, "%s: cannot dump options to stream\n", CMDLINE_PARSER_PACKAGE);
      return EXIT_FAILURE;
    }

  if (args_info->help_given)
    write_into_file(outfile, "help", 0, 0 );
  if (args_info->version_given)
    write_into_file(outfile, "version", 0, 0 );
  if (args_info->git_version_given)
    write_into_file(outfile, "git-version", 0, 0 );
  if (args_info->verbose_given)
    write_into_file(outfile, "verbose", 0, 0 );
  if (args_info->K_given)
    write_into_file(outfile, "K", args_info->K_orig, 0);
  if (args_info->directed_given)
    write_into_file(outfile, "directed", 0, 0 );
  if (args_info->weighted_given)
    write_into_file(outfile, "weighted", 0, 0 );
  if (args_info->selfloop_given)
    write_into_file(outfile, "selfloop", 0, 0 );
  if (args_info->seed_given)
    write_into_file(outfile, "seed", args_info->seed_orig, 0);
  if (args_info->GT_vector_given)
    write_into_file(outfile, "GT.vector", args_info->GT_vector_orig, 0);
  if (args_info->algo_metroK_given)
    write_into_file(outfile, "algo.metroK", args_info->algo_metroK_orig, 0);
  if (args_info->algo_1node_given)
    write_into_file(outfile, "algo.1node", args_info->algo_1node_orig, 0);
  if (args_info->algo_gibbs_given)
    write_into_file(outfile, "algo.gibbs", args_info->algo_gibbs_orig, 0);
  if (args_info->algo_m3_given)
    write_into_file(outfile, "algo.m3", args_info->algo_m3_orig, 0);
  if (args_info->algo_ejectabsorb_given)
    write_into_file(outfile, "algo.ejectabsorb", args_info->algo_ejectabsorb_orig, 0);
  if (args_info->iterations_given)
    write_into_file(outfile, "iterations", args_info->iterations_orig, 0);
  if (args_info->initGT_given)
    write_into_file(outfile, "initGT", 0, 0 );
  if (args_info->model_scf_given)
    write_into_file(outfile, "model.scf", 0, 0 );
  if (args_info->algo_sbm_cem_given)
    write_into_file(outfile, "algo.sbm.cem", 0, 0 );
  if (args_info->stringIDs_given)
    write_into_file(outfile, "stringIDs", 0, 0 );
  if (args_info->mega_given)
    write_into_file(outfile, "mega", 0, 0 );
  if (args_info->printEveryNIters_given)
    write_into_file(outfile, "printEveryNIters", args_info->printEveryNIters_orig, 0);
  if (args_info->assume_N_nodes_given)
    write_into_file(outfile, "assume_N_nodes", args_info->assume_N_nodes_orig, 0);
  if (args_info->alpha_given)
    write_into_file(outfile, "alpha", args_info->alpha_orig, 0);
  if (args_info->save_z_given)
    write_into_file(outfile, "save.z", args_info->save_z_orig, 0);
  if (args_info->gamma_s_given)
    write_into_file(outfile, "gamma.s", args_info->gamma_s_orig, 0);
  if (args_info->gamma_phi_given)
    write_into_file(outfile, "gamma.phi", args_info->gamma_phi_orig, 0);
  

  i = EXIT_SUCCESS;
  return i;
}

int
cmdline_parser_file_save(const char *filename, struct gengetopt_args_info *args_info)
{
  FILE *outfile;
  int i = 0;

  outfile = fopen(filename, "w");

  if (!outfile)
    {
      fprintf (stderr, "%s: cannot open file for writing: %s\n", CMDLINE_PARSER_PACKAGE, filename);
      return EXIT_FAILURE;
    }

  i = cmdline_parser_dump(outfile, args_info);
  fclose (outfile);

  return i;
}

void
cmdline_parser_free (struct gengetopt_args_info *args_info)
{
  cmdline_parser_release (args_info);
}

/** @brief replacement of strdup, which is not standard */
char *
gengetopt_strdup (const char *s)
{
  char *result = 0;
  if (!s)
    return result;

  result = (char*)malloc(strlen(s) + 1);
  if (result == (char*)0)
    return (char*)0;
  strcpy(result, s);
  return result;
}

int
cmdline_parser (int argc, char * const *argv, struct gengetopt_args_info *args_info)
{
  return cmdline_parser2 (argc, argv, args_info, 0, 1, 1);
}

int
cmdline_parser_ext (int argc, char * const *argv, struct gengetopt_args_info *args_info,
                   struct cmdline_parser_params *params)
{
  int result;
  result = cmdline_parser_internal (argc, argv, args_info, params, 0);

  if (result == EXIT_FAILURE)
    {
      cmdline_parser_free (args_info);
      exit (EXIT_FAILURE);
    }
  
  return result;
}

int
cmdline_parser2 (int argc, char * const *argv, struct gengetopt_args_info *args_info, int override, int initialize, int check_required)
{
  int result;
  struct cmdline_parser_params params;
  
  params.override = override;
  params.initialize = initialize;
  params.check_required = check_required;
  params.check_ambiguity = 0;
  params.print_errors = 1;

  result = cmdline_parser_internal (argc, argv, args_info, &params, 0);

  if (result == EXIT_FAILURE)
    {
      cmdline_parser_free (args_info);
      exit (EXIT_FAILURE);
    }
  
  return result;
}

int
cmdline_parser_required (struct gengetopt_args_info *args_info, const char *prog_name)
{
  FIX_UNUSED (args_info);
  FIX_UNUSED (prog_name);
  return EXIT_SUCCESS;
}


static char *package_name = 0;

/**
 * @brief updates an option
 * @param field the generic pointer to the field to update
 * @param orig_field the pointer to the orig field
 * @param field_given the pointer to the number of occurrence of this option
 * @param prev_given the pointer to the number of occurrence already seen
 * @param value the argument for this option (if null no arg was specified)
 * @param possible_values the possible values for this option (if specified)
 * @param default_value the default value (in case the option only accepts fixed values)
 * @param arg_type the type of this option
 * @param check_ambiguity @see cmdline_parser_params.check_ambiguity
 * @param override @see cmdline_parser_params.override
 * @param no_free whether to free a possible previous value
 * @param multiple_option whether this is a multiple option
 * @param long_opt the corresponding long option
 * @param short_opt the corresponding short option (or '-' if none)
 * @param additional_error possible further error specification
 */
static
int update_arg(void *field, char **orig_field,
               unsigned int *field_given, unsigned int *prev_given, 
               char *value, const char *possible_values[],
               const char *default_value,
               cmdline_parser_arg_type arg_type,
               int check_ambiguity, int override,
               int no_free, int multiple_option,
               const char *long_opt, char short_opt,
               const char *additional_error)
{
  FIX_UNUSED (field);
  char *stop_char = 0;
  const char *val = value;
  int found;
  char **string_field;

  stop_char = 0;
  found = 0;

  if (!multiple_option && prev_given && (*prev_given || (check_ambiguity && *field_given)))
    {
      if (short_opt != '-')
        fprintf (stderr, "%s: `--%s' (`-%c') option given more than once%s\n", 
               package_name, long_opt, short_opt,
               (additional_error ? additional_error : ""));
      else
        fprintf (stderr, "%s: `--%s' option given more than once%s\n", 
               package_name, long_opt,
               (additional_error ? additional_error : ""));
      return 1; /* failure */
    }

  FIX_UNUSED (default_value);
    
  if (field_given && *field_given && ! override)
    return 0;
  if (prev_given)
    (*prev_given)++;
  if (field_given)
    (*field_given)++;
  if (possible_values)
    val = possible_values[found];

  switch(arg_type) {
  case ARG_FLAG:
    *((int *)field) = !*((int *)field);
    break;
  case ARG_INT:
    if (val) *((int *)field) = strtol (val, &stop_char, 0);
    break;
  case ARG_FLOAT:
    if (val) *((float *)field) = (float)strtod (val, &stop_char);
    break;
  case ARG_STRING:
    if (val) {
      string_field = (char **)field;
      if (!no_free && *string_field)
        free (*string_field); /* free previous string */
      *string_field = gengetopt_strdup (val);
    }
    break;
  default:
    break;
  };

  /* check numeric conversion */
  switch(arg_type) {
  case ARG_INT:
  case ARG_FLOAT:
    if (val && !(stop_char && *stop_char == '\0')) {
      fprintf(stderr, "%s: invalid numeric value: %s\n", package_name, val);
      return 1; /* failure */
    }
    break;
  default:
    ;
  };

  /* store the original value */
  switch(arg_type) {
  case ARG_NO:
  case ARG_FLAG:
    break;
  default:
    if (value && orig_field) {
      if (no_free) {
        *orig_field = value;
      } else {
        if (*orig_field)
          free (*orig_field); /* free previous string */
        *orig_field = gengetopt_strdup (value);
      }
    }
  };

  return 0; /* OK */
}


int
cmdline_parser_internal (
  int argc, char * const *argv, struct gengetopt_args_info *args_info,
                        struct cmdline_parser_params *params, const char *additional_error)
{
  int c;	/* Character of the parsed option.  */

  int error = 0;
  struct gengetopt_args_info local_args_info;
  
  int override;
  int initialize;
  int check_required;
  int check_ambiguity;
  
  package_name = argv[0];
  
  override = params->override;
  initialize = params->initialize;
  check_required = params->check_required;
  check_ambiguity = params->check_ambiguity;

  if (initialize)
    cmdline_parser_init (args_info);

  cmdline_parser_init (&local_args_info);

  optarg = 0;
  optind = 0;
  opterr = params->print_errors;
  optopt = '?';

  while (1)
    {
      int option_index = 0;

      static struct option long_options[] = {
        { "help",	0, NULL, 'h' },
        { "version",	0, NULL, 'V' },
        { "git-version",	0, NULL, 0 },
        { "verbose",	0, NULL, 'v' },
        { "K",	1, NULL, 'K' },
        { "directed",	0, NULL, 'd' },
        { "weighted",	0, NULL, 'w' },
        { "selfloop",	0, NULL, 's' },
        { "seed",	1, NULL, 0 },
        { "GT.vector",	1, NULL, 0 },
        { "algo.metroK",	1, NULL, 0 },
        { "algo.1node",	1, NULL, 0 },
        { "algo.gibbs",	1, NULL, 0 },
        { "algo.m3",	1, NULL, 0 },
        { "algo.ejectabsorb",	1, NULL, 0 },
        { "iterations",	1, NULL, 'i' },
        { "initGT",	0, NULL, 0 },
        { "model.scf",	0, NULL, 0 },
        { "algo.sbm.cem",	0, NULL, 0 },
        { "stringIDs",	0, NULL, 0 },
        { "mega",	0, NULL, 0 },
        { "printEveryNIters",	1, NULL, 0 },
        { "assume_N_nodes",	1, NULL, 0 },
        { "alpha",	1, NULL, 'a' },
        { "save.z",	1, NULL, 'z' },
        { "gamma.s",	1, NULL, 0 },
        { "gamma.phi",	1, NULL, 0 },
        { 0,  0, 0, 0 }
      };

      c = getopt_long (argc, argv, "hVvK:dwsi:a:z:", long_options, &option_index);

      if (c == -1) break;	/* Exit from `while (1)' loop.  */

      switch (c)
        {
        case 'h':	/* Print help and exit.  */
          cmdline_parser_print_help ();
          cmdline_parser_free (&local_args_info);
          exit (EXIT_SUCCESS);

        case 'V':	/* Print version and exit.  */
          cmdline_parser_print_version ();
          cmdline_parser_free (&local_args_info);
          exit (EXIT_SUCCESS);

        case 'v':	/* detailed debugging.  */
        
        
          if (update_arg((void *)&(args_info->verbose_flag), 0, &(args_info->verbose_given),
              &(local_args_info.verbose_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "verbose", 'v',
              additional_error))
            goto failure;
        
          break;
        case 'K':	/* Number of clusters, K.  */
        
        
          if (update_arg( (void *)&(args_info->K_arg), 
               &(args_info->K_orig), &(args_info->K_given),
              &(local_args_info.K_given), optarg, 0, "-1", ARG_INT,
              check_ambiguity, override, 0, 0,
              "K", 'K',
              additional_error))
            goto failure;
        
          break;
        case 'd':	/* directed.  */
        
        
          if (update_arg((void *)&(args_info->directed_flag), 0, &(args_info->directed_given),
              &(local_args_info.directed_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "directed", 'd',
              additional_error))
            goto failure;
        
          break;
        case 'w':	/* weighted.  */
        
        
          if (update_arg((void *)&(args_info->weighted_flag), 0, &(args_info->weighted_given),
              &(local_args_info.weighted_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "weighted", 'w',
              additional_error))
            goto failure;
        
          break;
        case 's':	/* selfloops allowed.  */
        
        
          if (update_arg((void *)&(args_info->selfloop_flag), 0, &(args_info->selfloop_given),
              &(local_args_info.selfloop_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "selfloop", 's',
              additional_error))
            goto failure;
        
          break;
        case 'i':	/* How many iterations.  */
        
        
          if (update_arg( (void *)&(args_info->iterations_arg), 
               &(args_info->iterations_orig), &(args_info->iterations_given),
              &(local_args_info.iterations_given), optarg, 0, "120000", ARG_INT,
              check_ambiguity, override, 0, 0,
              "iterations", 'i',
              additional_error))
            goto failure;
        
          break;
        case 'a':	/* alpha. How uniform the cluster sizes.  */
        
        
          if (update_arg( (void *)&(args_info->alpha_arg), 
               &(args_info->alpha_orig), &(args_info->alpha_given),
              &(local_args_info.alpha_given), optarg, 0, "1", ARG_FLOAT,
              check_ambiguity, override, 0, 0,
              "alpha", 'a',
              additional_error))
            goto failure;
        
          break;
        case 'z':	/* save burnt-in z to this file.  */
        
        
          if (update_arg( (void *)&(args_info->save_z_arg), 
               &(args_info->save_z_orig), &(args_info->save_z_given),
              &(local_args_info.save_z_given), optarg, 0, "", ARG_STRING,
              check_ambiguity, override, 0, 0,
              "save.z", 'z',
              additional_error))
            goto failure;
        
          break;

        case 0:	/* Long option with no short option */
          /* detailed version description.  */
          if (strcmp (long_options[option_index].name, "git-version") == 0)
          {
          
          
            if (update_arg((void *)&(args_info->git_version_flag), 0, &(args_info->git_version_given),
                &(local_args_info.git_version_given), optarg, 0, 0, ARG_FLAG,
                check_ambiguity, override, 1, 0, "git-version", '-',
                additional_error))
              goto failure;
          
          }
          /* seed to drand48().  */
          else if (strcmp (long_options[option_index].name, "seed") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->seed_arg), 
                 &(args_info->seed_orig), &(args_info->seed_given),
                &(local_args_info.seed_given), optarg, 0, "0", ARG_INT,
                check_ambiguity, override, 0, 0,
                "seed", '-',
                additional_error))
              goto failure;
          
          }
          /* The ground truth. a file with N lines. Starts from ZERO..  */
          else if (strcmp (long_options[option_index].name, "GT.vector") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->GT_vector_arg), 
                 &(args_info->GT_vector_orig), &(args_info->GT_vector_given),
                &(local_args_info.GT_vector_given), optarg, 0, 0, ARG_STRING,
                check_ambiguity, override, 0, 0,
                "GT.vector", '-',
                additional_error))
              goto failure;
          
          }
          /* Use the simple Metropolis move on K.  */
          else if (strcmp (long_options[option_index].name, "algo.metroK") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->algo_metroK_arg), 
                 &(args_info->algo_metroK_orig), &(args_info->algo_metroK_given),
                &(local_args_info.algo_metroK_given), optarg, 0, "1", ARG_INT,
                check_ambiguity, override, 0, 0,
                "algo.metroK", '-',
                additional_error))
              goto failure;
          
          }
          /* Use the simple Metropolis on one node.  */
          else if (strcmp (long_options[option_index].name, "algo.1node") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->algo_1node_arg), 
                 &(args_info->algo_1node_orig), &(args_info->algo_1node_given),
                &(local_args_info.algo_1node_given), optarg, 0, "1", ARG_INT,
                check_ambiguity, override, 0, 0,
                "algo.1node", '-',
                additional_error))
              goto failure;
          
          }
          /* Use the simple Gibbs in the algorithm.  */
          else if (strcmp (long_options[option_index].name, "algo.gibbs") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->algo_gibbs_arg), 
                 &(args_info->algo_gibbs_orig), &(args_info->algo_gibbs_given),
                &(local_args_info.algo_gibbs_given), optarg, 0, "1", ARG_INT,
                check_ambiguity, override, 0, 0,
                "algo.gibbs", '-',
                additional_error))
              goto failure;
          
          }
          /* Use M3 in the algorithm.  */
          else if (strcmp (long_options[option_index].name, "algo.m3") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->algo_m3_arg), 
                 &(args_info->algo_m3_orig), &(args_info->algo_m3_given),
                &(local_args_info.algo_m3_given), optarg, 0, "1", ARG_INT,
                check_ambiguity, override, 0, 0,
                "algo.m3", '-',
                additional_error))
              goto failure;
          
          }
          /* Use N+F's eject-absorb move.  */
          else if (strcmp (long_options[option_index].name, "algo.ejectabsorb") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->algo_ejectabsorb_arg), 
                 &(args_info->algo_ejectabsorb_orig), &(args_info->algo_ejectabsorb_given),
                &(local_args_info.algo_ejectabsorb_given), optarg, 0, "1", ARG_INT,
                check_ambiguity, override, 0, 0,
                "algo.ejectabsorb", '-',
                additional_error))
              goto failure;
          
          }
          /* Initialize to the ground truth.  */
          else if (strcmp (long_options[option_index].name, "initGT") == 0)
          {
          
          
            if (update_arg((void *)&(args_info->initGT_flag), 0, &(args_info->initGT_given),
                &(local_args_info.initGT_given), optarg, 0, 0, ARG_FLAG,
                check_ambiguity, override, 1, 0, "initGT", '-',
                additional_error))
              goto failure;
          
          }
          /* Stochastic community finding.  */
          else if (strcmp (long_options[option_index].name, "model.scf") == 0)
          {
          
          
            if (update_arg((void *)&(args_info->model_scf_flag), 0, &(args_info->model_scf_given),
                &(local_args_info.model_scf_given), optarg, 0, 0, ARG_FLAG,
                check_ambiguity, override, 1, 0, "model.scf", '-',
                additional_error))
              goto failure;
          
          }
          /* Classification EM (CEM) for the SBM.  */
          else if (strcmp (long_options[option_index].name, "algo.sbm.cem") == 0)
          {
          
          
            if (update_arg((void *)&(args_info->algo_sbm_cem_flag), 0, &(args_info->algo_sbm_cem_given),
                &(local_args_info.algo_sbm_cem_given), optarg, 0, 0, ARG_FLAG,
                check_ambiguity, override, 1, 0, "algo.sbm.cem", '-',
                additional_error))
              goto failure;
          
          }
          /* string IDs in the input.  */
          else if (strcmp (long_options[option_index].name, "stringIDs") == 0)
          {
          
          
            if (update_arg((void *)&(args_info->stringIDs_flag), 0, &(args_info->stringIDs_given),
                &(local_args_info.stringIDs_given), optarg, 0, 0, ARG_FLAG,
                check_ambiguity, override, 1, 0, "stringIDs", '-',
                additional_error))
              goto failure;
          
          }
          /* dumb down the algorithm for *big* networks.  */
          else if (strcmp (long_options[option_index].name, "mega") == 0)
          {
          
          
            if (update_arg((void *)&(args_info->mega_flag), 0, &(args_info->mega_given),
                &(local_args_info.mega_given), optarg, 0, 0, ARG_FLAG,
                check_ambiguity, override, 1, 0, "mega", '-',
                additional_error))
              goto failure;
          
          }
          /* How often to print an update.  */
          else if (strcmp (long_options[option_index].name, "printEveryNIters") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->printEveryNIters_arg), 
                 &(args_info->printEveryNIters_orig), &(args_info->printEveryNIters_given),
                &(local_args_info.printEveryNIters_given), optarg, 0, "10", ARG_INT,
                check_ambiguity, override, 0, 0,
                "printEveryNIters", '-',
                additional_error))
              goto failure;
          
          }
          /* Pre-create N nodes (0 to N-1), which may be left with zero degree.  */
          else if (strcmp (long_options[option_index].name, "assume_N_nodes") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->assume_N_nodes_arg), 
                 &(args_info->assume_N_nodes_orig), &(args_info->assume_N_nodes_given),
                &(local_args_info.assume_N_nodes_given), optarg, 0, "0", ARG_INT,
                check_ambiguity, override, 0, 0,
                "assume_N_nodes", '-',
                additional_error))
              goto failure;
          
          }
          /* (for weighted only). Shape of Gamma prior.  */
          else if (strcmp (long_options[option_index].name, "gamma.s") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->gamma_s_arg), 
                 &(args_info->gamma_s_orig), &(args_info->gamma_s_given),
                &(local_args_info.gamma_s_given), optarg, 0, "1", ARG_FLOAT,
                check_ambiguity, override, 0, 0,
                "gamma.s", '-',
                additional_error))
              goto failure;
          
          }
          /* (for weighted only). Scale of Gamma prior.  */
          else if (strcmp (long_options[option_index].name, "gamma.phi") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->gamma_phi_arg), 
                 &(args_info->gamma_phi_orig), &(args_info->gamma_phi_given),
                &(local_args_info.gamma_phi_given), optarg, 0, "1", ARG_FLOAT,
                check_ambiguity, override, 0, 0,
                "gamma.phi", '-',
                additional_error))
              goto failure;
          
          }
          
          break;
        case '?':	/* Invalid option.  */
          /* `getopt_long' already printed an error message.  */
          goto failure;

        default:	/* bug: option not considered.  */
          fprintf (stderr, "%s: option unknown: %c%s\n", CMDLINE_PARSER_PACKAGE, c, (additional_error ? additional_error : ""));
          abort ();
        } /* switch */
    } /* while */




  cmdline_parser_release (&local_args_info);

  if ( error )
    return (EXIT_FAILURE);

  if (optind < argc)
    {
      int i = 0 ;
      int found_prog_name = 0;
      /* whether program name, i.e., argv[0], is in the remaining args
         (this may happen with some implementations of getopt,
          but surely not with the one included by gengetopt) */

      i = optind;
      while (i < argc)
        if (argv[i++] == argv[0]) {
          found_prog_name = 1;
          break;
        }
      i = 0;

      args_info->inputs_num = argc - optind - found_prog_name;
      args_info->inputs =
        (char **)(malloc ((args_info->inputs_num)*sizeof(char *))) ;
      while (optind < argc)
        if (argv[optind++] != argv[0])
          args_info->inputs[ i++ ] = gengetopt_strdup (argv[optind-1]) ;
    }

  return 0;

failure:
  
  cmdline_parser_release (&local_args_info);
  return (EXIT_FAILURE);
}
