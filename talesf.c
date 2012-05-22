/*

Copyright (c) 2011-2012, Daniel S. Standage <daniel.standage@gmail.com> and
Erin Doyle <edoyle@iastate.edu>. See README for license details.

For compiling, try this.

    gcc -Wall -O3 -o talesf *.c -lz -lm

*/

// System libraries
#include <getopt.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <zlib.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdarg.h>

// Include my libraries
#include "Array.h"
#include "Hashmap.h"

// Initialize the kseq library
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

typedef struct
{
  int strand;
  char *sequence;
  char *sequence_name;
  unsigned long index;
  double score;
} BindingSite;

/*
 * Utility
 */

// Allocate memory for an array of 4 doubles
double *double_array(double a, double c, double g, double t)
{
  double *array = malloc(sizeof(double)*4);
  array[0] = a;
  array[1] = c;
  array[2] = g;
  array[3] = t;
  return array;
}

// Allocate memory for an array of 4 integers
int *int_array(int a, int c, int g, int t)
{
  int *array = malloc(sizeof(int)*4);
  array[0] = a;
  array[1] = c;
  array[2] = g;
  array[3] = t;
  return array;
}

// Compute diresidue probabilities from observed counts
Hashmap *get_diresidue_probabilities(Array *rvdseq, double w)
{
  char **diresidues;
  Hashmap *diresidue_counts, *diresidue_probabilities;
  int i, num_diresidues;

  // First, use a hashmap to store counts
  diresidue_counts = hashmap_new(64);

  // Add known counts
  // The int_array provides counts in the order A C G T
  hashmap_add(diresidue_counts, "HD", int_array( 7, 99,  0,  1));
  hashmap_add(diresidue_counts, "NI", int_array(58,  6,  0,  0));
  hashmap_add(diresidue_counts, "NG", int_array( 6,  6,  1, 57));
  hashmap_add(diresidue_counts, "NN", int_array(21,  8, 26,  2));
  hashmap_add(diresidue_counts, "NS", int_array(20,  6,  4,  0));
  hashmap_add(diresidue_counts, "N*", int_array( 1, 11,  1,  7));
  hashmap_add(diresidue_counts, "HG", int_array( 1,  2,  0, 15));
  hashmap_add(diresidue_counts, "HA", int_array( 1,  4,  1,  0));
  hashmap_add(diresidue_counts, "ND", int_array( 0,  4,  0,  0));
  hashmap_add(diresidue_counts, "NK", int_array( 0,  0,  2,  0));
  hashmap_add(diresidue_counts, "HI", int_array( 0,  1,  0,  0));
  hashmap_add(diresidue_counts, "HN", int_array( 0,  0,  1,  0));
  hashmap_add(diresidue_counts, "NA", int_array( 0,  0,  1,  0));
  hashmap_add(diresidue_counts, "IG", int_array( 0,  0,  0,  1));
  hashmap_add(diresidue_counts, "H*", int_array( 0,  0,  0,  1));

  // Add unknown counts
  hashmap_add(diresidue_counts, "S*", int_array(1, 1, 1, 1));
  hashmap_add(diresidue_counts, "NH", int_array(1, 1, 1, 1));
  hashmap_add(diresidue_counts, "YG", int_array(1, 1, 1, 1));
  hashmap_add(diresidue_counts, "SN", int_array(1, 1, 1, 1));
  hashmap_add(diresidue_counts, "SS", int_array(1, 1, 1, 1));
  hashmap_add(diresidue_counts, "NC", int_array(1, 1, 1, 1));
  hashmap_add(diresidue_counts, "HH", int_array(1, 1, 1, 1));

  // Add any other RVDs in the sequence with equal weight
  for(i = 0; i < array_size(rvdseq); i++)
  {
    char *rvd = array_get(rvdseq, i);
    // hashmap_add doesn't do anything for existing key, only adds new ones
    hashmap_add(diresidue_counts, rvd, int_array(1, 1, 1, 1));
  }

  diresidue_probabilities = hashmap_new(64);
  diresidues = hashmap_keys(diresidue_counts);
  num_diresidues = hashmap_size(diresidue_counts);
  for(i = 0; i < num_diresidues; i++)
  {
    double p[4];
    int j;
    int *counts = hashmap_get(diresidue_counts, diresidues[i]);
    int sum = counts[0] + counts[1] + counts[2] + counts[3];

    for(j = 0; j < 4; j++)
    {
      p[j] = w * (double)counts[j] / (double)sum + (1.0 - w) / 4.0;
    }

    hashmap_add(diresidue_probabilities, diresidues[i], double_array(p[0], p[1], p[2], p[3]));
  }

  free(diresidues);
  hashmap_delete(diresidue_counts, free);
  return diresidue_probabilities;
}

// Convert diresidue probabilities into scores
Hashmap *convert_probabilities_to_scores(Hashmap *diresidue_probabilities)
{
  char **diresidues;
  int i, j, num_diresidues;
  Hashmap *diresidue_scores = hashmap_new(64);

  diresidues = hashmap_keys(diresidue_probabilities);
  num_diresidues = hashmap_size(diresidue_probabilities);
  for(i = 0; i < num_diresidues; i++)
  {
    double *probs = hashmap_get(diresidue_probabilities, diresidues[i]);
    for(j = 0; j < 4; j++)
      probs[j] = -1.0 * log(probs[j]);

    hashmap_add(diresidue_scores, diresidues[i], probs);
  }

  free(diresidues);
  return diresidue_scores;
}

// Given a sequence of RVDs, calculate the optimal score
double get_best_score(Array *rvdseq, Hashmap *rvdscores)
{
  double best_score = 0.0;
  int i, j;

  for(i = 0; i < array_size(rvdseq); i++)
  {
    char *rvd = array_get(rvdseq, i);
    double *scores = hashmap_get(rvdscores, rvd);
    double min_score = -1.0;
    if(scores == NULL)
      return -1.0;

    for(j = 0; j < 4; j++)
    {
      if(j == 0 || scores[j] < min_score)
        min_score = scores[j];
    }
    best_score += min_score;
  }

  return best_score;
}

void logger(FILE* log_file, char* message, ...) {

  time_t rawtime;
  va_list args;

  time(&rawtime);

  char* timestamp_message = calloc(strlen(message) + 256, sizeof(char));
  char* ctime_string = ctime(&rawtime);

  strcpy(timestamp_message, "[");
  strncat(timestamp_message, ctime_string, strlen(ctime_string) - 1);
  strcat(timestamp_message, "] ");
  strcat(timestamp_message, message);
  strcat(timestamp_message, "\n");

  va_start(args, message);
  vfprintf(log_file, timestamp_message, args);

  va_end(args);

  fflush(log_file);
  free(timestamp_message);

}

void create_options_string(char *options_str, char *rvd_str, Hashmap *prog_kwargs) {

  char cutoff_str[32];

  double cutoff = *((double *) hashmap_get(prog_kwargs, "cutoff"));
  int forward_only = *((int *) hashmap_get(prog_kwargs, "forward_only"));
  int c_upstream = *((int *) hashmap_get(prog_kwargs, "c_upstream"));

  strcat(options_str, "options_used:");

  if (!forward_only) {
    strcat(options_str, "search reverse complement, ");
  }

  strcat(options_str, "upstream_base = ");

  if (c_upstream != 1) {
    strcat(options_str, "T ");
  }

  if (c_upstream != 0) {
    strcat(options_str, "C ");
  }

  sprintf(cutoff_str, ", cutoff = %.2lf, ", cutoff);
  strcat(options_str, cutoff_str);

  sprintf(options_str, "RVDS = %s", rvd_str);

}

/*
 * Core
 */

int binding_site_compare_score(const void * a, const void * b)
{
  BindingSite *real_a = *((BindingSite **)a);
  BindingSite *real_b = *((BindingSite **)b);

  double real_a_score = floorf(real_a->score * 100 + 0.5) / 100;
  double real_b_score = floorf(real_b->score * 100 + 0.5) / 100;

  double score_diff = (real_a_score - real_b_score);

  if(score_diff < 0) {
    return -1;
  } else if(score_diff > 0) {
    return 1;
  } else {
    return 0;
  }

}

int print_results(Array *results, Array *rvd_seq, Hashmap *prog_kwargs, double best_score, char *output_filepath, FILE *log_file) {

  int forward_only = *((int *) hashmap_get(prog_kwargs, "forward_only"));
  char *organism_name = hashmap_get(prog_kwargs, "organism_name");

  char *source_str = "TALESF";
  char options_str[512];

  int num_rvds = array_size(rvd_seq);
  char *plus_strand_sequence;
  FILE *gff_out_file = NULL;
  FILE *tab_out_file = NULL;
  FILE *genome_browser_file = NULL;

  size_t output_filepath_length;
  char* temp_output_filepath;

  char *rvd_str = calloc(3 * num_rvds, sizeof(char));

  int is_genome = (*organism_name != '\0');
  int genome_using_gbrowse = (is_genome && (strcmp(organism_name, "oryza_sativa") == 0 || strcmp(organism_name, "arabidopsis_thaliana") == 0));

  int i;

  for(i = 0; i < num_rvds; i++)
  {
    char *rvd = array_get(rvd_seq, i);
    if(i == 0)
      strncpy(rvd_str, rvd, 2);
    else
      sprintf(rvd_str + (i*3)-1, "_%s", rvd);
  }
  rvd_str[3*num_rvds - 1] = '\0';

  create_options_string(options_str, rvd_str, prog_kwargs);

  if(output_filepath == NULL) {

    gff_out_file = stdout;

  } else {


    output_filepath_length = strlen(output_filepath) + 5;
    temp_output_filepath = calloc(output_filepath_length + 1, sizeof(char));

    sprintf(temp_output_filepath, "%s.txt", output_filepath);
    tab_out_file = fopen(temp_output_filepath, "w");
    memset(temp_output_filepath, '\0', output_filepath_length);
    sprintf(temp_output_filepath, "%s.gff3", output_filepath);
    gff_out_file = fopen(temp_output_filepath, "w");

    if(is_genome) {
      memset(temp_output_filepath, '\0', output_filepath_length);

      if(genome_using_gbrowse) {
        sprintf(temp_output_filepath, "%s.gff", output_filepath);
      } else {
        sprintf(temp_output_filepath, "%s.bed", output_filepath);
      }

      genome_browser_file = fopen(temp_output_filepath, "w");


    }

    free(temp_output_filepath);

  }

  if(!gff_out_file || !tab_out_file || (is_genome && !genome_browser_file))
  {
    fprintf(log_file, "Error: unable to open output file '%s'\n", output_filepath);
    return 1;
  }

  // Tab file header
  if (forward_only)
  {
    fprintf(tab_out_file, "table_ignores:Plus strand sequence\n");
  }

  fprintf(tab_out_file, options_str);

  fprintf(tab_out_file, "Best Possible Score:%.2lf\n", best_score);
  fprintf(tab_out_file, "Sequence Name\tStrand\tScore\tStart Position\tTarget Sequence\tPlus strand sequence\n");

  // GFF file header
  fprintf(gff_out_file, "##gff-version 3\n");

  if (forward_only)
  {
    fprintf(gff_out_file, "#table_display_tags:target_sequence\n");
  }
  else
  {
    fprintf(gff_out_file, "#table_display_tags:target_sequence,plus_strand_sequence\n");
  }

  fprintf(gff_out_file, "#%s", options_str);

  fprintf(gff_out_file, "#Best Possible Score:%.2lf\n", best_score);

  // bed file header

  if(genome_using_gbrowse) {
    fprintf(genome_browser_file, "##gff-version 3\n");
  }
  else if(is_genome) {
    fprintf(genome_browser_file, "track name=\"TAL Targets\" description=\"Targets for RVD sequence %s\" visibility=2 useScore=1\n", rvd_str);
  }

  for(i = 0; i < array_size(results); i++)
  {

    BindingSite *site = (BindingSite *)array_get(results, i);
    char *sequence = site->sequence;
    char strand = '+';
    char *tab_strand = "Plus";

    if(site->strand > 0)
      plus_strand_sequence = sequence;
    else
    {
      int j;
      int seq_len = num_rvds + 2;


      plus_strand_sequence = sequence;
      sequence = malloc(sizeof(char)*(seq_len+1));
      sequence[seq_len] = '\0';

      for(j = 0; j < seq_len; j++)
      {
        char base = site->sequence[seq_len - j - 1];
        if(base == 'A' || base == 'a')
          sequence[j] = 'T';
        else if(base == 'C' || base == 'c')
          sequence[j] = 'G';
        else if(base == 'G' || base == 'g')
          sequence[j] = 'C';
        else if(base == 'T' || base == 't')
          sequence[j] = 'A';
        else if(base == ' ')
          sequence[j] = ' ';
        else
        {
          fprintf(stderr, "Error: unexpected character '%c'\n", base);
          exit(1);
        }
      }
      strand = '-';
      tab_strand = "Minus";
    }

    fprintf( tab_out_file, "%s\t%s\t%.2lf\t%lu\t%s\t%s\n",
             site->sequence_name, tab_strand, site->score, site->index + 1, sequence, plus_strand_sequence);

    fprintf( gff_out_file, "%s\t%s\t%s\t%lu\t%lu\t%.2lf\t%c\t.\trvd_sequence=%s;target_sequence=%s;plus_strand_sequence=%s;\n",
             site->sequence_name, source_str, "TAL_effector_binding_site", site->index + 1,
             site->index + num_rvds, site->score, strand, rvd_str, sequence, plus_strand_sequence);

    if(is_genome && i < 10000) {

      if(genome_using_gbrowse) {

        fprintf( genome_browser_file, "chr%s\t%s\t%s\t%lu\t%lu\t%.2lf\t%c\t.\tName=site%d;\n",
                 site->sequence_name, source_str, "TAL_effector_binding_site", site->index + 1,
                 site->index + num_rvds, site->score, strand, i);

      } else {

        int bed_score = floorf((best_score / site->score * 1000) + 0.5);
        fprintf( genome_browser_file,"%s\t%lu\t%lu\tsite%d\t%d\t%c\n",
                 site->sequence_name, site->index, site->index + num_rvds - 1, i, bed_score, strand);

      }

    }

    if(plus_strand_sequence != sequence) {
      free(sequence);
    }

  }

  free(rvd_str);
  fclose(gff_out_file);
  fclose(tab_out_file);

  if(is_genome) {
    fclose(genome_browser_file);
  }

  return 0;

}

// Identify and print out TAL effector binding sites
void find_binding_sites(kseq_t *seq, Array *rvdseq, Hashmap *diresidue_scores, double cutoff, int forwardonly, int c_upstream, Array *results)
{
  unsigned long i, j;
  int num_rvds = array_size(rvdseq);
  int seq_name_len;

  if(num_rvds > seq->seq.l)
  {
    fprintf(stderr, "Warning: skipping sequence '%s' since it is shorter than the RVD sequence\n", seq->seq.s);
    return;
  }

  seq_name_len = strlen(seq->name.s);

  for(i = 1; i <= seq->seq.l - num_rvds; i++)
  {
    if((c_upstream != 0 && (seq->seq.s[i-1] == 'C' || seq->seq.s[i-1] == 'c')) || (c_upstream != 1 && (seq->seq.s[i-1] == 'T' || seq->seq.s[i-1] == 't')))
    {
      double cumscore = 0.0;
      for(j = 0; j < num_rvds; j++)
      {
        char *rvd = array_get(rvdseq, j);
        double *scores = hashmap_get(diresidue_scores, rvd);

        if(seq->seq.s[i+j] == 'A' || seq->seq.s[i+j] == 'a')
          cumscore += scores[0];
        else if(seq->seq.s[i+j] == 'C' || seq->seq.s[i+j] == 'c')
          cumscore += scores[1];
        else if(seq->seq.s[i+j] == 'G' || seq->seq.s[i+j] == 'g')
          cumscore += scores[2];
        else if(seq->seq.s[i+j] == 'T' || seq->seq.s[i+j] == 't')
          cumscore += scores[3];
        else
          cumscore += cutoff + 1;

        if(cumscore > cutoff)
          break;
      }

      if(cumscore <= cutoff)
      {

        BindingSite *site = malloc(sizeof(BindingSite));

        site->sequence = calloc(num_rvds + 2 + 1, sizeof(char));
        site->sequence_name = calloc(seq_name_len + 1, sizeof(char));
        site->sequence[num_rvds + 2] = '\0';
        site->sequence_name[seq_name_len] = '\0';

        site->strand = 1;
        site->index = i;

        strncpy(site->sequence, seq->seq.s + site->index - 1, 1);
        site->sequence[1] = ' ';
        strncpy(site->sequence + 2, seq->seq.s + site->index, num_rvds);
        strncpy(site->sequence_name, seq->name.s, seq_name_len);

        site->score = cumscore;

        #pragma omp critical (add_result)
        array_add(results, site);

      }
    }

    if(!forwardonly)
    {
      if((c_upstream != 0 && (seq->seq.s[i + num_rvds - 1] == 'G' || seq->seq.s[i + num_rvds - 1] == 'g')) || (c_upstream != 1 && (seq->seq.s[i + num_rvds - 1] == 'A' || seq->seq.s[i + num_rvds - 1] == 'a')))
      {
        double cumscore = 0.0;
        for(j = 0; j < num_rvds; j++)
        {
          char *rvd = array_get(rvdseq, num_rvds - j - 1);
          double *scores = hashmap_get(diresidue_scores, rvd);

          if(seq->seq.s[i+j-1] == 'A' || seq->seq.s[i+j-1] == 'a')
            cumscore += scores[3];
          else if(seq->seq.s[i+j-1] == 'C' || seq->seq.s[i+j-1] == 'c')
            cumscore += scores[2];
          else if(seq->seq.s[i+j-1] == 'G' || seq->seq.s[i+j-1] == 'g')
            cumscore += scores[1];
          else if(seq->seq.s[i+j-1] == 'T' || seq->seq.s[i+j-1] == 't')
            cumscore += scores[0];
          else
            cumscore += cutoff + 1;

          if(cumscore > cutoff)
            break;
        }

        if(cumscore <= cutoff)
        {

          BindingSite *site = malloc(sizeof(BindingSite));

          site->sequence = calloc(num_rvds + 2 + 1, sizeof(char));
          site->sequence_name = calloc(seq_name_len + 1, sizeof(char));
          site->sequence[num_rvds + 2] = '\0';
          site->sequence_name[seq_name_len] = '\0';

          site->strand = -1;
          site->index = i - 1;

          strncpy(site->sequence, seq->seq.s + site->index, num_rvds);
          site->sequence[num_rvds] = ' ';
          strncpy(site->sequence + num_rvds + 1, seq->seq.s + site->index + num_rvds, 1);
          strncpy(site->sequence_name, seq->name.s, seq_name_len);

          site->score = cumscore;

          #pragma omp critical (add_result)
          array_add(results, site);

        }
      }
    }
  }

}

int run_talesf_task(Hashmap *kwargs) {

  // Options
  char *seq_filename = hashmap_get(kwargs, "seq_filename");
  char *rvd_string = hashmap_get(kwargs, "rvd_string");
  char *output_filepath = hashmap_get(kwargs, "output_filepath");
  char *log_filepath = hashmap_get(kwargs, "log_filepath");

  double weight = *((double *) hashmap_get(kwargs, "weight"));
  double cutoff = *((double *) hashmap_get(kwargs, "cutoff"));

  int forward_only = *((int *) hashmap_get(kwargs, "forward_only"));
  int c_upstream = *((int *) hashmap_get(kwargs, "c_upstream"));
  int numprocs = *((int *) hashmap_get(kwargs, "num_procs"));

  // Program variable domain
  Array *rvd_seq;
  char *tok, cmd[256], line[32];
  gzFile seqfile;
  kseq_t *seq;
  int i, j, seq_num;

  FILE *log_file = stdout;

  if(log_filepath && strcmp(log_filepath, "NA") != 0) {
    log_file = fopen(log_filepath, "a");
  }

  if(!seq_filename || !rvd_string || !output_filepath) {
    logger(log_file, "Error: One or more arguments to run_talesf_task was null");
    return 1;
  }

  Array *results = array_new( sizeof(BindingSite *) );

  rvd_seq = array_new( sizeof(char *) );
  tok = strtok(rvd_string, " _");
  while(tok != NULL)
  {
    char *r = strdup(tok);
    array_add(rvd_seq, r);
    tok = strtok(NULL, " _");
  }

  // Get RVD/bp matching scores
  Hashmap *diresidue_probabilities = get_diresidue_probabilities(rvd_seq, weight);
  Hashmap *diresidue_scores = convert_probabilities_to_scores(diresidue_probabilities);
  hashmap_delete(diresidue_probabilities, NULL);

  // Compute optimal score for this RVD sequence
  double best_score = get_best_score(rvd_seq, diresidue_scores);

  // Determine number of sequences in file
  sprintf(cmd, "grep '^>' %s | wc -l", seq_filename);
  FILE *in = popen(cmd, "r");
  if(!in)
  {
    perror("Error: unable to check fasta file size\n");
    logger(log_file, "Error: unable to check fasta file size");
    return 1;
  }
  fgets(line, sizeof(line), in);
  pclose(in);
  seq_num = atoi(line);

  // Begin processing

  int abort = 0;

  omp_set_num_threads(numprocs);
  #pragma omp parallel private(i, j, seq, seqfile)
  {

    // Open genomic sequence file
    seqfile = gzopen(seq_filename, "r");
    if(!seqfile)
    {

      logger(log_file, "Error: unable to open sequence '%s'", seq_filename);
      abort = 1;

    } else {

      seq = kseq_init(seqfile);

      j = 0;
      #pragma omp for schedule(static)
      for(i = 0; i < seq_num; i++)
      {

        #pragma omp flush (abort)
        if (!abort) {
          while(j <= i)
          {
            int result = kseq_read(seq);
            if(result < 0)
            {
              logger(log_file, "Error: problem parsing data from '%s'", seq_filename);
              abort = 1;
            }
            j++;
          }

          logger(log_file, "Scanning %s for binding sites (length %ld)", seq->name.s, seq->seq.l);
          find_binding_sites(seq, rvd_seq, diresidue_scores, best_score * cutoff, forward_only, c_upstream, results);

        }

      }

      kseq_destroy(seq);
      gzclose(seqfile);

    }

  }

  if(!abort) {

    qsort(results->data, array_size(results), sizeof(BindingSite *), binding_site_compare_score);

    int print_results_result;

    print_results_result = print_results(results, rvd_seq, kwargs, best_score, output_filepath, log_file);

    logger(log_file, "Finished");

    if(print_results_result == 1) {
      abort = 1;
    }

  }

  // Free memory

  if(results) {

    for(i = 0; i < array_size(results); i++)
    {
      BindingSite *site = (BindingSite *)array_get(results, i);

      free(site->sequence);
      free(site->sequence_name);
      free(site);

    }

    array_delete(results);

  }


  if(rvd_seq) {

    for(i = 0; i < array_size(rvd_seq); i++)
    {
      char *temp = (char *)array_get(rvd_seq, i);
      free(temp);
    }


    array_delete(rvd_seq);

  }

  if(diresidue_scores) {
    hashmap_delete(diresidue_scores, free);
  }

  if(log_file != stdout) {
    fclose(log_file);
  }

  if(abort) {
    return 1;
  } else {
    return 0;
  }
}
