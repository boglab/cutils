/*

Copyright (c) 2011-2012, Daniel S. Standage <daniel.standage@gmail.com> and
Erin Doyle <edoyle@iastate.edu> with modifications by Nick Booher <njbooher@gmail.com>.
See README for license details.

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

typedef struct {
  char *sequence[2];
  char *sequence_name;
  unsigned long indexes[2];
  double scores[2];
  int spacer_length;
  int f_idx;
  int r_idx;
} BindingSite;

Hashmap *talesf_kwargs;

/*
 * Utility
 */

// Allocate memory for an array of 4 doubles
double *double_array(double a, double c, double g, double t) {
  double *array = malloc(sizeof(double)*4);
  array[0] = a;
  array[1] = c;
  array[2] = g;
  array[3] = t;
  return array;
}

// Allocate memory for an array of 4 integers
int *int_array(int a, int c, int g, int t) {
  int *array = malloc(sizeof(int)*4);
  array[0] = a;
  array[1] = c;
  array[2] = g;
  array[3] = t;
  return array;
}

// Compute diresidue probabilities from observed counts
Hashmap *get_diresidue_probabilities(Array *rvdseq, double w) {
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
  for (i = 0; i < array_size(rvdseq); i++)
  {
    char *rvd = array_get(rvdseq, i);
    if (hashmap_get(diresidue_counts, rvd) == NULL) {
      hashmap_add(diresidue_counts, rvd, int_array(1, 1, 1, 1));
    }
  }

  diresidue_probabilities = hashmap_new(64);
  diresidues = hashmap_keys(diresidue_counts);
  num_diresidues = hashmap_size(diresidue_counts);
  for (i = 0; i < num_diresidues; i++)
  {
    double p[4];
    int j;
    int *counts = hashmap_get(diresidue_counts, diresidues[i]);
    int sum = counts[0] + counts[1] + counts[2] + counts[3];

    for (j = 0; j < 4; j++)
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
Hashmap *convert_probabilities_to_scores(Hashmap *diresidue_probabilities) {
  char **diresidues;
  int i, j, num_diresidues;
  Hashmap *diresidue_scores = hashmap_new(64);

  diresidues = hashmap_keys(diresidue_probabilities);
  num_diresidues = hashmap_size(diresidue_probabilities);
  for (i = 0; i < num_diresidues; i++)
  {
    double *probs = hashmap_get(diresidue_probabilities, diresidues[i]);
    for (j = 0; j < 4; j++)
      probs[j] = -1.0 * log(probs[j]);

    hashmap_add(diresidue_scores, diresidues[i], probs);
  }

  free(diresidues);
  return diresidue_scores;
}

// Given a sequence of RVDs, calculate the optimal score
double get_best_score(Array *rvdseq, Hashmap *rvdscores) {
  double best_score = 0.0;
  int i, j;

  for (i = 0; i < array_size(rvdseq); i++)
  {
    char *rvd = array_get(rvdseq, i);
    double *scores = hashmap_get(rvdscores, rvd);
    double min_score = -1.0;
    if (scores == NULL)
      return -1.0;

    for (j = 0; j < 4; j++)
    {
      if (j == 0 || scores[j] < min_score)
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

void create_options_string(char *options_str, char *rvd_str) {

  char cutoff_str[32];
  char rvds_eq_str[64];

  double cutoff = *((double *) hashmap_get(talesf_kwargs, "cutoff"));
  int forward_only = *((int *) hashmap_get(talesf_kwargs, "forward_only"));
  int c_upstream = *((int *) hashmap_get(talesf_kwargs, "c_upstream"));

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

  sprintf(rvds_eq_str, "rvd_sequence = %s", rvd_str);
  strcat(options_str, rvds_eq_str);

  strcat(options_str, "\n");

}

Array *rvd_string_to_array(char *rvd_string) {

  char *tok;
  Array *rvd_seq;

  rvd_seq = array_new( sizeof(char *) );
  tok = strtok(rvd_string, " _");

  while (tok != NULL)
  {
    char *r = strdup(tok);
    array_add(rvd_seq, r);
    tok = strtok(NULL, " _");
  }

  return rvd_seq;

}

void str_replace(char *haystack, char *needle, char *replacement) {

  char *pos = strstr(haystack, needle);

  while (pos != NULL) {

    strcpy(pos, replacement);
    pos = strstr(haystack, needle);

  }

}


/*
 * Core
 */

int binding_site_compare_score(const void * a, const void * b) {
  BindingSite *real_a = *((BindingSite **)a);
  BindingSite *real_b = *((BindingSite **)b);

  int str_cmp_result = strcmp(real_a->sequence_name, real_b->sequence_name);

  if (str_cmp_result != 0) {

    return str_cmp_result;

  } else {

    double real_a_score = floorf(real_a->scores[0] * 100 + 0.5) / 100;
    double real_b_score = floorf(real_b->scores[0] * 100 + 0.5) / 100;

    double score_diff = (real_a_score - real_b_score);

    if (score_diff < 0) {
      return -1;
    } else if (score_diff > 0) {
      return 1;
    } else {

      double real_a_score2 = floorf(real_a->scores[1] * 100 + 0.5) / 100;
      double real_b_score2 = floorf(real_b->scores[1] * 100 + 0.5) / 100;
      double score_diff2 = (real_a_score2 - real_b_score2);

      if (score_diff2 < 0) {
        return -1;
      } else if (score_diff2 > 0) {
        return 1;
      } else {
        return 0;
      }

    }

  }

}

int binding_site_compare_pos(const void * a, const void * b) {

  BindingSite *real_a = *((BindingSite **)a);
  BindingSite *real_b = *((BindingSite **)b);

  int str_cmp_result = strcmp(real_a->sequence_name, real_b->sequence_name);
  long pos_diff = real_a->indexes[0] - real_b->indexes[0];
  long pos_diff2 = real_a->indexes[1] - real_b->indexes[1];

  if (str_cmp_result != 0) {

    return str_cmp_result;

  }
  else if (real_a->f_idx < real_b->f_idx) {

    return -1;

  } else if (real_a->r_idx < real_b->r_idx) {

    return -1;

  } else {

    if (pos_diff < 0) {
      return -1;
    } else if (pos_diff > 0) {
      return 1;
    } else {

      if (pos_diff2 < 0) {
        return -1;
      } else if (pos_diff2 > 0) {
        return 1;
      } else {
        return 0;
      }

    }

  }

}

int print_results(Array *results, Array **rvd_seqs, double best_score, double best_score2, FILE *log_file) {

  char *output_filepath = hashmap_get(talesf_kwargs, "output_filepath");

  char options_str[512];

  // strcat doesn't seem to work unless you do this
  *options_str = '\0';

  FILE *tab_out_file = NULL;

  size_t output_filepath_length;
  char* temp_output_filepath;

  char *rvd_string_printable = strdup(hashmap_get(talesf_kwargs, "rvd_string"));
  str_replace(rvd_string_printable, " ", "_");

  char *rvd_string2_printable = strdup(hashmap_get(talesf_kwargs, "rvd_string2"));
  str_replace(rvd_string2_printable, " ", "_");

  //create_options_string(options_str, rvd_str);

  output_filepath_length = strlen(output_filepath) + 5;
  temp_output_filepath = calloc(output_filepath_length + 1, sizeof(char));

  sprintf(temp_output_filepath, "%s.txt", output_filepath);
  tab_out_file = fopen(temp_output_filepath, "w");

  free(temp_output_filepath);

  if (!tab_out_file) {
    fprintf(log_file, "Error: unable to open output file '%s'\n", output_filepath);
    return 1;
  }

  //fprintf(tab_out_file, options_str);

  fprintf(tab_out_file, "TAL1 Best Possible Score:%.2lf\n", best_score);
  fprintf(tab_out_file, "TAL2 Best Possible Score:%.2lf\n", best_score2);
  fprintf(tab_out_file, "Sequence Name\tTAL 1\tTAL 2\tTAL 1 Score\tTAL 2 Score\tTAL 1 Start\tTAL 2 Start\tSpacer Length\tTAL 1 Target\tTAL 2 Target\n");

  int tal2_seq_len = array_size(rvd_seqs[1]) + 2;
  char *tal2_sequence = calloc(tal2_seq_len + 1, sizeof(char));

  for (unsigned long i = 0; i < array_size(results); i++)
  {

    BindingSite *site = (BindingSite *) array_get(results, i);

    char tal1_name[16];
    char tal2_name[16];

    sprintf(tal1_name, "TAL%d", site->f_idx + 1);
    sprintf(tal2_name, "TAL%d", site->r_idx + 1);

    for (int j = 0; j < tal2_seq_len; j++)
    {
      char base = site->sequence[1][tal2_seq_len - j - 1];
      if (base == 'A' || base == 'a')
        tal2_sequence[j] = 'T';
      else if (base == 'C' || base == 'c')
        tal2_sequence[j] = 'G';
      else if (base == 'G' || base == 'g')
        tal2_sequence[j] = 'C';
      else if (base == 'T' || base == 't')
        tal2_sequence[j] = 'A';
      else if (base == ' ')
        tal2_sequence[j] = ' ';
      else
      {
        fprintf(stderr, "Error: unexpected character '%c'\n", base);
        exit(1);
      }
    }

    fprintf( tab_out_file, "%s\t%s\t%s\t%.2lf\t%.2lf\t%lu\t%lu\t%d\t%s\t%s\n",
             site->sequence_name, tal1_name, tal2_name, site->scores[0], site->scores[1], site->indexes[0], site->indexes[1], site->spacer_length, site->sequence[0], tal2_sequence);


  }

  free(tal2_sequence);
  free(rvd_string_printable);
  free(rvd_string2_printable);
  fclose(tab_out_file);

  return 0;

}

double score_binding_site(kseq_t *seq, unsigned long i, Array *rvd_seq, Hashmap *diresidue_scores, double cutoff, int reverse) {

  double total_score = 0.0;

  if (!reverse) {

    for (unsigned long j = 0; j < array_size(rvd_seq); j++) {

      char *rvd = array_get(rvd_seq, j);
      double *scores = hashmap_get(diresidue_scores, rvd);

      if (seq->seq.s[i+j] == 'A' || seq->seq.s[i+j] == 'a')
        total_score += scores[0];
      else if (seq->seq.s[i+j] == 'C' || seq->seq.s[i+j] == 'c')
        total_score += scores[1];
      else if (seq->seq.s[i+j] == 'G' || seq->seq.s[i+j] == 'g')
        total_score += scores[2];
      else if (seq->seq.s[i+j] == 'T' || seq->seq.s[i+j] == 't')
        total_score += scores[3];
      else
        total_score += cutoff + 1;

      if (total_score > cutoff)
        break;

    }

  } else {

    for (unsigned long j = 0; j < array_size(rvd_seq); j++) {

      char *rvd = array_get(rvd_seq, array_size(rvd_seq) - j - 1);
      double *scores = hashmap_get(diresidue_scores, rvd);

      if (seq->seq.s[i+j] == 'A' || seq->seq.s[i+j-1] == 'a')
        total_score += scores[3];
      else if (seq->seq.s[i+j] == 'C' || seq->seq.s[i+j-1] == 'c')
        total_score += scores[2];
      else if (seq->seq.s[i+j] == 'G' || seq->seq.s[i+j-1] == 'g')
        total_score += scores[1];
      else if (seq->seq.s[i+j] == 'T' || seq->seq.s[i+j-1] == 't')
        total_score += scores[0];
      else
        total_score += cutoff + 1;

      if (total_score > cutoff)
        break;
    }

  }

  return total_score;

}

BindingSite *create_binding_site(kseq_t *seq, unsigned long i, unsigned long j, int num_forward_rvds, double forward_score, int num_reverse_rvds, double reverse_score, int spacer_size, int f_idx, int r_idx) {

  int seq_name_len = strlen(seq->name.s);

  BindingSite *site = malloc(sizeof(BindingSite));

  site->sequence_name = calloc(seq_name_len + 1, sizeof(char));
  site->sequence_name[seq_name_len] = '\0';
  strncpy(site->sequence_name, seq->name.s, seq_name_len);

  site->spacer_length = spacer_size;

  site->f_idx = f_idx;
  site->r_idx = r_idx;

  // Plus

  site->indexes[0] = i;

  site->sequence[0] = calloc(num_forward_rvds + 2 + 1, sizeof(char));
  site->sequence[0][num_forward_rvds + 2] = '\0';

  strncpy(site->sequence[0], seq->seq.s + i - 1, 1);

  // Upstream
  site->sequence[0][1] = ' ';
  strncpy(site->sequence[0] + 2, seq->seq.s + i, num_forward_rvds);

  site->scores[0] = forward_score;

  // Minus

  site->indexes[1] = j;

  site->sequence[1] = calloc(num_reverse_rvds + 2 + 1, sizeof(char));
  site->sequence[1][num_reverse_rvds + 2] = '\0';

  strncpy(site->sequence[1], seq->seq.s + j - num_reverse_rvds + 1, num_reverse_rvds);

  // Upstream
  site->sequence[1][num_reverse_rvds] = ' ';
  strncpy(site->sequence[1] + num_reverse_rvds + 1, seq->seq.s + j + 1, 1);

  site->scores[1] = reverse_score;

  return site;

}

// Identify and print out TAL effector binding sites
void find_binding_sites(FILE *log_file, kseq_t *seq, Array **rvd_seqs, Hashmap *diresidue_scores, double *cutoffs, Array *results) {

  int c_upstream = *((int *) hashmap_get(talesf_kwargs, "c_upstream"));
  int spacer_min = *((int *) hashmap_get(talesf_kwargs, "spacer_min"));
  int spacer_max = *((int *) hashmap_get(talesf_kwargs, "spacer_max"));

  if (array_size(rvd_seqs[0]) + array_size(rvd_seqs[0]) + spacer_min > seq->seq.l) {
    logger(log_file, "Warning: skipping sequence '%s' since it is shorter than the RVD sequence\n", seq->seq.s);
    return;
  }

  logger(log_file, "Scanning %s for binding sites (length %ld)", seq->name.s, seq->seq.l);

  for (int f_idx = 0; f_idx < 2; f_idx++) {
    
    for (int r_idx = 0; r_idx < 2; r_idx++) {
        
      Array *forward_rvd_seq = rvd_seqs[f_idx];
      Array *reverse_rvd_seq = rvd_seqs[r_idx];
      
      double forward_cutoff = cutoffs[f_idx];
      double reverse_cutoff = cutoffs[r_idx];
        
      int num_forward_rvds = array_size(forward_rvd_seq);
      int num_reverse_rvds = array_size(reverse_rvd_seq);

      for (unsigned long i = 1; i <= seq->seq.l - num_forward_rvds; i++) {
        
        char forward_upstream = seq->seq.s[i-1];
        
        if ((c_upstream != 0 && (forward_upstream == 'C' || forward_upstream == 'c')) || (c_upstream != 1 && (forward_upstream == 'T' || forward_upstream == 't'))) {
          
          double forward_score = score_binding_site(seq, i, forward_rvd_seq, diresidue_scores, forward_cutoff, 0);
          
          if (forward_score <= forward_cutoff) {
              
            for (int spacer_size = spacer_min; spacer_size < spacer_max + 1; spacer_size++) {
                
               unsigned long j = i + num_forward_rvds + spacer_size + num_reverse_rvds - 1;
               
               if (j >= seq->seq.l - 2) continue;
               
               char reverse_upstream = seq->seq.s[j + 1];
               
               if ((c_upstream != 0 && (reverse_upstream == 'G' || reverse_upstream == 'g')) || (c_upstream != 1 && (reverse_upstream == 'A' || reverse_upstream == 'a'))) {
                   
                 double reverse_score = score_binding_site(seq, j - num_reverse_rvds + 1, reverse_rvd_seq, diresidue_scores, reverse_cutoff, 1);
                 
                 if (reverse_score <= reverse_cutoff) {
                     
                   BindingSite *site = create_binding_site(seq, i, j, num_forward_rvds, forward_score, num_reverse_rvds, reverse_score, spacer_size, f_idx, r_idx);
                   
                   #pragma omp critical (add_result)
                   array_add(results, site);
                   
                 }
                   
               }
                
            }
            
          }
          
        }
        
      }
      
    }
    
  }

}

int run_talesf_task(Hashmap *kwargs) {

  talesf_kwargs = kwargs;

  // Options
  char *seq_filename = hashmap_get(kwargs, "seq_filename");
  char *rvd_string = hashmap_get(kwargs, "rvd_string");
  char *rvd_string2 = hashmap_get(kwargs, "rvd_string2");
  char *log_filepath = hashmap_get(kwargs, "log_filepath");

  double weight = *((double *) hashmap_get(kwargs, "weight"));
  double cutoff = *((double *) hashmap_get(kwargs, "cutoff"));

  int numprocs = *((int *) hashmap_get(kwargs, "num_procs"));

  // Setup the logger

  FILE *log_file = stdout;

  if (log_filepath && strcmp(log_filepath, "NA") != 0) {
    log_file = fopen(log_filepath, "a");
  }

  // Determine number of sequences in the input file

  int seq_num;
  char cmd[256], line[32];

  sprintf(cmd, "grep '^>' %s | wc -l", seq_filename);
  FILE *fasta_filesize_in = popen(cmd, "r");

  if (!fasta_filesize_in) {
    perror("Error: unable to check fasta file size\n");
    logger(log_file, "Error: unable to check fasta file size");
    if (log_file != stdout) fclose(log_file);
    return 1;
  }

  fgets(line, sizeof(line), fasta_filesize_in);
  pclose(fasta_filesize_in);
  seq_num = atoi(line);

  // Process RVD sequences

  Array *results = array_new( sizeof(BindingSite *) );

  Array *rvd_seq = rvd_string_to_array(rvd_string);
  Array *rvd_seq2 = rvd_string_to_array(rvd_string2);
  
  Array *rvd_seqs[2];
  
  rvd_seqs[0] = rvd_seq;
  rvd_seqs[1] = rvd_seq2;
  
  Array *joined_rvd_seq = array_concat(rvd_seq, rvd_seq2);

  // Get RVD/bp matching scores

  Hashmap *diresidue_probabilities = get_diresidue_probabilities(joined_rvd_seq, weight);
  Hashmap *diresidue_scores = convert_probabilities_to_scores(diresidue_probabilities);
  hashmap_delete(diresidue_probabilities, NULL);

  // Compute optimal score for the RVD sequences

  double best_score = get_best_score(rvd_seq, diresidue_scores);
  double best_score2 = get_best_score(rvd_seq2, diresidue_scores);
  
  // Define score cutoffs for match sites

  double cutoffs[2];
  
  cutoffs[0] = cutoff * best_score;
  cutoffs[1] = cutoff * best_score2;

  // Begin processing

  int abort = 0;

  omp_set_num_threads(numprocs);

  #pragma omp parallel
  {

    // Open sequence file
    gzFile seqfile = gzopen(seq_filename, "r");

    if (!seqfile) {

      logger(log_file, "Error: unable to open sequence '%s'", seq_filename);
      abort = 1;

    } else {

      kseq_t *seq = kseq_init(seqfile);

      int j = 0;

      #pragma omp for schedule(static)
      for (int i = 0; i < seq_num; i++) {

        #pragma omp flush (abort)
        if (!abort) {

          while (j <= i) {

            int result = kseq_read(seq);

            if (result < 0) {
              logger(log_file, "Error: problem parsing data from '%s'", seq_filename);
              abort = 1;
            }

            j++;

          }

          if (!abort) {
            find_binding_sites(log_file, seq, rvd_seqs, diresidue_scores, cutoffs, results);
          }

        }

      }

      kseq_destroy(seq);
      gzclose(seqfile);

    }

  }

  if (!abort) {

    qsort(results->data, array_size(results), sizeof(BindingSite *), binding_site_compare_pos);

    abort = print_results(results, rvd_seqs, best_score, best_score2, log_file);

    logger(log_file, "Finished");

  }

  // Free memory

  if (results) {

    for (int i = 0; i < array_size(results); i++) {

      BindingSite *site = (BindingSite *) array_get(results, i);

      free(site->sequence[0]);
      free(site->sequence[1]);
      free(site->sequence_name);
      free(site);

    }

    array_delete(results, NULL);

  }


  if (rvd_seq) {
    array_delete(rvd_seq, free);
  }

  if (rvd_seq2) {
    array_delete(rvd_seq2, free);
  }

  if (joined_rvd_seq) {
    array_delete(joined_rvd_seq, NULL);
  }

  if (diresidue_scores) {
    hashmap_delete(diresidue_scores, free);
  }

  if (log_file != stdout) {
    fclose(log_file);
  }

  return abort;

}
