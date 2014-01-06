#include <stdio.h>
#include <stdarg.h>
#include <time.h>
#include <math.h>

#include "bcutils.h"

double *double_array(double a, double c, double g, double t, double dummy) {
  double *array = malloc(sizeof(double)*5);
  array[0] = a;
  array[1] = c;
  array[2] = g;
  array[3] = t;
  array[4] = dummy;
  return array;
}

int *int_array(int a, int c, int g, int t) {
  int *array = malloc(sizeof(int)*4);
  array[0] = a;
  array[1] = c;
  array[2] = g;
  array[3] = t;
  return array;
}

Hashmap *get_diresidue_probabilities(Array *rvdseq, double w) {
  char **diresidues;
  Hashmap *diresidue_counts, *diresidue_probabilities;
  int i, num_diresidues;

  // First, use a hashmap to store counts
  diresidue_counts = hashmap_new(128);

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
  hashmap_add(diresidue_counts, "NH", int_array( 0,  0,  1,  0));
  
  // Add unknown counts
  hashmap_add(diresidue_counts, "S*", int_array(1, 1, 1, 1));
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

  diresidue_probabilities = hashmap_new(128);
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

    hashmap_add(diresidue_probabilities, diresidues[i], double_array(p[0], p[1], p[2], p[3], 0));
  }

  free(diresidues);
  hashmap_delete(diresidue_counts, free);
  return diresidue_probabilities;
}

Hashmap *convert_probabilities_to_scores(Hashmap *diresidue_probabilities) {
  char **diresidues;
  int i, j, num_diresidues;
  Hashmap *diresidue_scores = hashmap_new(128);

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

Array *rvd_string_to_array(char *rvd_string) {

  char *tok;
  Array *rvd_seq;

  // strtok screws up the rvd string
  char *mutable_rvd_string = strdup(rvd_string);

  rvd_seq = array_new( sizeof(char *) );
  tok = strtok(mutable_rvd_string, " _");

  while (tok != NULL)
  {
    char *r = strdup(tok);
    array_add(rvd_seq, r);
    tok = strtok(NULL, " _");
  }

  free(mutable_rvd_string);

  return rvd_seq;

}

void str_replace(char *haystack, char *needle, char *replacement) {

  char *pos = strstr(haystack, needle);

  while (pos != NULL) {

    strcpy(pos, replacement);
    pos = strstr(haystack, needle);

  }

}
