#ifndef BCUTILS_BCUTILS
#define BCUTILS_BCUTILS

#include <stdio.h>
#include <stdarg.h>
#include <time.h>
#include <math.h>

#include "Array.h"
#include "Hashmap.h"

// Allocate memory for an array of 5 doubles
double *double_array(double a, double c, double g, double t, double dummy);

// Allocate memory for an array of 4 integers
int *int_array(int a, int c, int g, int t);

// Compute diresidue probabilities from observed counts
Hashmap *get_diresidue_probabilities(Array *rvdseq, double w);

// Convert diresidue probabilities into scores
Hashmap *convert_probabilities_to_scores(Hashmap *diresidue_probabilities);

// Given a sequence of RVDs, calculate the optimal score
double get_best_score(Array *rvdseq, Hashmap *rvdscores);

void logger(FILE* log_file, char* message, ...);

Array *rvd_string_to_array(char *rvd_string);

void str_replace(char *haystack, char *needle, char *replacement);

#endif
