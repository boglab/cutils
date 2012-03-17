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

int binding_site_compare_score(const void * a, const void * b)
{
  BindingSite *real_a = *((BindingSite **)a);
  BindingSite *real_b = *((BindingSite **)b);

  double score_diff = ((real_a->score) - (real_b->score));

  if(score_diff < -0.01) {
    return -1;
  } else if(score_diff > 0.01) {
    return 1;
  } else {
    return 0;
  }

}

void print_results(Array *results, char *sourcestr, Array* rvdseq, double best_score, int forwardonly, char* output_filepath, int create_tabfile) {

    int num_rvds = array_size(rvdseq);
    char strand = '+';
    char *tab_strand = "Plus";
    char *plus_strand_sequence;
    FILE *gff_out_file, *tab_out_file;

    size_t output_filepath_length;
    char* temp_output_filepath;

    char *rvdstring = calloc(3 * num_rvds, sizeof(char));

    int i;

    for(i = 0; i < num_rvds; i++)
    {
      char *rvd = array_get(rvdseq, i);
      if(i == 0)
        strncpy(rvdstring, rvd, 2);
      else
        sprintf(rvdstring + (i*3)-1, "_%s", rvd);
    }
    rvdstring[3*num_rvds - 1] = '\0';

    if(create_tabfile)
    {

        output_filepath_length = strlen(output_filepath) + 5;
        temp_output_filepath = calloc(output_filepath_length + 1, sizeof(char));

        sprintf(temp_output_filepath, "%s.txt", output_filepath);
        tab_out_file = fopen(temp_output_filepath, "w");
        memset(temp_output_filepath, '\0', output_filepath_length);
        sprintf(temp_output_filepath, "%s.gff3", output_filepath);
        gff_out_file = fopen(temp_output_filepath, "w");
        free(temp_output_filepath);

    }
    else
    {
        gff_out_file = fopen(output_filepath, "w");
    }

    if(!gff_out_file || (create_tabfile && !tab_out_file))
    {
        fprintf(stderr, "Error: unable to open output file '%s'\n", output_filepath);
        return;
    }

    if(create_tabfile)
    {
        if (forwardonly)
        {
            fprintf(tab_out_file, "table_ignores:Plus strand sequence\n");
        }

        fprintf(tab_out_file, "Best Possible Score:%.2lf\n", best_score);
        fprintf(tab_out_file, "Genome Coordinates\tStrand\tScore\tTarget Sequence\tPlus strand sequence\n");

    }


    fprintf(gff_out_file, "##gff-version 3\n");

    if (forwardonly)
    {
        fprintf(gff_out_file, "#table_display_tags:target_sequence\n");
    }
    else
    {
        fprintf(gff_out_file, "#table_display_tags:target_sequence,plus_strand_sequence\n");
    }

    fprintf(gff_out_file, "#Best Possible Score:%.2lf\n", best_score);

    for(i = 0; i < array_size(results); i++)
    {

      BindingSite *site = (BindingSite *)array_get(results, i);
      char *sequence = site->sequence;

      if(site->strand > 0)
        plus_strand_sequence = sequence;
      else
      {
        int j;

        plus_strand_sequence = sequence;
        sequence = malloc(sizeof(char)*(num_rvds+1));
        sequence[num_rvds] = '\0';

        for(j = 0; j < num_rvds; j++)
        {
          char base = site->sequence[num_rvds - j - 1];
          if(base == 'A' || base == 'a')
            sequence[j] = 'T';
          else if(base == 'C' || base == 'c')
            sequence[j] = 'G';
          else if(base == 'G' || base == 'g')
            sequence[j] = 'C';
          else if(base == 'T' || base == 't')
            sequence[j] = 'A';
          else
          {
            fprintf(stderr, "Error: unexpected character '%c'\n", base);
            exit(1);
          }
        }
        strand = '-';
        tab_strand = "Minus";
      }

      if(create_tabfile)
      {
          fprintf( tab_out_file, "%s\t%s\t%.2lf\t%lu\t%s\t%s\n",
                   site->sequence_name, tab_strand, site->score, site->index + 1, sequence, plus_strand_sequence);
      }

      fprintf( gff_out_file, "%s\t%s\t%s\t%lu\t%lu\t%.2lf\t%c\t.\trvd_sequence=%s;target_sequence=%s;\n",
               site->sequence_name, sourcestr, "TAL_effector_binding_site", site->index + 1,
               site->index + num_rvds, site->score, strand, rvdstring, sequence);

      if(plus_strand_sequence != sequence) {
          free(sequence);
      }

    }

    free(rvdstring);
    fclose(gff_out_file);

    if(create_tabfile)
    {
        fclose(tab_out_file);
    }

}

// Print usage statement
void print_usage(FILE *outstream, char *progname)
{
  fprintf( outstream, "\nUsage: %s [options] genomeseq \"rvdseq\"\n"
           "  Options:\n"
           "    -t|--tabfile          generate a tab-delimited output file in addition to gff3\n"
           "    -f|--forwardonly      only search the forward strand of the genomic sequence\n"
           "    -h|--help             print this help message and exit\n"
           "    -n|--numprocs         the number of processors to use; default is 1\n"
           "    -o|--outfile          file to which output will be written; default is terminal\n"
           "                          (STDOUT)\n"
           "    -s|--source           text for the third column of the GFF3 output; default is\n"
           "                          TALESF\n"
           "    -w|--weight           user-defined weight; default is 0.9\n"
           "    -x|--cutoffmult       multiple of best score at which potential sites will be\n"
           "                          filtered; default is 3.0\n\n", progname );
}

// Identify and print out TAL effector binding sites
void find_binding_sites(kseq_t *seq, Array *rvdseq, Hashmap *diresidue_scores, double cutoff, int forwardonly, Array *results)
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
    if(seq->seq.s[i-1] == 'T' || seq->seq.s[i-1] == 't')
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

        site->sequence = calloc(num_rvds + 1, sizeof(char));
        site->sequence_name = calloc(seq_name_len + 1, sizeof(char));
        site->sequence[num_rvds] = '\0';
        site->sequence_name[seq_name_len] = '\0';

        site->strand = 1;
        site->index = i;

        strncpy(site->sequence, seq->seq.s + site->index, num_rvds);
        strncpy(site->sequence_name, seq->name.s, seq_name_len);

        site->score = cumscore;

        #pragma omp critical (add_result)
        array_add(results, site);

      }
    }

    if(!forwardonly)
    {
      if(seq->seq.s[i + num_rvds - 1] == 'A' || seq->seq.s[i + num_rvds - 1] == 'a')
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

          site->sequence = calloc(num_rvds + 1, sizeof(char));
          site->sequence_name = calloc(seq_name_len + 1, sizeof(char));
          site->sequence[num_rvds] = '\0';
          site->sequence_name[seq_name_len] = '\0';

          site->strand = -1;
          site->index = i - 1;

          strncpy(site->sequence, seq->seq.s + site->index, num_rvds);
          strncpy(site->sequence_name, seq->name.s, seq_name_len);

          site->score = cumscore;

          #pragma omp critical (add_result)
          array_add(results, site);

        }
      }
    }
  }

}

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
Hashmap *get_diresidue_probabilities(double w)
{
  char **diresidues;
  Hashmap *diresidue_counts, *diresidue_probabilities;
  int i, num_diresidues;

  // First, use a hashmap to store counts
  diresidue_counts = hashmap_new(64);

  // Add known counts
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

// Main method
int main(int argc, char **argv)
{
  // Program arguments and options
  char *progname, *seqfilename, *rvdstring, *output_filepath = NULL, sourcestr[128];
  int forwardonly, numprocs;
  double w, x;

  // Program variable domain
  Array *rvdseq;
  char *tok, cmd[256], line[32];
  gzFile seqfile;
  kseq_t *seq;
  int i, j, seqnum, rank, create_tabfile;

  Array *results = array_new( sizeof(BindingSite *) );

  // Set option defaults
  create_tabfile = 0;
  forwardonly = 0;
  strcpy(sourcestr, "TALESF");
  numprocs = 1;
  w = 0.9;
  x = 3.0;

  // Parse options
  progname = argv[0];
  int opt, optindex;
  const char *optstr = "tfhn:o:s:w:x:";
  const struct option otsf_options[] =
  {
    { "tabfile", no_argument, NULL, 't' },
    { "forwardonly", no_argument, NULL, 'f' },
    { "help", no_argument, NULL, 'h' },
    { "numprocs", required_argument, NULL, 'n' },
    { "outfile", required_argument, NULL, 'o' },
    { "source", required_argument, NULL, 's' },
    { "weight", required_argument, NULL, 'w' },
    { "cutoffmult", required_argument, NULL, 'x' },
    { NULL, no_argument, NULL, 0 },
  };

  for( opt = getopt_long(argc, argv + 0, optstr, otsf_options, &optindex);
       opt != -1;
       opt = getopt_long(argc, argv + 0, optstr, otsf_options, &optindex) )
  {
    switch(opt)
    {
      case 'f':
        forwardonly = 1;
        break;

      case 't':
        create_tabfile = 1;
        break;

      case 'h':
        print_usage(stdout, progname);
        return 0;

      case 'n':
        if( sscanf(optarg, "%d", &numprocs) == EOF )
        {
          fprintf(stderr, "Error: unable to convert numprocs '%s' to an integer\n", optarg);
          return 1;
        }
        break;

      case 'o':

        if(output_filepath != NULL) {
            free(output_filepath);
            output_filepath = NULL;
        }
        output_filepath = calloc(strlen(optarg) + 1, sizeof(char));
        strcpy(output_filepath, optarg);

        break;

      case 's':
        strcpy(sourcestr, optarg);
        break;

     case 'w':
       if( sscanf(optarg, "%lf", &w) == EOF )
       {
         fprintf(stderr, "Error: unable to convert weight '%s' to a double\n", optarg);
         return 1;
       }
       break;

     case 'x':
       if( sscanf(optarg, "%lf", &x) == EOF )
       {
         fprintf(stderr, "Error: unable to convert cutoff multiple '%s' to a double\n", optarg);
         return 1;
       }
       break;
    }
  }

  // Parse arguments
  if(argc - optind != 2)
  {
    fputs("Error: must provide genomic sequence (file) and RVD sequence (string)...no more, no less\n", stderr);
    print_usage(stderr, progname);
    return 1;
  }
  seqfilename = argv[optind];
  rvdstring = argv[optind + 1];

  // Print verbose output
  fprintf(stderr, "\n%-20s '%s'\n", "Sequence data:", seqfilename);
  fprintf(stderr, "%-20s '%s'\n", "RVD sequence:", rvdstring);

  rvdseq = array_new( sizeof(char *) );
  tok = strtok(rvdstring, " _");
  while(tok != NULL)
  {
    char *r = strdup(tok);
    array_add(rvdseq, r);
    tok = strtok(NULL, " _");
  }

  // Get RVD/bp matching scores
  Hashmap *diresidue_probabilities = get_diresidue_probabilities(w);
  Hashmap *diresidue_scores = convert_probabilities_to_scores(diresidue_probabilities);
  hashmap_delete(diresidue_probabilities, NULL);

  // Compute optimal score for this RVD sequence
  double best_score = get_best_score(rvdseq, diresidue_scores);
  fprintf(stderr, "%-20s %.4lf\n\n", "RVD best score:", best_score);

  // Determine number of sequences in file
  sprintf(cmd, "grep '^>' %s | wc -l", seqfilename);
  FILE *in = popen(cmd, "r");
  if(!in)
  {
    perror("Error: unable to check fasta file size\n");
    return 1;
  }
  fgets(line, sizeof(line), in);
  pclose(in);
  seqnum = atoi(line);

  // Begin processing

  omp_set_num_threads(numprocs);
  #pragma omp parallel private(i, j, seq, seqfile, rank)
  {
    rank = omp_get_thread_num();

    // Open genomic sequence file
    seqfile = gzopen(seqfilename, "r");
    if(!seqfile)
    {
      fprintf(stderr, "Error: unable to open genomic sequence '%s'\n", seqfilename);
      exit(1);
    }
    seq = kseq_init(seqfile);

    j = 0;
    #pragma omp for schedule(static)
    for(i = 0; i < seqnum; i++)
    {
      while(j <= i)
      {
        int result = kseq_read(seq);
        if(result < 0)
        {
          fprintf(stderr, "Error: problem parsing data from '%s'\n", seqfilename);
          exit(1);
        }
        j++;
      }

      fprintf(stderr, "  Processor %d working on sequence '%s' (length %ld)\n", rank, seq->name.s, seq->seq.l);
      find_binding_sites(seq, rvdseq, diresidue_scores, best_score * x, forwardonly, results);

    }

    kseq_destroy(seq);
    gzclose(seqfile);

  }

  qsort(results->data, array_size(results), sizeof(BindingSite *), binding_site_compare_score);

  print_results(results, sourcestr, rvdseq, best_score, forwardonly, output_filepath, create_tabfile);

  // Free memory

  free(output_filepath);
  for(i = 0; i < array_size(results); i++)
  {
      BindingSite *site = (BindingSite *)array_get(results, i);

      free(site->sequence);
      free(site->sequence_name);
      free(site);

  }

  array_delete(results);


  for(i = 0; i < array_size(rvdseq); i++)
  {
    char *temp = (char *)array_get(rvdseq, i);
    free(temp);
  }


  array_delete(rvdseq);
  hashmap_delete(diresidue_scores, free);

  return 0;
}
