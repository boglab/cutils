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

// Include my libraries
#include "Array.h"
#include "Hashmap.h"

// Initialize the kseq library
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

typedef struct
{
  int strand;
  kseq_t *seq;
  Array *rvdseq;
  char *rvdstring;
  char *sourcestr;
  unsigned long index;
  double score;
} BindingSite;

// Print the binding site in GFF3 format to the given output stream
void binding_site_print(BindingSite *site, FILE *outstream)
{
  int num_rvds = array_size(site->rvdseq);
  char strand = '+';
  char *sequence = malloc(sizeof(char)*(num_rvds+1));
  sequence[num_rvds] = '\0';
  if(site->strand > 0)
    strncpy(sequence, site->seq->seq.s + site->index, num_rvds);
  else
  {
    int i;
    for(i = 0; i < num_rvds; i++)
    {
      char base = site->seq->seq.s[site->index + num_rvds - i - 1];
      if(base == 'A' || base == 'a')
        sequence[i] = 'T';
      else if(base == 'C' || base == 'c')
        sequence[i] = 'G';
      else if(base == 'G' || base == 'g')
        sequence[i] = 'C';
      else if(base == 'T' || base == 't')
        sequence[i] = 'A';
      else
      {
        fprintf(stderr, "Error: unexpected character '%c'\n", base);
        exit(1);
      }
    }
    strand = '-';
  }

  fprintf( outstream, "%s\t%s\t%s\t%lu\t%lu\t%.2lf\t%c\t.\trvd_sequence=%s;target_sequence=%s;\n",
           site->seq->name.s, site->sourcestr, "TAL_effector_binding_site", site->index + 1,
           site->index + num_rvds, site->score, strand, site->rvdstring, sequence);
  free(sequence);
}

// Print usage statement
void print_usage(FILE *outstream, char *progname)
{
  fprintf( outstream, "\nUsage: %s [options] genomeseq \"rvdseq\"\n"
           "  Options:\n"
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
void find_binding_sites(kseq_t *seq, Array *rvdseq, Hashmap *diresidue_scores, double cutoff, char *sourcestr, int forwardonly, FILE *outstream)
{
  unsigned long i, j;
  char *rvdstring;

  if(array_size(rvdseq) > seq->seq.l)
  {
    fprintf(stderr, "Warning: skipping sequence '%s' since it is shorter than the RVD sequence\n", seq->seq.s);
    return;
  }

  rvdstring = malloc(sizeof(char)*3*array_size(rvdseq));
  for(i = 0; i < array_size(rvdseq); i++)
  {
    char *rvd = array_get(rvdseq, i);
    if(i == 0)
      strncpy(rvdstring, rvd, 2);
    else
      sprintf(rvdstring + (i*3)-1, "_%s", rvd);
  }
  rvdstring[3*array_size(rvdseq)] = '\0';

  for(i = 1; i <= seq->seq.l - array_size(rvdseq); i++)
  {
    if(seq->seq.s[i-1] == 'T' || seq->seq.s[i-1] == 't')
    {
      double cumscore = 0.0;
      for(j = 0; j < array_size(rvdseq); j++)
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
        BindingSite site = { 1, seq, rvdseq, rvdstring, sourcestr, i, cumscore };
        binding_site_print(&site, outstream);
      }
    }

    if(!forwardonly)
    {
      if(seq->seq.s[i + array_size(rvdseq) - 1] == 'A' || seq->seq.s[i + array_size(rvdseq) - 1] == 'a')
      {
        double cumscore = 0.0;
        for(j = 0; j < array_size(rvdseq); j++)
        {
          char *rvd = array_get(rvdseq, array_size(rvdseq) - j - 1);
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
          BindingSite site = { -1, seq, rvdseq, rvdstring, sourcestr, i-1, cumscore };
          binding_site_print(&site, outstream);
        }
      }
    }
  }

  free(rvdstring);
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
  char *progname, *seqfilename, *rvdstring, sourcestr[128];
  int forwardonly, numprocs;
  FILE *outstream;
  double w, x;

  // Program variable domain
  Array *rvdseq;
  char *tok, cmd[256], line[32];
  gzFile seqfile;
  kseq_t *seq;
  int i, j, seqnum, rank;

  // Set option defaults
  forwardonly = 0;
  outstream = stdout;
  strcpy(sourcestr, "TALESF");
  numprocs = 1;
  w = 0.9;
  x = 3.0;

  // Parse options
  progname = argv[0];
  int opt, optindex;
  const char *optstr = "fhn:o:s:w:x:";
  const struct option otsf_options[] =
  {
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
        outstream = fopen(optarg, "w");
        if(!outstream)
        {
          fprintf(stderr, "Error: unable to open output file '%s'\n", optarg);
          return 1;
        }
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
  fprintf(outstream, "##gff-version 3\n");

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
      
    //while((i = kseq_read(seq)) >= 0)
    //{
      fprintf(stderr, "  Processor %d working on sequence '%s' (length %ld)\n", rank, seq->name.s, seq->seq.l);
      find_binding_sites(seq, rvdseq, diresidue_scores, best_score * x, sourcestr, forwardonly, outstream);
    }
    kseq_destroy(seq);
    gzclose(seqfile);  
  }

  // Free memory
  for(i = 0; i < array_size(rvdseq); i++)
  {
    char *temp = (char *)array_get(rvdseq, i);
    free(temp);
  }
  array_delete(rvdseq);
  hashmap_delete(diresidue_scores, free);

  return 0;
}
