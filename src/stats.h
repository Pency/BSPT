#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <zlib.h>

#include "seq.h"

typedef struct {
	int min;
	int max;
	int count;
	unsigned long long sum;
	int bases_quality_count[50];
	int nucleotide_count[5];
} cycle;

typedef struct nlist
{
	int id;
	int count;
	struct nlist *next;
} int_list;

typedef struct clist
{
	int id;
	char* text;
	struct clist *next;
} char_list;

typedef struct {
//statistic variables of sequenc legth
	int max_length;
	int min_length;
	float mean_length;
	long int total_reads;
	long int total_length;
//statistic varialbes of quality
	int max_qual;
	int min_qual;
	long int qual_hist[MAX_QUALITY] ;
	long int gc_content[20] ;
	int_list *length_hist;

	cycle cycles[MAX_READ_LENGTH];
} stats;

void singleUsage();

void init_values(stats *seq_stats);

void insert_node();
void update_count();

int quartile_value(int n, int array[]);

float stats_gc_content(SEQ_QUAL* read);

int stats_fastq(char* file, int base, FILE* output);

#ifdef HAVE_BZLIB
int stats_fastq_bz2(char* file);
#endif

int stats_fasta(char* file, FILE* output);
