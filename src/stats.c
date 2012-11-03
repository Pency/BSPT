#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <zlib.h>

#include "seq.h"
#include "html_utils.h"


void singleUsage(){
	printf("Program:  statistic function of NGS data preprocess toolkits\n");
	printf("Version:  %s by Yanpc\n\n",VERSION);
	printf("Usage:    pt stats <file>\n\n");
}

void init_values(stats* seq_stats){
        int i,j;
//statistic variables of sequenc legth
	seq_stats->max_length = 0;
	seq_stats->min_length = 999999;
	seq_stats->mean_length = 0.0;
	seq_stats->total_reads = 0;
	seq_stats->total_length = 0;
//statistic varialbes of quality
	seq_stats->max_qual = 0;
	seq_stats->min_qual = 1000;

        for (i=0;i<MAX_READ_LENGTH;i++) {
                for (j=0;j<5;++j) {
                        seq_stats->cycles[i].min = 100 ;
                        seq_stats->cycles[i].max = -100;
                }
        }
	for(i=0;i<MAX_QUALITY;i++)
		seq_stats->qual_hist[i]=0;
	for(i=0;i<20;i++)
		seq_stats->gc_content[i]=0;

	seq_stats->length_hist = malloc(sizeof(int_list)); 
}

int quartile_value(int n, int array[]){

        if (n<=0) {
                error_msg("internal error at quartile_value ()\n");
                return -1;
        }

	int i = 0 ;
	long int accumulated = 0;
	for(i=0;i<=MAX_QUALITY;i++){
		if(array[i]==0)
			continue;
		accumulated += array[i];
                if (accumulated >= n)
                        break;
        }

        if (n>=accumulated) {
                error_msg("internal error at quartile_value ()\n");
                return -1;
        }

        return i;
}

float stats_gc_content(SEQ_QUAL* read){
	float gc = 0.0;
	int i = 0;
	
	for(i=0; i<read->length; i++){
		if(read->seq[i] == 'G' || read->seq[i] == 'C')
			gc++;
		if(read->seq[i] == 'g' || read->seq[i] == 'c')
			gc++;
	}
	return (gc * 100 / read->length);
}

void update_count(int_list *head, int id, int count){
	int_list *temp, *new;
	int state = 0;
	if(head->next == NULL){
		new = malloc(sizeof(int_list));
		new->count = count;
		new->id = id;
		head->next = new;
		state = 1;
	}else{
		temp = head;
		do{
			temp = temp->next;
			if(temp->id == id){
				temp->count += count;
				state = 1;
			} 
			if(temp->next == NULL || temp->next->id > id){
				new = malloc(sizeof(int_list));
				new->count = count;
				new->id = id;
				new->next = temp->next;
				temp->next = new;
				state = 1;
			}
		}while(state == 0);
	}

}

int stats_fastq(char* file, int base, FILE* output){
	int i = 0;
	int baseQual;
	int nuc_idx=0;
	int length[MAX_READ_LENGTH];
	for(i=0;i<MAX_READ_LENGTH;i++)
		length[i]=0;

	char nucleotide_index[255];
	nucleotide_index['A'] = A ;
	nucleotide_index['a'] = A ;
	nucleotide_index['C'] = C ;
	nucleotide_index['c'] = C ;
	nucleotide_index['G'] = G ;
	nucleotide_index['g'] = G ;
	nucleotide_index['T'] = T ;
	nucleotide_index['t'] = T ;
	nucleotide_index['N'] = N ;
	nucleotide_index['n'] = N ;

	stats seq_stats;
	init_values(&seq_stats);
	SEQ_QUAL item = init_read();

	int index = 0;

	gzFile zfp = gzopen_report(file,"r");

	while(read_fastq(zfp, &item, index++) > 0){
		int tmp_gc = (int) (stats_gc_content(&item) /20.0);
		seq_stats.gc_content[tmp_gc]++;
		seq_stats.total_reads++;
		for(i=0; i<item.length; i++){
			baseQual = item.qual[i] - base;
			nuc_idx = nucleotide_index[(int)item.seq[i]];
			seq_stats.cycles[i].max = MAX(seq_stats.cycles[i].max,baseQual);
			seq_stats.cycles[i].min = MIN(seq_stats.cycles[i].min,baseQual);
			seq_stats.cycles[i].sum += baseQual;
			seq_stats.cycles[i].count++;
			seq_stats.cycles[i].nucleotide_count[nuc_idx]++;
			seq_stats.cycles[i].bases_quality_count[baseQual]++;
			seq_stats.qual_hist[baseQual]++;
		}
		seq_stats.total_length += (long int) item.length;
		length[item.length]++;
	}

	for(i=0;i<MAX_READ_LENGTH;i++){
		if(length[i] == 0)
			continue;
		update_count(seq_stats.length_hist, i, length[i]);
		seq_stats.max_length = MAX(seq_stats.max_length, i);
		seq_stats.min_length = MIN(seq_stats.min_length, i);
	}

	seq_stats.mean_length = (float) (seq_stats.total_length/seq_stats.total_reads);

	print_stats(output, &seq_stats);
	gzclose(zfp);
	free_read(&item);
	return 0;
}

int stats_fastq_bz2(char* file){

	return 0;
}

int stats_fasta(char* file, FILE* output){
	int index = 0;
//statistic variables of sequenc legth
	stats seq_stats;
	init_values(&seq_stats);
	SEQ_QUAL item=init_read();
	gzFile zfp = gzopen_report(file,"r");
	while(read_fasta(zfp, &item, index++) > 0){
		seq_stats.total_reads++;
		seq_stats.total_length += item.length;
		update_count(seq_stats.length_hist, item.length, 1);
	}
	seq_stats.mean_length = (float) (seq_stats.total_length/seq_stats.total_reads);

	print_stats(output, &seq_stats);

	gzclose(zfp);
	free_read(&item);
	return 0;
}

int stats_main(int argc, char* argv[]){

	if(argc < 2){
		singleUsage();
		exit(1);
	}
	char *file = argv[1];
	//strcpy(file, argv[1]);
	int f1=0;

	f1 = detect_filetype(file);

	if(f1 & FILE_UNKN){
		error_msg("Unknown sequence format, not supported in this version\n");
		return 1;
	}

	if(f1 & FILE_BZ2){
		error_msg("Not support BZIP format compressed file in this version\n");
		return 1;
	}

	int base = f1 & FILE_PHRED33 ? BASE_SANGER : BASE_SOLEXA;
	int nucleotide_index[255];

	FILE* result=fopen_report("result.html","w+");
	head(result);

//start processing GZIP or regular FASTQ/FASTA sequence
	if(f1 & FILE_FASTA)
		stats_fasta(file, result);

	if(f1 & FILE_FASTQ)
		stats_fastq(file, base, result);

	foot(result);

	return 0;
}
