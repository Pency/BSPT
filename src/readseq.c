#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <zlib.h>

#include "seq.h"

int read_fasta(gzFile zfp, SEQ_QUAL *item, int id)
/*
default reading file use gzopen, because gzopen can
 handle ether gzip format file or not gzip format file, 
*/
{

	char c='\n';
	int i=0;
	int max=BUFFER_LENGTH;

	if((c=gzgetc(zfp)) != '>'){
		if(gzeof(zfp)) return -1;
		error_msg("sequence %d is not a FASTA format(%c)\n", id,c);
		return -1;
	}
	
	while(!gzeof(zfp) && (c=gzgetc(zfp)) != '\n'){
		if(c != '\n'){
			if(i+1 >= max){
				max += BUFFER_LENGTH;
				item->name = realloc(item->name, sizeof(char) * (max + BUFFER_LENGTH));
			}
			item->name[i++] = c;
		}
	}
	i=0;
	while(!gzeof(zfp) && (c=gzgetc(zfp)) != '>'){
		if(isalpha(c)){
			if(i+1 >= max){
				max += BUFFER_LENGTH;
				item->seq = realloc(item->seq, sizeof(char) * (max + BUFFER_LENGTH));
			}
			item->seq[i++] = c;
		}
	}

	if(!gzeof(zfp))
		gzseek(zfp,0,SEEK_CUR-1);

	item->length = strlen(item->seq);
	item->id=id;
	return item->length;
}

int read_fastq(gzFile zfp, SEQ_QUAL *item, int id)
/*
default reading file use gzopen, because gzopen can
 handle ether gzip format file or not gzip format file, 
*/
{

	char line[MAX_READ_LENGTH];

	if(gzgetc(zfp) != '@'){
		if(gzeof(zfp) != 0) return -1;
		error_msg("Item %d is not a FASTQ format, \"@\" missed\n", id);
		return -1;
	}
	
	gzgets(zfp,item->name,MAX_READ_LENGTH);
	gzgets(zfp,item->seq,MAX_READ_LENGTH);

	gzgets(zfp,line,MAX_READ_LENGTH);
	if(line[0] != '+'){
		error_msg("Item %d is not a FASTQ format, \"+\" missed\n", id);
		return -1;
	}

	gzgets(zfp,item->qual,MAX_READ_LENGTH);

	item->id = id;
	item->length = strlen(item->seq) - 1;	//exlude '\n'
	item->start = 0;
	item->end = item->length - 1;
	item->name[strlen(item->name)-1] = '\0';
	item->seq[item->length] = '\0';
	item->qual[strlen(item->qual)-1] = '\0';
	
	return item->length;
}

int read_fastaq(gzFile zfps, gzFile zfpq, SEQ_QUAL *item, int id){

	char c='\n';
	char qual[64];
	int i=0,j=0;
	int max=BUFFER_LENGTH;

	if(read_fasta(zfps, item, id) < 0)	return -1;

	if(gzgetc(zfpq) != '>'){
		if(gzeof(zfpq)) return -1;
		error_msg("sequence %d has no a FASTA format quality", id);
		return -1;
	}
	
	while(!gzeof(zfpq) && (c=gzgetc(zfpq)) != '\n');

	while(!gzeof(zfpq) && (c=gzgetc(zfpq)) != '>'){
		if(c != '\n'){
			if(c != ' '){
				qual[j++] = c;
			}else{
				if(i+1 >= max){
					max += BUFFER_LENGTH;
					item->qual = realloc(item->qual, sizeof(char) * BUFFER_LENGTH);
				}
				qual[j]='\0';
				item->qual[i++] = atoi(qual);
				j=0;
			}
		}
	}

	if(!gzeof(zfpq))
		gzseek(zfpq,0,SEEK_CUR-1);

	item->id = id;
	item->start = 0;
	item->end = item->length - 1;
	return item->length;
}

#ifdef HAVA_BZLIB
int read_line_bz2(BZFILE* fp, char *line){

	char c;
	int bzerror;
	int i=0;
	do{
		BZ2_bzRead(&bzerror, fp, &c, sizeof(char));
		if(bzerror ==BZ_OK && c != '\n')
			line[i++]=c;
	}while(bzerror == BZ_OK && c != '\n');
	line[i]='\0';
	return bzerror;
}

int read_fasta_bz2(BZFILE* fp, SEQ_QUAL *item, int id){

	char c='\n';
	int i=0;
	int max=BUFFER_LENGTH;

	if(fgetc(fp) != '>'){
		if(feof(fp)) return -1;
		fprintf(stderr, "Warning: sequence %d is not a FASTA format\n", id);
		return -1;
	}
	
	while(!feof(fp) && (c=fgetc(fp)) != '\n'){
		if(c != '\n'){
			if(i+1 >= max){
				max += BUFFER_LENGTH;
				item->name = realloc(item->name, sizeof(char) * (max + BUFFER_LENGTH));
			}
			item->name[i++] = c;
		}
	}
	i=0;
	while(!feof(fp) && (c=fgetc(fp)) != '>'){
		if(isalpha(c)){
			if(i+1 >= max){
				max += BUFFER_LENGTH;
				item->seq = realloc(item->seq, sizeof(char) * (max + BUFFER_LENGTH));
			}
			item->seq[i++] = c;
		}
	}

	if(!feof(fp))	ungetc(c,fp);

	item->length = strlen(item->seq);
	item->id=id;
	return item->length;
}

int read_fastq_bz2(BZFILE* fp, SEQ_QUAL *item, int id){

	char line[MAX_READ_LENGTH];
	char c;
	int bzerror;
	int i=0;
	bzerror = read_line_bz2(fp, &line);
	if(line[0] != '@'){
		error_msg("Item %d is not a FASTQ format, missing character '@'\n", id);
		error_msg("read in : %s\n",line);
		return -1;
	}

	for(i=1;i<strlen(line);i++)
		item->name[i-1]=line[i];
	item->name[i-1]='\0';

	bzerror = read_line_bz2(fp,item->seq);
	bzerror = read_line_bz2(fp, &line);

	if(line[0] != '+'){
		error_msg("Item %d is not a FASTQ format, missing character '+'\n", id);
		error_msg("read in : %s\n",line);
		return -1;
	}

	bzerror = read_line_bz2(fp,item->qual);

	if(bzerror != BZ_OK) return -1;

	item->id = id;
	item->length = strlen(item->seq) - 1;	//exlude '\n'
	item->start = 0;
	item->end = item->length - 1;
	item->name[strlen(item->name)-1] = '\0';
	item->seq[item->length] = '\0';
	item->qual[strlen(item->qual)-1] = '\0';
	
	return item->length;
}
#endif

void output_fastq(FILE* out, SEQ_QUAL *item){
	int i = 0;
	fprintf(out,"@%s\n",item->name);
	for(i=item->start;i<=item->end;i++)
		fprintf(out,"%c",item->seq[i]);
	fprintf(out,"\n+\n");
	for(i=item->start;i<=item->end;i++)
		fprintf(out,"%c",item->qual[i]);

	fprintf(out,"\n");
}

void output_fasta(FILE* seq, FILE* qual, SEQ_QUAL *item){
	int i = 0;
	fprintf(seq,">%s\n",item->name);
	for(i=item->start;i<=item->end;i++)
		fprintf(seq,"%c",item->seq[i]);
	fprintf(seq,"\n");
	fprintf(qual,">%s\n",item->name);
	for(i=item->start;i<=item->end;i++)
		fprintf(qual,"%d ",item->qual[i]);
	fprintf(qual,"\b\n");
}

SEQ_QUAL init_read(void){
	SEQ_QUAL item;
	item.name   = malloc(sizeof(char) * MAX_READ_LENGTH);
	item.seq    = malloc(sizeof(char) * MAX_READ_LENGTH);
	item.qual   = malloc(sizeof(char) * MAX_READ_LENGTH);
	item.length = 0;
	item.start  = 0;
	item.end    = 0;
	item.id     = 0;
	return item;
}

int check_read(SEQ_QUAL *item, int mode){
	if(mode == 2){
		if(strlen(item->qual) < item->length){
			warning_msg("sequence %d : quality string is shorter than sequence string, quality length will be used as sequence length\n",item->id);
			item->length = strlen(item->qual);
			item->seq[item->length]='\0';
		}else if(strlen(item->qual) > item->length)
			warning_msg("sequence %d : quality string is longer than sequence string, the longer part will be discarded\n",item->id);
			item->qual[item->length]='\0';
	}

	if(mode == 3){
		if(strlen(item->qual) == 0){
			warning_msg("sequence %d : no quality string, quality will be assiated automatically\n",item->id);
			item->qual = realloc(item->qual, sizeof(char) * (item->length + 1));
			memset(item->qual, 'H', item->length);
			//item->qual[0]='\0';
		}
	}

	if(strlen(item->name) == 0){
		warning_msg("sequence %d is not complete: no sequence name specified, automatically associated\n", item->id);
		sprintf(item->name, "BSPT%6d", item->id);
	}
	if(strlen(item->seq) == 0){
		error_msg("sequence %d is not complete: no bases\n", item->id);
		return -1;
	}
	return 1;
}

void free_read(SEQ_QUAL *item){

	free(item->name);
	free(item->seq);
	free(item->qual);

}
