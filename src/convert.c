#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>

#include "seq.h"

void convert_usage(){

	char* usage=
"Decsription: convert is one tool of bspt NGS data preprocess toolkits\n" \
"It detect file format, sequence format and sequencing platform\n" \
" automatically. Coversion is preformanced according detected format.\n" \
"When it processing compressed file, it read content directly from file\n" \
" and decompressed result will be outputed\n\n";

	printf("%s",usage);
	printf("Version:  %s by Yanpc\n\n",VERSION);
	printf("Usage:    bspt convert <file>\n\n");
}

int to_fasta(char* fastq){
	char tmp[128];
	int index = 0;
	char fname[128];
	file_name(&fname,fastq);
	gzFile zfp=gzopen_report(fastq,"r");
	if(!zfp)	return -1;
	SEQ_QUAL item = init_read();
	sprintf(tmp, "%s.fasta",fname);
	FILE *fs=fopen_report(tmp,"w");
	if(!fs)	return -1;
	sprintf(tmp, "%s.qual",fname);
	FILE *fq=fopen_report(tmp,"w");
	if(!fq)	return -1;
	while(read_fastq(zfp, &item, index++) >= 0){
		check_read(&item,2);
		output_fasta(fs, fq, &item);
	}
	free_read(&item);
	gzclose(zfp);
	fclose(fs);
	fclose(fq);
	return 0;
}

#ifdef HAVE_BZLIB
int to_fasta_bz2(char* fastq){
	char tmp[128];
	int index = 0;
	char fname[128];
	file_name(&fname,fastq);
	BZFILE *zfp=bzopen_report(fastq,"r");
	if(!zfp)	return -1; 
	sprintf(tmp, "%s.fasta",fname);
	FILE *fs=fopen_report(tmp,"w+");
	if(!fs)	return -1;
	sprintf(tmp, "%s.qual",fname);
	FILE *fq=fopen_report(tmp,"w+");
	if(!fq)	return -1;

	SEQ_QUAL item = init_read();
	while(read_fastq_bz2(zfp, &item, index++) >= 0){
		check_read(&item,2);
		output_fasta(fs, fq, &item);
	}
	free_read(&item);
	BZ2_bzclose(zfp);
	fclose(fs);
	fclose(fq);
	return 0;
}
#endif

int to_fastq(char* fasta, char* qual){
	SEQ_QUAL item = init_read();
	char tmp[128];
	int index = 0;
	FILE *fps=fopen_report(fasta,"r");
	if(!fps)	return -1;
	FILE *fpq=fopen_report(qual,"r");
	if(!fpq)	return -1;
	FILE *fq=fopen_report(strcat(fasta, ".fq"),"w+");
	if(!fq)	return -1;
	if(!fpq){
		warning_msg("can not find file %s with %s, trate as simple FASTA file and convert to FASTQ format\n", tmp,fasta);
		while(read_fasta(fps, &item, index++) >= 0){
			check_read(&item,3);
			output_fastq(fq, &item);
		}
	}else{
		while(read_fastaq(fps, fpq, &item, index++) >= 0){
			check_read(&item,2);
			output_fastq(fq, &item);
		}
		fclose(fpq);
	}
	free_read(&item);
	fclose(fps);
	fclose(fq);
	return 0;
}

int convert(int argc, char* argv[])
{

	if(argc < 2){
		convert_usage();
		exit(1);
	}

	char file[256];

	strncpy(file, argv[1],256);
	int type = detect_filetype(file);

	if(type & FILE_UNKN){
		error_msg("Unknown sequence format, not supported in this version\n");
		return 1;
	}

//start processing GZIP or regular FASTQ/FASTA format sequence

	if(type & FILE_FASTQ){
		notice_msg("\tFASTQ format convert to FASTA format in two files\n");
		if(type & FILE_BZ2)
#ifdef HAVE_BZLIB
			to_fasta_bz2(argv[1]);
#else
			error_msg("Do not support bzip2 file");
#endif
		else
			to_fasta(argv[1]);
			
	}else if(type & FILE_FASTAQ){
		notice_msg("FASTA format convert to FASTQ\n");
		to_fastq(argv[1],argv[2]);
	}

	return 0;
}
