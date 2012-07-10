#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>

#include "seq.h"

typedef struct {
	int qual_cutoff;
	int len_cutoff;
	char r1[128];
	char r2[128];
	char output[128];
	int no_head;
	int mode;
} TRIM_OPTS;

void trim_usage(){
	printf("Program:  trim function of NGS data preprocess toolkits\n");
	printf("Version:  %s by Yanpc\n\n","0.1.3");
	printf("Usage:    pt trim [options]\n\n");
	printf("Options:  -l\tmin length of read after trimming\n");
	printf("\t  -q\tquality threshold\n\t  -1\tsingle end or one end of pair-end\n");
	printf("\t  -2\tother end of pair-end\n\t  -o\toutput dir\n");
	printf("\t  -t\ttrim mode[1,2,3], default is 1\n");
	//printf("\t  -v\twhether to trim first base\n");
	printf("\t  -h\tfor this help\n");
}

TRIM_OPTS get_trim_options(int argc, char* argv[]){

	if(argc < 2){
		trim_usage();
		exit(1);
	}

	int opt;
	TRIM_OPTS mopts;
	mopts.qual_cutoff = 15;
	mopts.len_cutoff = 25;
	mopts.no_head = 0;
	mopts.mode = 1;
	strcpy(mopts.r1, "none");
	strcpy(mopts.r2, "none");
	strcpy(mopts.output, ".");

	while ((opt = getopt(argc, argv, "o:1:2:hq:l:t:v")) != -1) {
		switch (opt) {
			case 'o':	strncpy(mopts.output, optarg, 128);
				break;
			case 'q':	mopts.qual_cutoff = atoi(optarg);
				break;
			case '1':	strncpy(mopts.r1, optarg, 128);
				break;
			case '2':	strncpy(mopts.r2, optarg, 128);
				break;
			case 'v':	mopts.no_head = 1;
				break;
			case 't':	mopts.mode = atoi(optarg);
				break;
			case 'l':	mopts.len_cutoff = atoi(optarg);
				break;
			case 'h':	trim_usage();
				exit(1);
			default:	trim_usage();
				exit(1);
		}
	}

	if(strcmp(mopts.r1,"none") == 0){
		printf("Option -1 must be set a input file\n");
		exit(-1);
	}
	if(mopts.mode < 0 || mopts.mode > 3){
		printf("Option -t 1, 2 or 3\n");
		exit(-1);
	}
	if(mopts.len_cutoff < 0){
		printf("Option -l must be a positive number\n");
		exit(-1);
	}
	if(mopts.qual_cutoff < 0 || mopts.qual_cutoff > 47){
		printf("Option -q must between 0 and 47\n");
		exit(-1);
	}
	return mopts;
}

int QATrim(SEQ_QUAL *item, int threshold)
// trim mode 3, find out longest sequencial high quality sequence,
// rest parts of sequence are trimmed.
{
	int left_len=0;
	int i=0,start=0,tstart=0;
	int tend=0;
	int flag = 0;
	int sl=0;
	for(i=0;i<item->length;i++){
		if(flag == 0){
			if(item->qual[i] >= threshold){
				tstart=i;
				flag = 1;
			}
		}else{
			if(item->qual[i] < threshold){
				tend = i - 1;
				flag = 0;
				if(sl < (tend - tstart + 1)){
					sl = tend - tstart;
					start = tstart;
				}
			}
		}
	}
	if(flag == 1){
		tend = i -1;
		if(sl < (tend - tstart + 1)){
			sl = tend - tstart;
			start = tstart;
		}
	}
	item->start = start;
	item->end = start + sl;
	return (start + item->length - left_len);
}

int BWATrim(SEQ_QUAL *item, int threshold)
// Trim mode 1, BWA algorithm, trimming from right side
{
	int left_len = (item->length - 1);
	int i=0,j=0;
	int argmax=0;
	int arg=0;
	for(i=0;i<strlen(item->qual);i++){
		arg = 0;
		for(j=i+1;j<strlen(item->qual);j++)
			arg += threshold - item->qual[j];
		if(argmax < arg){
			left_len = i;
			argmax = arg;
		}
	}
	item->start = 0;
	item->end = left_len;
	return 0;
}

int YPCTrim(SEQ_QUAL *item, int threshold)
// Trim mode 2, adjusted BWA algorithm, trimming from both sides
{
	int left_len = (item->length - 1);
	int start = 0;
	int i=0,j=0;
	int argmax=0,argmin=999999;
	int arg=0;
	for(i=0;i<strlen(item->qual);i++){
		arg = 0;
		for(j=i+1;j<strlen(item->qual);j++)
			arg += threshold - item->qual[j];
		if(argmax < arg){
			left_len = i;
			argmax = arg;
		}
		if(argmin > arg){
			start = i;
			argmin = arg;
		}
	}
	item->start = start;
	item->end = left_len;
	if(argmin < 0)
		item->start++;
	return start;
}

int trim_pe_fastq(TRIM_OPTS *opt){

	int index = 1;
	int stat_single1 = 0;
	int stat_single2 = 0;
	int stat_paired = 0;
	char fn[128];
	char outfile[128];

	gzFile fp1=gzopen_report(opt->r1,"r");
	if(!fp1)	return -1;
	gzFile fp2=gzopen_report(opt->r2,"r");
	if(!fp2)	return -1;
	file_name(outfile,opt->r1);
	sprintf(fn,"%s/%s.trm",opt->output,outfile);
	FILE *fo1=fopen_report(fn,"w+");
	if(!fo1)	return -1;
	file_name(outfile,opt->r2);
	sprintf(fn,"%s/%s.trm",opt->output,outfile);
	FILE *fo2=fopen_report(fn,"w+");
	if(!fo2)	return -1;
	sprintf(fn,"%s/%s.trm.s",opt->output,outfile);
	FILE *fos=fopen_report(fn,"w+");
	if(!fos)	return -1;
   	SEQ_QUAL item1=init_read();
	SEQ_QUAL item2=init_read();

	while(read_fastq(fp1,&item1,index) > 0 && read_fastq(fp2,&item2,index) > 0){
		if(opt->mode > 2){
			if(opt->mode == 2){
				YPCTrim(&item1, opt->qual_cutoff);
				YPCTrim(&item2, opt->qual_cutoff);
			}else{
				QATrim(&item1, opt->qual_cutoff);
				QATrim(&item2, opt->qual_cutoff);
			}
		}else{
			if(opt->mode == 1){
				BWATrim(&item1, opt->qual_cutoff);
				BWATrim(&item2, opt->qual_cutoff);
			}
		}

		if((item1.end - item1.start +1) >= opt->len_cutoff)
			if((item2.end - item2.start +1) >= opt->len_cutoff){
				stat_paired++;
				output_fastq(fo1, &item1);
				output_fastq(fo2, &item2);
			}else{
				stat_single1++;
				output_fastq(fos, &item1);
			}
		else
			if((item2.end - item2.start +1) >= opt->len_cutoff){
				stat_single2++;
				output_fastq(fos, &item2);
			}

		index++;
	}
	printf("Totally %d reads were processed\n",(index-1)*2);
	printf("  file [ %s ]: %d reads were left (%.2f%)\n",opt->r1,stat_paired+stat_single1,(float) (stat_paired+stat_single1)*100/(index-1));
	printf("  file [ %s ]: %d reads were left (%.2f%)\n",opt->r2,stat_paired+stat_single2,(float) (stat_paired+stat_single2)*100/(index-1));
	printf("After trimming %d reads are paired in each file (%.2f%)\n",stat_paired,(float) stat_paired*100/(index-1));
	printf("  file [ %s ]: %d reads were left as single end\n",opt->r1,stat_single1);
	printf("  file [ %s ]: %d reads were left as single end\n",opt->r2,stat_single2);
	free_read(&item1);
	free_read(&item2);
	gzclose(fp1);
	gzclose(fp2);
	fclose(fo1);
	fclose(fo2);
	fclose(fos);

	return 0;
}

#ifdef HAVE_BZLIB
int trim_pe_fastq_bz2(TRIM_OPTS *opt){

	int index = 1;
	int stat_single1 = 0;
	int stat_single2 = 0;
	int stat_paired = 0;
	char fn[128];
	char outfile[128];

	BZFILE *fp1=bzopen_report(opt->r1,"r");
	if(!fp1)	return -1;
	BZFILE *fp2=bzopen_report(opt->r2,"r");
	if(!fp2)	return -1;
	file_name(outfile,opt->r1);
	sprintf(fn,"%s/%s.trm",opt->output,outfile);
	FILE *fo1=fopen_report(fn,"w+");
	if(!fo1)	return -1;
	file_name(outfile,opt->r2);
	sprintf(fn,"%s/%s.trm",opt->output,outfile);
	FILE *fo2=fopen_report(fn,"w+");
	if(!fo2)	return -1;
	sprintf(fn,"%s/%s.trm.s",opt->output,outfile);
	FILE *fos=fopen_report(fn,"w+");
	if(!fos)	return -1;
   	SEQ_QUAL item1=init_read();
	SEQ_QUAL item2=init_read();

	while(read_fastq_bz2(fp1,&item1,index) > 0 && read_fastq_bz2(fp2,&item2,index) > 0){
		if(opt->mode > 2){
			if(opt->mode == 2){
				YPCTrim(&item1, opt->qual_cutoff);
				YPCTrim(&item2, opt->qual_cutoff);
			}else{
				QATrim(&item1, opt->qual_cutoff);
				QATrim(&item2, opt->qual_cutoff);
			}
		}else{
			if(opt->mode == 1){
				BWATrim(&item1, opt->qual_cutoff);
				BWATrim(&item2, opt->qual_cutoff);
			}
		}

		if((item1.end - item1.start +1) >= opt->len_cutoff)
			if((item2.end - item2.start +1) >= opt->len_cutoff){
				stat_paired++;
				output_fastq(fo1, &item1);
				output_fastq(fo2, &item2);
			}else{
				stat_single1++;
				output_fastq(fos, &item1);
			}
		else
			if((item2.end - item2.start +1) >= opt->len_cutoff){
				stat_single2++;
				output_fastq(fos, &item2);
			}

		index++;
	}
	printf("Totally %d reads were processed\n",(index-1)*2);
	printf("  file [ %s ]: %d reads were left (%.2f%)\n",opt->r1,stat_paired+stat_single1,(float) (stat_paired+stat_single1)*100/(index-1));
	printf("  file [ %s ]: %d reads were left (%.2f%)\n",opt->r2,stat_paired+stat_single2,(float) (stat_paired+stat_single2)*100/(index-1));
	printf("After trimming %d reads are paired in each file (%.2f%)\n",stat_paired,(float) stat_paired*100/(index-1));
	printf("  file [ %s ]: %d reads were left as single end\n",opt->r1,stat_single1);
	printf("  file [ %s ]: %d reads were left as single end\n",opt->r2,stat_single2);
	free_read(&item1);
	free_read(&item2);
	BZ2_bzclose(fp1);
	BZ2_bzclose(fp2);
	fclose(fo1);
	fclose(fo2);
	fclose(fos);

	return 0;
}
#endif

int trim_se_fastq(TRIM_OPTS *opt){
	int index = 1;
	int stat_left = 0;
	char fn[128];
	char outfile[128];
	file_name(outfile, opt->r1);
	sprintf(fn,"%s/%s.trm",opt->output,outfile);
	FILE *fo=fopen_report(fn,"w+");
	gzFile fp=gzopen_report(opt->r1,"r");
	SEQ_QUAL item=init_read();

	while(read_fastq(fp,&item, index) > 0){
		check_read(&item,2);
		if(opt->mode > 2)
			if(opt->mode == 2)
				YPCTrim(&item, opt->qual_cutoff);
			else
				QATrim(&item, opt->qual_cutoff);
		else
			if(opt->mode == 1)
				BWATrim(&item, opt->qual_cutoff);
		if((item.end - item.start +1) >= opt->len_cutoff){
			stat_left++;
			output_fastq(fo, &item);
		}
		index++;
	}
	printf("Totally %d reads were processed\n",index-1);
	printf("  file [ %s ]: %d reads were left (%.2f%)\n",opt->r1,stat_left,(float) stat_left*100/(index-1));
	free_read(&item);
	gzclose(fp);
	fclose(fo);
	
	return 0;
}	

#ifdef HAVE_BZLIB
int trim_se_fastq_bz2(TRIM_OPTS *opt){
	int index = 1;
	int stat_left = 0;
	char fn[128];
	char outfile[128];
	file_name(outfile, opt->r1);
	sprintf(fn,"%s/%s.trm",opt->output,outfile);
	FILE *fo=fopen_report(fn,"w+");
	if(!fo)	return -1;
	BZFILE *fp=bzopen_report(opt->r1,"r");
	if(!fp)	return -1;
	SEQ_QUAL item=init_read();

	while(read_fastq_bz2(fp,&item, index) > 0){
		check_read(&item,2);
		if(opt->mode > 2)
			if(opt->mode == 2)
				YPCTrim(&item, opt->qual_cutoff);
			else
				QATrim(&item, opt->qual_cutoff);
		else
			if(opt->mode == 1)
				BWATrim(&item, opt->qual_cutoff);
		if((item.end - item.start +1) >= opt->len_cutoff){
			stat_left++;
			output_fastq(fo, &item);
		}
		index++;
	}
	printf("Totally %d reads were processed\n",index-1);
	printf("  file [ %s ]: %d reads were left (%.2f%)\n",opt->r1,stat_left,(float) stat_left*100/(index-1));
	free_read(&item);
	BZ2_bzclose(fp);
	fclose(fo);
	
	return 0;
}	
#endif

int trim_se_fastaq(TRIM_OPTS *opt){
	SEQ_QUAL item = init_read();
	int index = 1;
	int stat_left=0;
	gzFile fs=gzopen_report(opt->r1,"r");
	if(!fs)	return -1;
	FILE *fso=fopen(strcat(opt->output,".trm"),"w+");
	if(!fso)	return -1;
	gzFile fq=gzopen_report(opt->r1,"r");
	if(!fq)	return -1;
	FILE *fqo=fopen(strcat(opt->output,".trm"),"w+");
	if(!fqo)	return -1;

	while(read_fastaq(fs,fq,&item, index) > 0){
		check_read(&item,2);
		if(opt->mode > 2){
			if(opt->mode == 2){
				YPCTrim(&item, opt->qual_cutoff);
			}else{
				QATrim(&item, opt->qual_cutoff);
			}
		}else{
			if(opt->mode == 1)
				BWATrim(&item, opt->qual_cutoff);
		}
		if(item.length >= opt->len_cutoff)
			output_fastq(fso, &item);
		index++;
	}
	printf("Totally %d reads were processed\n",index-1);
	printf("  file [ %s ]: %d reads were left (%.2f%)\n",opt->r1,stat_left,(float) stat_left*100/(index-1));
	free_read(&item);
	gzclose(fs);
	gzclose(fq);
	fclose(fso);
	fclose(fqo);
	
	return 0;
}	

int trim_se_sff(TRIM_OPTS *opt){
	warning_msg("SFF format file is not supported in this version\n");
	return 1;
}

int trim_main(int argc, char* argv[]){

	TRIM_OPTS opt = get_trim_options(argc, argv);
	
	//fprintf(stderr, "Auto-detect file format:\n");
	int flag_a = detect_filetype(opt.r1);

	if(flag_a == -1)	return 1;

	if(flag_a & FILE_NONE || flag_a & FILE_UNKN || flag_a & FILE_FASTA)
		return 1;

	if(mkdir_report(opt.output) != 0)	return 1;

	opt.qual_cutoff += flag_a & FILE_PHRED33 ? BASE_SANGER : BASE_SOLEXA;

	if(strcmp(opt.r2,"none") != 0){
		int flag_b = detect_filetype(opt.r2);

		if(flag_a != flag_b){
			fprintf(stderr,"This two files:%s,%s are not same type\n",opt.r1,opt.r2);
			return 1;
		}

		if( flag_a & FILE_FASTQ ) {
			if( flag_a & FILE_BZ2 )
#ifdef HAVE_BZLIB
				trim_pe_fastq_bz2(&opt);
#else
				error_msg("Do not support bzip2 file");
#endif
			else
				trim_pe_fastq(&opt);
		}else{
			error_msg("unkonwn file type");
		}

	}else{

		if( flag_a & FILE_FASTQ ) {
			if( flag_a & FILE_BZ2 )
#ifdef HAVE_BZLIB
				trim_se_fastq_bz2(&opt);
#else
				error_msg("Do not support bzip2 file");
#endif
			else
				trim_se_fastq(&opt);
		}else if( flag_a & FILE_SFF ){
			trim_se_sff(&opt);
		}else if( flag_a & FILE_FASTAQ ){
			trim_se_fastaq(&opt);
		}else{
			error_msg("unkonwn file type");
		}

	}
	
	return 0;
}	

#ifdef TEST
int main(int argc, char* argv[]){
	trim(argc,argv);
}
#endif
