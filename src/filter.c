#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <getopt.h>

#include "seq.h"

typedef struct {
	int qualityCutoff;
	int meanQuality;
	int lengthCutoff;
	int nNumber;
	char r1[128];
	char r2[128];
	char output[128];
	int percent;
} FLT_OPTS;

void filter_usage(){
	printf("\n");
	printf("Program:  filter function of NGS data preprocess toolkits\n");
	printf("Version:  %s by Yanpc\n\n", VERSION);
	printf("Usage:    pt filter [options]\n\n");
	printf("Options:  -q\tquality threshold [15]\n");
	printf("\t  -p\tmin percentage of base quality above -q, -q forced [0.3]\n\t  -m\tmean quality of read (supoior to -q and -p) [20]\n");
	printf("\t  -1\tsingle end or one end of pair-end\n");
	printf("\t  -2\tother end of pair-end\n\t  -o\tprefix of output\n");
	printf("\t  -h\tfor this help\n");
}

FLT_OPTS getFilterOptions(int argc, char* argv[]){
	char opts[]="o:1:2:hq:m:p:l:n:";
	int opt;

	if(argc < 2){
		filter_usage();
		exit(1);
	}
	FLT_OPTS mopts;
	mopts.qualityCutoff = 0;
	mopts.lengthCutoff = 0;
	mopts.nNumber = 99;
	mopts.meanQuality = 20;
	mopts.percent = 0;
	strcpy(mopts.r1, "none");
	strcpy(mopts.r2, "none");
	strcpy(mopts.output, "out");

	while ((opt = getopt(argc, argv, opts)) != -1) {
		switch (opt) {
			case 'o':	strncpy(mopts.output, optarg, 128);
				break;
			case 'q':	mopts.qualityCutoff = atoi(optarg);
				break;
			case '1':	strncpy(mopts.r1, optarg, 128);
				break;
			case '2':	strncpy(mopts.r2, optarg, 128);
				break;
			case 'p':	mopts.percent = atoi(optarg);
				break;
			case 'm':	mopts.meanQuality = atoi(optarg);
				break;
			case 'n':	mopts.nNumber = atoi(optarg);
				break;
			case 'l':	mopts.lengthCutoff = atoi(optarg);
				break;
			case 'h':	filter_usage();
				exit(1);
			//case '?':	printf("unknown option: %s \n", opt);
			//	break;
			default:	filter_usage();
				exit(1);
		}
	}

	if(strcmp(mopts.r1,"none") == 0){
		printf("must put a input file -1\n");
		exit(-1);
	}
	return mopts;
}

int filter_all(SEQ_QUAL *item, FLT_OPTS *flt){
	int left=1;
	int i=0;
	int argH=0;
	int argQ=0;
	int argN=0;
	//filter by length
	if(item->length < flt->lengthCutoff) return -1;
	for(i=0;i<item->length;i++){
		if(flt->qualityCutoff <= item->qual[i])
			argH++;
		argQ += item->qual[i];
		if(item->seq[i] == 'N' || item->seq[i] == 'n')
			argN++;
	}

	//filter by N number
	if(argN > flt->nNumber) return -1;
	//filter by mean quality
	if((argQ/item->length) < flt->meanQuality) return -1;
	//filter by high quality percent
	if((argH/item->length) < flt->percent) return -1;
	return left;
}

int filter_pe_fastq(FLT_OPTS *opts){
	int left1=0,left2=0;
	int stat_single1 = 0;
	int stat_single2 = 0;
	int stat_paired = 0;
	int index=1;
	char fn[128];
	char outfile[128];
	SEQ_QUAL item1=init_read();
	SEQ_QUAL item2=init_read();

	gzFile fp1=gzopen_report(opts->r1,"r");
	if(!fp1)	return -1;
	gzFile fp2=gzopen_report(opts->r2,"r");
	if(!fp2)	return -1;
	file_name(outfile,opts->r1);
	sprintf(fn,"%s/%s.flt",opts->output,outfile);
	FILE *fo1=fopen_report(fn,"w+");
	if(!fo1)	return -1;
	file_name(outfile,opts->r2);
	sprintf(fn,"%s/%s.flt",opts->output,outfile);
	FILE *fo2=fopen_report(fn,"w+");
	if(!fo2)	return -1;
	sprintf(fn,"%s/%s.flt.s",opts->output,outfile);
	FILE *fos=fopen_report(fn,"w+");
	if(!fos)	return -1;

	while(read_fastq(fp1,&item1,index) >= 0 && read_fastq(fp2,&item2,index) >= 0)
	{
		left1=filter_all(&item1, opts);
		left2=filter_all(&item2, opts);

		if(left1 == 1 && left2 == 1){
			output_fastq(fo1, &item1);
			output_fastq(fo2, &item2);
			stat_paired++;
		}else{
			if(left1 == 1){
				output_fastq(fos, &item1);
				stat_single1++;
			}
			if(left2 == 1){
				output_fastq(fos, &item2);
				stat_single2++;
			}
		}
		index++;
	}

	printf("Totally %d reads were processed\n",(index-1)*2);
	printf("  file [ %s ]: %d reads were left (%.2f%)\n",opts->r1,stat_paired+stat_single1,(float) (stat_paired+stat_single1)*100/(index-1));
	printf("  file [ %s ]: %d reads were left (%.2f%)\n",opts->r2,stat_paired+stat_single2,(float) (stat_paired+stat_single2)*100/(index-1));
	printf("After filtering %d reads are paired in each file (%.2f%)\n",stat_paired,(float) stat_paired*100/(index-1));
	printf("  file [ %s ]: %d reads were left as single end\n",opts->r1,stat_single1);
	printf("  file [ %s ]: %d reads were left as single end\n",opts->r2,stat_single2);
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
int filter_pe_fastq_bz2(FLT_OPTS *opts){
	int left1=0,left2=0;
	int stat_single1 = 0;
	int stat_single2 = 0;
	int stat_paired = 0;
	int index=1;
	char fn[128];
	char outfile[128];
	SEQ_QUAL item1=init_read();
	SEQ_QUAL item2=init_read();

	BZFILE *fp1=bzopen_report(opts->r1,"r");
	BZFILE *fp2=bzopen_report(opts->r2,"r");
	file_name(outfile,opt->r1);
	sprintf(fn,"%s/%s.flt",opts->output,outfile);
	FILE *fo1=fopen_report(fn,"w+");
	if(!fo1)	return -1;
	file_name(outfile,opts->r2);
	sprintf(fn,"%s/%s.flt",opts->output,outfile);
	FILE *fo2=fopen_report(fn,"w+");
	if(!fo2)	return -1;
	sprintf(fn,"%s/%s.flt.s",opts->output,outfile);
	FILE *fos=fopen_report(fn,"w+");
	if(!fos)	return -1;

	while(read_fastq_bz2(fp1,&item1,index) >= 0 && read_fastq_bz2(fp2,&item2,index) >= 0)
	{
		left1=filter_all(&item1, opts);
		left2=filter_all(&item2, opts);

		if(left1 == 1 && left2 == 1){
			output_fastq(fo1, &item1);
			output_fastq(fo2, &item2);
			stat_single1++;
			stat_single2++;
			stat_paired++;
		}else{
			if(left1 == 1){
				output_fastq(fos, &item1);
				stat_single1++;
			}
			if(left2 == 1){
				output_fastq(fos, &item2);
				stat_single2++;
			}
		}
		index++;
	}

	printf("Totally %d reads were processed\n",(index-1)*2);
	printf("  file [ %s ]: %d reads were left (%.2f%)\n",opts->r1,stat_paired+stat_single1,(float) (stat_paired+stat_single1)*100/(index-1));
	printf("  file [ %s ]: %d reads were left (%.2f%)\n",opts->r2,stat_paired+stat_single2,(float) (stat_paired+stat_single2)*100/(index-1));
	printf("After filtering %d reads are paired in each file (%.2f%)\n",stat_paired,(float) stat_paired*100/(index-1));
	printf("  file [ %s ]: %d reads were left as single end\n",opts->r1,stat_single1);
	printf("  file [ %s ]: %d reads were left as single end\n",opts->r2,stat_single2);
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

int filter_se_fastq(FLT_OPTS *opts){
	int index=1;
	int stat_left = 0;

	SEQ_QUAL item=init_read();
	gzFile fp=gzopen_report(opts->r1,"r");
	FILE *fo=fopen_report(strcat(opts->output,".flt"),"w+");
	int left = 0;

	while(read_fastq(fp,&item,index++) > 0){
		check_read(&item,2);
		left=filter_all(&item, opts);

		if(left == 1){
			output_fastq(fo, &item);
			stat_left++;
		}
	}

	printf("Totally %d reads were processed\n",index-1);
	printf("  file [ %s ]: %d reads were left (%.2f%)\n",opts->r1,stat_left,(float) stat_left*100/(index-1));
	free_read(&item);
	fclose(fp);
	fclose(fo);
	
	return 0;
}	

#ifdef HAVA_BZLIB
int filter_se_fastq_bz2(FLT_OPTS *opts){
	int index=1;
	int stat_left = 0;

	SEQ_QUAL item=init_read();
	gzFile fp=gzopen_report(opts->r1,"r");
	FILE *fo=fopen_report(strcat(opts->output,".flt"),"w+");
	int left = 0;

	while(read_fastq(fp,&item,index++) > 0){
		check_read(&item,2);
		left=filter_all(&item, opts);

		if(left == 1){
			output_fastq(fo, &item);
			stat_left++;
		}
	}

	printf("Totally %d reads were processed\n",index-1);
	printf("  file [ %s ]: %d reads were left (%.2f%)\n",opts->r1,stat_left,(float) stat_left*100/(index-1));
	free_read(&item);
	fclose(fp);
	fclose(fo);
	
	return 0;
}	
#endif

int filter_se_fastaq(FLT_OPTS *opts){
	warning_msg("SFF format file is not supported in this version\n");
	return 1;
}

int filter_se_sff(FLT_OPTS *opts){
	warning_msg("SFF format file is not supported in this version\n");
	return 1;
}

int filter_main(int argc, char* argv[]){
	FLT_OPTS mopts=getFilterOptions(argc,argv);
	int flag_a = detect_filetype(mopts.r1);

	if(flag_a & FILE_NONE || flag_a & FILE_UNKN || flag_a & FILE_FASTA)
		return 1;

	mopts.qualityCutoff += flag_a & FILE_PHRED33 ? BASE_SANGER : BASE_SOLEXA;
	mopts.meanQuality += flag_a & FILE_PHRED33 ? BASE_SANGER : BASE_SOLEXA;
	if(mkdir_report(mopts.output) != 0)	return 1;

	if(strcmp(mopts.r2,"none") != 0){
		int flag_b = detect_filetype(mopts.r2);
		if(flag_a != flag_b){
			fprintf(stderr,"This two files:%s,%s are not same type\n",mopts.r1,mopts.r2);
			return 1;
		}
		if( flag_a & FILE_FASTQ ) {
			if( flag_a & FILE_BZ2 )
#ifdef HAVE_BZLIB
				filter_pe_fastq_bz2(&mopts);
#else
				error_msg("Do not support bzip2 file");
#endif
			else
				filter_pe_fastq(&mopts);
		}else{
			error_msg("unkonwn file type");
		}

	}else{

		if( flag_a & FILE_FASTQ ) {
			if( flag_a & FILE_BZ2 )
#ifdef HAVE_BZLIB
				filter_se_fastq_bz2(&mopts);
#else
				error_msg("Do not support bzip2 file");
#endif
			else
				filter_se_fastq(&mopts);
		}else if( flag_a & FILE_SFF )
			filter_se_sff(&mopts);
		else if( flag_a & FILE_FASTAQ )
			filter_se_fastaq(&mopts);
		else
			error_msg("unkonwn file type");
	}
	
	return 0;
}	
