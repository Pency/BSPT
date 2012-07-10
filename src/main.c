#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>

#include "main.h"
#include "seq.h"

void Usage(){
	printf("Program:  NGS data preprocess toolkits\n");
	printf("Version:  %s by Yanpc\n\n",VERSION);
	printf("Usage:    bspt <command> [options]\n\n");
	printf("Command:  convert\tSequence format converter\n");
	printf("\t  filter\tquality filter\n");
	printf("\t  stats \tstatistic data\n");
	printf("\t  trim  \tquality trim \n");
}

int main(int argc, char* argv[]){
    
	if(argc < 2){
		Usage();
		return -1;
	}

	if (strcmp(argv[1], "convert") == 0) return convert(argc-1,argv+1);
	else if (strcmp(argv[1], "filter") == 0) return filter_main(argc-1, argv+1);
	else if (strcmp(argv[1], "stats") == 0) return stats_main(argc-1,argv+1);
	else if (strcmp(argv[1], "trim") == 0) return trim_main(argc-1,argv+1);
	else{
			Usage();
			return 1;
	}
}	
