#include <zlib.h>

#include <bzlib.h>

#ifndef SEQ_H
#define SEQ_H

#define BASE_SANGER 33
#define BASE_SOLEXA 64
#define BUFFER_LENGTH 4096
#define MAX_READ_LENGTH 500
#define MAX_QUALITY 50
#define VERSION "1.5.2"

#ifndef MAX
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#endif

typedef struct {
	int id;			//read index
	int length;		//total length of sequence
	int start;		//start position of good quality part
	int end;		//end position of good quality part
	char *seq;		//all sequence base
	char *qual;		//all quality base
	char *name;		//sequence name
} SEQ_QUAL;

typedef enum {A,C,G,T,N} NI ;

enum fileType {
	FILE_NONE	 = 0x00000000,	//empty file or not readable file
	FILE_FASTQ	 = 0x00000001,	//FASTQ file
	FILE_FASTAQ	 = 0x00000002,	//FASTA + QUAL files
	FILE_FASTA	 = 0x00000004,	//FASTA file with no QUAL
	FILE_SFF	 = 0x00000008,	//sff file type
	FILE_PHRED33	 = 0x00000010,	//Sanger format (phred+33)
	FILE_PHRED64	 = 0x00000020,	//Solexa 1.3+ like format (phred+64)
	FILE_SOLEXA	 = 0x00000040,	//Solexa 1.2 like format
	FILE_UNQUAL	 = 0x00000080,	//unkown quality format
	FILE_GZ		 = 0x00000100,	//GZIP compressed file
	FILE_BZ2	 = 0x00000200,	//BZIP compressed file
	FILE_UNKN	 = 0x00000800,	//unkonwn file type
};

/* Basic sequence read function for different sequencing platform and 
file format.*/
// read FASTQ file using GZIP api 
int read_fastq(gzFile zfp, SEQ_QUAL *item, int id);

// read FASTA file using GZIP api 
int read_fasta(gzFile zfp, SEQ_QUAL *item, int id);

// read FASTAQ file using GZIP api 
int read_fastaq(gzFile zfps, gzFile zfpq, SEQ_QUAL *item, int id);

// read SFF file using GZIP api 
int read_sff(gzFile zfps, gzFile zfpq, SEQ_QUAL *item, int id);

#ifdef HAVE_BZLIB
// read FASTQ file using BZIP2 zpi 
int read_fastq_bz2(BZFILE* bzf, SEQ_QUAL *item, int id);

// read FASTA file using BZIP2 api 
int read_fasta_bz2(BZFILE* bzfps, SEQ_QUAL *item, int id);
#endif

//int readLine(FILE* fp, char *line);

void output_fastq(FILE* out, SEQ_QUAL *item);
void output_fasta(FILE* seq, FILE* qual, SEQ_QUAL *item);

int detect_filetype(char *file);

int notice_msg(const char *format, ...);
int warning_msg(const char *format, ...);
int error_msg(const char *format, ...);

int file_exist(char *file);
int dir_exist(char *path);
int file_name(char *file_path, char *file_name);

int mkdir_report(char *path);
FILE* fopen_report(char *file, char *mode);
gzFile gzopen_report(char *file, char *mode);
BZFILE* bzopen_report(char *file, char *mode);
//int type;	/*0->UNKONWN, 1->FASTQ, 2->FASTA*/

SEQ_QUAL init_read(void);
int check_read(SEQ_QUAL *item, int mode);

void free_read(SEQ_QUAL *item);
#endif
