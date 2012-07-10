#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "seq.h"

int file_exist(char *file)
{
	struct stat buf;
	if(stat(file,&buf) != 0 || S_ISDIR(buf.st_mode))
		return -1;

	return 0;
}

int dir_exist(char *path)
{
	struct stat buf;
	if(stat(path,&buf) != 0 || !S_ISDIR(buf.st_mode))
		return -1;

	return 0;
}

int file_name(char* file_name, char* file_path)	//parse full path string and get file name
{
	struct stat buf;
	if(stat(file_path, &buf) == 0 && S_ISDIR(buf.st_mode)){
		warning_msg("This path [%s] is not contain file name.\n",file_path);
		return -1;
	}

	int i = 0;
	int lgth = strlen(file_path);
	if(lgth <= 0 )	return -1;

	for(i=lgth; i>=0; i--)
		if(file_path[i] == '/')
			break;

	int start = i+1;
	for(i=0; (start+i)<=lgth; i++)
		file_name[i] = file_path[start+i];

	file_name[i]='\0';
	return 0;
}

int mkdir_report(char *path)
{
	struct stat buf;
	if(stat(path, &buf) != 0){
		warning_msg("There is no directory [%s], will be created!\n",path);
		if(mkdir(path, 10705) != 0){
			error_msg("Directory [%s] can not be created!\n",path);
			return -1;
		}else{
			notice_msg("Directory has been created successfully!\n");
		}
	}else if(!S_ISDIR(buf.st_mode)){
			error_msg("%s is not dir\n",path);
			return -1;
	}
	return 0;
}

FILE* fopen_report(char *file, char *mode)
{

	if(strcmp(mode,"r") == 0){
		if(file_exist(file) != 0){
			error_msg("file %s is not exist\n", file);
			return (FILE*) NULL;
		}
	}

	if(strcmp(mode,"w") == 0){
		if(file_exist(file) == 0){
			error_msg("file %s is already exist\n", file);
			return (FILE*) NULL;
		}
	}

	FILE *file_handle=fopen(file,mode);
	if(!file_handle){
		error_msg("Can not open the file %s\n", file);
	}
	return file_handle;
}

gzFile gzopen_report(char *file, char *mode)
{

	if(strcmp(mode,"r") == 0){
		if(file_exist(file) != 0){
			error_msg("file %s is not exist\n", file);
			return (gzFile) NULL;
		}
	}

	if(strcmp(mode,"w") == 0){
		if(file_exist(file) == 0){
			error_msg("file %s is already exist\n", file);
			return (gzFile) NULL;
		}
	}

	gzFile file_handle=gzopen(file,mode);
	if(!file_handle){
		error_msg("Can not open the file %s\n", file);
	}
	return file_handle;
}

#ifdef HAVE_BZLIB
BZFILE* bzopen_report(char *file, char *mode)
{

	if(strcmp(mode,"r") == 0){
		if(file_exist(file) != 0){
			error_msg("file %s is not exist\n", file);
			return (BZFILE*) NULL;
		}
	}

	if(strcmp(mode,"w") == 0){
		if(file_exist(file) == 0){
			error_msg("file %s is already exist\n", file);
			return (BZFILE*) NULL;
		}
	}

	BZFILE *file_handle=BZ2_bzopen(file,mode);
	if(!file_handle){
		error_msg("Can not open the file %s\n", file);
	}
	return file_handle;
}
#endif

int is_gz(char *file)
/* according to (RFC 1952)
GZIP file format specification version 4.3, gzip header include 
ID1 (IDentification 1) and ID2 (IDentification 2).
These have the fixed values ID1 = 31 (0x1f, \037), ID2 = 139
 (0x8b, \213), to identify the file as being in gzip format.
*/
{
	unsigned int ID1=0;
	unsigned int ID2=0;
	FILE *fp=fopen_report(file,"r");
	fread(&ID1, 1, 1, fp);
	fread(&ID2, 1, 1,fp);
	fclose(fp);
	if(ID1 == 31 && ID2 == 139)
		return 1;
	return 0;
}

int is_bz2(char *file){
	char magic[3];
	FILE *fp=fopen_report(file,"r");
	fread(&magic, 1, 16, fp);
	magic[2]='\0';
	fclose(fp);
	if(strcmp(magic,"BZ") == 0)
		return 1;
	return 0;
}

void report_filetype(int type)
{
	char file_format[20] = "";
	char seq_format[50] = "";
	if(type & FILE_FASTQ){
		if(type & FILE_PHRED33)		strcpy(seq_format, "FASTQ with phred+33 quality score");
		if(type & FILE_PHRED64)		strcpy(seq_format, "FASTQ with phred+64 quality score");
	}else if(type & FILE_FASTA){
		strcpy(seq_format, "FASTA with no quality score");
		if(type & FILE_PHRED33)		strcpy(seq_format, "FASTA with phred+33 quality score");
		if(type & FILE_PHRED64)		strcpy(seq_format, "FASTA with phred+64 qualtiy score");
	}
	if(type & FILE_GZ)		sprintf(file_format, "compressed in GZIP");
	if(type & FILE_BZ2)		sprintf(file_format, "compressed in BZIP");

	notice_msg("%s %s\n", seq_format,file_format);

}

int detect_datatype(char *file){

	int dataType=0;
	int i = 0;
	int max = 0;
	int min = 999;
	int sample = 100;

	SEQ_QUAL item = init_read();
	gzFile zfp = gzopen_report(file,"r");

	if(gzgetc(zfp) == '>'){
		dataType |= FILE_FASTA;
	}else{
		gzseek(zfp, 0L, SEEK_SET);
		if(read_fastq(zfp,&item,i) >= 0){
			dataType |= FILE_FASTQ;
			do{
				for(i=0;i+1<strlen(item.qual);i++){
                                        min = MIN(min, item.qual[i]);
                                        max = MAX(max, item.qual[i]);
                                }
                                if((sample--) == 0){
                                        if(max >= 75)
						dataType|=FILE_PHRED64;
                                        else{
						dataType|=FILE_PHRED33;
						if(min > 58)
							warning_msg("Can not identified quality score type in 100 read samples, assume phred+33\n");
					}
					break;
				}
			}while(read_fastq(zfp,&item, i));
		}else
			dataType |= FILE_UNKN;
	}

	gzclose(zfp);
	free_read(&item);
	return dataType;
}

#ifdef HAVE_BZLIB
int detect_datatype_bz2(char *file){

	int dataType=0;
	int i = 0;
	int max = 0;
	int min = 999;
	int sample = 100;

	SEQ_QUAL item = init_read();
	BZFILE *zfp = bzopen_report(file,"r");

	if(read_fastq_bz2(zfp,&item,i) >= 0){
		dataType |= FILE_FASTQ;
		do{
			for(i=0;i+1<strlen(item.qual);i++){
				min = MIN(min, item.qual[i]);
				max = MAX(max, item.qual[i]);
			}
			if((sample--) == 0){
				if(max >= 75)
					dataType|=FILE_PHRED64;
				else{
					dataType|=FILE_PHRED33;
					if(min > 58)
						warning_msg("Can not identified quality score type in 100 read samples, assume phred+33\n");
				}
				break;
			}
		}while(read_fastq_bz2(zfp,&item, i));

	}else{
		dataType |= FILE_UNKN;
	}
	BZ2_bzclose(zfp);
	free_read(&item);
	return dataType;
}
#endif

int detect_filetype(char *file){

	int fileType=0;

	if(is_bz2(file) > 0){
#ifdef HAVE_BZLIB
		fileType = detect_datatype_bz2(file);
#endif
		fileType|= FILE_BZ2;
	}else{
		fileType = detect_datatype(file);
		if(is_gz(file) > 0)
			fileType|= FILE_GZ;
	}

	report_filetype(fileType);

	return fileType;
}
