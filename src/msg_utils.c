#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <errno.h>
#include <sys/time.h>

void notice_msg(const char *format, ...)
{

	va_list args;
	va_start(args, format);
	fprintf(stderr, "Notice: ");
	vfprintf(stderr, format, args);
	va_end(args);

}

void warning_msg(const char *format, ...)
{
  va_list args;

  /*if (!timeIsSet)
  {
    gettimeofday(&tvStart, NULL);
    timeIsSet = true;
  }

  gettimeofday(&tvNow, NULL);
  timersub(&tvNow, &tvStart, &tvDiff);
	*/
  //fprintf("[%ld.%06ld] ", (long) tvDiff.tv_sec, (long) tvDiff.tv_usec);
	fprintf(stderr,"Warning: ");
	va_start(args, format);
	vfprintf(stderr,format, args);
	va_end(args);

}

void error_msg(const char * format, ...) 
{
	va_list args;

	va_start(args, format);
	fprintf(stderr,"Error: ");
	vfprintf(stderr, format, args);
	va_end(args);
}
