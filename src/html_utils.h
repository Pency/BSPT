#include <stdio.h>

#include "stats.h"

void flotr2_js(FILE* output);
void barplot(char* figure, char* data, FILE* output);
void boxplot(char* figure, char* data, FILE* output);
void plot_cycle(FILE* output, cycle *cycles);
void print_stats(FILE* output, stats *seq_stats);
void head(FILE* output);
void foot(FILE* output);
