#pragma once

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <libgen.h>
#include <math.h>
#include <time.h>

#include "define.h"

#define ERROR_MSG(msg) {                                     \
	fprintf(stderr, "[Error] %s, function %s, line %i: %s\n",\
					__FILE__, __func__, __LINE__, msg);      \
	exit(EXIT_FAILURE);                                      \
};

#define CHECK_PTR(ptr) {                                     \
	if (ptr == NULL) {                                       \
		ERROR_MSG("Null pointer");                           \
	}                                                        \
};

static inline int lm_index(int l, int m) {
	return l * (l + 1) + m;
}

static inline void myfree(void* ptr) 
{
	if (ptr != NULL) {
		free(ptr);
		ptr = NULL;
	}
}

static inline FILE* myfopen(char* path, char* mode)
{
	FILE* f = fopen(path, mode);
	
	if (f == NULL) {
	    fprintf(stderr, "Unable to open %s\n", path );
		exit(EXIT_FAILURE);
	}
	
	return f;
}
