#ifndef _parser_h
#define _parser_h

#include "bsconfig.h"

double 
to_double (const char* s);

long
to_long (const char* s);

int
to_int (const char* s);

void
parse_parameters (double* S, 
		  double* E, 
		  double* r, 
		  double* sigma, 
		  double* T, 
		  long* M,
		  const char* filename);

void parse_parameters (BSConfig* config, const char* filename);

#endif /* _parser_h */
