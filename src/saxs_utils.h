/** @file \brief asdf
 */
#pragma once

#include "common.h"

#include "mol2/atom_group.h"
#include "mol2/sasa.h"
#include "mol2/transform.h"

#ifndef min
#define min(a,b)				\
	({ __typeof__ (a) _a = (a);		\
		__typeof__ (b) _b = (b);	\
		_a < _b ? _a : _b; })
#endif

#ifndef max
#define max(a,b)				\
	({ __typeof__ (a) _a = (a);		\
		__typeof__ (b) _b = (b);	\
		_a > _b ? _a : _b; })
#endif

typedef struct dec_ftresult_t
{
	struct ftresult *res;
	double rgyr;
	double chi;
	int rank;
} dec_ftresult;

double mol_atom_group_max_dist(const struct mol_atom_group *ag);

double mol_atom_group_average_radius(const struct mol_atom_group *ag);

double* mkarray(double begin, double end, int qnum);

void faccs(
	double *fractional_sa,
	const struct mol_atom_group *ag,
	double r_solv);
	
void upcase(char *s);

static inline double sinc(double x)
{
	return x == 0.0 ? 1 : sin(x)/x;
}

static inline double square(double x)
{
	return x*x;
}
