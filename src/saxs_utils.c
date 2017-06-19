#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "saxs_utils.h"

double mol_atom_group_max_dist(const struct mol_atom_group *ag)
{
	double max_d = 0.0;
	double d;

	for (size_t i = 0; i < ag->natoms; i++) {
		for (size_t j = i + 1; j < ag->natoms; j++) {
			d = MOL_VEC_EUCLIDEAN_DIST_SQ(ag->coords[i], ag->coords[j]);

			max_d = fmax(max_d, d);
		}
	}

	return sqrt(max_d);
}

double mol_atom_group_average_radius(const struct mol_atom_group *ag)
{
	double av_r = 0.0;

	for (size_t i = 0; i < ag->natoms; i++) {
		av_r += ag->vdw_radius[i];
	}

	return av_r / ag->natoms;
}

void faccs(
	double *fractional_sa,
	const struct mol_atom_group *ag,
	double r_solv)
{
	static const short cont_acc = 1;
	static const double fourpi = 4 * M_PI;

	accs(fractional_sa, ag, r_solv, cont_acc);

	for (size_t i = 0; i < ag->natoms; ++i) {
		// sa expected by saxs is the fraction of accessible surface area
		double ri = ag->vdw_radius[i];
		if (ri == 0 || isnan(fractional_sa[i]))
			fractional_sa[i] = 0.0;
		else
			fractional_sa[i] = fractional_sa[i] / (fourpi * ri*ri);
	}
}

void upcase(char *s)
{
	while(*s) {
		*s = toupper(*s);
		s++;
	}
}

double* mkarray(double begin, double end, int qnum)
{
	if (begin > end || begin < 0 || end < 0)
	{
		fprintf(stderr, "WRONG Q VALUES WHEN CREATING ARRAY");
		exit(EXIT_FAILURE);
	}
	double step = (end - begin)/(int)(qnum + 1);
	double* q = (double*)calloc(qnum, sizeof(double));
	q[0] = begin;
	for (int i = 1; i < qnum; i++)
		q[i] = q[i-1] + step;
	return q;
}

