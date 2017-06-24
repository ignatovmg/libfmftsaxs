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
	double step = (end - begin)/ (qnum - 1);
	double* q = (double*)calloc(qnum, sizeof(double));
	q[0] = begin;
	for (int i = 1; i < qnum; i++)
		q[i] = q[i-1] + step;
	return q;
}

void sxs_fill_active_rotation_matrix(struct mol_matrix3 *rm, double alpha, double beta, double gamma)
{
	rm->m11 = cos(gamma)*cos(beta)*cos(alpha) - sin(gamma)*sin(alpha);
	rm->m21 = cos(gamma)*cos(beta)*sin(alpha) + sin(gamma)*cos(alpha);
	rm->m31 = -cos(gamma)*sin(beta);
	rm->m12 = -sin(gamma)*cos(beta)*cos(alpha) - cos(gamma)*sin(alpha);
	rm->m22 = -sin(gamma)*cos(beta)*sin(alpha) + cos(gamma)*cos(alpha);
	rm->m32 = sin(gamma)*sin(beta);
	rm->m13 = sin(beta)*cos(alpha);
	rm->m23 = sin(beta)*sin(alpha);
	rm->m33 = cos(beta);
//	{ cos(g)*cos(b)*cos(a) - sin(g)*sin(a),  cos(g)*cos(b)*sin(a) + sin(g)*cos(a), -cos(g)*sin(b)},
//	{-sin(g)*cos(b)*cos(a) - cos(g)*sin(a), -sin(g)*cos(b)*sin(a) + cos(g)*cos(a),  sin(g)*sin(b)},
//	{ sin(b)*cos(a)                       ,  sin(b)*sin(a)                       ,  cos(b)       }
}

void sxs_mult_rot_mats(struct mol_matrix3 *c, struct mol_matrix3 *a, struct mol_matrix3 *b)
{
	c->m11 = a->m11 * b->m11 + a->m12 * b->m21 + a->m13 * b->m31;
	c->m12 = a->m11 * b->m12 + a->m12 * b->m22 + a->m13 * b->m32;
	c->m13 = a->m11 * b->m13 + a->m12 * b->m23 + a->m13 * b->m33;

	c->m21 = a->m21 * b->m11 + a->m22 * b->m21 + a->m23 * b->m31;
	c->m22 = a->m21 * b->m12 + a->m22 * b->m22 + a->m23 * b->m32;
	c->m23 = a->m21 * b->m13 + a->m22 * b->m23 + a->m23 * b->m33;

	c->m31 = a->m31 * b->m11 + a->m32 * b->m21 + a->m33 * b->m31;
	c->m32 = a->m31 * b->m12 + a->m32 * b->m22 + a->m33 * b->m32;
	c->m33 = a->m31 * b->m13 + a->m32 * b->m23 + a->m33 * b->m33;
}

