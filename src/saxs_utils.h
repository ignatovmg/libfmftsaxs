/** @file 
 * \brief Some useful functions. Most of them are copied from
 * libsaxs.
 */
 
/**
 * \addtogroup common
 * @{
 */
#pragma once

#include "common.h"

#include "mol2/atom_group.h"
#include "mol2/sasa.h"
#include "mol2/transform.h"

double mol_atom_group_max_dist(const struct mol_atom_group *ag);

double mol_atom_group_average_radius(const struct mol_atom_group *ag);

/** 
 * Create array uniformly sampled array
 * from `begin` to `end`.
 * @param begin Start array.
 * @param end   End array.
 * @param qnum  Array length.
 * @return Array.
 */
double* sxs_mkarray(double begin, double end, int qnum);

/**
 * Create rotation matrix from Euler angles \f$(\alpha, \beta, \gamma)\f$.
 * @param rm Rotation matrix.
 * @param alpha Alpha angle.
 * @param beta Beta angle.
 * @param gamma Gamma angle.
 */
void sxs_fill_active_rotation_matrix(struct mol_matrix3 *rm, double alpha, double beta, double gamma);

/**
 * Multiply rotation matrices c = a*b.
 */
void sxs_mult_rot_mats(struct mol_matrix3 *c, struct mol_matrix3 *a, struct mol_matrix3 *b);

/**
 * Fill solvent accessibility for each atom.
 */
void sxs_faccs(
	double *fractional_sa,
	const struct mol_atom_group *ag,
	double r_solv);
	
static void upcase(char *s)
{
	while(*s) {
		*s = toupper(*s);
		s++;
	}
}

/*static inline double sinc(double x)
{
	return x == 0.0 ? 1 : sin(x)/x;
}*/

/*static inline double square(double x)
{
	return x*x;
}*/

/** @} */
