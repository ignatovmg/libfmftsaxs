/** @file
 * \brief These functions are borrowed from libfmft. 
 * They compute things related to the spherical space.
 */
 
/** \addtogroup sph_interface
 * @{
 */
#pragma once

#include "common.h"

#include <gmp.h>

// libfmft/polynomials.h
/**
 * Standard LM indexer.
 * \param N Highest polynomial order.
 * \param l In range [0, N-1]
 * \param m In range [-l, l]
 */
static inline int index_LM(int N, int l, int m)
{
	return N*(m + N-1) + l;
}

/**
	Generates an N*N normalization
	array for a set of spherical harmonics
	up to order N. (Y(l,m) = norm_l,m * P_lm(cos(theta)) * exp(i*m*phi))
	\param N Highest polynomial order of spherical harmonics.
	\return Pointer to the generated array. Individual elements of this array can be accessed with index_LM() indexer.
 */
double* generate_spherical_norm(int N);


/**
	Fills a pre-allocated array P
	with associated Legendre polynomials' values
	calculated for specific parameter value x.
	\param N Highest polynomial order of spherical harmonics.
	\param x Associated Legendre polinomials' parameter.
	\param P Pointer to a pre-allocated array of size (N*(2*N - 1)). Individual elements of this array can be accessed with index_LM() indexer.
 */
void fill_assoc_Legendre_array(int N, double x, double* P);




// libfmft/translational_matrix.h
/**
 * Allocates an array of arbitrary precision floating point elements
 * of size n+1
 * and fills it with values of descending (i.e. normal) factorial
 * of each i in range [0, n].
 * \param n Highest number, for which the factorial will be calculated.
 * \return Pointer to the filled array.
 */
mpf_t* generate_desc_fact_array_arb(const int n);

/**
 * Calculates wigner 3-j symbol
 * \param desc Pointer to an array of descending factorials generated by generate_desc_fact_array_arb()
 */
void wigner_3j_symbol_arb(mpf_t res, const int j1, const int j2, const int j3, const int m1, const int m2, const int m3, const mpf_t* desc);





// libfmft/rotational_matrix.h
struct d_array
{
    int L; /**< Highest polynomial order, for which matrix elements are stored. */
    double *data; /**< array of size (L+1)*(2*L + 1)*(2*L + 1). */
};

static inline int index_d_array(int L, int l, int m, int m1)
{
	return (2*L + 1)*((2*L + 1)*l + (L + m)) + (L + m1);
}

/**
 * Allocates an instance of d_array struct
 * that is intended to store Wigner d matrix elements up to polynomial order L.
 * \param L Highest polynomial order, for which matrix elements will be stored.
 * \return Pointer to the allocated struct.
 */
struct d_array* allocate_d_array(const int L);

/**
 * Deallocates an instance of d_array struct.
 * \param d Pointer to an instance of d_array struct.
 */
void deallocate_d_array(struct d_array* d);

/**
 * Allocates an instance of d_array struct and
 * fills it with wigner d matrix elements calculated for a
 * single value of Euler beta angle up to polynomial order L.
 * \return Pointer to the filled struct.
 */
struct d_array* generate_d_array(const int L, const double beta);

/** @} */


