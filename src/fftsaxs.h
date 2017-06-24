/** @file
 * \brief This file contains the main routine for FFT-based SAXS calculation.
 *
 * Conformation of a dimer is defined by translation of a ligand along \f$z\f$
 * and 5 Euler angles (z-y-z). The receptor is rotated by \f$(0, \beta_A, \gamma_A)\f$
 * and the ligand is rotated and translated by \f$(\alpha_B, \beta_B, \gamma_B)\f$.
 * Then 
 * \f[
 * I(q) = |A^\prime(\vec{q}) + B^\prime(\vec{q})|_{\langle\Omega\rangle} = 
 * I_A(q) + I_B(q) + F(q)
 * \f]
 * \f$I_A(q)\f$, \f$I_B(q)\f$ are intensities from each of the subunints separately,
 * are computed only once. Here \f$F(q)\f$ is computed via FFT. 
 * \f{eqnarray*}{
 * F(q) = \sum_{m = -L}^{L} \sum_{m_1 = -L}^{L} \sum_{m_2 = -L}^{L} e^{i m\alpha_B} 
 * e^{i m_1(-\gamma_A)} e^{i m_2\gamma_B}
 * \sum_{l = max(|m|, |m_1|)}^{L} \sum_{l_1 = max(|m|, |m_2|)}^{L}
 * d_{m m_1}^{l}(\beta_A) d_{m m_2}^{l_1}(\beta_B)  A_{l m_1}(q) 
 * {[T^{|m|}_{l l_1}(qz) B_{l_1m_2}(q)]}^\star 
 * \f}
 * T is a translation matrix.
 * \f[ 
 * T^{|m|}_{ll_1}(qz) = \sum_{p = |l_1 - l|}^{l_1 + l} (-1)^m j_p(q z) d_{l m}(l_1, p)
 * \f]
 */
#pragma once

#include "common.h"

#ifdef _MPI_
#include <mpi.h>
#endif

#include <complex.h>
#include <fftw3.h>
#include <gmp.h>
#include <time.h>

#include "sfbessel.h"
#include "pdb2spf.h"
#include "min_saxs.h"
#include "profile.h"
#include "index.h"

/**
 * Computes SAXS scores for each grid point present in `index_list`. 
 * How to compute index of a particular grid point:
 *
 *
 * \f[ index = ((((z_n(L+1) + \beta^A_n)(L+1) + \beta^B_n)(2L+1) + \alpha^B_n)(2L+1) + \gamma^A_n)(2L+1) + \gamma^B_n \f]
 *
 * \f$z_n\f$ is the index of a translation step in `zvals`.
 * \f$\beta\f$ angles are sampled uniformly from 0.0 to 3.14: 
 *
 * Sampled values of \f$\beta^A\f$, \f$\beta^B\f$ are
 * \f$ [0.0, \frac{\pi}{L+1}, 2\frac{\pi}{L+1}, ..., \pi] \f$. 
 *
 * Sampled values of \f$\alpha^B\f$, \f$\gamma^B\f$ are
 * \f$[0.0, \frac{2\pi}{2L+1}, 2\frac{2\pi}{2L+1}, ..., 2L\frac{2\pi}{2L+1}]\f$.
 *
 * Sampled values of \f$\gamma^A\f$ are negatives of the latter
 * \f$[0.0, -\frac{2\pi}{2L+1}, -2\frac{2\pi}{2L+1}, ..., -2L\frac{2\pi}{2L+1}]\f$
 *
 * So if you pick \f$z_n = 1\f$, \f$\beta^A_n = 1\f$, \f$\beta^B_n = 2\f$, \f$\alpha^B_n = 20\f$, 
 * \f$\gamma^A_n = 10\f$, \f$\gamma^B_n = 10\f$, given \f$L=15\f$, your 
 * A is rotated by \f$ (0, \frac{180}{16}, -10*\frac{360}{31}) \f$ and B is rotated by
 * \f$ (20*\frac{360}{31}, 2*\frac{180}{16}, 10*\frac{360}{31}) \f$, and its index on the 
 * grid is \f$z_n*16^2*31^3 + 555778 = 801794\f$.
 * Euler rotation axes are z-y-z.
 *
 * All the indices, for which you want score to be computed are written 
 * in `index_list`.
 *
 * @param[in] A SPF coefficients of the receptor (see pdb2spf.h).
 * @param[in] B SPF coefficients of the receptor (see pdb2spf.h).
 * @param[in] qvals Scattering angles.
 * @param[in] qnum  Length of `qvals`.
 * @param[in] params Optimization params (see min_saxs.h).
 * @param[in] zvals Array of translation steps.
 * @param[in] znum Number of translation steps.
 * @param[out] scores_list Array of length `nout`, where \f$\chi\f$-scores are written
 * @param[in] index_list  Array of length `nout` with grid point indices, for which \f$\chi\f$-scores
 * are computed.
 * @param[out] c1_list Array of length `nout`, where optimized parameters are written
 * @param[out] c2_list Array of length `nout`, where optimized parameters are written
 * @param[in] nout Length of `c1_list`, `c2_list`, `index_list`, `scores_list`.
 * @param[in] L Expansion depth.
 * @param[in] skip If `skip=1`, then only those scores are computed,
 * which are provided in the `index_list`.
 */ 
void compute_saxs_scores(
	double* scores_list,
	double*     c1_list,
	double*     c2_list,
	int* index_list, int nout,
	struct sxs_spf_full* A,
	struct sxs_spf_full* B,
	struct sxs_opt_params* params,
	double* qvals, int qnum, 
	double* zvals, int znum, 
	int L, int skip);

// for basic clustering
/*typedef struct compare_class_
{
	public:
	double* arr;
	bool operator()(int i, int j)
	{
		return (arr[i] < arr[j]);
	}
}compare_class;

struct neighbors
{
	int num;
	int *alpha_id;
	int *beta_id;
	int *gamma_id;
	int howmany;
};

struct trig_tab
{
	int num;
	double *cos;
	double *sin;
	double *beta_cos;
	double *beta_sin;
};*/


