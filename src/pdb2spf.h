/** @file
 * \brief Contains structures and functions for computing 
 * SPF coefficients from atom groups.
 *
 * Atom group SAXS intencity can be computed as \f$ 
 * I(q) = \sum_{l=0}^L \sum_{m=-l}^l {|A_{lm}^{v}(q) - 
 * G(c_1,q)A_{lm}^{d}(q) + c_2 A_{lm}^{w}(q)|}^2
 * \f$. This file provides tools to compute and manipulate 
 * \f$A_{lm}\f$ coefficients.
 * 
 */
#pragma once

#include "common.h"

#include "borrowed.h"
#include "saxs_utils.h"
#include "form_factor_table.h"
#include "sfbessel.h"

/**
 * \brief Stores SPF coefficients \f$A_{lm}\f$ for a single q value.
 * Contains \f$(L+1)^2\f$ coefficients in total of maximum order \f$L\f$
 */
struct sxs_spf_sing 
{
	int L;		  /**< Maximum depth */
	double *re;	  /**< Real values */
	double *im;	  /**< Imaginary values */
};

/**
 * \brief Stores SPF coefficients \f$A_{lm}^{v}(q), A_{lm}^{d}(q), A_{lm}^{w}(q)\f$ 
 * for the entire scattering curve.
 */
struct sxs_spf_full 
{
	int L;            /**< Maximum depth */
	int qnum;         /**< Number of q values */
	double rm;        /**< Molecule mean radius */
	
	struct sxs_spf_sing** V;   /**< Vacuum part: \f$A_{lm}^{v}(q) = 4\pi i^{l} 
	                                \sum_{i} f_{i}^{v}(q) j_l(qr_i) Y_{lm}(\omega)\f$ */
	struct sxs_spf_sing** D;   /**< Dummy part:  \f$A_{lm}^{d}(q) = 4\pi i^{l} 
	                                \sum_{i} f_{i}^{d}(q) j_l(qr_i) Y_{lm}(\omega)\f$ */
	struct sxs_spf_sing** W;   /**< Water part:  \f$A_{lm}^{w}(q) = 4\pi i^{l} 
	                                \sum_{i} SASA_i f^{w}(q) j_l(qr_i) Y_{lm}(\omega)\f$ */
};

/**
 * Compute SPF coefficients for atom group.
 * @param[in] ag       Input molecule.
 * @param[in] ff_table Form-factor table.
 * @param[in] qvals    Array with q values.
 * @param     qnum     `qvals` length.
 * @param     L        Expansion depth.
 * @param     water    If equal to 1, sxs_spf_full::W is computed, 
 * otherwise is set to zeros.
 * @return SPF coefficients for all q.
 */
struct sxs_spf_full* atom_grp2spf(
	struct mol_atom_group* ag,
	struct saxs_form_factor_table* ff_table, 
	double* qvals,
	int qnum,
	int L,
	int water);
	
/**
 * Compute SPF coefficients for atom group and fill in `spf_coefs`.
 * @param[out] spf_coefs  Preallocated structure.
 * @param[in]  ag         Input molecule.
 * @param[in]  ff_table   Form-factor table.
 * @param[in]  qvals      Array with q values.
 * @param      qnum       `qvals` length.
 * @param      L          Expansion depth.
 * @param[in]  saxs_sa    Solvent accessibility of each atom.
 */
void atom_grp2spf_inplace(
	struct sxs_spf_full* spf_coefs,
	struct mol_atom_group* ag,
	struct saxs_form_factor_table* ff_table, 
	double* qvals,
	int qnum,
	int L,
	double* saxs_sa);

/**
 * Allocates an instance of sxs_spf_sing.
 * @param L Expansion depth.
 * @return An instance of sxs_spf_sing.
 */
struct sxs_spf_sing* sxs_spf_sing_create(int L);

/**
 * Allocates internal fields of sxs_spf_sing.
 * @param L Expansion depth.
 * @param spf An instance of sxs_spf_sing.
 */
void sxs_spf_sing_init(struct sxs_spf_sing* spf, int L);

/**
 * Frees internal fields of sxs_spf_sing.
 * @param s Non-empty instance of sxs_spf_sing.
 */
void sxs_spf_sing_destroy(struct sxs_spf_sing* s);

/**
 * Frees memory of a sxs_spf_sing.
 * @param s Non-empty instance of sxs_spf_sing.
 */
void sxs_spf_sing_free(struct sxs_spf_sing* s);

/**
 * Allocates an instance of sxs_spf_full.
 * @param L Expansion depth.
 * @param qnum Number of sampled scattering angles \f$q\f$.
 * @return An instance of sxs_spf_full.
 */
struct sxs_spf_full* sxs_spf_full_create(int L, int qnum);

/**
 * Allocates internal fields of sxs_spf_full.
 * @param L Expansion depth.
 * @param qnum Number of sampled scattering angles \f$q\f$.
 * @param spf An instance of sxs_spf_full.
 */
void sxs_spf_full_init(struct sxs_spf_full* spf, int L, int qnum);

/**
 * Frees memory of a sxs_spf_full.
 * @param s Non-empty instance of sxs_spf_full.
 */
void sxs_spf_full_free(struct sxs_spf_full* s);

/**
 * Frees internal memory of sxs_spf_full.
 * @param s Non-empty instance of sxs_spf_full.
 */
void sxs_spf_full_destroy(struct sxs_spf_full* s);

/**
 * Resets sxs_spf_full::V, sxs_spf_full::D, sxs_spf_full::W.
 * @param s An instance of sxs_spf_full.
 */
void sxs_spf_full_reset(struct sxs_spf_full* s);

/**
 * Print coefficients of s into path.
 * @param path Path to destination file.
 * @param s    An instance of sxs_spf_full.
 */
void sxs_spf_full_write(char* path, struct sxs_spf_full* s);

/**
 * Read SPF coefficients from path.
 * @param path Path to source file.
 * @return An instance of sxs_spf_full.
 */
struct sxs_spf_full* sxs_spf_full_read(char* path);

/**
 * Read SPF coefficients from stream.
 * @param f Source file.
 * @return An instance of sxs_spf_full.
 */
struct sxs_spf_full* sxs_spf_full_fread(FILE* f);
