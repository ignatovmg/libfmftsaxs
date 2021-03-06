/** @file 
 * \brief Contains basic functions to operate with SAXS profile.
 */
 
/** \defgroup profile_interface SAXS profile interface
 *  \brief Functions and structures working with SAXS profiles.
 *  @{
 */
#pragma once

#include "common.h"

#include "pdb2spf.h"

/** 
 * \brief This structure stores SAXS curve. In the case of a dimer
 * it can contain Euler coordinates of subunits, together with 
 * fitting parameters to the experimental curve. If dimer profile is
 * computed using FFT, 6 partial cross-terms are filled first and 
 * then intencity is compiled, after fitting the parameters 
 * \f$c_1\f$ and \f$c_2\f$ (see massha2.h for details).
 */
struct sxs_profile 
{
	int qnum;         /**< Number of scattering angles \f$q\f$, at which intencity is sampled.*/
	double rerr;      /**< Relative error (defined by #REL_ERR).*/
	
	double* in;       /**< Scattering intensity. */
	double* err;      /**< Intensity error. */
	double* qvals;    /**< Scattering angle values. */
	
	double* VV;       /**< \f$ A_v(q)B_v(q)\vert_{\langle\Omega\rangle}\f$ */
	double* VD;       /**< \f$ A_v(q)B_d(q) + A_d(q)B_v(q)\vert_{\langle\Omega\rangle}\f$ */
	double* VW;       /**< \f$ A_v(q)B_w(q) + A_w(q)B_v(q)\vert_{\langle\Omega\rangle}\f$ */
	double* DD;       /**< \f$ A_d(q)B_d(q)\vert_{\langle\Omega\rangle}\f$ */
	double* DW;       /**< \f$ A_d(q)B_w(q) + A_w(q)B_d(q)\vert_{\langle\Omega\rangle}\f$ */
	double* WW;       /**< \f$ A_w(q)B_w(q)\vert_{\langle\Omega\rangle}\f$ */
	
	double score;     /**< Profile \f$\chi\f$-score. */
	double scale;     /**< Profile scale.*/
	double c1;        /**< Parameter to adjust mean atomic radius. \f$0.96 < c_1 < 1.04\f$ */
	double c2;        /**< Parameter to adjust solvation shell density. \f$-2 < c_1 < 4\f$ */
	
	double b1;
	double g1;
	double a2;
	double b2;
	double g2;
};

/**
 * Allocates an instance of sxs_profile.
 * @param qvals Array of scattering angles.
 * @param qnum  Number of scattering angles.
 * @param cross_terms_flag If equals 1, then the 
 * cross-terms of sxs_profile::VV are allocated as well.
 * @return An instance of sxs_profile. Returns NULL if case of error.
 */
struct sxs_profile* sxs_profile_create(double* qvals, int qnum, int cross_terms_flag);

/**
 * Allocates internal fields of sxs_profile.
 * @param profile An instance of sxs_profile.
 * @param qvals Array of scattering angles.
 * @param qnum  Number of scattering angles.
 * @param cross_terms_flag If equals 1, then the 
 * cross-terms of sxs_profile::VV are allocated as well.
 */
void sxs_profile_init(struct sxs_profile* profile, double* qvals, int qnum, int cross_terms_flag);

/**
 * Allocates cross-terms of sxs_profile.
 * @param profile A non-empty instance of sxs_profile.
 * @param qnum    Number of scattering angles.
 */
void sxs_profile_alloc_cross_terms(struct sxs_profile* profile, int qnum);

/** 
 * Frees profile internal fields memory, except for sxs_profile::qvals.
 * sxs_profile::qvals is set to NULL, so make sure you saved the pointer
 * before hand.
 * @param profile An instance of sxs_profile.
 */
void sxs_profile_destroy(struct sxs_profile* profile);

/** 
 * Frees profile memory. Make sure you saved the pointer to sxs_profile::qvals
 * because it not freed, but simply set to zero.
 * @param profile An instance of sxs_profile.
 */
void sxs_profile_free(struct sxs_profile* profile);

/**
 * Print profile to file.
 * @param path Destination path.
 * @param profile An instance of sxs_profile.
 */
void sxs_profile_write(char* path, struct sxs_profile* profile);

/**
 * Read profile from path.
 * @param path Source path.
 */
struct sxs_profile* sxs_profile_read(char* path);

/**
 * Read profile from file.
 * @param f Source stream.
 */
struct sxs_profile* sxs_profile_fread(FILE* f);

/**
 * Fill profile intensity using coeffients in `spf`.
 * @param profile A non-empty instance of sxs_profile.
 * @param spf     A non-empty instance of sxs_spf_full, containing 
 * SPF expansion coefficients computed through atom_grp2spf().
 * @param c1 See sxs_profile::c1 (Default 1.0).
 * @param c2 See sxs_profile::c2 (Default 0.0).
 */
void sxs_profile_from_spf(struct sxs_profile* profile, struct sxs_spf_full* spf, double c1, double c2);

/** @}*/

