/** @file 
 * \brief This file contains structures and functions 
 * for \f$c_1\f$ and \f$c_2\f$ parameter optimization using L-BFGS-B
 * library. 
 * 
 * Minimized scoring function is called \f$\chi\f$-score
 * \f$ \chi = \sqrt{\frac{1}{N} \sum_{i = 0}^{N} \frac{(I_{exp}(i) - kI(i))^2}{\sigma^2_i}}\f$.
 * This file also includes scaling function and function for contructing 
 * profile from cross-terms \f$A_xB_x\f$ from given pair of parameters.
 * To do optimization with respect to experimental profile, experiment
 * is fed to scoring_helper(), which transforms it into compressed form,
 * which together with other optimization params goes to sxs_lbfgs_fitting().
 */
#pragma once

#include "common.h"

#include "profile.h"
#include "pdb2spf.h"

#include "lbfgsb.h"

/**
 * \brief This structure stores the parameters for L-BFGS-B library
 * and other params.
 */
struct sxs_opt_params
{
	double rm;            /**< Mean atomic radius. */
	double mult;          /**< \f$ \left(\frac{4\pi}{3}\right)^\frac{3}{2} r_m^2 \f$. */
	double peak;          /**< Peak intencity value of experimental profile. Can be used for initial 
		                       scaling. */

	double* a;            /**< Compressed form of experimental profile. */

	double  wa[43251];    /**< L-BFGS-B library internal workspace. */   
	integer iwa[3072];    /**< L-BFGS-B library internal workspace. */   
	double  dsave[29];    /**< L-BFGS-B library internal workspace. */   
	integer isave[44];    /**< L-BFGS-B library internal workspace. */   
	logical lsave[4];     /**< L-BFGS-B library internal workspace. */   
};

/**
 * Allocates an instance of sxs_opt_params.
 *
 * @param exp Experimental profile.
 * @param qvals Array of scattering angles.
 * @param qnum  Number of scattering angles.
 * @param rm    Mean atomic radius of reference structure.
 * @return An intance of sxs_opt_params. NULL in case of error.
 */
struct sxs_opt_params* sxs_opt_params_create(struct sxs_profile* exp, double* qvals, int qnum, double rm);

/**
 * Allocates internal fields of sxs_opt_params.
 *
 * @param params An instance of sxs_opt_params.
 * @param exp Experimental profile.
 * @param qvals Array of scattering angles.
 * @param qnum  Number of scattering angles.
 * @param rm    Mean atomic radius of reference structure.
 */
void sxs_opt_params_init(struct sxs_opt_params* params, struct sxs_profile* exp, double* qvals, int qnum, double rm);

/** 
 * Frees internal fields of sxs_opt_params.
 *
 * @param params An instance of sxs_opt_params.
 */
void sxs_opt_params_destroy(struct sxs_opt_params* params);

/** 
 * Frees sxs_opt_params.
 * @param params An instance of sxs_opt_params.
 */
void sxs_opt_params_free(struct sxs_opt_params* params);

/**
 * Fits every entry profile in `profiles` using L-BFGS-B algorithm. Skips profile 
 * number i if `mask[i]` != 1. sxs_profile::qvals, sxs_profile::in, sxs_profile::err 
 * must be non-empty and all cross-terms sxs_profile::VV ... must be precomputed
 * (see ).
 *
 * @param profiles Array of sxs_profile with filled cross-terms sxs_profile::VV ...
 * @param params   Optimization parameters.
 * @param mask     Mask for skipping profiles.
 * @param n        Profile number.
 */
void sxs_fit_params(struct sxs_profile** profiles, struct sxs_opt_params* params, int* mask, int n);

/**
 * Fits `profile` using L-BFGS-B algorithm. sxs_profile::qvals, sxs_profile::in, sxs_profile::err 
 * must be non-empty and all cross-terms sxs_profile::VV ... must be precomputed
 * (see massha2.h).
 *
 * @param profile  Profile with filled cross-terms sxs_profile::VV ...
 * @param params   Optimization parameters.
 */
void sxs_lbfgs_fitting(struct sxs_profile* profile, struct sxs_opt_params* params);

/** 
 * Compresses `exp` experimental profile to 6*qnum array of `double`.
 * Is used by sxs_opt_params_init.
 * @param exp Experimental profile.
 * @param qnum Number of elements in qvals.
 * @param qvals Array of scattering angles.
 * @return 6*qnum array of type `double`.
 */
double* scoring_helper(struct sxs_profile* exp, int qnum, double* qvals);

/**
 * Builds full intensity from the cross-terms in `profile`, hence sxs_profile::VV ... must be filled.
 * \f{eqnarray*}{
 * I(q) = A^v{B^v}^\star(q) + G^2(c_1,q) A^d{B^d}^\star(q) + c^2_2 A^w{B^w}^\star(q) \\ 
 *   - G(c_1,q)(A^v{B^d}^\star(q) + A^d{B^v}^\star(q)) \\ 
 *   - c_2 G(c1,q) (A^d{B^w}^\star(q) + A^w{B^d}^\star(q)) \\ 
 *   + c_2 (A^v{B^w}^\star(q) + A^w{B^v}^\star(q)) \vert_{\langle\Omega\rangle}
 * \f}
 * Here `rm` and \f$c_1\f$ (\f$0.94 < c_1 < 1.04\f$) are used to compute \f$G\f$
 * \f[
 * G(c_1, q) = c_1^3 \exp \left[ \frac{-{\left(\frac{4\pi}{3}\right)}^\frac{3}{2} 
 * r_m^2 (c_1^2-1) q^2}{4\pi} \right]
 * \f]
 * \f$c_2\f$ accounts for solvation shell density (\f$-2 < c_2 < 4\f$).
 *
 * @param profile Profile to which the intensity is written.
 * @param rm      Mean atomic radius.
 * @param c1      First parameter (\f$0.94 < c_1 < 1.04\f$).
 * @param c2      Second parameter (\f$-2 < c_2 < 4\f$).
 */
void compile_intensity(struct sxs_profile* profile, double rm, double c1, double c2);

/**
 * Finds the best scaling factor for the `profile` given `c1` and `c2`. Uses cross-terms sxs_profile::VV ...
 * not sxs_profile::in.
 *
 * @param profile Profile with filled cross-terms sxs_profile::VV ...
 * @param params  Parameters for optimization.
 * @param c1      First parameter (\f$0.94 < c_1 < 1.04\f$).
 * @param c2      Second parameter (\f$-2 < c_2 < 4\f$).
 * @return Best scaling parameter for given experimental curve and c1, c2.
 */
double best_scale(struct sxs_profile* profile, struct sxs_opt_params* params, double c1, double c2);

/**
 * Finds \f$chi\f$-score of the `profile` given `c1` and `c2`, without opitimzation. 
 * Uses cross-terms sxs_profile::VV ...
 * not sxs_profile::in.
 *
 * @param profile Profile with filled cross-terms sxs_profile::VV ...
 * @param params  Parameters for optimization.
 * @param c1      First parameter (\f$0.94 < c_1 < 1.04\f$).
 * @param c2      Second parameter (\f$-2 < c_2 < 4\f$).
 * @return \f$\chi\f$-score of the `profile` for given experimental curve and c1, c2.
 */
double point_score(struct sxs_profile* profile, struct sxs_opt_params* params, double c1, double c2);

/** 
 * Fills cross-terms sxs_profile::VV ... given SPF coefficients in `s`.
 *
 * @param profile An instance of sxs_profile.
 * @param s       Structure filled with SPF coefficients (see atom_grp2spf).
 */
void spf2cross_terms(struct sxs_profile* profile, struct sxs_spf_full* s);

/** 
 * Fills cross-terms sxs_profile::VV ... given SPF coefficients in `s` and fits 
 * `profile` to the experiment using `params`.
 *
 * @param profile An instance of sxs_profile.
 * @param s       Structure filled with SPF coefficients (see atom_grp2spf).
 * @param params  Optimization parameters.
 */
void spf2fitted_profile(struct sxs_profile* profile, struct sxs_spf_full* s, struct sxs_opt_params* params);
