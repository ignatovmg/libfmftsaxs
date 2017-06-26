/** @file
 * \brief Function to compute spherical Bessel function.
 */

/** \addtogroup sph_interface
 * @{
 */

#pragma once

#include "common.h"

/** 
 * Compute spherical Bessel function.
 * @param order Order.
 * @param x     X value.
 * @return Function value.
 */
double sxs_sbessel(int order, double x);

/** @} */
