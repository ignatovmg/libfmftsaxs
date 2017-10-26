/** @file
 * \brief Constants and flags.
 */
 
/**
 * \addtogroup common
 * @{
 */
#pragma once

#define SQRTPI 1.77245385091

#ifndef M_PI
	#define M_PI 3.14159265358
#endif

#define QMAX 0.5
#define QNUM 50

/**@def REL_ERR
 * Default intencity error.
 */
#define REL_ERR 0.05

#define MINIMIZER_ITERMAX 1000
#define VERBOSE_LBFGS -1

#define C1_LOWER  0.96
#define C1_UPPER  1.04
#define C2_LOWER  0.00
#define C2_UPPER  2.00

#define C1_DEFAULT 1.0
#define C2_DEFAULT 0.5

/** @} */
