#ifndef _DEFINE_H_
#define _DEFINE_H_

#define SQRTPI 	1.77245385091

#ifndef M_PI
	#define M_PI    3.14159265358
#endif

#define QMAX 0.5
#define QNUM 50
#define REL_ERR 0.05

#define MINIMIZER_ITERMAX 1000
#define VERBOSE_LBFGS -1

#define C1_LOWER  0.96
#define C1_UPPER  1.04
#define C2_LOWER -2.00
#define C2_UPPER  4.00

#define C1_DEFAULT 1.0
#define C2_DEFAULT 0.0

#define SKIP_PROFILES 1

#define RECORD_MANIFOLD_RESULTS 0

#define SAVE_CURVES 0
#define DIRECT_CURVES_PATH "testing/sandbox/simple"
#define STORED_CURVES_PATH "testing/sandbox/massha"
#define CURVES_PARAMS_PATH "testing/sandbox/params"

#endif /* _DEFINE_H_ */
