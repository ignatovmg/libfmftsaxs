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
#include "mpapi2.h"
#include "index.h"

void compute_saxs_scores(
	struct sxs_spf_full* A,
	struct sxs_spf_full* B,
	double* qvals,
	int qnum, 
	struct sxs_opt_params* params,
	double z_beg, 
	double z_step,
	int znum, 
	double* scores_out,
	int* index_out,
	double* c1_out,
	double* c2_out,
	int nout,
	int L);

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


