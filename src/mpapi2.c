#include "mpapi2.h"

struct sxs_profile* sxs_profile_create(double* qvals, int qnum, int cross_terms_flag)
{
	struct sxs_profile* profile = calloc(1, sizeof(struct sxs_profile));
	sxs_profile_init(profile, qvals, qnum, cross_terms_flag);
	
	return profile;
}

void sxs_profile_init(struct sxs_profile* profile, double* qvals, int qnum, int cross_terms_flag)
{
	if (profile != NULL) {
		profile->qnum  = qnum;
		profile->rerr  = REL_ERR;
		profile->qvals = qvals;
	
		profile->in  = calloc(qnum, sizeof(double));
		profile->err = calloc(qnum, sizeof(double));
	
		profile->score = INFINITY;
		profile->c1 = C1_DEFAULT;
		profile->c2 = C2_DEFAULT;
		profile->scale = 1.0;
	
		if (cross_terms_flag == 1) {
			sxs_profile_alloc_cross_terms(profile, qnum);
		} else {
			profile->VV = NULL;
			profile->VD = NULL;
			profile->VW = NULL;
			profile->DD = NULL;
			profile->DW = NULL;
			profile->WW = NULL;
		}
	}
}

void sxs_profile_alloc_cross_terms(struct sxs_profile* profile, int qnum)
{
	profile->VV = calloc(qnum, sizeof(double));
	profile->VD = calloc(qnum, sizeof(double));
	profile->VW = calloc(qnum, sizeof(double));
	profile->DD = calloc(qnum, sizeof(double));
	profile->DW = calloc(qnum, sizeof(double));
	profile->WW = calloc(qnum, sizeof(double));
}

void sxs_profile_destroy(struct sxs_profile* profile)
{
	profile->qvals = NULL;
	myfree(profile->in);
	myfree(profile->err);

	myfree(profile->VV);
	myfree(profile->VD);
	myfree(profile->VW);
	myfree(profile->DD);
	myfree(profile->DW);
	myfree(profile->WW);
}

void sxs_profile_free(struct sxs_profile* profile)
{
	sxs_profile_destroy(profile);
	myfree(profile);
}

void sxs_profile_write(char* path, struct sxs_profile* profile)
{
	FILE* f = myfopen(path, "w");
	
	double err;
	for (int i = 0; i < profile->qnum; i++) {
		if (profile->err[i] > 0.0) {
			err = profile->err[i];
		} else {	
			err = profile->in[i] * profile->rerr;
		}
		
		fprintf(f, "%.4f %.4f %.4f\n", profile->qvals[i], profile->in[i], err);
	}
	
	fclose(f);
}

struct sxs_profile* sxs_profile_read(char* path)
{
	FILE* f = fopen(path, "r");
	
	struct sxs_profile* profile = sxs_profile_fread(f);
	
	fclose(f);

	return profile;
}

struct sxs_profile* sxs_profile_fread(FILE* f)
{
	if (f != NULL) {
		int cnum;
		int qnum = 0;
		double buf1, buf2, buf3;
		while ((cnum = fscanf(f, "%lf %lf %lf", &buf1, &buf2, &buf3)) != EOF) {
			if (cnum != 3) {
				ERROR_MSG("Wrong input file format.");
			}
		
			qnum++;
		}
		rewind(f);
	
		double* qvals = calloc(qnum, sizeof(double));
		struct sxs_profile* profile = sxs_profile_create(qvals, qnum, 0);
	
		for (int q = 0; q < qnum; q++) {
			fscanf(f, "%lf %lf %lf\n", &(qvals[q]), &(profile->in[q]), &(profile->err[q]));
		}
		
		return profile;
	}
	
	return NULL;
}

void scale_profiles(struct sxs_profile** profiles, int n, double scale, int* mask) 
{
	struct sxs_profile* profile;
	int qnum, q; 
	
	for (int i = 0; i < n; i++) {
	
		if (mask == NULL || mask[i] == 1) {
			profile = profiles[i];
			qnum = profile->qnum;
			
			for (q = 0; q < qnum; q++) {
				profile->in[q]   *= scale;
				profile->err[q]  *= scale;
				profile->VV[q] *= scale;
				profile->DD[q] *= scale;
				profile->WW[q] *= scale;
				profile->VW[q] *= scale;
				profile->VD[q] *= scale;
				profile->DW[q] *= scale;
			}
		}
	}
}

void scale_profiles_by_peak(struct sxs_profile** profiles, int n, double peak, int* mask) 
{
	struct sxs_profile* profile;
	double scale;
	int qnum, q;
	
	for (int i = 0; i < n; i++) {
	
		if (mask == NULL || mask[i] == 1) {
			profile = profiles[i];
			qnum = profile->qnum;
			
			scale = profile->VV[0] +
					profile->DD[0] +
					profile->WW[0] +
					profile->VW[0] -
					profile->VD[0] -
					profile->DW[0];
			
			scale = peak / scale;
			
			for (q = 0; q < qnum; q++) {
				profile->in[q]   *= scale;
				profile->err[q]  *= scale;
				profile->VV[q] *= scale;
				profile->DD[q] *= scale;
				profile->WW[q] *= scale;
				profile->VW[q] *= scale;
				profile->VD[q] *= scale;
				profile->DW[q] *= scale;
			}
		}
	}
}

void sxs_profile_from_spf(struct sxs_profile* profile, struct sxs_spf_full* s, double c1, double c2)
{
	CHECK_PTR(profile);
	CHECK_PTR(s);
	
	assert(profile->qnum == s->qnum);
	
	int L = s->L;
	int qnum = profile->qnum;
	double* qvals = profile->qvals;
	
	double rm   = s->rm;
	double rerr = profile->rerr;
	double mult = pow(4.0 * M_PI / 3.0, 3.0 / 2.0) * rm * rm / (16.0 * M_PI);
	double corr = - mult * (c1 * c1 - 1.0);
	
	int idx;
	double G, in, a_re, a_im;
	double *v_re, *v_im, *d_re, *d_im, *w_re, *w_im;
	for (int q = 0; q < qnum; q++) {
		
		v_re = s->V[q]->re;
		v_im = s->V[q]->im;
		d_re = s->D[q]->re;
		d_im = s->D[q]->im;
		w_re = s->W[q]->re;
		w_im = s->W[q]->im;
		
		G = c1 * c1 * c1 * exp(corr * qvals[0] * qvals[0]);
		
		in = 0.0;
		for (int l = 0; l <= L; l++) {
			for (int m = -l; m <= l; m++) {
				
				idx = lm_index(l, m);
				a_re = v_re[idx] - G * d_re[idx] + c2 * w_re[idx];
				a_im = v_im[idx] - G * d_im[idx] + c2 * w_im[idx];
				in += a_re * a_re + a_im * a_im;
			}
		}
		
		profile->in[q]  = in;
		profile->err[q] = in * rerr;
	}	
	
	profile->c1 = c1;
	profile->c2 = c2;
}


