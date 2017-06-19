#include "min_saxs.h"

static void gradient(double* grad, double* score_val, double* point, struct sxs_profile* profile, struct sxs_opt_params* params)
{
	int  qnum = profile->qnum;
	double* q = profile->qvals;
	double c1 = point[0];
	double c2 = point[1];
		   
	double* a = params->a;
	double k = profile->scale;
	
	double grad0 = 0.0;
	double grad1 = 0.0;
	
	double mult = params->mult;   // pow(4.0 * M_PI / 3.0, 3.0 / 2.0) * rm * rm / (16.0 * M_PI);
	double corr = - mult * (c1 * c1 - 1.0);
	double G = c1 * c1 * c1 * exp(corr * q[0] * q[0]);
	double G_der = G*(3.0 / c1 - 2.0 * c1 * mult * q[0] * q[0]);
	
	double* AvBv = profile->AvBv;
	double* AvBd = profile->AvBd;
	double* AvBw = profile->AvBw;
	double* AdBv = profile->AdBv;
	double* AdBd = profile->AdBd;
	double* AdBw = profile->AdBw;
	double* AwBv = profile->AwBv;
	double* AwBd = profile->AwBd;
	double* AwBw = profile->AwBw;
	
	double in_prev = AvBv[0] - 
					 G  * (AvBd[0] + AdBv[0]) +
					 c2 * (AvBw[0] + AwBv[0]) +
					 G * G * AdBd[0] -
					 G  * c2 * (AdBw[0] + AwBd[0]) +
					 c2 * c2 *  AwBw[0];
				
	double in_der_c1_prev = - G_der * (AvBd[0] + AdBv[0]) +
					          2.0 * G * G_der * AdBd[0] -
					          G_der * c2 * (AdBw[0] + AwBd[0]);
					   
	double in_der_c2_prev = AvBw[0] + AwBv[0] -
					        G * (AdBw[0] + AwBd[0]) +
					        2.0 * c2 * AwBw[0];
					   
	double tan, tan_c1_der, tan_c2_der, buf;
	double in , in_der_c1,  in_der_c2;
	double q_prev = -1.0;
	double q_cur;
	double score = 0.0;
	double c1_cube = c1 * c1 * c1;
	
	for (int i = 0; i < qnum; i++) {
		q_cur = q[i];
		G = c1_cube * exp(corr * q_cur * q_cur);
		G_der = G * (3.0 / c1 - 2.0 * c1 * mult * q_cur * q_cur);
	
		in = AvBv[i] -
		     G  * (AvBd[i] + AdBv[i]) +
		     c2 * (AvBw[i] + AwBv[i]) +
			 G * G * AdBd[i] -
			 G  * c2 * (AdBw[i] + AwBd[i]) +
			 c2 * c2 * AwBw[i];
				
		in_der_c1 = G_der * ( -(AvBd[i] + AdBv[i]) +
					            2.0 * G * AdBd[i] -
					            c2 * (AdBw[i] + AwBd[i]) );
					
		in_der_c2 = AvBw[i] + AwBv[i] -
				    G * (AdBw[i] +  AwBd[i]) +
				    2.0 * c2 * AwBw[i];
				
		buf = 1.0 / (q_cur - q_prev);
		tan = (in - in_prev) * buf;
		tan_c1_der = (in_der_c1 - in_der_c1_prev) * buf;
		tan_c2_der = (in_der_c2 - in_der_c2_prev) * buf;
		
		grad0 += 2.0 * k * ( -(in_der_c1 - tan_c1_der * q_cur) * a[i*6+1] -
				               tan_c1_der * a[i*6+2] +
				               k * ( (in - tan * q_cur) * (in_der_c1 - tan_c1_der * q_cur) * a[i*6+3] +
				                     (in * tan_c1_der + in_der_c1 * tan - 2.0 * tan * tan_c1_der * q_cur) * a[i*6+4] +
				                     tan * tan_c1_der * a[i*6+5] ) );
				    
		grad1 += 2.0 * k * ( -(in_der_c2 - tan_c2_der * q_cur) * a[i*6+1] -
				              tan_c2_der * a[i*6+2] +
				              k * ( (in - tan * q_cur) * (in_der_c2 - tan_c2_der * q_cur) * a[i*6+3]  +
				                    (in * tan_c2_der + in_der_c2 * tan - 2.0 * tan * tan_c2_der * q_cur) * a[i*6+4] +
				                    tan * tan_c2_der * a[i*6+5] ) ); 
				    
		buf = in - tan * q_cur;
		score += a[i*6]
		         + k * (- 2.0 * buf * a[i*6+1]
		                - 2.0 * tan * a[i*6+2]
		                + k * (   buf * buf * a[i*6+3]
		                        + 2.0 * buf * tan * a[i*6+4]
		                        + tan * tan * a[i*6+5]) );  
		
		in_prev = in;
		in_der_c1_prev = in_der_c1;
		in_der_c2_prev = in_der_c2;
		q_prev = q_cur;
	}
	
	grad[0] = grad0;
	grad[1] = grad1;

	*score_val = score;
}


struct sxs_opt_params* sxs_opt_params_create(struct sxs_profile* exp, double* qvals, int qnum, double rm)
{
	struct sxs_opt_params* params = calloc(1, sizeof(struct sxs_opt_params));
	sxs_opt_params_init(params, exp, qvals, qnum, rm);
	
	return params;
}

void sxs_opt_params_init(struct sxs_opt_params* params, struct sxs_profile* exp, double* qvals, int qnum, double rm)
{
	if (params != NULL) {
		params->a = scoring_helper(exp, qnum, qvals);
		params->rm = rm;
		params->mult = pow(4.0 * M_PI / 3.0, 1.5) * rm * rm / (16.0 * M_PI);
		params->peak = exp->in[0]; 
	}
}

void sxs_opt_params_destroy(struct sxs_opt_params* params)
{
	if (params != NULL) {
		myfree(params->a);
		params->a = NULL;
	
		double*  wa    = params->wa;
		integer* iwa   = params->iwa;
		 double* dsave = params->dsave;
		integer* isave = params->isave;
		logical* lsave = params->lsave;
		
		memset(wa,    0, 43251*sizeof( double));
		memset(iwa,   0,  3072*sizeof(integer));
		memset(dsave, 0,    29*sizeof( double));
		memset(isave, 0,    44*sizeof(integer));
		memset(lsave, 0,     4*sizeof(logical));
	}
}

void sxs_opt_params_free(struct sxs_opt_params* params)
{
	sxs_opt_params_destroy(params);
	myfree(params);
}


void sxs_fit_params(struct sxs_profile** profiles, struct sxs_opt_params* params, int* mask, int n)
{
	//int iters = 0;
	int count = 0;

	int qnum, q;
	struct sxs_profile* profile;
	
	double scale;
	double peak = params->peak;
	
	for (int i = 0; i < n; i++) {
	
		if (mask[i] == 1) {
			profile = profiles[i];
			qnum = profiles[i]->qnum;
			
			scale = profile->AvBv[0] +
					profile->AdBd[0] +
					profile->AwBw[0] +
					profile->AvBw[0] + 
					profile->AwBv[0] - 
					profile->AvBd[0] -
					profile->AdBv[0] -
					profile->AwBd[0] -
					profile->AdBw[0];
			
			//printf("peak = %e, scale %e\n", peak, scale);
			
			scale = peak / scale;
			
			for (q = 0; q < qnum; q++) {
				profile->AvBv[q] *= scale;
				profile->AdBd[q] *= scale;
				profile->AwBw[q] *= scale;
				profile->AvBw[q] *= scale;
				profile->AwBv[q] *= scale;
				profile->AvBd[q] *= scale;
				profile->AdBv[q] *= scale;
				profile->AwBd[q] *= scale;
				profile->AdBw[q] *= scale;
			}
			
			sxs_lbfgs_fitting(profile, params);
			count ++;
		}
	}
}

void sxs_lbfgs_fitting(struct sxs_profile* profile, struct sxs_opt_params* params)
{
	//clock_t time = clock();
	//double t1 = 0.0, t2 = 0.0;
	
    double score = 0.0;
    double g[2] = {0.0, 0.0};                    // gradient
    double l[2] = {C1_LOWER,   C2_LOWER};        // lower bounds
    double u[2] = {C1_UPPER,   C2_UPPER};        // upper bounds 
    double x[2] = {C1_DEFAULT, C2_DEFAULT};      // current point
    
    static integer m = 3; 
    static integer n = 2;              // problem dimension
    static integer nbd[2] = {2, 2};    // 2 means the variable has both lower and upper bounds
    
    double*  wa    = params->wa;
    integer* iwa   = params->iwa;
     double* dsave = params->dsave;
    integer* isave = params->isave;
    logical* lsave = params->lsave;
    
    memset(wa,    0, 43251*sizeof( double));
    memset(iwa,   0,  3072*sizeof(integer));
    memset(dsave, 0,    29*sizeof( double));
    memset(isave, 0,    44*sizeof(integer));
    memset(lsave, 0,     4*sizeof(logical));

    integer taskValue = 0;
    integer *task=&taskValue;
    integer csaveValue = 0;
    integer *csave=&csaveValue;
    
    double factr = 1e+7;
    double pgtol = 1e-5;
    integer iprint = VERBOSE_LBFGS;
    
    *task = (integer)START;

	L111:
		setulb(&n, &m, x, l, u, nbd, &score, g, &factr, &pgtol, wa, iwa, task, 
		       &iprint, csave, lsave, isave, dsave);

		if ( IS_FG(*task) ) {

		    profile->scale = best_scale(profile, params, x[0], x[1]);
		    gradient(g, &score, x, profile, params);
		    
		    goto L111;
		}

		if ( *task==NEW_X ) {
		    goto L111;
		}
	
	profile->score = sqrt(score);
	profile->scale = best_scale(profile, params, x[0], x[1]);
	profile->c1 = x[0];
	profile->c2 = x[1];
	
	compile_intensity(profile, params->rm, x[0], x[1]);

	//getchar();
	//return params->isave[30];
}

double best_scale(struct sxs_profile* profile, struct sxs_opt_params* params, double c1, double c2)
{
	int  qnum = profile->qnum;
	double* q = profile->qvals;
	
	double* a = params->a;
	
	double mult = params->mult;
	double corr = - mult * (c1 * c1 - 1.0);
	double G = c1 * c1 * c1 * exp(corr * q[0] * q[0]);
	
	double* AvBv = profile->AvBv;
	double* AvBd = profile->AvBd;
	double* AvBw = profile->AvBw;
	double* AdBv = profile->AdBv;
	double* AdBd = profile->AdBd;
	double* AdBw = profile->AdBw;
	double* AwBv = profile->AwBv;
	double* AwBd = profile->AwBd;
	double* AwBw = profile->AwBw;
	
	double in_prev = AvBv[0] -
				     G  * (AvBd[0] + AdBv[0]) +
				     c2 * (AvBw[0] + AwBv[0]) +
				     G  * G  *  AdBd[0] -
				     G  * c2 * (AdBw[0] + AwBd[0]) +
				     c2 * c2 *  AwBw[0];

	double c1_cube = c1 * c1 * c1;
	double buf, tan, in;
	double q_prev = -1.0, q_cur;
	double up = 0.0, down = 0.0;
	
	for (int i = 0; i < qnum; i++) {
		q_cur = q[i];
		G = c1_cube * exp(corr * q_cur * q_cur);
		
		in = AvBv[i] -
		     G *  (AvBd[i] + AdBv[i]) +
		     c2 * (AvBw[i] + AwBv[i]) +
			 G * G * AdBd[i] -
			 G * c2 * (AdBw[i] + AwBd[i]) +
			 c2 * c2 * AwBw[i];
			 
		tan = (in - in_prev) / (q_cur - q_prev);
		buf = in - tan * q_cur;
		
		up += buf * a[i*6+1] + 
		      tan * a[i*6+2];
		      
		down += buf * buf * a[i*6+3] + 
		        2.0 * tan * buf * a[i*6+4] + 
		        tan * tan * a[i*6+5];
	
		in_prev = in;
		q_prev = q_cur;
	}
	
	double k = up / down;
	//printf("k = %e\n", k);
	return k;
}

void compile_intensity(struct sxs_profile* profile, double rm, double c1, double c2)
{
	int  qnum = profile->qnum;
	double* q = profile->qvals;
	
	double* AvBv = profile->AvBv;
	double* AvBd = profile->AvBd;
	double* AvBw = profile->AvBw;
	double* AdBv = profile->AdBv;
	double* AdBd = profile->AdBd;
	double* AdBw = profile->AdBw;
	double* AwBv = profile->AwBv;
	double* AwBd = profile->AwBd;
	double* AwBw = profile->AwBw;
	
	double mult = pow(4.0 * M_PI / 3.0, 3.0 / 2.0) * rm * rm / (16.0 * M_PI);
	double corr = - mult * (c1 * c1 - 1.0);
	double rerr = profile->rerr;
	double G;
	
	for (int i = 0; i < qnum; i++) {
		G = c1 * c1 * c1 * exp(corr * q[i] * q[i]);
		
		profile->in[i] = AvBv[i] -
						 G * (AvBd[i] + AdBv[i]) +
						 c2 * (AvBw[i] + AwBv[i]) +
						 G * G * AdBd[i] -
						 G * c2 * (AdBw[i] + AwBd[i]) +
						 c2 * c2 * AwBw[i];
						 
		profile->err[i] = rerr * profile->in[i];
	}
}


double* scoring_helper(struct sxs_profile* exp, int qnum, double* qvals)
{
	double* a = calloc(6 * qnum, sizeof(double));
	
	memset(a, 0, 6 * qnum * sizeof(double));
	int j = 0;
	
	double lower, upper;
	for (int i = 0; i < qnum; i++) {
		lower = qvals[0];
		upper = qvals[i];
		
		if (i > 0) {
			lower = qvals[i-1];
		}

		while(exp->qvals[j] >= lower && exp->qvals[j] <= upper) {
			a[i*6+0] += exp->in[j] * exp->in[j] / (exp->err[j] * exp->err[j]);
			a[i*6+1] += exp->in[j] / (exp->err[j] * exp->err[j]);
			a[i*6+2] += exp->qvals[j] * exp->in[j] / (exp->err[j] * exp->err[j]);
			a[i*6+3] += 1.0 / (exp->err[j] * exp->err[j]);
			a[i*6+4] += exp->qvals[j] / (exp->err[j] * exp->err[j]);
			a[i*6+5] += exp->qvals[j] * exp->qvals[j] / (exp->err[j] * exp->err[j]);
			j++;
		}
		
		//printf("a[%i]:\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", i, a[i*6+0], a[i*6+1], a[i*6+2], a[i*6+3], a[i*6+4], a[i*6+5]);
	}
	
	double norm = (double)j;
	
	for (int i = 0; i < 6*qnum; i++) {
		a[i] /= norm;
	}
	
	return a;
}

double point_score(struct sxs_profile* profile, struct sxs_opt_params* params, double c1, double c2)
{
	double* AvBv = profile->AvBv;
	double* AvBd = profile->AvBd;
	double* AvBw = profile->AvBw;
	double* AdBv = profile->AdBv;
	double* AdBd = profile->AdBd;
	double* AdBw = profile->AdBw;
	double* AwBv = profile->AwBv;
	double* AwBd = profile->AwBd;
	double* AwBw = profile->AwBw;

	int qnum  = profile->qnum;
	double* q = profile->qvals;
	double* a = params->a;
	double rm = params->rm;
	double corr = - params->mult * (c1 * c1 - 1.0);   
	double G = c1 * c1 * c1 * exp(corr * q[0] * q[0]);
	
	double in_prev = AvBv[0] - 
					 G  * (AvBd[0] + AdBv[0]) +
					 c2 * (AvBw[0] + AwBv[0]) +
					 G * G * AdBd[0] -
					 G  * c2 * (AdBw[0] + AwBd[0]) +
					 c2 * c2 *  AwBw[0];
				   
	double k = best_scale(profile, params, c1, c2);
	
	double buf, tan, in;
	double q_prev = -1.0, q_cur;
	double score = 0.0;
	double c1_cube = c1 * c1 * c1;
	
	for (int i = 0; i < qnum; i++) {
		q_cur = q[i];
		G = c1_cube * exp(corr*q_cur*q_cur);
			 
		in = AvBv[i] -
		     G  * (AvBd[i] + AdBv[i]) +
		     c2 * (AvBw[i] + AwBv[i]) +
			 G * G * AdBd[i] -
			 G  * c2 * (AdBw[i] + AwBd[i]) +
			 c2 * c2 * AwBw[i];
			 	
		tan = (in-in_prev)/(q_cur-q_prev);
		buf = in - tan * q_cur;
		score += a[i*6]
		         + k * (- 2.0 * buf * a[i*6+1]
		                - 2.0 * tan * a[i*6+2]
		                + k * (   buf * buf * a[i*6+3]
		                        + 2.0 * buf * tan * a[i*6+4]
		                        + tan * tan * a[i*6+5]) );  
		
		in_prev = in;
		q_prev = q_cur;
	}
	
	//printf("score = %+6.3e\n\n", score);
	
	return score;
}

void spf2cross_terms(struct sxs_profile* profile, struct sxs_spf_full* s)
{
	int qnum = profile->qnum;
	int L = s->L;

	double* AvBv = profile->AvBv;
	double* AvBd = profile->AvBd;
	double* AvBw = profile->AvBw;
	double* AdBv = profile->AdBv;
	double* AdBd = profile->AdBd;
	double* AdBw = profile->AdBw;
	double* AwBv = profile->AwBv;
	double* AwBd = profile->AwBd;
	double* AwBw = profile->AwBw;
	
	struct sxs_spf_sing** V = s->V;
	struct sxs_spf_sing** D = s->D;
	struct sxs_spf_sing** W = s->W;
	
	for (int q = 0; q < qnum; q++) {
	
		AvBv[q] = 0.0;
		AvBd[q] = 0.0;
		AvBw[q] = 0.0;
		AdBv[q] = 0.0;
		AdBd[q] = 0.0;
		AdBw[q] = 0.0;
		AwBv[q] = 0.0;
		AwBd[q] = 0.0;
		AwBw[q] = 0.0;
		
		int id = 0;
		double val;
		for (int l = 0; l <= L; l++) {
			for (int m = -l; m <= l; m++) {
				val = V[q]->re[id] * V[q]->re[id] + V[q]->im[id] * V[q]->im[id];
				AvBv[q] += val;
				
				val = D[q]->re[id] * D[q]->re[id] + D[q]->im[id] * D[q]->im[id];
				AdBd[q] += val;
				
				val = W[q]->re[id] * W[q]->re[id] + W[q]->im[id] * W[q]->im[id];
				AwBw[q] += val;
				
				val = V[q]->re[id] * D[q]->re[id] + V[q]->im[id] * D[q]->im[id];
				AvBd[q] += val;
				
				val = V[q]->re[id] * W[q]->re[id] + V[q]->im[id] * W[q]->im[id];
				AvBw[q] += val;
				
				val = D[q]->re[id] * W[q]->re[id] + D[q]->im[id] * W[q]->im[id];
				AdBw[q] += val;
				
				id++;
			}
		}
		
		AdBv[q] = AvBd[q];
		AwBv[q] = AvBw[q];
		AwBd[q] = AdBw[q];
	}
}

void spf2fitted_profile(struct sxs_profile* profile, struct sxs_spf_full* s, struct sxs_opt_params* params)
{
	spf2cross_terms(profile, s);
	sxs_lbfgs_fitting(profile, params);
}


