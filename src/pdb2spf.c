#include "pdb2spf.h"

typedef struct sph_coords{
	double r;
	double theta;
	double fi;
}s_coords;

s_coords cart2sph(double x, double y, double z)
{
	s_coords coords;
	coords.r = sqrt(x * x + y * y + z * z);
	coords.theta = acos(z / coords.r);
	
	if (y > 0.0) {
		coords.fi =  acos(x / sqrt(x * x + y * y));
	} else {
		coords.fi = -acos(x / sqrt(x * x + y * y)) + 2.0 * M_PI;
	}

	return coords;
}

void atom_grp2spf_inplace(
	struct sxs_spf_full* spf_coefs,
	struct mol_atom_group* ag,
	struct saxs_form_factor_table* ff_table, 
	double* qvals,
	int qnum,
	int L,
	double* saxs_sa)
{
	sxs_spf_full_reset(spf_coefs);
	spf_coefs->rm = mol_atom_group_average_radius(ag);
	
	int natoms = ag->natoms;
	int spf_size = (L+1)*(L+1);
	int leg_size = (L+1)*(2*L+1);
	
	double h2o_ff;
	if (saxs_sa != NULL) {
		h2o_ff = ff_table->factors[s_OH2].zero_ff;
	}

	double* y_norm = generate_spherical_norm(L + 1);
	double* sph_harm_re = calloc(spf_size, sizeof(double));
	double* sph_harm_im = calloc(spf_size, sizeof(double));
	double* assoc_legendre = calloc(leg_size, sizeof(double));
	
	int idx, index_L_lm;
	double hrm_val, bes_val, val_re, val_im;
	double vacuum_ff, dummy_ff, h2o_ff_cur;
	double *v_re, *v_im, *d_re, *d_im, *w_re, *w_im;
	
	double four_pi = 4.0 * M_PI;
	double four_pi_i_pow_l_re[] = {four_pi, 0.0, -four_pi, 0.0};
	double four_pi_i_pow_l_im[] = {0.0, four_pi, 0.0, -four_pi};
	
	s_coords sph;
	struct mol_vector3 *atom;
	const struct saxs_form_factor *ff;
	
	for(int aid = 0; aid < natoms; aid++) {
		atom = &(ag->coords[aid]);
		sph  = cart2sph(atom->X, atom->Y, atom->Z);
		ff   = get_ff(ff_table, ag, aid);

		memset(assoc_legendre,  0, leg_size * sizeof(double));
		fill_assoc_Legendre_array(L + 1, cos(sph.theta), assoc_legendre);

		h2o_ff_cur = 0.0;
		if (saxs_sa != NULL) {
			h2o_ff_cur = h2o_ff * saxs_sa[aid];
		}
		
		for (int l = 0; l <= L; l++) {
			for (int m = -l; m <= l; m++) {
				idx = lm_index(l, m);
				
				hrm_val = y_norm[index_LM(L + 1, l,  m)] * assoc_legendre[index_LM(L + 1, l, abs(m))];
				sph_harm_re[idx] = hrm_val * cos(  sph.fi * m);
				sph_harm_im[idx] = hrm_val * sin(- sph.fi * m); // conjugate
			}
		}
		
		vacuum_ff = ff->vacuum_ff;
		dummy_ff  = ff->dummy_ff;
		
		for(int q = 0; q < qnum; q++) {
			
			v_re = spf_coefs->V[q]->re;
			v_im = spf_coefs->V[q]->im;
			d_re = spf_coefs->D[q]->re;
			d_im = spf_coefs->D[q]->im;
			w_re = spf_coefs->W[q]->re;
			w_im = spf_coefs->W[q]->im;

			for (int l = 0; l <= L; l++) {
				bes_val = sf_bessel(l, qvals[q] * sph.r);
		
				for (int m = -l; m <= l; m++) {
					idx = lm_index(l, m);
					
					val_re = bes_val * sph_harm_re[idx];
					val_im = bes_val * sph_harm_im[idx];
					
					v_re[idx] += vacuum_ff * val_re;
					v_im[idx] += vacuum_ff * val_im;
					d_re[idx] += dummy_ff * val_re;
					d_im[idx] += dummy_ff * val_im;
					w_re[idx] += h2o_ff_cur * val_re;
					w_im[idx] += h2o_ff_cur * val_im;
				}
			}
		}
	}
	
	double coef_re, coef_im;
	for(int q = 0; q < qnum; q++) {
		v_re = spf_coefs->V[q]->re;
		v_im = spf_coefs->V[q]->im;
		d_re = spf_coefs->D[q]->re;
		d_im = spf_coefs->D[q]->im;
		w_re = spf_coefs->W[q]->re;
		w_im = spf_coefs->W[q]->im;
	
		for (int l = 0; l <= L; l++) {
			coef_re = four_pi_i_pow_l_re[l % 4];
			coef_im = four_pi_i_pow_l_im[l % 4];
			
			for (int m = -l; m <= l; m++) {
				idx = lm_index(l, m);
				
				val_re = v_re[idx]; val_im = v_im[idx];
				v_re[idx] = coef_re * val_re - coef_im * val_im;
				v_im[idx] = coef_re * val_im + coef_im * val_re;
				
				val_re = d_re[idx]; val_im = d_im[idx];
				d_re[idx] = coef_re * val_re - coef_im * val_im;
				d_im[idx] = coef_re * val_im + coef_im * val_re;
				
				val_re = w_re[idx]; val_im = w_im[idx];
				w_re[idx] = coef_re * val_re - coef_im * val_im;
				w_im[idx] = coef_re * val_im + coef_im * val_re;
			}
		}
	}
	free(y_norm);
	free(assoc_legendre);
	free(sph_harm_re);
	free(sph_harm_im);
}

struct sxs_spf_full* atom_grp2spf(
	struct mol_atom_group* ag,
	struct saxs_form_factor_table* ff_table, 
	double* qvals,
	int qnum,
	int L,
	int water)
{
	struct sxs_spf_full* spf_coefs = sxs_spf_full_create(L, qnum);
	
	double* saxs_sa = NULL;
	if (water == 1) {
		saxs_sa = calloc(ag->natoms, sizeof(double));
		faccs(saxs_sa, ag, 1.4);
	}
	
	atom_grp2spf_inplace(spf_coefs, ag, ff_table, qvals, qnum, L, saxs_sa);
	
	free(saxs_sa);
	
	return spf_coefs;
}


struct sxs_spf_sing* sxs_spf_sing_create(int L)
{
	struct sxs_spf_sing* s = calloc(1, sizeof(struct sxs_spf_sing));
	sxs_spf_sing_init(s, L);
	
	return s;
}

void sxs_spf_sing_init(struct sxs_spf_sing* s, int L)
{
	if (s != NULL) {
		s->L = L;
		s->re = calloc((L+1)*(L+1), sizeof(double));
		s->im = calloc((L+1)*(L+1), sizeof(double));
	}
}

void sxs_spf_sing_free(struct sxs_spf_sing* s)
{
	sxs_spf_sing_destroy(s);
	myfree(s);
}

void sxs_spf_sing_destroy(struct sxs_spf_sing* s)
{
	if (s != NULL) {
		myfree(s->re);
		myfree(s->im);
	}
}

struct sxs_spf_full* sxs_spf_full_create(int L, int qnum) 
{
	struct sxs_spf_full* s = malloc(sizeof(struct sxs_spf_full));
	sxs_spf_full_init(s, L, qnum);
	
	return s;
}

void sxs_spf_full_init(struct sxs_spf_full* s, int L, int qnum) 
{
	if (s != NULL) {
		s->L = L;
		s->qnum = qnum;
	
		s->V = calloc(qnum, sizeof(struct sxs_spf_sing*));
		s->D = calloc(qnum, sizeof(struct sxs_spf_sing*));
		s->W = calloc(qnum, sizeof(struct sxs_spf_sing*));
	
		for (int i = 0; i < qnum; i++) {
			s->V[i] = sxs_spf_sing_create(L);
			s->D[i] = sxs_spf_sing_create(L);
			s->W[i] = sxs_spf_sing_create(L);
		}
	}
}

void sxs_spf_full_free(struct sxs_spf_full* s) 
{
	sxs_spf_full_destroy(s);
	myfree(s);
}

void sxs_spf_full_destroy(struct sxs_spf_full* s) 
{
	if (s != NULL) {
		for (int i = 0; i < s->qnum; i++) {
			sxs_spf_sing_free(s->V[i]);
			sxs_spf_sing_free(s->D[i]);
			sxs_spf_sing_free(s->W[i]);
		}
	}
}

void sxs_spf_full_reset(struct sxs_spf_full* s) 
{
	int L = s->L;
	int size = (L+1)*(L+1);
	
	for (int i = 0; i < s->qnum; i++) {
		memset(s->V[i]->re, 0, size * sizeof(double));
		memset(s->V[i]->im, 0, size * sizeof(double));
		memset(s->D[i]->re, 0, size * sizeof(double));
		memset(s->D[i]->im, 0, size * sizeof(double));
		memset(s->W[i]->re, 0, size * sizeof(double));
		memset(s->W[i]->im, 0, size * sizeof(double));
	}
}

void sxs_spf_full_write(char* path, struct sxs_spf_full* s)
{
	FILE* f = fopen(path, "w");
	fprintf(f, "% i % i % .4f\n", s->L, s->qnum, s->rm);
	
	int idx;
	struct sxs_spf_sing *V, *D, *W;
	for (int q = 0; q < s->qnum; q++) {
	
		V = s->V[q];
		D = s->D[q];
		W = s->W[q];
		
		for (int l = 0; l <= s->L; l++) {
			for (int m = -l; m <= l; m++) {
			
				idx = lm_index(l, m);
				fprintf(f, "% .4e % .4e % .4e % .4e % .4e % .4e\n", V->re[idx], 
				                                              V->im[idx],
				                                              D->re[idx],
				                                              D->im[idx],
				                                              W->re[idx],
				                                              W->im[idx]);
			}
		}
	}
	
	fclose(f);
}

struct sxs_spf_full* sxs_spf_full_fread(FILE* f)
{
	if (f != NULL) {
		int L, qnum, nread;
		double rm;

		nread = fscanf(f, "%i %i %lf\n", &L, &qnum, &rm);
	
		if (nread != 3) {
			ERROR_MSG("Wrong input file format.");
		}
	
		struct sxs_spf_full* s = sxs_spf_full_create(L, qnum);
		s->rm = rm;
	
		int idx;
		struct sxs_spf_sing *V, *D, *W;
		for (int q = 0; q < qnum; q++) {
	
			V = s->V[q];
			D = s->D[q];
			W = s->W[q];
			for (int l = 0; l <= L; l++) {
				for (int m = -l; m <= l; m++) {
			
					idx = lm_index(l, m);
					nread = fscanf(f, "%lf %lf %lf %lf %lf %lf\n", &(V->re[idx]), 
						                             &(V->im[idx]),
						                             &(D->re[idx]),
						                             &(D->im[idx]),
						                             &(W->re[idx]),
						                             &(W->im[idx]));
				
					if (nread != 6) {
						ERROR_MSG("Wrong input file format.");
					}
				}
			}
		}
		
		return s;
	}
	
	return NULL;
}

struct sxs_spf_full* sxs_spf_full_read(char* path)
{
	FILE* f = fopen(path, "r");
	struct sxs_spf_full* s = sxs_spf_full_fread(f);
	fclose(f);
	
	return s;
}





