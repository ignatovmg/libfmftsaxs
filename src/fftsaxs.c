#include "fftsaxs.h"

static inline int d_symb_index(int l, int m, int l1, int p, int L) 
{
	return (lm_index(l, m) * (L+1) + l1) * (2*L+1) + p;
}

static inline int t_index(int q, int m, int l, int l1, int L) 
{
	return ((q * (L+1) + abs(m)) * (L+1) + l) * (L+1) + l1;
}

static inline int sum1_index(int q, int m, int m2, int l, int L) 
{
	return ((q * (2*L+1) + m + L) * (2*L+1) + m2 + L) * (L+1) + l;
}

static inline int sum2_index(int m, int m1, int m2, int L) 
{
	m  = (m  < 0 ? (2*L+1+m ) : m );
	m1 = (m1 < 0 ? (2*L+1+m1) : m1);
	m2 = (m2 < 0 ? (2*L+1+m2) : m2);
	
	return (m * (2*L+1) + m1) * (2*L+1) + m2;
}

static double* comp_const_int(struct sxs_spf_sing** A1, 
                              struct sxs_spf_sing** A2, 
                              struct sxs_spf_sing** B1, 
                              struct sxs_spf_sing** B2, 
                              int qnum)
{
	int index;
	double* in = (double*)calloc(qnum, sizeof(double));
	
	for (int q = 0; q < qnum; q++) {
		in[q] = 0.0;
		for (int l = 0; l <= A1[q]->L; l++) 
			for(int m = -l; m <= l; m++) {
				index = lm_index(l, m);
				in[q] += A1[q]->re[index]*A2[q]->re[index]+
						 A1[q]->im[index]*A2[q]->im[index];
						 
				in[q] += B1[q]->re[index]*B2[q]->re[index]+
						 B1[q]->im[index]*B2[q]->im[index];
			}
	}
	
	return in;
}

static void fill_const(struct sxs_profile** profiles, 
					double* const_int_VV, 
					double* const_int_VD, 
					double* const_int_VW, 
					double* const_int_DD, 
					double* const_int_DW,
					double* const_int_WW, 
					int* mask,
					int L)
{
	double *VV, *VD, *VW, *DD, *DW, *WW;
	struct sxs_profile* profile;
	
	for (int i = 0; i < (2*L+1)*(2*L+1)*(2*L+1); i++) {
		if (mask[i] == 1) {
			profile = profiles[i];
			VV = profile->VV;
			VD = profile->VD;
			VW = profile->VW;
			DD = profile->DD;
			DW = profile->DW;
			WW = profile->WW;
		
			for (int q = 0; q < profile->qnum; q++) {
				VV[q] = const_int_VV[q];
				VD[q] = const_int_VD[q] * 2.0;
				VW[q] = const_int_VW[q] * 2.0;
				DD[q] = const_int_DD[q];
				DW[q] = const_int_DW[q] * 2.0;
				WW[q] = const_int_WW[q];
			}
		}
	}
}

static void fill_var(struct sxs_profile** profiles, fftw_complex* fft, int* mask, int q, int L)
{
	int offset[6];
	int fft_size = (2*L+1)*(2*L+1)*(2*L+1);
	struct sxs_profile* profile;
	
	for (int i = 0; i < 6; i++) {
		offset[i] = i * fft_size;
	}
	
	for (int i = 0; i < fft_size; i++) {
		if (mask[i] == 1) {
			profile = profiles[i];
			profile->VV[q] += 2.0 * creal(fft[offset[0] + i]);
			profile->VD[q] += 2.0 * creal(fft[offset[1] + i]);
			profile->VW[q] += 2.0 * creal(fft[offset[2] + i]);
			profile->DD[q] += 2.0 * creal(fft[offset[3] + i]);
			profile->DW[q] += 2.0 * creal(fft[offset[4] + i]);
			profile->WW[q] += 2.0 * creal(fft[offset[5] + i]);
		}
	};
}

static void besselj_n(double* bessel, double zval, double* qvals, int qnum, int L)
{
	CHECK_PTR(bessel);

	double qval;
	for (int q = 0; q < qnum; q++) {
		qval = qvals[q];
		for (int p = 0; p < 2*L+1; p++) {
			bessel[q * (2*L+1) + p] = sf_bessel(p, zval * qval);
		}
	}
}

/**
 * Compute d-symbol d_lm(l1,p).
 *
 *  / l  p  l1  \     / l  p  k \
 *                *              * (2p+1)*((2l1+1)(2l+1))^0.5 = d_lm(l1,p)
 *  \ 0  0   0  /     \-m  0  m /
 *
 * p : [0, 2*L]
 * l : [0, L]
 * l1: [0, L]
 * m : [-l, l]
 *
 * Dimensions:
 * [l, m, l1, p]
 */
static void compute_d_symb(double* d_symb, int L, mpf_t* fact)
{
	CHECK_PTR(d_symb);

	double k1, k2, k3, k4;
	int l, m, l1, p, idx1;
	
	mpf_t w3j1, w3j2;
	mpf_init(w3j1);
	mpf_init(w3j2);

	for (l = 0; l < L+1; l++) {
		k1 = 2*l+1;
		
		for (l1 = 0; l1 < L+1; l1++) {
			k2 = sqrt((2*l1+1) * k1);
			
			for (p = abs(l - l1); p <= l + l1; p++) {
				//k3 = (2*p+1) * k2 * wigner_3j_symbol_arb(l, p, l1, 0, 0, 0, fact).get_d();
				wigner_3j_symbol_arb(w3j1, l, p, l1, 0, 0, 0, fact);
				k3 = (2*p+1) * k2 * mpf_get_d(w3j1);
				
				for (m = -min(l,l1); m <= min(l,l1); m++) {
					//k4 = wigner_3j_symbol_arb(l, p, l1, -m, 0, m, fact).get_d();
					wigner_3j_symbol_arb(w3j2, l, p, l1, -m, 0, m, fact);
					k4 = mpf_get_d(w3j2);
					
					idx1 = d_symb_index(l, m, l1, p, L);
					d_symb[idx1] = k3 * k4;
				}
			}
		}
	}
	
	mpf_clear(w3j1);
	mpf_clear(w3j2);
}

/** 
 * Computes translation matrix T^{|m|}_{ll_1}(qz)
 *
 * Dimensions in row-major order:
 * [q, abs(m), l, l1]
 */
static void fill_t_matrix(
	double* t_matrix_re, 
	double* t_matrix_im, 
	double* d_symb, 
	double* bessel, 
	double* qvals, 
	int qnum, int L)
{
	CHECK_PTR(t_matrix_re);
	CHECK_PTR(t_matrix_im);
	CHECK_PTR(d_symb);
	CHECK_PTR(bessel);
	CHECK_PTR(qvals);

	double i_pow_p_re[] = {1.0, 0.0, -1.0, 0.0};
	double i_pow_p_im[] = {0.0, 1.0,  0.0,-1.0};
	double bin[] = {1.0, -1.0};
	
	double min_one_pow_m;
	double re_val;
	double im_val;
	double val;
	
	int m, l, l1, p, q, bsidx, didx, tidx1, tidx2;
	for (q = 0; q < qnum; q++) {
		bsidx = q*(2*L+1);
		
		for (m = 0; m <= L; m++) { // T(-m) = T(m)
			min_one_pow_m = bin[m % 2];
			tidx1 = (q * (L+1) + m) * (L+1);
			
			for (l = m; l <= L; l++) {
				for (l1 = l; l1 <= L; l1++) { // T(l, l1) = T(l1, l) -> (l1 = l -> L)
					re_val = 0.0;
					im_val = 0.0;
				
					didx = d_symb_index(l, m, l1, 0, L);
					for (p = abs(l - l1); p <= l + l1; p++) {
						val = min_one_pow_m * d_symb[didx + p] * bessel[bsidx + p];
						
						//double dif1 = d_symb[d_symb_index(l, m, l1, p, L)];
						//double dif2 = d_symb[d_symb_index(l, -m, l1, p, L)];
						//if (dif1 - dif2 > 0.0000001) { printf("%.1e %.1e %.1e\n", dif1 - dif2, dif1, dif2); }
						
						re_val += i_pow_p_re[p % 4] * val;
						im_val += i_pow_p_im[p % 4] * val;
					}
					
					tidx2 = (tidx1 + l) * (L+1) + l1; // TODO: change
					//tidx2 = t_index(q, m, l, l1, L);
					t_matrix_re[tidx2] = re_val;
					t_matrix_im[tidx2] = im_val;
					
					tidx2 = (tidx1 + l1) * (L+1) + l; // TODO: change
					//tidx2 = t_index(q, m, l1, l, L);
					t_matrix_re[tidx2] = re_val;
					t_matrix_im[tidx2] = im_val;
				}
			}
		}
	}
}

/**
 * Compute first-sum
 *
 * Output dimensions:
 * [q, m, m2, l]
 */
static void compute_sum1(
	double* sum_re, 
	double* sum_im,
	double* t_matrix_re,
	double* t_matrix_im,
	struct d_array* d_wigner,
	struct sxs_spf_sing** B, int qnum, int L)
{
	CHECK_PTR(sum_re);
	CHECK_PTR(sum_im);

	double* d_data = d_wigner->data;
	double* B_re, *B_im, *t_loc_re, *t_loc_im;
	
	double B_re_val, 
			B_im_val, 
			val_re,
			val_im,
			d_val;
	
	int square = (L+1)*(L+1);
	int cube = square*(L+1);
	
	int q, l, m2, l1, m;
	int sidx1, sidx2, sidx3, tidx1, tidx2, tidx3, bidx;
	
	for (q = 0; q < qnum; q++) {
		B_re = B[q]->re;
		B_im = B[q]->im;
		t_loc_re = t_matrix_re + q*cube;
		t_loc_im = t_matrix_im + q*cube;
		
		sidx1 = q*(2*L+1);
		
		for (m = -L; m <= L; m++) {
			tidx1 = abs(m) * (L+1);
			sidx2 = (sidx1 + m + L) * (2*L+1);
			
			for (m2 = -L; m2 <= L; m2++) {
				sidx3 = (sidx2 + m2 + L) * (L+1);
				
				for (l = abs(m); l <= L; l++) {
					tidx2 = (tidx1 + l) * (L+1);
					
					val_re = 0.0;
					val_im = 0.0;
					for (l1 = max(abs(m2), abs(m)); l1 <= L; l1++) {
						
						d_val = d_data[index_d_array(L, l1, m, m2)];
						
						bidx = lm_index(l1, m2);
						B_re_val = B_re[bidx];
						B_im_val = B_im[bidx];
						
						tidx3 = tidx2 + l1;
						val_re += d_val * (B_re_val * t_loc_re[tidx3] - B_im_val * t_loc_im[tidx3]); // B, T - conjugate
						val_im -= d_val * (B_re_val * t_loc_im[tidx3] + B_im_val * t_loc_re[tidx3]);
								  
						/*tidx3 = t_index(q, abs(m), l, l1, L);
						val_re += d_val * B_re_val * t_matrix_re[tidx3] - 
								  d_val * B_im_val * t_matrix_im[tidx3]; // B, T - conjugate
								  
						val_im -= d_val * B_re_val * t_matrix_im[tidx3] + 
								  d_val * B_im_val * t_matrix_re[tidx3];*/		  
						//if (q > 3)		  
						//printf("\t(%i, %i, %i, %i): d_val = %e, A_re = %e, t = %e\n", m, m2, l, l1, d_val, B_re_val,  t_loc_re[tidx3]);
					}	
					
					
					sum_re[sidx3 + l] = val_re;
					sum_im[sidx3 + l] = val_im;
					//sum_re[sum1_index(q, m, m2, l, L)] = val_re;
					//sum_im[sum1_index(q, m, m2, l, L)] = val_im;
					
					//if (q > 3)		
					//printf("%e %e\n", val_re, val_im);
				}
				//if (q > 3)		
				//getchar();
			}
		}
	}	
}

/**
 * Compute FFT coefficients.
 *
 * Output dimensions:
 * [m, m1, m2]
 */
static void compute_sum2_old(
	fftw_complex* sum2,
	double* sum1_re,
	double* sum1_im,
	struct d_array* d_wigner,
	struct sxs_spf_sing* A, int q, int L)
{
	CHECK_PTR(sum2);
	
	double* d_data = d_wigner->data;
	double* A_re, *A_im;
	
	double A_re_val, 
			A_im_val, 
			val_re, 
			val_im,
			d_val;
	
	int l, m2, m1, m;
	int sum1idx1, sum1idx2, sum1idx3, sum1idx4;
	int sum2idx1, sum2idx2, sum2idx3, aidx;
	
	A_re = A->re;
	A_im = A->im;

	sum1idx1 = q * (2*L+1);
	for (m = -L; m <= L; m++) {
		sum1idx2 = (sum1idx1 + m + L) * (2*L+1);
		sum2idx1 = (m < 0 ? (2*L+1+m) : m) * (2*L+1);
	
		for (m1 = -L; m1 <= L; m1++) {
			sum2idx2 = ((m1 < 0 ? (2*L+1+m1) : m1) + sum2idx1) * (2*L+1);
		
			for (m2 = -L; m2 <= L; m2++) {
				sum1idx3 = (sum1idx2 + m2 + L) * (L+1);
				sum2idx3 = (m2 < 0 ? (2*L+1+m2) : m2) + sum2idx2;

				val_re = 0.0;
				val_im = 0.0;
				for (l = max(abs(m1), abs(m)); l < L+1; l++) {
					d_val = d_data[index_d_array(L, l, m, m1)];

					aidx = lm_index(l, m1);
					A_re_val = A_re[aidx];
					A_im_val = A_im[aidx];
					
					sum1idx4 = sum1idx3 + l;
					//sum1idx4 = sum1_index(q, m, m2, l, L);
					val_re += d_val * A_re_val * sum1_re[sum1idx4] - 
							  d_val * A_im_val * sum1_im[sum1idx4];
							  
					val_im += d_val * A_re_val * sum1_im[sum1idx4] +
							  d_val * A_im_val * sum1_re[sum1idx4];
					
					//if (q > 3)		  
					//printf("\t(%i, %i, %i, %i): d_val = %e, A_re = %e, %e\n", m, m1, m2, l, d_val, A_re_val, sum1_re[sum1idx4]);
				}	
				
				//sum2idx3 = sum2_index(m, m1, m2, L);
				sum2[sum2idx3] = val_re + I * val_im;
				//if (q > 3)		
				//printf("%i %e %e\n", sum2idx3, val_re, val_im);
			}
			//if (q > 3)		
			//getchar();
		}
	}
}

/**
 * Compute FFT coefficients.
 *
 * Output dimensions:
 * 6 * [m, m1, m2]
 */
static void compute_sum2(
	fftw_complex* sum2,
	double* sum_v_re, double* sum_v_im,
	double* sum_d_re, double* sum_d_im,
	double* sum_w_re, double* sum_w_im,
	struct d_array* d_wigner,
	struct sxs_spf_full* A, int q, int L)
{
	CHECK_PTR(sum2);
	
	double* d_data = d_wigner->data;
	int fft_rank = (2*L+1);
	int d_stride = fft_rank * fft_rank;
			
	int offset[6];
	int fft_size = fft_rank * fft_rank * fft_rank;
	
	for (int i = 0; i < 6; i++) {
		offset[i] = i * fft_size;
	}
	
	double val_re[6], val_im[6];
	size_t val_size = 6 * sizeof(double);
	
	double *V_re, *V_im, *D_re, *D_im, *W_re, *W_im;
	V_re = A->V[q]->re;
	V_im = A->V[q]->im;
	D_re = A->D[q]->re;
	D_im = A->D[q]->im;
	W_re = A->W[q]->re;
	W_im = A->W[q]->im;
	
	double Vval_re, Vval_im, 
	       Dval_re, Dval_im, 
	       Wval_re, Wval_im, 
	       d_val;
	       
	int l, m2, m1, m;
	int sum1idx1, sum1idx2, sum1idx3, sum1idx4;
	int sum2idx1, sum2idx2, sum2idx3, aidx, didx;	

	sum1idx1 = q * (2*L+1);
	for (m = -L; m <= L; m++) {
		sum1idx2 = (sum1idx1 + m + L) * (2*L+1);
		sum2idx1 = (m < 0 ? (2*L+1+m) : m) * (2*L+1);
	
		for (m1 = -L; m1 <= L; m1++) {
			sum2idx2 = ((m1 < 0 ? (2*L+1+m1) : m1) + sum2idx1) * (2*L+1);
		
			for (m2 = -L; m2 <= L; m2++) {
				sum1idx3 = (sum1idx2 + m2 + L) * (L+1);
				sum2idx3 = (m2 < 0 ? (2*L+1+m2) : m2) + sum2idx2;
				didx = (L + m) * fft_rank + (L + m1);

				memset(val_re, 0, val_size);
				memset(val_im, 0, val_size);
				
				for (l = max(abs(m1), abs(m)); l < L+1; l++) {
					d_val = d_data[didx + l * d_stride];
					
					aidx = lm_index(l, m1);
					Vval_re = V_re[aidx];
					Vval_im = V_im[aidx];
					Dval_re = D_re[aidx];
					Dval_im = D_im[aidx];
					Wval_re = W_re[aidx];
					Wval_im = W_im[aidx];
					
					sum1idx4 = sum1idx3 + l;
					val_re[0] += d_val * (  Vval_re * sum_v_re[sum1idx4] - Vval_im * sum_v_im[sum1idx4]);            
					val_im[0] += d_val * (  Vval_re * sum_v_im[sum1idx4] + Vval_im * sum_v_re[sum1idx4]);
					
					val_re[1] += d_val * (  Vval_re * sum_d_re[sum1idx4] + Dval_re * sum_v_re[sum1idx4]
					                      - Vval_im * sum_d_im[sum1idx4] - Dval_im * sum_v_im[sum1idx4]);	  
					val_im[1] += d_val * (  Vval_re * sum_d_im[sum1idx4] + Dval_re * sum_v_im[sum1idx4]
					                      + Vval_im * sum_d_re[sum1idx4] + Dval_im * sum_v_re[sum1idx4]);

					val_re[2] += d_val * (  Vval_re * sum_w_re[sum1idx4] + Wval_re * sum_v_re[sum1idx4]
					                      - Vval_im * sum_w_im[sum1idx4] - Wval_im * sum_v_im[sum1idx4]);	  
					val_im[2] += d_val * (  Vval_re * sum_w_im[sum1idx4] + Wval_re * sum_v_im[sum1idx4]
					                      + Vval_im * sum_w_re[sum1idx4] + Wval_im * sum_v_re[sum1idx4]);
					                      
					val_re[3] += d_val * (  Dval_re * sum_d_re[sum1idx4] - Dval_im * sum_d_im[sum1idx4]);
					val_im[3] += d_val * (  Dval_re * sum_d_im[sum1idx4] + Dval_im * sum_d_re[sum1idx4]);

					val_re[4] += d_val * (  Dval_re * sum_w_re[sum1idx4] + Wval_re * sum_d_re[sum1idx4] 
					                      - Dval_im * sum_w_im[sum1idx4] - Wval_im * sum_d_im[sum1idx4]);	  
					val_im[4] += d_val * (  Dval_re * sum_w_im[sum1idx4] + Wval_re * sum_d_im[sum1idx4]
					                      + Dval_im * sum_w_re[sum1idx4] + Wval_im * sum_d_re[sum1idx4]);
					                      
					val_re[5] += d_val * (Wval_re * sum_w_re[sum1idx4] - Wval_im * sum_w_im[sum1idx4]);	  
					val_im[5] += d_val * (Wval_re * sum_w_im[sum1idx4] + Wval_im * sum_w_re[sum1idx4]);
				}	
				
				sum2[sum2idx3 + offset[0]] = val_re[0] + I * val_im[0];
				sum2[sum2idx3 + offset[1]] = val_re[1] + I * val_im[1];
				sum2[sum2idx3 + offset[2]] = val_re[2] + I * val_im[2];
				sum2[sum2idx3 + offset[3]] = val_re[3] + I * val_im[3];
				sum2[sum2idx3 + offset[4]] = val_re[4] + I * val_im[4];
				sum2[sum2idx3 + offset[5]] = val_re[5] + I * val_im[5];
			}
		}
	}
}

/**
 * Nested loop 3D Fourier transform implementation.
 */
static void simple_ft(
    fftw_complex* fft,
    fftw_complex* sum2,
    int* mask, int L)
{
	int ft_rank = 2*L+1;
	int ft_size = ft_rank * ft_rank * ft_rank;
	double step = 2*M_PI / ft_rank;
	
	int offset[6];
	for (int i = 0; i < 6; i++) {
		offset[i] = i * ft_size;
	}
	
	double* exp_re = calloc(ft_rank, sizeof(double));
	double* exp_im = calloc(ft_rank, sizeof(double));
	for (int i = 0; i < ft_rank; i++) {
		exp_re[i] =  cos(i * step);
		exp_im[i] = -sin(i * step);
	}

	int tmp, mm1m2, id;
	int a2, g1, g2, m, m1, m2, exp_id;
	double exp_re_val, exp_im_val;
	double val[6];
	
	for (int i = 0; i < ft_size; i++) {
	
		if (mask[i] == 1) {
			g2  = i % ft_rank;
			tmp = i / ft_rank;
			g1  = tmp % ft_rank;
			a2  = tmp / ft_rank;
			
			memset(val, 0, 6 * sizeof(double));
			
			exp_id = 0;
			mm1m2  = 0;
			for (m = 0; m < ft_rank; m++) {
				for (m1 = 0; m1 < ft_rank; m1++) {
					for (m2 = 0; m2 < ft_rank; m2++) {
					
						exp_re_val = exp_re[exp_id % ft_rank];
						exp_im_val = exp_im[exp_id % ft_rank];
						
						val[0] += exp_re_val * creal(sum2[offset[0] + mm1m2]) - 
						          exp_im_val * cimag(sum2[offset[0] + mm1m2]);
						val[1] += exp_re_val * creal(sum2[offset[1] + mm1m2]) - 
						          exp_im_val * cimag(sum2[offset[1] + mm1m2]);
						val[2] += exp_re_val * creal(sum2[offset[2] + mm1m2]) - 
						          exp_im_val * cimag(sum2[offset[2] + mm1m2]);
						val[3] += exp_re_val * creal(sum2[offset[3] + mm1m2]) - 
						          exp_im_val * cimag(sum2[offset[3] + mm1m2]);
						val[4] += exp_re_val * creal(sum2[offset[4] + mm1m2]) - 
						          exp_im_val * cimag(sum2[offset[4] + mm1m2]);
						val[5] += exp_re_val * creal(sum2[offset[5] + mm1m2]) - 
						          exp_im_val * cimag(sum2[offset[5] + mm1m2]);
						
						exp_id += g2;
						mm1m2++;
					}
					exp_id += g1;
				}
				exp_id += a2;
			}
			
			fft[offset[0] + i] = val[0];
			fft[offset[1] + i] = val[1];
			fft[offset[2] + i] = val[2];
			fft[offset[3] + i] = val[3];
			fft[offset[4] + i] = val[4];
			fft[offset[5] + i] = val[5];
		}
	}
	
	free(exp_re);
	free(exp_im);
}

void compute_saxs_scores(
	double* scores_list,
	double*     c1_list,
	double*     c2_list,
	int* index_list, int nout,
	struct sxs_spf_full* A,
	struct sxs_spf_full* B,
	struct sxs_opt_params* params,
	double* qvals, int qnum, 
	double* zvals, int znum, 
	int L, int skip) 
{
#ifdef _MPI_	
	int myrank = 0;
	int mpi_flag;
	if (MPI_Initialized(&mpi_flag) != MPI_SUCCESS) {
		MPI_Abort (MPI_COMM_WORLD, 1);
	}
	if (mpi_flag == 1) {
		if (MPI_Comm_rank(MPI_COMM_WORLD, &myrank) != MPI_SUCCESS) {
			MPI_Abort (MPI_COMM_WORLD, 1);
		}
	}
#endif			
	int nbeta = L + 1;
	
	double beta_bgn = 0.0;
	double beta_end = M_PI;
	double beta_step = (beta_end - beta_bgn) / (nbeta - 1);
	
	double* const_int_VV = comp_const_int(A->V, A->V, B->V, B->V, qnum);
	double* const_int_VD = comp_const_int(A->V, A->D, B->V, B->D, qnum);
	double* const_int_VW = comp_const_int(A->V, A->W, B->V, B->W, qnum);
	double* const_int_DD = comp_const_int(A->D, A->D, B->D, B->D, qnum);
	double* const_int_DW = comp_const_int(A->D, A->W, B->D, B->W, qnum);
	double* const_int_WW = comp_const_int(A->W, A->W, B->W, B->W, qnum);

	double* bessel = calloc(qnum*(2*L+1), sizeof(double));
	
	double* d_symb = calloc((L+1)*(L+1)*(L+1)*(2*L+1), sizeof(double));
	
	// small Wigner d-matrices for each beta in [0, PI]
	struct d_array** d_wigner = calloc(nbeta, sizeof(struct d_array*));
	
	// translation matrix
	double* t_matrix_re = calloc(qnum*(L+1)*(L+1)*(L+1), sizeof(double));
	double* t_matrix_im = calloc(qnum*(L+1)*(L+1)*(L+1), sizeof(double));
	
	// sums over B coefficients
	int sum_stride = qnum*(2*L+1)*(2*L+1)*(L+1);
	double** sum_v_re = calloc(nbeta, sizeof(double*));
	double** sum_v_im = calloc(nbeta, sizeof(double*));
	double** sum_d_re = calloc(nbeta, sizeof(double*));
	double** sum_d_im = calloc(nbeta, sizeof(double*));
	double** sum_w_re = calloc(nbeta, sizeof(double*));
	double** sum_w_im = calloc(nbeta, sizeof(double*));

	for (int i = 0; i < nbeta; i++) {
		sum_v_re[i] = calloc(sum_stride, sizeof(double));
		sum_v_im[i] = calloc(sum_stride, sizeof(double));
		sum_d_re[i] = calloc(sum_stride, sizeof(double));
		sum_d_im[i] = calloc(sum_stride, sizeof(double));
		sum_w_re[i] = calloc(sum_stride, sizeof(double));
		sum_w_im[i] = calloc(sum_stride, sizeof(double));
	}
	
	// fft setup
	int fft_size = (2*L+1)*(2*L+1)*(2*L+1);
	int fft_many_size[] = {2*L+1, 2*L+1, 2*L+1};

	fftw_complex* sum2 = (fftw_complex*) fftw_malloc(6 * fft_size * sizeof(fftw_complex));
	fftw_complex* fft  = (fftw_complex*) fftw_malloc(6 * fft_size * sizeof(fftw_complex));
	fftw_plan plan = fftw_plan_many_dft(3, fft_many_size, 6, 
										sum2, NULL, 1, fft_size, 
										fft,  NULL, 1, fft_size, 
										FFTW_FORWARD, FFTW_PATIENT);
	// allocate profiles
	struct sxs_profile** profiles = calloc(fft_size, sizeof(struct sxs_profile*)); 
	for (int i = 0; i < fft_size; i++) {
		profiles[i] = sxs_profile_create(qvals, qnum, 1);
	}
	
	// allocate storage
	int size5d = fft_size * nbeta * nbeta;
	int*      index_5d = calloc(size5d, sizeof(int));
	double*  scores_5d = calloc(size5d, sizeof(double));
	double*      c1_5d = calloc(size5d, sizeof(double));
	double*      c2_5d = calloc(size5d, sizeof(double));
	
	// mask to skip minimization of redundant profiles
	int* mask = calloc(fft_size, sizeof(int));
	
	// compute d_lm(l1, p)
	mpf_t* fact = generate_desc_fact_array_arb(4*L+1);

	compute_d_symb(d_symb, L, fact);
	
	// compute d-matrices
	for (int i = 0; i < nbeta; i++) {
		d_wigner[i] = generate_d_array(L, beta_bgn + i * beta_step);
	}
	
	// skip conformations not present in the list
	int* z_skip = NULL;
	int* zb1_skip = NULL;
 	int* zb2_skip = NULL;
 	int* zb1b2_skip = NULL;
	
	if (skip == 1) {
		z_skip  = calloc(znum,   sizeof(int));
		zb1_skip = calloc(znum * nbeta, sizeof(int));
		zb2_skip = calloc(znum * nbeta, sizeof(int));
		zb1b2_skip  = calloc(znum * nbeta * nbeta, sizeof(int));
		
		int list_id;
		struct sxs_index saxs_id;
		for (int i = 0; i < nout; i++) {
		
			sxs_disassemble_index (&saxs_id, index_list[i], nbeta, L);
			
		 	z_skip[saxs_id.z] = 1;
		 	zb1_skip[saxs_id.z * nbeta + saxs_id.b1] = 1;
		 	zb2_skip[saxs_id.z * nbeta + saxs_id.b2] = 1;
		 	zb1b2_skip[(saxs_id.z * nbeta + saxs_id.b1) * nbeta + saxs_id.b2] = 1;
		}
	}
	
	// z loop
	int counter = 0;
	for (int z_idx = 0; z_idx < znum; z_idx++) { 
	
		double zval = zvals[z_idx];
		
		if (skip == 1) {
			if (z_skip[z_idx] != 1) {
#ifdef _MPI_
				SXS_PRINTF("process %i: skipping z = %.2f\n", myrank, zval);	
#else			
				SXS_PRINTF("skipping z: %.2f\n", zval);
#endif						
				continue;
			}
		}
		
#ifdef _MPI_	
		SXS_PRINTF("process %i: z = %.2f\n", myrank, zval);
#else
		SXS_PRINTF("z = %.2f\n", zval);
#endif		
		
		
		besselj_n(bessel, zval, qvals, qnum, L);
		
		fill_t_matrix(
			t_matrix_re, 
			t_matrix_im, 
			d_symb, bessel, 
			qvals, qnum, L);
			
		// precompute sums for ligand beta angles
		for (int beta2 = 0; beta2 < nbeta; beta2++) {
		
			if (skip == 1) {
				if (zb2_skip[z_idx * nbeta + beta2] != 1) {
					continue;
				}
			}
		
			compute_sum1(
				sum_v_re[beta2], 
				sum_v_im[beta2],
				t_matrix_re, 
				t_matrix_im,
				d_wigner[beta2],
				B->V, qnum, L);
				
			compute_sum1(
				sum_d_re[beta2], 
				sum_d_im[beta2],
				t_matrix_re, 
				t_matrix_im,
				d_wigner[beta2],
				B->D, qnum, L);
				
			compute_sum1(
				sum_w_re[beta2], 
				sum_w_im[beta2],
				t_matrix_re, 
				t_matrix_im,
				d_wigner[beta2],
				B->W, qnum, L);
		}
		
		// main loop
		for (int beta1 = 0; beta1 < nbeta; beta1++) { 
			double b1_val = beta_bgn + beta1 * beta_step;
		
			if (skip == 1) {
				if (zb1_skip[z_idx * nbeta + beta1] != 1) {
#ifdef _MPI_
					//SXS_PRINTFF("process %i: skipping z = %.2f\tbeta1 = %.2f\t\n", myrank, zval, b1_val);	
#else			
					//SXS_PRINTFF("skipping beta1: %.2f\n", b1_val);
#endif						
					continue;
				}
			}
			
#ifdef _MPI_
			//SXS_PRINTFF("process %i: z = %.2f\tbeta1 = %.2f\n", myrank, zval, b1_val);	
#else
			//SXS_PRINTFF("beta1: %.2f\n", b1_val);
#endif
		
			struct d_array* d_beta1 = d_wigner[beta1];
		
			for (int beta2 = 0; beta2 < nbeta; beta2++) {
			
				double b2_val = beta_bgn + beta2 * beta_step;
				
				//if (beta1 != 5 || beta2 != 10) { continue; }
#ifdef _MPI_
			//SXS_PRINTFF("process %i: z = %.2f\tbeta1 = %.2fbeta2 = %.2f\n", myrank, zval, b1_val, b2_val);	
#else
			//SXS_PRINTFF("beta2: %.2f\n", b2_val);
#endif

				if (skip == 1) {
					if (zb1b2_skip[(z_idx * nbeta + beta1) * nbeta + beta2] != 1) {
#ifdef _MPI_
						//SXS_PRINTFF("process %i: skipping z = %.2f\tbeta1 = %.2f\tbeta2 = %.2f\n", myrank, zval, b1_val, b2_val);	
#else			
						//SXS_PRINTFF("skipping beta2: %.2f\n", b2_val);
#endif						
						continue;
					}
				}
				
				clock_t time1 = clock();
				clock_t time2 = clock();

				// fill the mask for minimizer
				int nprofile = 0;
				if (skip == 1) {
					memset(mask, 0, fft_size * sizeof(int));
					
					int tmp;
					struct sxs_index saxs_id;
					
					for (int i = 0; i < nout; i++) {
						tmp = index_list[i];

						
						sxs_disassemble_index(&saxs_id, tmp, nbeta, L);
						if (saxs_id.z == z_idx && saxs_id.b1 == beta1 && saxs_id.b2 == beta2) {
							mask[(saxs_id.a2*(2*L+1) + saxs_id.g1)*(2*L+1) + saxs_id.g2] = 1;
							nprofile++;
						}
					}
				} else {
					for (int i = 0; i < fft_size; i++) {
						mask[i] = 1;
						nprofile++;
					}
				}
				
				//printf("%i profiles\n\t%.5f: (fill mask)\n", nprofile, (double)(clock()-time1) / CLOCKS_PER_SEC);
				//time1 = clock();
				
				fill_const(profiles, 
					const_int_VV, 
					const_int_VD, 
					const_int_VW, 
					const_int_DD, 
					const_int_DW,
					const_int_WW, 
					mask,
					L);
					
				//printf("\t%.5f: (fill_const)\n", (double)(clock()-time1) / CLOCKS_PER_SEC);
				//time1 = clock();

				for (int q = 0; q < qnum; q++) {
					//if (q == 20)
					//	printf("q = %i\nnprofile = %i\n", q, nprofile);
						
					clock_t time = clock();

					compute_sum2(sum2, 
					             sum_v_re[beta2], 
					             sum_v_im[beta2], 
					             sum_d_re[beta2], 
					             sum_d_im[beta2], 
					             sum_w_re[beta2], 
					             sum_w_im[beta2], 
					             d_beta1, A, q, L);
					
					//if (q == 20) {
					//	printf("\t%.5f: (sum2)\n", (double)(clock()-time) / CLOCKS_PER_SEC);
					//}
					//time = clock();
					
					//if (q == 20)
					//	printf("%.5f\n", (double)(clock()-time) / CLOCKS_PER_SEC);
					//time = clock();
					
					//time = clock();
					if (nprofile >= 15) {
						fftw_execute(plan);
					
						/*if (q == 20) {
							printf("\t%.5f: (fast ft)\n", (double)(clock()-time) / CLOCKS_PER_SEC);
						}
						time = clock();*/
					} else {
						simple_ft(fft, sum2, mask, L);
					
						/*if (q == 20) {
							printf("\t%.5f: (smpl ft)\n", (double)(clock()-time) / CLOCKS_PER_SEC);
						}*/
					}
					
					fill_var(profiles, fft, mask, q, L);
					
					//if (q == 20) {
					//	printf("\t%.5f: (fill var)\n", (double)(clock()-time) / CLOCKS_PER_SEC);
					//}
				}
				
				//printf("\t%.5f: (fill FFT)\n", (double)(clock()-time1) / CLOCKS_PER_SEC);
				//time1 = clock();
				
				// minimization over c1, c2
				sxs_fit_params(profiles, params, mask, fft_size);
				
				// keep the result
				int offset = (beta1 * nbeta + beta2) * fft_size;
				for (int i = 0; i < fft_size; i++) {
					scores_5d[offset + i] = profiles[i]->score;
					    c1_5d[offset + i] = profiles[i]->c1;
					    c2_5d[offset + i] = profiles[i]->c2;
				}
				
				//printf("(beta1, beta2) = (%.3f, %.3f): %.5f\n\n", b1_val, b2_val, (double)(clock()-time2) / CLOCKS_PER_SEC);
				//time2 = clock();
			}
		}
		
		// fill the output
		if (nout > 0)
		{
			int saxs_id;
			for (int i = 0; i < nout; i++)
			{
				saxs_id = index_list[i];
				if (saxs_id / size5d == z_idx)
				{
					saxs_id %= size5d;
					scores_list[i] = scores_5d[saxs_id];
					    c1_list[i] = c1_5d[saxs_id];
					    c2_list[i] = c2_5d[saxs_id];
				}
			}
		}
	}
	
	// let them free
	for (int i = 0; i < fft_size; i++) {
		sxs_profile_free(profiles[i]);
	}
	sxs_myfree(profiles);
	
	for (int i = 0; i < nbeta; i++) {
		deallocate_d_array(d_wigner[i]);
	}
	sxs_myfree(d_wigner);
	
	for (int i = 0; i < nbeta; i++) {
		sxs_myfree(sum_v_re[i]);
		sxs_myfree(sum_v_im[i]);
		sxs_myfree(sum_d_re[i]);
		sxs_myfree(sum_d_im[i]);
		sxs_myfree(sum_w_re[i]);
		sxs_myfree(sum_w_im[i]);
	}
	sxs_myfree(sum_v_re);
	sxs_myfree(sum_v_im);
	sxs_myfree(sum_d_re);
	sxs_myfree(sum_d_im);
	sxs_myfree(sum_w_re);
	sxs_myfree(sum_w_im);
	sxs_myfree(const_int_VV);
	sxs_myfree(const_int_VD);
	sxs_myfree(const_int_VW);
	sxs_myfree(const_int_DD);
	sxs_myfree(const_int_DW);
	sxs_myfree(const_int_WW);
	sxs_myfree(bessel);
	sxs_myfree(d_symb);
	sxs_myfree(t_matrix_re);
	sxs_myfree(t_matrix_im);
	fftw_free(sum2);
	fftw_free(fft);
	fftw_destroy_plan(plan);
	sxs_myfree(index_5d);
	sxs_myfree(scores_5d);
	sxs_myfree(c1_5d);
	sxs_myfree(c2_5d);
	sxs_myfree(mask);
	sxs_myfree(z_skip);
	sxs_myfree(zb1_skip);
	sxs_myfree(zb2_skip);
	sxs_myfree(zb1b2_skip);
}

