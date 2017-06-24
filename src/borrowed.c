#include "borrowed.h"

static int min(int i, int j)
{
	if(i<j)
		return i;
	else
		return j;
}

static int max(int i, int j)
{
	if(i>j)
		return i;
	else
		return j;
}

double* generate_spherical_norm(int N)
{
	double *Y_norm;
	Y_norm = (double*) calloc (N*(2*N - 1), sizeof(double));

    for(int l = 0; l < N; ++l)
    {
        double fact = 1.;
        double tmp0 = (2 * l + 1.)/(4. * M_PI);
        for(int m = 0; m <= l; ++m)
        {
        	if (m > 0) fact *= (l - m + 1) * (l + m);
            double tmp = sqrt(tmp0/fact);
            double sgn = (m % 2 == 0) ? (1) : (-1);
            Y_norm[index_LM(N, l,  m)] =         tmp;
            Y_norm[index_LM(N, l, -m)] =   sgn * tmp;
        }
    }
    return Y_norm;
}

void fill_assoc_Legendre_array(int N, double x, double* P)
{
    double y = sqrt(1.-x*x);

    P[index_LM(N, 0, 0)] = 1;//initialization
    if (N == 1) return;
    P[index_LM(N, 1, 0)] = x;
    for (int l = 2; l < N; ++l)  //First recursion for m=0 values.
    {
            P[index_LM(N,l,0)] = (2.*l-1.)/l*x*P[index_LM(N,(l-1),0)] - (l-1.)/l*P[index_LM(N,(l-2),0)];
    }

    for(int m = 1; m < N - 1; ++m)
    {
            P[index_LM(N, m, m)] = -1.*(2.*m-1.)*y*P[index_LM(N,(m-1),(m-1))]; //The P_m^m quantities
            P[index_LM(N,(m+1), m)] = (2.*m+1.)*x*P[index_LM(N,m,m)];//The P_m+1^m quantities.

            for (int l = m + 2; l < N; ++l)  //Recursion for the rest of P_l^m
            {
                    P[index_LM(N,l,m)] = (2.*l-1.)/(l-m)*x*P[index_LM(N,(l-1),m)] - (l+m-1.)/(l-m)*P[index_LM(N,(l-2),m)];
            }
    }

    P[index_LM(N,(N-1),(N-1))] = -1.*(2.*(N-1.)-1.)*y*P[index_LM(N,(N-2),(N-2))]; //The P_bw-1^(+/-)bw-1 quantities
    return;
}







mpf_t* generate_desc_fact_array_arb(const int n)
{
    mpf_t *desc_fact = calloc(n+1, sizeof(mpf_t));
    mpf_init_set_d(desc_fact[0], 1.0);
	
    for (unsigned int i = 1; i <= n; ++i) {
    	mpf_init(desc_fact[i]);
    	mpf_mul_ui(desc_fact[i], desc_fact[i-1], i);
    }

    return desc_fact;
}

void wigner_3j_symbol_arb(mpf_t res, const int j1, const int j2, const int j3, const int m1, const int m2, const int m3, const mpf_t* desc)
{
//======================================================================
// After Wigner3j.m by David Terr, Raytheon, 6-17-04
//
// Compute the Wigner 3j symbol using the Racah formula.
//
//  / j1 j2 j3 \    '
//  |          |    '
//  \ m1 m2 m3 /    '
//
// Reference: Wigner 3j-Symbol entry of Eric Weinstein's Mathworld:
// http://mathworld.wolfram.com/Wigner3j-Symbol.html
//======================================================================

    // Error checking
    if ( ( 2*j1 != floor(2*j1) ) || ( 2*j2 != floor(2*j2) ) || ( 2*j3 != floor(2*j3) ) || ( 2*m1 != floor(2*m1) ) || ( 2*m2 != floor(2*m2) ) || ( 2*m3 != floor(2*m3) ) )
    {
        printf("All arguments must be integers or half-integers.\n");
        exit(0);
    }

    // Additional check if the sum of the second row equals zero
    if ( m1+m2+m3 != 0 )
    {
    	printf("3j-Symbol unphysical\n");
    }

    if ( j1 - m1 != floor ( j1 - m1 ) )
    {
    	printf("2*j1 and 2*m1 must have the same parity\n");
    }

    if ( j2 - m2 != floor ( j2 - m2 ) )
    {
    	printf("2*j2 and 2*m2 must have the same parity\n");
    }

    if ( j3 - m3 != floor ( j3 - m3 ) )
    {
    	printf("2*j3 and 2*m3 must have the same parity\n");
    }

    if ( (j3 > j1 + j2)  || ( j3 < abs(j1 - j2)) )
    {
    	printf("j3 is out of bounds.\n");
    }

    if (abs(m1) > j1)
    {
    	printf("m1 is out of bounds.\n");
    }

    if (abs(m2) > j2)
    {
    	printf("m2 is out of bounds.\n");
    }

    if (abs(m3) > j3)
    {
    	printf("m3 is out of bounds.\n");
    }

    int t1 = j2 - m1 - j3;
    int t2 = j1 + m2 - j3;
    int t3 = j1 + j2 - j3;
    int t4 = j1 - m1;
    int t5 = j2 + m2;

    int tmin = max( 0,  max( t1, t2 ) );
    int tmax = min( t3, min( t4, t5 ) );

	mpf_t wigner; 
	mpf_init_set_d(wigner, 0.0);
	
	mpf_t one; 
	mpf_init(one);
	
	mpf_t bufer; 
	mpf_init(bufer);
	
    for(int t = tmin; t <= tmax; ++t) {
    
    	if((t)%2 == 0) {
    		mpf_set_d(one,  1.0);
    	} else {
    		mpf_set_d(one, -1.0);
    	}

        mpf_mul(bufer, desc[t4-t], desc[t5-t]);
        mpf_mul(bufer, desc[t3-t], bufer);
        mpf_mul(bufer, desc[t-t2], bufer);
        mpf_mul(bufer, desc[t-t1], bufer);
        mpf_mul(bufer, desc[t], bufer);
        
        mpf_div(bufer, one, bufer);
        mpf_add(wigner, wigner, bufer);
        
    	/*wigner += one / ( get_desc_fact_arb(t) *\
    						get_desc_fact_arb(t-t1) *\
    						get_desc_fact_arb(t-t2) *\
    						get_desc_fact_arb(t3-t) *\
    						get_desc_fact_arb(desc,t4-t) *\
    						get_desc_fact_arb(desc,t5-t) );*/
    }
    
	if((j1 - j2 - m3)%2==0) {
		mpf_set_d(one,  1.0);
	} else {
		mpf_set_d(one, -1.0);
	}

	mpf_mul(bufer, desc[j3+m3], desc[j3-m3]);
	mpf_mul(bufer, desc[j2-m2], bufer);
	mpf_mul(bufer, desc[j2+m2], bufer);
	mpf_mul(bufer, desc[j1-m1], bufer);
	mpf_mul(bufer, desc[j1+m1], bufer);
	mpf_div(bufer, bufer, desc[j1+j2+j3+1]);
	mpf_mul(bufer, desc[-j1+j2+j3], bufer);
	mpf_mul(bufer, desc[j1-j2+j3], bufer);
	mpf_mul(bufer, desc[j1+j2-j3], bufer);
	mpf_sqrt(bufer, bufer);
	mpf_mul(bufer, one, bufer);
	mpf_mul(wigner, wigner, bufer);

    /*wigner *= one * sqrt( get_desc_fact_arb(desc, j1+j2-j3) * get_desc_fact_arb(desc, j1-j2+j3) * get_desc_fact_arb(desc, -j1+j2+j3) /\
    						get_desc_fact_arb(desc, j1+j2+j3+1) * get_desc_fact_arb(desc, j1+m1) * get_desc_fact_arb(desc, j1-m1) *\
    						get_desc_fact_arb(desc, j2+m2)      * get_desc_fact_arb(desc, j2-m2) * get_desc_fact_arb(desc, j3+m3) *\
    						get_desc_fact_arb(desc,j3-m3)\
    		            );*/

	mpf_clear(bufer);
	mpf_clear(one);
	
	mpf_set(res, wigner);
	mpf_clear(wigner);
}





struct d_array* allocate_d_array(const int L)
{
	struct d_array* d;
	d = (struct d_array*) calloc (1, sizeof(struct d_array));
	d->L = L;
	d->data = (double*) calloc ((L+1)*(2*L + 1)*(2*L + 1), sizeof(double));
	return d;
}

void deallocate_d_array(struct d_array* d)
{
	free(d->data);
	free(d);
}

struct d_array* generate_d_array(const int L, const double beta)
{
	struct d_array* d;
	d = allocate_d_array(L);
	/* Wolfram Mathematica*/

    d->data[index_d_array(L, 0,  0,  0)] =  1.0;
    d->data[index_d_array(L, 1, -1, -1)] =  1.0 * (1.0 + cos(beta)) / 2.0;
    d->data[index_d_array(L, 1, -1,  0)] =  1.0 * sin(beta) / sqrt(2.0);
    d->data[index_d_array(L, 1, -1,  1)] =  1.0 * (1.0 - cos(beta)) / 2.0;
    d->data[index_d_array(L, 1,  0, -1)] = -1.0 * sin(beta) / sqrt(2.0);
    d->data[index_d_array(L, 1,  0,  0)] =  1.0 * cos(beta);
    d->data[index_d_array(L, 1,  0,  1)] =  1.0 * sin(beta) / sqrt(2.0);
    d->data[index_d_array(L, 1,  1, -1)] =  1.0 * (1.0 - cos(beta)) / 2.0;
    d->data[index_d_array(L, 1,  1,  0)] = -1.0 * sin(beta) / sqrt(2.0);
    d->data[index_d_array(L, 1,  1,  1)] =  1.0 * (1.0 + cos(beta)) / 2.0;

    for(int l = 2; l <= L; ++l)
    {
    	for(int m = -l; m <= l; ++m)
    	{
    		double fact, fact1 = 1.0;
            //fact = sqrt( (2 * l + 1) / 2.);
	        fact=1;
            double d1 = fact;
            double d2 = fact;
            double d3 = fact;
            double d4 = fact;
            /* Eq. 26 Kostelec 2003 */
            for (int k = 1; k <= 2 * l; ++k)
            {
                    fact = sqrt(k / (((k <= l + m) ? k : 1.0) * ((k <= l - m) ? k : 1.0)));
                    fact1 *= fact;
                    d1 *= fact * ((k <= l + m) ? cos(beta / 2.0) : 1.0)* ((k <= l - m) ? -sin(beta / 2.0) : 1.0);
                    d2 *= fact * ((k <= l - m) ? cos(beta / 2.0) : 1.0)* ((k <= l + m) ?  sin(beta / 2.0) : 1.0);
                    d3 *= fact * ((k <= l + m) ? cos(beta / 2.0) : 1.0)* ((k <= l - m) ?  sin(beta / 2.0) : 1.0);
                    d4 *= fact * ((k <= l - m) ? cos(beta / 2.0) : 1.0)* ((k <= l + m) ? -sin(beta / 2.0) : 1.0);
            }

            d->data[index_d_array(L, l,  l,  m)] = d1;
            d->data[index_d_array(L, l, -l,  m)] = d2;
            d->data[index_d_array(L, l,  m,  l)] = d3;
            d->data[index_d_array(L, l,  m, -l)] = d4;

            for(int m1 = -l; m1 <=l; ++m1)
            {

                int j = l - 1;
                double val, val1;
                /* Eq. 28 Kostelec 2003 */
                if ((m1 > -l) && (m1 < l) && (m > -l) && (m < l))
                {
					val  = 1.0;
                    val *= (j + 1) * (2.0 * j + 1) / sqrt(((j + 1) * (j + 1) - m * m) * ((j + 1) * (j + 1) - m1 * m1));
                    val *= (cos(beta) - (double)(m * m1) / (double)(j * (j + 1))) * d->data[index_d_array(L, j, m, m1)];

					val1  = -1.0;
                    val1 *= sqrt((j * j - m * m) * (j * j - m1 * m1));
                    val1 *= (double)(j + 1.0) * d->data[index_d_array(L, j - 1, m, m1)] / (double)j;
                    val1 /= sqrt(((j + 1.0) * (j + 1.0) - (double)(m * m)) * ((j + 1.0) * (j + 1.0) - (double)(m1*m1)));

                    val += val1;
                    d->data[index_d_array(L, l, m, m1)] = val;
                }

            }
        }
    }

    return d;
}



