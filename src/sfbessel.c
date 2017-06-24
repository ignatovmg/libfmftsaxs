#include "sfbessel.h"

double fact(int x)
{	
	int i;
	if (x == 0) {
		return 1.0;
	} else {
		if (x > 0) {
			double res = 1.0;
			for(i = 1; i <= x; i++)
				res *= (double)i;
			return res;
		}
		else {
			ERROR_MSG("FACTORIAL ERROR\n");
		}
	}		
}

double doublefact(int x)
{	
	int i;
	if (x == 0) {
		return 1.0;
	} else {
		if (x > 0) {
			double res = 1.0;
			for(i = 1; i <= x; i += 2)
				res *= (double)i;
			return res;
		} else {
			ERROR_MSG("FACTORIAL ERROR\n");
			exit(-1);
		}
	}		
}


double sf_bessel(int l, double x)
{
	if (x > 0.0) {
		double two_l_plus_one = 2.0 * (double)l + 1.0;
		double summand = 1.0 / doublefact(two_l_plus_one);
		double summ = summand;
		for(int k = 1; fabs(summand / summ) > 0.00001; k++) {
			summand *= (-1.0) * x * x / (2.0 * (double)k * (2.0 * (double)k + two_l_plus_one));
			summ += summand;
		}
		
		return summ * pow(x, l);
	}
	else
		if (l == 0)
			return 1.0;
		else 
			if(l > 0)
				return 0.0;
		
}
