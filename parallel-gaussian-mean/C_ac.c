/***************************************************************************
	C_ac.c  -  Calculate critical value threshold for given alpha and n
	-------------------
	email                : hawyang@princeton.edu

	See:
	L. Watkins and H. Yang, J Phys Chem B (2005), 109 1):617-28.
		doi:10.1021/jp0467548
	K.J. Worsley, Biometrika (1986), 73(1):91-104. doi: 10.1093/biomet/73.1.91

 ***************************************************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_erf.h>
#include "changepoint.h"

#define MAX_X 10.0
#define DX    0.01
#define TOL   1.0e-8

// Calculate critical value threshold for given alpha and n
double C_ac(double alpha, int n) {
	double h,l,T,thr,x,x_min,x1,x2,p,p_max;

	h = pow(log((double)n),1.5) / (double) n;
	l = h;
	T = log((1.0-h)/h * (1.0-l)/l);
	for (x_min=0.0,p_max=-DBL_MAX,x=0.01;x<=MAX_X;x+=DX) {
		p = x*exp(-x*x/2.0)/sqrt(2.0)/gsl_sf_gamma(0.5)*(T-T/x/x+4.0/x/x);
		if ( p > p_max ) {
			x_min = x;
			p_max = p;
		}
	}
	x1 = x_min - DX;
	if ( x1 < 0.0 ) 
		x1 = 0.0;
	x2 = x_min + DX;
	thr = zbrent( &Vostrikova, x_min, MAX_X, TOL, T, alpha );
	return thr*thr;
}

