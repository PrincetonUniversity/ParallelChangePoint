/***************************************************************************
    Vostrikova.c  -  Helper function to find critical values
    -------------------
	email                : hawyang@princeton.edu
	
    See:
    L. Watkins and H. Yang, J Phys Chem B (2005), 109 1):617-28.
		doi:10.1021/jp0467548
	Marc Noe, Ann. Math. Statist (1972), 43(1):58-64.
 ***************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_erf.h>
#include "changepoint.h"

double Vostrikova(double x, double T, double alpha) {
	return x*exp(-x*x/2.0)/sqrt(2.0)/gsl_sf_gamma(0.5)*(T-T/x/x+4.0/x/x)-alpha;
}
