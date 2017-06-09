/*************************************************************************
FindCP.c  -  Change Point Detection Method
-------------------
last modified   : Mon Dec 19 2016 (NS)
email           : hawyang@princeton.edu, nsong@princeton.edu
*************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_erf.h>
#include "changepoint.h"

#define MAX_X 10.0
#define DX    0.01
#define TOL   1.0e-8

// Functions for calculation of critical values
double C_ac();
double Vostrikova();
struct changepoint* AddCPNode();

// Find a change point in traj with bounds cpl (left) and cpr (right), exclusive. 
// Type-I error of alpha and confidence interval of beta. Found change point added to 
// tree pointed to by cp. ca points to confidence region critical values. Na is length
// of traj, rc indicates recursive tracing
int FindCP(struct changepoint **cp, double *traj, int cpl, int cpr, double alpha, 
	double beta, int * Ncp, double *ca, int Na, int rc ) {
	int n,i,k,k_max;				//indices to iterate through traj	
	int LB, RB;                 	//left and right bounds of confidence region
	int cp2=0, cp1=0, cp_max;   	//largest change point in absolute index

	//Relevant variables for gaussian mean change point detection
	double kL, nL;					// Double analogs in integers
	double x_old, dx, mean_dx, var_dx; // Used for calculation of log-likelihood
	double critical_region;			// Critical value for specified alpha and N
	double llrt_max=0.0;			// Maximum log-likelihood
	double *cumsum, *llrt;			// Arrays to store values
	enum bool dummy_bool = false;	// If applicable, denotes if found CP has been inserted
	
	n = cpr - cpl;
	LB = cpl;
	RB = cpr;
	cp_max = 0;
	if ( n > 1 ) {
		nL = (double) n;
		cumsum = (double *) malloc( (n+1) * sizeof(double) );
    	llrt = (double *) malloc( (n+1) * sizeof(double) );

		// Calculate the critical region using Horvath's approximation
		critical_region = C_ac(alpha,n);
		
		if(cpl == 0)
			x_old = traj[1];
		else
			x_old = traj[cpl];
		cumsum[0] = 0;

	
		for (i=1, mean_dx=var_dx=0.0; i<n; i++) {
			// estimate Gaussian means at different k
			cumsum[i] = cumsum[i-1] + traj[cpl+i];
			// estimate sample variance using gaussian difference
			dx = traj[cpl+i] - x_old;
			mean_dx += dx;
			var_dx += dx*dx;
			x_old = traj[cpl+i];
		}
		var_dx = fabs(var_dx - mean_dx*mean_dx)/(nL-1.0);
		
		for ( k=1,k_max=0; k<n; k++ ) {
			kL = (double) k;

			// traditional Generalized Log-likelihood ratio test
			llrt[k] = cumsum[k] * cumsum[k] / kL + (cumsum[n-1]-cumsum[k]) * (cumsum[n-1]-cumsum[k]) / 
				(nL-kL-1.0) - cumsum[n-1] * cumsum[n-1] / (nL-1.0);

			llrt[k] = 2.0*llrt[k] / var_dx;
			if ( llrt[k] > llrt_max ) {
			// find the maximum of the likelihood functions
				llrt_max = llrt[k];
				k_max = k;
			}
		}

		// Compare maximum log likelihood ratio to critical region
		if ( llrt_max > critical_region ) {
			// Find the left-hand bound of critical region
			k = k_max;
			while ( (k>1) && (llrt_max - llrt[k-1] < ca[n]) )
				k--;
			LB = k;

			// Find the right-hand bound of critical region
			k = k_max;
			while ( (k < n-1) && (llrt[k+1]+ca[n]-llrt_max > 0) ) 
				k++;
			RB = k;
			cp_max = cpl + RB;

			// Add change point to current tree
			*cp = AddCPNode(k_max+cpl, LB+cpl, RB+cpl, *cp, Ncp, &dummy_bool);
			// No recursion for checking change points
			if(rc==3){
				/*********************
				Free up used workspace
				*********************/
				// Return index of change point found
				return k_max+cpl;
			}
			// Recursively find change points
			if((rc==1)||(rc==0)){
				// Go to the left branch
				cp1 = FindCP(cp, traj, cpl, LB+cpl, alpha, beta, Ncp, ca, Na, 1);
				// Go to the right branch
				cp1 = FindCP(cp, traj, RB+cpl, cpr, alpha, beta, Ncp, ca, Na, 1);
				if (cp1 > cp_max ) 
					cp_max = cp1;
			}
		}
		free(llrt);
		free(cumsum);
	}
	// Return maximum index of change point found
	return cp_max;
}
