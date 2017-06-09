/*************************************************************************
CheckCP.c - Check change points sequentially
            Remove and/or replace spurious change points
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_math.h>
#include "changepoint.h"

struct changepoint** DelCPNode();

// Check change points in tree with root node at *cp_root for time series traj
// Change points were found with Type-I error of alpha and confidence interval 
// of beta. Ncpdlt is the number of change points that were deleted. ca describes 
// the critical value for the confidence interval. Na is the length of time series
void CheckCP(struct changepoint** cp_root, double* traj, double alpha, double beta, 
	int* Ncpdlt, int trace, double* ca, int N) {
	int            *cps=NULL, *cpsl=NULL, *cpsr=NULL, cp1;
	int            Ncp=0, Ny=0, Ny2=0;
	int            i, istart=1, iend=0;
	enum bool      deletion_adjustment=true, h=false;
	*Ncpdlt=0;

	// While there are still change points to check
	while(deletion_adjustment){
		deletion_adjustment=false;
		// Make arrays to store change points and their confidence intervals
		MakeCPArray(*cp_root, &cps, &cpsl, &cpsr, &Ncp);

		// Since the LB and RB are exclusive, want to start at 0 and N+1
		cps = (int *) realloc(cps, (Ncp+2)*sizeof(int) );
		cps[0] = -1;
		cps[Ncp+1] = N+2;
		cpsl = (int *) realloc( cpsl, (Ncp+2)*sizeof(int) );
		cpsl[0] = -1;
		cpsl[Ncp+1] = N+2;
		cpsr = (int *) realloc( cpsr, (Ncp+2)*sizeof(int) );
		cpsr[0] = -1;
		cpsr[Ncp+1] = N+2;

		// Pick up deletion adjustment where it last left off or at 1
		if(iend > 1) istart=iend;
		else istart=1;

		// Iterate through change points in array
		// Check by deleting change point and finding it again
		for(i=istart;i<Ncp+1;i++){
			// Delete change point
			*cp_root = DelCPNode(cps[i], (*cp_root), &h);
			Ncp--;
			// If change point is less than right bound of previous, it is invalid
			if(cps[i]<cpsr[i-1]){
				deletion_adjustment=true;
				(*Ncpdlt)++;
				fprintf(stderr, "Illegal p at %d: less than right bound of previous cp\n", 
					cps[i]);
				break;
			}
			// Find change point
			cp1=FindCP(cp_root,traj,cps[i-1]+1,cps[i+1]-1,alpha,beta,&Ncp,ca,N,3);

			// No change point found
			if (cp1==0){
				fprintf(stderr, "%d/%d Checking node %d %d %d...", i, Ncp, cps[i-1]+1,cps[i], 
					cps[i+1]-1);
				fprintf(stderr, "not found! RESTARTING\n", cp1);
				deletion_adjustment=true;
				(*Ncpdlt)++;
				iend=i;
				break;
			}
			// Change point found at different index
			else if(cp1!=cps[i]){
				fprintf(stderr, "%d/%i Checking node %d %d %d...", i, Ncp, cps[i-1]+1, cps[i],
					cps[i+1]-1);
				fprintf(stderr,"found at %d.\n",  cp1);
				deletion_adjustment=true;
				iend=i;
				break;
			}
		}
		free(cps);
		free(cpsr);
		free(cpsl);
		cps=NULL;
		cpsr=NULL;
		cpsl=NULL;
		Ncp = 0;
	}
}