/*************************************************************************
main.c:  Parallelized change point detection
         Originally written using OpenMPI v1.4.2
         Arguments are filename, alpha, beta

filename:   File with Time Series Data
alpha:      Type-I error
1-beta:     Confidence Interval is chosen among {0.99, 0.95, 0.9, 0.69}
*************************************************************************/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <gsl/gsl_version.h>
#include "changepoint.h"
#include "critical_values.h"
#include <mpi.h>
#include <math.h>

#define TOL 1.0e-10

int main(int argc, char *argv[])
{
	FILE *fpin, *fpout;		// Filenames and output names
	char in_name[255],  filename[255], *endptr[255], out_name[255];	// Store strings
	char buffer[65536];     // Input buffer
	double alpha = 0.075;   // Type-I error, mis-specify transition
	double beta = 0.075;    // Confidence interval, mis-specify transition
	double delta_t = 1;     // Time unit between measurements
	int L = 0;              // Total number of data points
	int N_ca;               // Total number of ca
	int ui;                 // Dummy index
	int cpl;                // Left bound 
	int cpr;                // Right bound
	int cp1;
	int Ncp = 0;            // Change points found in each process
	int Ncpdlt;             // Change points removed in each process
	double tmp;             // Temporary variable for input
	double *traj;           // Trajectory or time series
	struct changepoint *cp_root=NULL;        // Change Point Binary Tree
	clock_t time0, time1;   // Variables for timing
	int h,i,j,k;            // Dummy indices
	int trace=0;            // Debug flag
	double *ca=NULL;        // Critical region threshold for confidence interval
	enum bool done = false;
	int NA_BASE;
	int NA_OVERLAP = (int) 200;

	// Variables related to parallelization 
	int AVG_L;                            // Average length passed to one process
	int offset = 0;                       // Offset of time series
	double * split;                       // Temporary storage for time series segment
	int * sNcp;                           // Number of change points from each process
	int** cpArray = malloc(sizeof(int*)); // Store change points from each process
	int** lbArray = malloc(sizeof(int*));
	int** rbArray = malloc(sizeof(int*));
	int rNcp, pCP = 0;                    // Total number of change points found, pooled
	int end, start, first, last, middle;  // Variables used in pooling
	int lb1, lb2, rb1, rb2;					

	// MPI related variables
	int my_id, root_process, ierr, num_procs, an_id, sender;
	int send_data_tag = 2;
	MPI_Request request;
	MPI_Status status;
	int flag;

	// Get command line arguments
	if (argc != 4) {
		// No file name entered in the command line
		printf("\nchangepoint %s%s build %s (GSL version %s)\n",        
						CHANGEPOINT_VERSION, PLATFORM, COMPILE_DATE, GSL_VERSION);
		printf("Syntax : changepoint.exe filename alpha beta\n");
		printf("Example: changepoint myfile 0.05 0.95\n");
		printf("         produces myfile.cp with type-I error (alpha) of 5%% and \n");
		printf("         confidence interval of 0.95\n");
		printf("BUG    : Please send emails to hawyang-at-princeton.edu\n\n");
		exit(1);
	}
	else if (argc == 4) {
		// Set the output filename
		strcpy(in_name, argv[1]);
		strcpy(out_name, argv[1]);
		alpha = strtod(argv[2], endptr);
		beta = strtod(argv[3], endptr);
	}

	// Load confidence interval (beta) critical region data
	if (beta == 0.99) {
		N_ca = N_ca99;
		ca = (double *) realloc(ca,N_ca*sizeof(double));
		for (ui=0; ui<N_ca; ui++)
			ca[ui] = ca99[ui];
		}
	else if (beta == 0.95) {
		N_ca = N_ca95;
		ca = (double *) realloc(ca,N_ca*sizeof(double));
		for (ui=0; ui<N_ca; ui++)
			ca[ui] = ca95[ui];
		}
	else if (beta == 0.9) {
		N_ca = N_ca90;
		ca = (double *) realloc(ca,N_ca*sizeof(double));
		for (ui=0; ui<N_ca; ui++)
			ca[ui] = ca90[ui];
		}
	else if (beta == 0.69) {
		N_ca = N_ca69;
		ca = (double *) realloc(ca,N_ca*sizeof(double));
		for (ui=0; ui<N_ca; ui++)
			ca[ui] = ca69[ui];
		}
	else {
		printf("Critical region threshold for confidence region of %.2f is not supplied.\n",
				beta);
		printf("Available are: 0.01, 0.05, 0.1, and 0.31.\n");
		exit(1);
	}

	// Initiate MPI Processes
	NA_BASE = N_ca;
	root_process = 0;

	ierr = MPI_Init(&argc, &argv);
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	// Read in time series
	if(my_id == root_process) {
		if ((fpin=fopen(in_name,"r")) == NULL) {
			// file does not exist, exit with error
			printf("File [%s] does not exist. Exiting ...\n", in_name);
			exit(1);
		}
		else {
			// File exists. start reading in data - only one set though
			traj = (double*) malloc(sizeof(double*));
			traj[0] = 0.0;

			while (!feof(fpin)) {
				
				fscanf(fpin, "%le", &tmp);
				L++;
				traj = (double *) realloc(traj,(L+1)*sizeof(double));
				traj[L] = tmp;
			}
			L--;
			fclose(fpin);
			printf("  Time Series is of length %lu.\n", (long unsigned int) L);
		}
		time0 = clock();

		// Average batch length
		AVG_L = L/num_procs;
		
		/******************************************************
		* Can calculate NA_Overlap in this section *
		******************************************************/
	
		// Send data to children processes
		for (an_id = 1; an_id < num_procs; an_id++) {
			// Calculate offset for each process
			offset = an_id*AVG_L - NA_OVERLAP;
			// Make sure analyze entire time series
			if(an_id == num_procs-1)
				AVG_L += (L % num_procs);
			// Send offset, length of time series segment, actual data in time series,
			ierr = MPI_Send(&offset, 1 , MPI_INT, an_id, send_data_tag, MPI_COMM_WORLD);
			ierr = MPI_Send(&AVG_L, 1 , MPI_INT, an_id, send_data_tag, MPI_COMM_WORLD);
			ierr = MPI_Send(&traj[offset], AVG_L + NA_OVERLAP, MPI_DOUBLE, an_id, 
				send_data_tag, MPI_COMM_WORLD);
		}

		// Multiple change point detection
		if (NA_BASE > AVG_L)
			cp1 = FindCP(&cp_root, traj, 0, AVG_L, alpha, beta, &Ncp, ca, AVG_L, 1);
		else {
			cpl = 0;
			cpr = NA_BASE;
			done = false;
			while (!done) {
				cp1 = FindCP(&cp_root, traj, cpl, cpr, alpha, beta, &Ncp, ca, AVG_L, 1);
				if (cpr < AVG_L) {
					if(cp1 == 0)  cpl = cpr - NA_OVERLAP;
					else  cpl = cp1;
					cpr = cpl + NA_BASE;
					if(cpr > AVG_L)  cpr = AVG_L;
					done = false;
				}
				else  done = true;
			}
		}
		// Check for spurious change points
		CheckCP(&cp_root, traj, alpha, beta, &Ncpdlt, trace, ca, AVG_L);

		// Store changepoints as array for each process
		cpArray = realloc(cpArray, num_procs * sizeof(int*));
		lbArray = realloc(lbArray, num_procs * sizeof(int*));
		rbArray = realloc(rbArray, num_procs * sizeof(int*));
		Ncp = 0;

		MakeCPArray(cp_root, &cpArray[0], &lbArray[0], &rbArray[0], &Ncp);
		sNcp = malloc(num_procs * sizeof(int));
		sNcp[0] = Ncp;

		// Receive detected change points from other processes
		for(i = 1; i < num_procs; i++) {

			ierr = MPI_Recv(&rNcp, 1, MPI_INT, MPI_ANY_SOURCE, send_data_tag, 
				MPI_COMM_WORLD, &status);
			sender = status.MPI_SOURCE;

			sNcp[sender] = rNcp - 1;
			Ncp += (rNcp - 1);

			cpArray[sender] = (int *) malloc(rNcp*sizeof(int));
			lbArray[sender] = (int *) malloc(rNcp*sizeof(int));
			rbArray[sender] = (int *) malloc(rNcp*sizeof(int));

			ierr = MPI_Recv(cpArray[sender], rNcp, MPI_INT, sender, send_data_tag, 
				MPI_COMM_WORLD, &status);
			ierr = MPI_Recv(lbArray[sender], rNcp, MPI_INT, sender, send_data_tag, 
				MPI_COMM_WORLD, &status);
			ierr = MPI_Recv(rbArray[sender], rNcp, MPI_INT, sender, send_data_tag, 
				MPI_COMM_WORLD, &status);
		}
		time1 = clock();
		printf("  %i change points found. [%.0f s]\n", Ncp, 
			(double)(time1-time0)/CLOCKS_PER_SEC);

		// Compare overlap regions for duplicate changepoints
		printf("] Pooling change points together from %d processes...\n", num_procs);
		time0 = clock();
		for(i = 1; i < num_procs; i++) {
			k = i - 1;
			offset = i*AVG_L - NA_OVERLAP;

			h = 1;
			j = 1;

			first = sNcp[k] - NA_OVERLAP;
			if(first < 1)
				first = 1;

			last = sNcp[k];
			middle = (first+last)/2;
			
			// Binary search for starting indices
			while (first <= last) {			
				if (cpArray[k][middle] < offset)
					first = middle + 1;    
				else if (cpArray[k][middle] == offset) {
					j = middle;
					break;
				}
				else
					last = middle - 1;
				
				middle = (first + last)/2;
			}
			if (first > last)
				j = last;

			if(j == 0)
				j = 1;

			// Iterate through change points detected in overlap region
			while(j <= sNcp[k] && h <= sNcp[i] && lbArray[i][h] <= offset+NA_OVERLAP) {

				// Take into consideration overlapping confidence intervals
				if(cpArray[i][h] <= rbArray[k][j] && cpArray[i][h] >= lbArray[k][j] ||
					cpArray[k][j] <= rbArray[i][h] && cpArray[k][j] >= lbArray[i][h]) {
					
					done = false;

					// Use confidence intervals as weights
					lb2 = (double) cpArray[i][h] - lbArray[i][h];
					lb1 = (double) cpArray[k][j] - lbArray[k][j];
					rb2 = (double) rbArray[i][h] - cpArray[i][h];
					rb1 = (double) rbArray[k][j] - cpArray[k][j];

					// Pool change points together
					if(lb1!=0 && lb2!=0 && rb1!=0 && rb2!=0) {
						cpArray[i][h] = (int)(((1.0/(lb1*lb1) + 1.0/(rb1*rb1)) * cpArray[k][j] + 
							(1.0/(lb2 * lb2) + 1.0/(rb2*rb2)) * cpArray[i][h])/(1.0/(lb1*lb1) + 
							1.0/(rb1*rb1) + 1.0/(lb2*lb2) + 1.0/(rb2*rb2)));
						lbArray[i][h] = (int)((1.0/(lb1*lb1) * lbArray[k][j] + 1.0/(lb2 * lb2) * 
							lbArray[i][h])/(1.0/(lb1*lb1) + 1.0/(lb2*lb2)));
						rbArray[i][h] = (int)((1.0/(rb1*rb1) * rbArray[k][j] + 1.0/(rb2 * rb2) * 
							rbArray[i][h])/(1.0/(rb1*rb1) + 1.0/(rb2*rb2)) + 0.5);
					}
					else {
						if(lb1==0 && lb2 == 0) {
							cpArray[i][h] = (int)((cpArray[k][j]+cpArray[i][h])/2.0);
							lbArray[i][h] = (int)((cpArray[k][j]+cpArray[i][h])/2.0);
						}
						else if(lb1==0) {
							cpArray[i][h] = cpArray[k][j];
							lbArray[i][h] = cpArray[k][j];
						}
						else if(lb1 != 0 && lb2 != 0) {
							lbArray[i][h] = (int)((1.0/(lb1*lb1) * lbArray[k][j] + 1.0/(lb2 * lb2) * 
								lbArray[i][h])/(1.0/(lb1*lb1) + 1.0/(lb2*lb2)));
						}	  	
						if(rb1==0 && rb2 == 0) {
							cpArray[i][h] = (int)((cpArray[k][j]+cpArray[i][h])/2.0+0.5);
							rbArray[i][h] = (int)((cpArray[k][j]+cpArray[i][h])/2.0+0.5);
						}
						else if(rb1 ==0) {
							cpArray[i][h] = cpArray[k][j];
							rbArray[i][h] = cpArray[k][j];
						}
						else if(rb1 != 0 && rb2 != 0) {
							rbArray[i][h] = (int)((1.0/(rb1*rb1) * rbArray[k][j] + 1.0/(rb2 * rb2) * 
								rbArray[i][h])/(1.0/(rb1*rb1) + 1.0/(rb2*rb2)) + 0.5);
						}  
					}

					// Merge change points. Delete duplicates in array k (less number of elements to copy)
					for(start = j + 1; start <= sNcp[k]; start++) {
						cpArray[k][start-1] = cpArray[k][start];
						lbArray[k][start-1] = lbArray[k][start];
						rbArray[k][start-1] = rbArray[k][start];
					}
					sNcp[k]--;
					h++;
					pCP++;
				}
				else if (cpArray[i][h] > cpArray[k][j]) {
					j++;
				}
				else {
					h++;
				}
			}
		}
		time1 = clock();
		printf("  %i change points pooled. [%.0f s]\n", pCP, 
			(double)(time1-time0)/CLOCKS_PER_SEC);

		// Save results...
		strcpy(filename, out_name);
		strcat(filename, ".cp");
		printf("  saving file: %s \n", filename);
		fpout = fopen(filename, "w");
		for(i = 0; i < num_procs; i++) {
			for (k = 1; k <= sNcp[i]; k++)
				fprintf(fpout, "%d %d %d\n", cpArray[i][k], lbArray[i][k], 
					rbArray[i][k]);
		}
		fclose(fpout);

		// Free Workspace
		free(traj);
		free(ca);
		for(i=0; i<num_procs; i++) {
			free(cpArray[i]);
			free(lbArray[i]);
			free(rbArray[i]);
		}
		free(cpArray);
		free(lbArray);
		free(rbArray);
		free(sNcp);

		ierr = MPI_Finalize();
	
	}

	// Find change points in child processes
	else {

		// Receive offset
		ierr = MPI_Recv(&offset, 1, MPI_INT, root_process, 
			send_data_tag, MPI_COMM_WORLD, &status);

		// Receive length of time series to be processed
		ierr = MPI_Recv(&AVG_L, 1, MPI_INT, root_process, send_data_tag, 
			MPI_COMM_WORLD, &status);

		// Store received time series data in split
		split = (double *) malloc((AVG_L + NA_OVERLAP) * sizeof(double));

		ierr = MPI_Recv(split, AVG_L+NA_OVERLAP, MPI_DOUBLE, root_process, 
			send_data_tag, MPI_COMM_WORLD, &status);

		// Multiple change point detection on data in split
		if (NA_BASE > AVG_L)
			cp1 = FindCP(&cp_root, split, 0, AVG_L, alpha, beta, &Ncp, ca, 
				AVG_L + NA_OVERLAP-1, 1);
		else {
			cpl = 0;
			cpr = NA_BASE;
			done = false;
			while (!done) {
				cp1 = FindCP(&cp_root, split, cpl, cpr, alpha, beta, &Ncp, ca, 
					AVG_L + NA_OVERLAP-1, 1);
				if (cpr < AVG_L) {
					if (cp1 == 0)
						cpl = cpr - NA_OVERLAP; // no change points found
					else
						cpl = cp1;
					cpr = cpl + NA_BASE;
					if (cpr > AVG_L)
						cpr = AVG_L;
					done = false;
				}
				else
					done = true;
			}
		}

		// Check change points sequentially and remove spurious ones 
		CheckCP(&cp_root, split, alpha, beta, &Ncpdlt, trace, 
			ca, AVG_L + NA_OVERLAP);
		Ncp = 0;

		// dummy value allocation
		if(Ncp <= 1) {
			cpArray[0] = malloc(sizeof(int));
			lbArray[0] = malloc(sizeof(int));
			rbArray[0] = malloc(sizeof(int));
		}

		// Add offsets to detected change points
		MakeCPArray(cp_root, cpArray, lbArray, rbArray, &Ncp);
		for (j = 1; j <= Ncp; j++) {
			(*cpArray)[j] = (*cpArray)[j] + offset;
			(*lbArray)[j] = (*lbArray)[j] + offset;
			(*rbArray)[j] = (*rbArray)[j] + offset;
		}
		Ncp++;

		// Send detected change points and confidence intervals to master process
		ierr = MPI_Send(&Ncp, 1, MPI_INT, root_process, send_data_tag, 
			MPI_COMM_WORLD);
		ierr = MPI_Send(*cpArray, Ncp, MPI_INT, root_process, send_data_tag, 
			MPI_COMM_WORLD);
		ierr = MPI_Send(*lbArray, Ncp, MPI_INT, root_process, send_data_tag, 
			MPI_COMM_WORLD);
		ierr = MPI_Send(*rbArray, Ncp, MPI_INT, root_process, send_data_tag, 
			MPI_COMM_WORLD);

		free(split);
		free(cpArray[0]);
		free(rbArray[0]);
		free(lbArray[0]);
		free(cpArray);
		free(rbArray);
		free(lbArray);

		// End parallelization
		ierr = MPI_Finalize();
	}
	return EXIT_SUCCESS;
}

