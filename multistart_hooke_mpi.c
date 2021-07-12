#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <mpi.h>

#define MAXVARS		(250)	/* max # of variables	     */
#define RHO_BEGIN	(0.5)	/* stepsize geometric shrink */
#define EPSMIN		(1E-6)	/* ending value of stepsize  */
#define IMAX		(5000)	/* max # of iterations	     */

/* global variables */
unsigned long funevals = 0;


/* Rosenbrocks classic parabolic valley ("banana") function */
double f(double *x, int n)
{
    double fv;
    int i;

	funevals++;
    fv = 0.0;
    for (i=0; i<n-1; i++)   /* rosenbrock */
        fv = fv + 100.0*pow((x[i+1]-x[i]*x[i]),2) + pow((x[i]-1.0),2);

    return fv;
}

/* given a point, look for a better one nearby, one coord at a time */
double best_nearby(double delta[MAXVARS], double point[MAXVARS], double prevbest, int nvars)
{
	double z[MAXVARS];
	double minf, ftmp;
	int i;
	minf = prevbest;
	for (i = 0; i < nvars; i++)
		z[i] = point[i];
	for (i = 0; i < nvars; i++) {
		z[i] = point[i] + delta[i];
		ftmp = f(z, nvars);
		if (ftmp < minf)
			minf = ftmp;
		else {
			delta[i] = 0.0 - delta[i];
			z[i] = point[i] + delta[i];
			ftmp = f(z, nvars);
			if (ftmp < minf)
				minf = ftmp;
			else
				z[i] = point[i];
		}
	}
	for (i = 0; i < nvars; i++)
		point[i] = z[i];

	return (minf);
}


int hooke(int nvars, double startpt[MAXVARS], double endpt[MAXVARS], double rho, double epsilon, int itermax)
{
	double delta[MAXVARS];
	double newf, fbefore, steplength, tmp;
	double xbefore[MAXVARS], newx[MAXVARS];
	int i, j, keep;
	int iters, iadj;

	for (i = 0; i < nvars; i++) {
		newx[i] = xbefore[i] = startpt[i];
		delta[i] = fabs(startpt[i] * rho);
		if (delta[i] == 0.0)
			delta[i] = rho;
	}
	iadj = 0;
	steplength = rho;
	iters = 0;
	fbefore = f(newx, nvars);
	newf = fbefore;
	while ((iters < itermax) && (steplength > epsilon)) {
		iters++;
		iadj++;
#if DEBUG
		printf("\nAfter %5d funevals, f(x) =  %.4le at\n", funevals, fbefore);
		for (j = 0; j < nvars; j++)
			printf("   x[%2d] = %.4le\n", j, xbefore[j]);
#endif
		/* find best new point, one coord at a time */
		for (i = 0; i < nvars; i++) {
			newx[i] = xbefore[i];
		}
		newf = best_nearby(delta, newx, fbefore, nvars);
		/* if we made some improvements, pursue that direction */
		keep = 1;
		while ((newf < fbefore) && (keep == 1)) {
			iadj = 0;
			for (i = 0; i < nvars; i++) {
				/* firstly, arrange the sign of delta[] */
				if (newx[i] <= xbefore[i])
					delta[i] = 0.0 - fabs(delta[i]);
				else
					delta[i] = fabs(delta[i]);
				/* now, move further in this direction */
				tmp = xbefore[i];
				xbefore[i] = newx[i];
				newx[i] = newx[i] + newx[i] - tmp;
			}
			fbefore = newf;
			newf = best_nearby(delta, newx, fbefore, nvars);
			/* if the further (optimistic) move was bad.... */
			if (newf >= fbefore)
				break;

			/* make sure that the differences between the new */
			/* and the old points are due to actual */
			/* displacements; beware of roundoff errors that */
			/* might cause newf < fbefore */
			keep = 0;
			for (i = 0; i < nvars; i++) {
				keep = 1;
				if (fabs(newx[i] - xbefore[i]) > (0.5 * fabs(delta[i])))
					break;
				else
					keep = 0;
			}
		}
		if ((steplength >= epsilon) && (newf >= fbefore)) {
			steplength = steplength * rho;
			for (i = 0; i < nvars; i++) {
				delta[i] *= rho;
			}
		}
	}
	for (i = 0; i < nvars; i++)
		endpt[i] = xbefore[i];

	return (iters);
}


double get_wtime(void)
{
    struct timeval t;

    gettimeofday(&t, NULL);

    return (double)t.tv_sec + (double)t.tv_usec*1.0e-6;
}


struct local_min{
	double local_fx;
	double local_pt[MAXVARS];
	int local_trial;
	int local_jj;
	unsigned long local_funevals;
};


int main(int argc, char *argv[])
{
	//Initialize MPI
	MPI_Init(&argc, &argv);

	//Get the essentials
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	int rank ;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//All these are local to each MPI_Process
	double startpt[MAXVARS], endpt[MAXVARS];
	int itermax = IMAX;
	double rho = RHO_BEGIN;
	double epsilon = EPSMIN;
	//counters and temporary variables
	int nvars;
	int trial, ntrials;
	double fx;
	int i, j, jj;
	double t0, t1;

	//these need to be reduced
	struct local_min best;
	struct local_min *gather_best;
	best.local_fx = 1e10;
	best.local_trial = -1;
	best.local_jj = -1;
	best.local_funevals = 0;
	for (i = 0; i < MAXVARS; i++) best.local_pt[i] = 0.0;
	
	ntrials = 1024*128;	/* number of trials */
	nvars = 16;		/* number of variables (problem dimension) */

	//Get the starting and ending iteration of the current MPI process
	int to_exec = ntrials/size;
	int start_trial = rank*(to_exec);
	int end_trial = start_trial + to_exec - 1;
	if(rank == size - 1){
		end_trial = ntrials;
	}

	//erand buffer
	long init_seed = time(0); 
	unsigned short ebuf[3];
	ebuf[0] = 0;
	ebuf[1] = 0;
	ebuf[2] = rank + init_seed;

	if (rank==0) {t0 = get_wtime();}
	for (trial = start_trial; trial < end_trial; trial++) {
		/* starting guess for rosenbrock test function, search space in [-4, 4) */
		for (i = 0; i < nvars; i++) {
			startpt[i] = 4.0*erand48(ebuf)-4.0;
		}
		jj = hooke(nvars, startpt, endpt, rho, epsilon, itermax);
		fx = f(endpt, nvars);

		//update local minimum
		if (fx < best.local_fx) {
			best.local_trial = trial;
			best.local_jj = jj;
			best.local_fx = fx;
			for (i = 0; i < nvars; i++)
				best.local_pt[i] = endpt[i];
		}
	}
	best.local_funevals = funevals;
	//Zero gather all local minimums
	if (rank == 0){
		gather_best = (struct local_min *)malloc(sizeof(struct local_min)*size);
	}
	MPI_Gather(&best, sizeof(struct local_min), MPI_BYTE, gather_best, sizeof(struct local_min), MPI_BYTE, 0, MPI_COMM_WORLD);
	if (rank == 0){
		//This information is in the array that it gathered and thus I need to zero out the counter so
		//I get the correct amount of iterations that were executed
		best.local_funevals = 0;
		for(i = 0; i < size; i++){
			//Add all local_funevals together (could be implemented with a MPI_Reduce instead)
			best.local_funevals += gather_best[i].local_funevals;
			
			//set as the root minimum, the minimum fx of all MPI_Processes
			if(best.local_fx > gather_best[i].local_fx){
				best.local_trial = gather_best[i].local_trial;
				best.local_jj = gather_best[i].local_jj;
				best.local_fx = gather_best[i].local_fx;
				for (j = 0; j < nvars; j++)
					best.local_pt[j] = gather_best[i].local_pt[j];
			}
		}
	}
	if (rank==0) {t1 = get_wtime();}
	//only root prints
	if (rank == 0){
		printf("\n\nFINAL RESULTS:\n");
		printf("Elapsed time = %.3lf s\n", t1-t0);
		printf("Total number of trials = %d\n", ntrials);
		printf("Total number of function evaluations = %ld\n", best.local_funevals);
		printf("Best result at trial %d used %d iterations, and returned\n", best.local_trial, best.local_jj);
		for (j = 0; j < nvars; j++) {
			printf("x[%3d] = %15.7le \n", j, best.local_pt[j]);
		}
		printf("f(x) = %15.7le\n", best.local_fx);


		//save the results in a txt file, so I may later plot them
		FILE *fp;
		fp=fopen("Results/res_mpi.txt", "a");
		fprintf(fp, "%.3lf size %d\n", t1-t0, size);
		fclose(fp);

		free(gather_best);
	}

	//Finalize MPI
	MPI_Finalize();

	return 0;
}



//IDEA:
//In case I want to minimize the amount of messages being sent over the network,
//in order to minimize traffic (will require one extra RTT (Round Trip Time) though)
//Based on:
//Each process (root process included) sends the contents of its send buffer to the root process.
//The root process receives the messages and stores them in rank order.
//ref: https://www.open-mpi.org/doc/v3.1/man3/MPI_Gather.3.php

//Initialize int flag = 0;
//Gather only best_fx
//Locate the minimum and the rank that sent it
//Ssend a message (ie an integer flag, setting it to 1) only to that one task (request for more information (best_pt[] etc.))
//The rest of the ranks do an Irecv and block temporarily on a Barrier. Once The message is received root rank will join the rest
//on the barrier ("free them")
//Now the only rank that has flag==1 is the one with the best_fx.
//Using an if statement send the rest of the info to the root rank.
//Root rank should expect a message after the barrier.

//If funevals was not a thing then I could also use that last if statement in order to print the result, insead of
//sending everything to the root.