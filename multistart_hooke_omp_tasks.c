#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>

#define MAXVARS		(250)	/* max # of variables	     */
#define RHO_BEGIN	(0.5)	/* stepsize geometric shrink */
#define EPSMIN		(1E-6)	/* ending value of stepsize  */
#define IMAX		(5000)	/* max # of iterations	     */

/* global variables */
unsigned long funevals = 0;


/* Rosenbrocks classic parabolic valley ("banana") function */
double f(double *x, int n, unsigned long *local_func_evals)
{
    double fv;
    int i;

	(*local_func_evals)++;
    fv = 0.0;
    for (i=0; i<n-1; i++)   /* rosenbrock */
        fv = fv + 100.0*pow((x[i+1]-x[i]*x[i]),2) + pow((x[i]-1.0),2);

    return fv;
}

/* given a point, look for a better one nearby, one coord at a time */
double best_nearby(double delta[MAXVARS], double point[MAXVARS], double prevbest, int nvars,
					unsigned long *local_func_evals)
{
	double z[MAXVARS];
	double minf, ftmp;
	int i;
	minf = prevbest;
	for (i = 0; i < nvars; i++)
		z[i] = point[i];
	for (i = 0; i < nvars; i++) {
		z[i] = point[i] + delta[i];
		ftmp = f(z, nvars, local_func_evals);
		if (ftmp < minf)
			minf = ftmp;
		else {
			delta[i] = 0.0 - delta[i];
			z[i] = point[i] + delta[i];
			ftmp = f(z, nvars, local_func_evals);
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


int hooke(int nvars, double startpt[MAXVARS], double endpt[MAXVARS], double rho, double epsilon, int itermax,
			unsigned long *local_func_evals){
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
	fbefore = f(newx, nvars, local_func_evals);
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
		newf = best_nearby(delta, newx, fbefore, nvars, local_func_evals);
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
			newf = best_nearby(delta, newx, fbefore, nvars, local_func_evals);
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

int main(int argc, char *argv[])
{
	//The default number of threads is 8
	int thr = 8;
	//If there is an argument, use it to redefine the number of threads
	argc > 1 ? thr = atoi(argv[1]) : printf("Defaulting to %d threads!\n", thr);
	
	double startpt[MAXVARS], endpt[MAXVARS];
	int itermax = IMAX;
	double rho = RHO_BEGIN;
	double epsilon = EPSMIN;
	int nvars;
	int trial, ntrials;
    double fx;
	int i, jj;
	double t0, t1;

	double best_fx = 1e10;
	double best_pt[MAXVARS];
	int best_trial = -1;
	int best_jj = -1;

	for (i = 0; i < MAXVARS; i++) best_pt[i] = 0.0;

	ntrials = 128*1024;	/* number of trials */
	nvars = 16;		/* number of variables (problem dimension) */
    //initial seed will be used as part of the 3rd Xi (for erand48)
	long init_Seed = time(0);

	t0 = get_wtime();
    #pragma omp parallel private(startpt, endpt, jj, fx) num_threads(thr)
    {
		//a single thread will create all tasks (the rest of the threads will run them)
		#pragma omp single nowait
		{
			//only need to init single's buffer, since that is the only one available to the task (I think)
			unsigned short seedBuf[3];
        	seedBuf[0] = 0;
        	seedBuf[1] = 0;
        	seedBuf[2] = omp_get_thread_num() + init_Seed;
			unsigned long local_func_evals = 0;

			//printf("Single's seedBuf[2] = %d\n", seedBuf[2]);
			for (trial = 0; trial < ntrials; trial++) {
				/* starting guess for rosenbrock test function, search space in [-4, 4) */
				for (i = 0; i < nvars; i++) {
					startpt[i] = 4.0*erand48(seedBuf)-4.0;
				}
				#pragma omp task
				{
					//printf("seedBuf[2] = %d\t->\tthread: %d\n", seedBuf[2], omp_get_thread_num());
					jj = hooke(nvars, startpt, endpt, rho, epsilon, itermax, &local_func_evals);
					fx = f(endpt, nvars, &local_func_evals);

//The if statement needs to be evaluated for each iteration,
//thus I thought maybe it is a good idea not to allow every iteration of
//each thread to access the critical region, especially if (best_fx > fx) returns false.
//It is also important to note that the second if statement is necessary, since two threads, with values
//fx1 and fx2 (with fx1 < fx2) may simultaneously get through the first if statement.
//In the case that fx1 gets first into the critical section, if the if statement (in critical section)
//was not present, best_fx = fx2
					if (best_fx > fx){
						#pragma omp critical
						{
							if(best_fx > fx){
								best_trial = trial;
								best_jj = jj;
								best_fx = fx;
								for (i = 0; i < nvars; i++)
									best_pt[i] = endpt[i];
							}
						}
					}
					//lock once every task is completed. If I were to lock on every function evaluation it would be wastefull
					//integrate with critical? Is atomic worth locking and unlocking twice?
					#pragma omp atomic
					funevals += local_func_evals;
				}
			}
		}
		//taskwait is required, since I removed the implicit barrier at the end of the single region
		//could have chose to keep the implicit barrier and not to call taskwait
		#pragma omp taskwait
    }
	t1 = get_wtime();

	printf("\n\nFINAL RESULTS: threads=%d\n", thr);
	printf("Elapsed time = %.3lf s\n", t1-t0);
	printf("Total number of trials = %d\n", ntrials);
	printf("Total number of function evaluations = %ld\n", funevals);
	printf("Best result at trial %d used %d iterations, and returned\n", best_trial, best_jj);
	for (i = 0; i < nvars; i++) {
		printf("x[%3d] = %15.7le \n", i, best_pt[i]);
	}
	printf("f(x) = %15.7le\n", best_fx);

	//save the results in a txt file, so I may later plot them
	FILE *fp;
	fp=fopen("Results/res_omp_tasks.txt", "a");
	fprintf(fp, "%.3lf threads %d\n", t1-t0, thr);
	fclose(fp);
	return 0;
}