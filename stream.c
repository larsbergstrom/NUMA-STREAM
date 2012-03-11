/*-----------------------------------------------------------------------*/
/* Program: Stream                                                       */
/* Revision: $Id: stream.c,v 5.9 2009/08/20 21:14:27 mccalpin Exp mccalpin $ */
/* Original code developed by John D. McCalpin                           */
/* Programmers: John D. McCalpin                                         */
/*              Joe R. Zagar                                             */
/*                                                                       */
/* This program measures memory transfer rates in MB/s for simple        */
/* computational kernels coded in C.                                     */
/*-----------------------------------------------------------------------*/
/* Copyright 1991-2005: John D. McCalpin                                 */
/*-----------------------------------------------------------------------*/
/* License:                                                              */
/*  1. You are free to use this program and/or to redistribute           */
/*     this program.                                                     */
/*  2. You are free to modify this program for your own use,             */
/*     including commercial use, subject to the publication              */
/*     restrictions in item 3.                                           */
/*  3. You are free to publish results obtained from running this        */
/*     program, or from works that you derive from this program,         */
/*     with the following limitations:                                   */
/*     3a. In order to be referred to as "STREAM benchmark results",     */
/*         published results must be in conformance to the STREAM        */
/*         Run Rules, (briefly reviewed below) published at              */
/*         http://www.cs.virginia.edu/stream/ref.html                    */
/*         and incorporated herein by reference.                         */
/*         As the copyright holder, John McCalpin retains the            */
/*         right to determine conformity with the Run Rules.             */
/*     3b. Results based on modified source code or on runs not in       */
/*         accordance with the STREAM Run Rules must be clearly          */
/*         labelled whenever they are published.  Examples of            */
/*         proper labelling include:                                     */
/*         "tuned STREAM benchmark results"                              */
/*         "based on a variant of the STREAM benchmark code"             */
/*         Other comparable, clear and reasonable labelling is           */
/*         acceptable.                                                   */
/*     3c. Submission of results to the STREAM benchmark web site        */
/*         is encouraged, but not required.                              */
/*  4. Use of this program or creation of derived works based on this    */
/*     program constitutes acceptance of these licensing restrictions.   */
/*  5. Absolutely no warranty is expressed or implied.                   */
/*-----------------------------------------------------------------------*/
/*
 * 2011 - Lars Bergstrom (larsberg@cs.uchicago.edu)
 *
 * This version of stream.c contains some small modifications to the original stream benchmark.
 *
 * First, you will need to install libnuma on your machine. This library is not suppored on
 * OSX, but is available via any package manager on Linux or can be locally installed from:
 * http://oss.sgi.com/projects/libnuma/
 *
 * Then, set the number of threads equal to the number of cores on your machine:
 * export OMP_NUM_THREADS=48
 *
 * Finally, compile using GCC:
 * gcc -O3 -std=c99 -fopenmp -lnuma -DN=80000000 -DNTIMES=100 stream.c -o stream-gcc
 *
 * Variants:
 * - To simulate non-NUMA aware access, define the constant: -DNON_NUMA
 * This constant will cause the memory to be touched by a thread on one node and then
 * for all of the stream operations to happen from another node.
 * - To show what happens when you perform only accesses that do not coincide with
 * a prior cache block, compile with -DSTRIDE=8 (assuming 64-byte cache lines and 8-byte
 * doubles).
 *
 * Original tech report:
 * Measuring NUMA effects with the STREAM benchmark. Lars Bergstrom.
 * University of Chicago Computer Science Department Technical Report TR-2011-02, May 2011.
 * http://www.cs.uchicago.edu/research/publications/techreports/TR-2011-02
 *
 */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <float.h>
# include <limits.h>
# include <sys/time.h>
# include <numa.h>

/* INSTRUCTIONS:
 *
 *	1) Stream requires a good bit of memory to run.  Adjust the
 *          value of 'N' (below) to give a 'timing calibration' of 
 *          at least 20 clock-ticks.  This will provide rate estimates
 *          that should be good to about 5% precision.
 */

#ifndef N
#   define N   20000000
#endif
#ifndef NTIMES
#   define NTIMES   10
#endif
#ifndef STRIDE
#define STRIDE 1
#endif


/*
 *	3) Compile the code with full optimization.  Many compilers
 *	   generate unreasonably bad code before the optimizer tightens
 *	   things up.  If the results are unreasonably good, on the
 *	   other hand, the optimizer might be too smart for me!
 *
 *         Try compiling with:
 *               cc -O stream_omp.c -o stream_omp
 *
 *         This is known to work on Cray, SGI, IBM, and Sun machines.
 *
 *
 *	4) Mail the results to mccalpin@cs.virginia.edu
 *	   Be sure to include:
 *		a) computer hardware model number and software revision
 *		b) the compiler flags
 *		c) all of the output from the test case.
 * Thanks!
 *
 */

# define HLINE "-------------------------------------------------------------\n"

# ifndef MIN
# define MIN(x,y) ((x)<(y)?(x):(y))
# endif
# ifndef MAX
# define MAX(x,y) ((x)>(y)?(x):(y))
# endif

double *a, *b, *c;

static double	avgtime[4] = {0}, maxtime[4] = {0},
		mintime[4] = {FLT_MAX,FLT_MAX,FLT_MAX,FLT_MAX};

static char	*label[4] = {"Copy:      ", "Scale:     ",
    "Add:       ", "Triad:     "};

static double	bytes[4] = {
    2 * sizeof(double) * N / STRIDE,
    2 * sizeof(double) * N / STRIDE,
    3 * sizeof(double) * N / STRIDE,
    3 * sizeof(double) * N / STRIDE
    };

extern double mysecond();
extern void checkSTREAMresults (long length, double *a, double *b, double *c);
#ifdef _OPENMP
extern int omp_get_num_threads();
#endif
int
main()
    {
    int			quantum, checktick();
    int			BytesPerWord;

    register int	i, j, k;
    double		scalar, t, times[4][NTIMES];
    
/* --- SETUP --- determine precision and check timing --- */

    printf(HLINE);
    printf("STREAM version $Revision: 5.9 $\n");
    printf(HLINE);
    BytesPerWord = sizeof(double);
    printf("This system uses %d bytes per DOUBLE PRECISION word.\n",
	BytesPerWord);

    if (numa_available() < 0) {
        printf ("System does not support NUMA.\n");
        exit (1);
    }

	a = malloc(sizeof(double)*(N));
	b = malloc(sizeof(double)*(N));
	c = malloc(sizeof(double)*(N));

	if (a==0 | b==0 | c==0) {
		printf("ERROR: one of the malloc's failed\n");
		printf("     a = %p, b = %p, c = %p\n", a,b,c);
		exit(1);
	}

    printf(HLINE);
    printf("Array size = %d\n" , N);
    printf("Total memory required = %.1f MB.\n",
	(3.0 * BytesPerWord) * ( (double) N / 1048576.0));
    printf("Each test is run %d times, but only\n", NTIMES);
    printf("the *best* time for each is used.\n");

#ifdef _OPENMP
    printf(HLINE);
#pragma omp parallel 
    {
#pragma omp master
	{
	    k = omp_get_num_threads();
	    printf ("Number of Threads requested = %i\n",k);


        }
    }
#endif

    int num_nodes = numa_max_node()+1;
    printf ("Number of available nodes = %i\n", num_nodes);

#pragma omp parallel for private(i)
	for (j=0; j<k; j++) {
		i = j%num_nodes;
        numa_run_on_node(i);
	}


    /* Get initial value for system clock. */
#pragma omp parallel for
    for (j=0; j<N; j++) {
		a[j] = 1.0;
		b[j] = 2.0;
		c[j] = 0.0;
	}

#ifdef NON_NUMA
#pragma omp parallel for private(i)
	for (j=0; j<k; j++) {
		i = (j+1)%num_nodes;
        numa_run_on_node(i);
	}
    printf("Execution will be non-NUMA aware.\n");
#endif
    printf(HLINE);

    if  ( (quantum = checktick()) >= 1) 
	printf("Your clock granularity/precision appears to be "
	    "%d microseconds.\n", quantum);
    else {
	printf("Your clock granularity appears to be "
	    "less than one microsecond.\n");
	quantum = 1;
    }

    t = mysecond();
#pragma omp parallel for
    for (j = 0; j < N; j++)
		a[j] = 2.0E0 * a[j];
    t = 1.0E6 * (mysecond() - t);

    printf("Each test below will take on the order"
	" of %d microseconds.\n", (int) t  );
    printf("   (= %d clock ticks)\n", (int) (t/quantum) );
    printf("Increase the size of the arrays if this shows that\n");
    printf("you are not getting at least 20 clock ticks per test.\n");

    printf(HLINE);

    printf("WARNING -- The above is only a rough guideline.\n");
    printf("For best results, please be sure you know the\n");
    printf("precision of your system timer.\n");
    printf(HLINE);
    
    /*	--- MAIN LOOP --- repeat test cases NTIMES times --- */

    scalar = 3.0;
    for (k=0; k<NTIMES; k++)
	{
	times[0][k] = mysecond();
#pragma omp parallel for
	for (j=0; j<N; j+=STRIDE)
	    c[j] = a[j];
	times[0][k] = mysecond() - times[0][k];
	
	times[1][k] = mysecond();
#pragma omp parallel for
	for (j=0; j<N; j+=STRIDE)
	    b[j] = scalar*c[j];
	times[1][k] = mysecond() - times[1][k];
	
	times[2][k] = mysecond();
#pragma omp parallel for
	for (j=0; j<N; j+=STRIDE)
	    c[j] = a[j]+b[j];
	times[2][k] = mysecond() - times[2][k];
	
	times[3][k] = mysecond();
#pragma omp parallel for
	for (j=0; j<N; j+=STRIDE)
	    a[j] = b[j]+scalar*c[j];
	times[3][k] = mysecond() - times[3][k];
	}

    /*	--- SUMMARY --- */

    for (k=1; k<NTIMES; k++) /* note -- skip first iteration */
	{
	for (j=0; j<4; j++)
	    {
	    avgtime[j] = avgtime[j] + times[j][k];
	    mintime[j] = MIN(mintime[j], times[j][k]);
	    maxtime[j] = MAX(maxtime[j], times[j][k]);
	    }
	}
    
    printf("Function      Rate (MB/s)   Avg time     Min time     Max time\n");
    for (j=0; j<4; j++) {
	avgtime[j] = avgtime[j]/(double)(NTIMES-1);

	printf("%s%11.4f  %11.4f  %11.4f  %11.4f\n", label[j],
	       1.0E-06 * bytes[j]/mintime[j],
	       avgtime[j],
	       mintime[j],
	       maxtime[j]);
    }
    printf(HLINE);

    /* --- Check Results --- */
    checkSTREAMresults((long) N, a, b, c);
    printf(HLINE);

    return 0;
}

# define	M	20

int
checktick()
    {
    int		i, minDelta, Delta;
    double	t1, t2, timesfound[M];

/*  Collect a sequence of M unique time values from the system. */

    for (i = 0; i < M; i++) {
	t1 = mysecond();
	while( ((t2=mysecond()) - t1) < 1.0E-6 )
	    ;
	timesfound[i] = t1 = t2;
	}

/*
 * Determine the minimum difference between these M values.
 * This result will be our estimate (in microseconds) for the
 * clock granularity.
 */

    minDelta = 1000000;
    for (i = 1; i < M; i++) {
	Delta = (int)( 1.0E6 * (timesfound[i]-timesfound[i-1]));
	minDelta = MIN(minDelta, MAX(Delta,0));
	}

   return(minDelta);
    }



/* A gettimeofday routine to give access to the wall
   clock timer on most UNIX-like systems.  */

#include <sys/time.h>

double mysecond()
{
        struct timeval tp;
        int i;

        i = gettimeofday(&tp,NULL);
        return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}


/* ----------------------------------------------- 
   Check the results to make sure all the loops
   have actually been run.   
   This revised version (in 5.9 and above) sums the
   absolute errors across the arrays, rather than 
   summing the values in the arrays and comparing
   with the expected sum.  This version is much 
   less sensitive to accumulation of roundoff error.
-------------------------------------------------- */
void checkSTREAMresults (long length, double *a, double *b, double *c)
{
	double aj,bj,cj,scalar;
	double asum,bsum,csum;
	double epsilon;
	int	j,k,fail=0;

    /* reproduce initialization */
	aj = 1.0;
	bj = 2.0;
	cj = 0.0;
    /* a[] is modified during timing check */
	aj = 2.0E0 * aj;
    /* now execute timing loop */
	scalar = 3.0;
	for (k=0; k<NTIMES; k++)
        {
            cj = aj;
            bj = scalar*cj;
            cj = aj+bj;
            aj = bj+scalar*cj;
        }
    /* now aj, bj, and cj have values that should match each element */
    /* of arrays a[], b[], and c[] -- unless I modified the code to */
    /* fiddle with some entries to confuse optimizers -- watch for this */

#ifdef VERBOSE
	printf ("Comparison of specific values at midpoint of arrays: \n");
	printf ("        Expected  : %f %f %f \n",aj,bj,cj);
	printf ("        Observed  : %f %f %f \n",a[length/2],b[length/2],c[length/2]);
#endif

#ifndef abs
#define abs(a) ((a) >= 0 ? (a) : -(a))
#endif
	asum = 0.0;
	bsum = 0.0;
	csum = 0.0;
	for (j=0; j<length; j+=STRIDE) {
		asum += abs(a[j] - aj);
		bsum += abs(b[j] - bj);
		csum += abs(c[j] - cj);
	}
	asum = asum / (double) (length);
	csum = bsum / (double) (length);
	csum = csum / (double) (length);
#ifdef VERBOSE
	printf ("Average Absolute Error : \n");
	printf ("    arrays: a, b, c  : %f %f %f \n",asum,bsum,csum);
#endif

	epsilon = 1.e-8;

	if (asum > epsilon) {
		printf ("Failed Validation on array a[]\n");
		printf ("        Max Allowable Error  : %f \n",epsilon);
		printf ("        Observed Error       : %f \n",asum);
		fail = 1;
	}
	if (bsum > epsilon) {
		printf ("Failed Validation on array b[]\n");
		printf ("        Max Allowable Error  : %f \n",epsilon);
		printf ("        Observed Error       : %f \n",bsum);
		fail = 1;
	}
	if (csum > epsilon) {
		printf ("Failed Validation on array c[]\n");
		printf ("        Max Allowable Error  : %f \n",epsilon);
		printf ("        Observed Error       : %f \n",csum);
		fail = 1;
	}
	if (fail == 0) {
		printf ("Solution Validates\n");
	}
}
