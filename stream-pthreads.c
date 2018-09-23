/*-----------------------------------------------------------------------*/
/* Program: Stream                                                       */
/* Revision: $Id: stream.c,v 5.9 2009/04/11 16:35:00 mccalpin Exp $ */
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
#define _GNU_SOURCE

# include <stdio.h>
# include <math.h>
# include <float.h>
# include <limits.h>
# include <sys/time.h>
# include <malloc.h>
# include <pthread.h>
# include <numa.h>

/* INSTRUCTIONS:
 *
 *	1) Stream requires a good bit of memory to run.  Adjust the
 *          value of 'N' (below) to give a 'timing calibration' of 
 *          at least 20 clock-ticks.  This will provide rate estimates
 *          that should be good to about 5% precision.
 */

#ifndef NUM_THREADS
#define NUM_THREADS 48
#endif

#ifndef MEM_OFF
#define MEM_OFF 0
#endif

#ifndef STRIDE
#define STRIDE (64/sizeof(double))
#endif


/*
node   0   1   2   3   4   5   6   7 
0:  10  16  16  22  16  22  16  22 
1:  16  10  16  22  22  16  22  16 
2:  16  16  10  16  16  16  16  22 
3:  22  22  16  10  16  16  22  16 
4:  16  22  16  16  10  16  16  16 
5:  22  16  16  16  16  10  22  22 
6:  16  22  16  22  16  22  10  16 
7:  22  16  22  16  16  22  16  10 
*/

#ifndef N
/* #   define N	16000000 */
#   define N	4800000
#endif
#ifndef NTIMES
#   define NTIMES	10
#endif
#define REPEAT (5*STRIDE)
#ifndef OFFSET
#   define OFFSET	0
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

static double	**a, **b, **c;

static double	avgtime[5] = {0}, maxtime[5] = {0},
    mintime[5] = {FLT_MAX,FLT_MAX,FLT_MAX,FLT_MAX,FLT_MAX};

static char	*label[5] = {"Read:      ","Copy:      ", "Scale:     ",
    "Add:       ", "Triad:     "};

static double	bytes[5] = {
    1 * sizeof(double) * N * NUM_THREADS * REPEAT / STRIDE,
    2 * sizeof(double) * N * NUM_THREADS * REPEAT / STRIDE,
    2 * sizeof(double) * N * NUM_THREADS * REPEAT / STRIDE,
    3 * sizeof(double) * N * NUM_THREADS * REPEAT / STRIDE,
    3 * sizeof(double) * N * NUM_THREADS * REPEAT / STRIDE
    };

extern double mysecond();
extern void checkSTREAMresults();

typedef struct Arg_T  {
    int proc;
    int allocNode;
} arg_t;

double results[NUM_THREADS];
    
void *readProc(void *arg) {
    arg_t *pArg = (arg_t*)arg;
    int me = pArg->proc;
    if (numa_run_on_node (me) == -1) {
        printf("unable to set affinity to processor %d\n", me);
    }

    double *a2 = a[me];
    double total = 0.0;

    for (int i = 0; i < REPEAT; i++) {
        for (int j = 0; j < N; j+=STRIDE) {
            total += a2[j];
        }
    }

    /* This construct is just to prevent the compiler from
       eliminating the reads */
    results[me] = total;
    return NULL;
}

void *copyProc(void *arg) {
    arg_t *pArg = (arg_t*)arg;
    int me = pArg->proc;
    if (numa_run_on_node (me) == -1) {
        printf("unable to set affinity to processor %d\n", me);
    }

    double *a2 = a[me];
    double *c2 = c[me];

    for (int i = 0; i < REPEAT; i++) {
        for (int j = 0; j < N; j+=STRIDE) {
            c2[j] = a2[j];
        }
    }
    return NULL;
}

void *scaleProc(void *arg) {
    arg_t *pArg = (arg_t*)arg;
    int me = pArg->proc;
    if (numa_run_on_node (me) == -1) {
        printf("unable to set affinity to processor %d\n", me);
    }

    double *b2 = b[me];
    double *c2 = c[me];

    for (int i = 0; i < REPEAT; i++) {
        for (int j = 0; j < N; j+=STRIDE) {
            b2[j] = 3.0*c2[j];
        }
    }
}

void *addProc(void *arg) {
    arg_t *pArg = (arg_t*)arg;
    int me = pArg->proc;
    if (numa_run_on_node (me) == -1) {
        printf("unable to set affinity to processor %d\n", me);
    }

    double *a2 = a[me];
    double *b2 = b[me];
    double *c2 = c[me];

    for (int i = 0; i < REPEAT; i++) {
        for (int j = 0; j < N; j+=STRIDE) {
            c2[j] = a2[j]+b2[j];
        }
    }
    return NULL;
}

void *triadProc(void *arg) {
    arg_t *pArg = (arg_t*)arg;
    int me = pArg->proc;
    if (numa_run_on_node (me) == -1) {
        printf("unable to set affinity to processor %d\n", me);
    }

    double *a2 = a[me];
    double *b2 = b[me];
    double *c2 = c[me];

    for (int i = 0; i < REPEAT; i++) {
        for (int j = 0; j < N; j+=STRIDE) {
            a2[j] = b2[j]+3.0*c2[j];
        }
    }
    return NULL;
}

int
main(int argc, char *argv[])
    {
    int			quantum, checktick();
    int			BytesPerWord;
    register int	j, k;
    double		t, times[5][NTIMES];

    /* --- SETUP --- determine precision and check timing --- */

    printf(HLINE);
    printf("STREAM version $Revision: 5.9 $\n");
    printf(HLINE);
    BytesPerWord = sizeof(double);
    printf("This system uses %d bytes per DOUBLE PRECISION word.\n",
	BytesPerWord);

    printf(HLINE);
#ifdef NO_LONG_LONG
    printf("Array size = %d, Offset = %d\n" , N, OFFSET);
#else
    printf("Array size = %llu, Offset = %d\n", (unsigned long long) N, OFFSET);
#endif

    printf("Total memory required = %.1f MB.\n",
	(3.0 * BytesPerWord) * ( (double) N / 1048576.0) * NUM_THREADS);
    printf("Each test is run %d times, but only\n", NTIMES);
    printf("the *best* time for each is used.\n");

    printf(HLINE);
    printf ("Number of Threads requested = %i\n", NUM_THREADS);

    printf(HLINE);

    if (numa_available() < 0) {
        printf ("System does not support NUMA.\n");
        exit (1);
    }

    int num_nodes = numa_max_node() + 1;
    printf ("Number of available nodes = %i\n", num_nodes);

    pthread_t* threads=(pthread_t *)malloc(NUM_THREADS*sizeof(pthread_t));
    pthread_attr_t pthread_custom_attr;
    pthread_attr_init(&pthread_custom_attr);
    
    arg_t *p=(arg_t*)malloc(sizeof(arg_t)*NUM_THREADS);
    for (int i=0; i<NUM_THREADS; i++)
    {
        p[i].proc=i%num_nodes;
        p[i].allocNode = ((i+MEM_OFF)%num_nodes);
    }
    
    /* Allocate memory for the threads */
    a = (double**)malloc(sizeof(double*)*NUM_THREADS);
    b = (double**)malloc(sizeof(double*)*NUM_THREADS);
    c = (double**)malloc(sizeof(double*)*NUM_THREADS);

    for (int i = 0; i < NUM_THREADS; i++) {
        a[i] = numa_alloc_onnode(sizeof(double)*N, p[i].allocNode);
        b[i] = numa_alloc_onnode(sizeof(double)*N, p[i].allocNode);
        c[i] = numa_alloc_onnode(sizeof(double)*N, p[i].allocNode);
        if ((a[i] == NULL) || (b[i] == NULL) || (c[i] == NULL)) {
            printf ("Failed to allocate %d bytes on node %d\n", (int)sizeof(double)*N, p[i].allocNode);
            exit (-1);
        }

        for (int j = 0; j < N; j++) {
            a[i][j] = 1.0;
            b[i][j] = 2.0;
            c[i][j] = 0.0;
        }
    }

    printf(HLINE);

    if  ( (quantum = checktick()) >= 1) 
	printf("Your clock granularity/precision appears to be "
	    "%d microseconds.\n", quantum);
    else {
	printf("Your clock granularity appears to be "
	    "less than one microsecond.\n");
	quantum = 1;
    }

    printf(HLINE);
    
    /*	--- MAIN LOOP --- repeat test cases NTIMES times --- */
    for (k=0; k<NTIMES; k++)
	{
        /* READ */
        times[0][k] = mysecond();
        for (int i=0; i<NUM_THREADS; i++)
        {
            pthread_create(&threads[i], &pthread_custom_attr, readProc, (void *)(p+i));
        }
        for (int i=0; i<NUM_THREADS; i++)
        {
            pthread_join(threads[i], NULL);
        }
        times[0][k] = mysecond() - times[0][k];

        printf (".");
        /* COPY */
        times[1][k] = mysecond();
        for (int i=0; i<NUM_THREADS; i++)
        {
            pthread_create(&threads[i], &pthread_custom_attr, copyProc, (void *)(p+i));
        }
        for (int i=0; i<NUM_THREADS; i++)
        {
            pthread_join(threads[i], NULL);
        }
        times[1][k] = mysecond() - times[1][k];

        /* SCALE */
        times[2][k] = mysecond();
        for (int i=0; i<NUM_THREADS; i++)
        {
            pthread_create(&threads[i], &pthread_custom_attr, scaleProc, (void *)(p+i));
        }
        for (int i=0; i<NUM_THREADS; i++)
        {
            pthread_join(threads[i], NULL);
        }
        times[2][k] = mysecond() - times[2][k];

        /* ADD */
        times[3][k] = mysecond();
        for (int i=0; i<NUM_THREADS; i++)
        {
            pthread_create(&threads[i], &pthread_custom_attr, addProc, (void *)(p+i));
        }
        for (int i=0; i<NUM_THREADS; i++)
        {
            pthread_join(threads[i], NULL);
        }
        times[3][k] = mysecond() - times[3][k];

        /* TRIAD */
        times[4][k] = mysecond();
        for (int i=0; i<NUM_THREADS; i++)
        {
            pthread_create(&threads[i], &pthread_custom_attr, triadProc, (void *)(p+i));
        }
        for (int i=0; i<NUM_THREADS; i++)
        {
            pthread_join(threads[i], NULL);
        }
        times[4][k] = mysecond() - times[4][k];
    }
    printf ("\n");


    
    /*	--- SUMMARY --- */

    for (k=1; k<NTIMES; k++) /* note -- skip first iteration */
	{
	for (j=0; j<5; j++)
	    {
	    avgtime[j] = avgtime[j] + times[j][k];
	    mintime[j] = MIN(mintime[j], times[j][k]);
	    maxtime[j] = MAX(maxtime[j], times[j][k]);
	    }
	}
    
    printf("Function      Rate (MB/s)   Latency(ns)   Avg time     Min time     Max time\n");
    for (j=0; j<5; j++) {
	avgtime[j] = avgtime[j]/(double)(NTIMES-1);

	printf("%s%11.4f  %11.4f  %11.4f  %11.4f  %11.4f\n", label[j],
	       1.0E-06 * bytes[j]/avgtime[j],
           (avgtime[j]/(bytes[j]/(sizeof(double)*NUM_THREADS)))*1.0E9,
	       avgtime[j],
	       mintime[j],
	       maxtime[j]);
    }
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
