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
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <string.h>
#include <unistd.h>
#include "KUtils.h"
//#include <sys/time.h>


/* INSTRUCTIONS:
 *
 *	1) Stream requires a good bit of memory to run.  Adjust the
 *          value of 'N' (below) to give a 'timing calibration' of
 *          at least 20 clock-ticks.  This will provide rate estimates
 *          that should be good to about 5% precision.
 */

//#define TUNED

#ifndef N
#   define N 600000
#endif
#pragma message ("Size of array(bytes)=" STRINGIFY(N))
/*#ifndef NTIMES
#   define NTIMES	10
#endif*/
int NTIMES = 0;

#ifndef OFFSET
#   define OFFSET	0
#endif
#ifndef CLOCKS_PER_SEC
//#   define CLOCKS_PER_SEC 1000000000UL
#   define CLOCKS_PER_SEC 1UL
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

clock_t clock(void);      //forward declaration of sting function

static double	a[N+OFFSET],
		b[N+OFFSET],
		c[N+OFFSET];

static double	avgtime[4] = {0, 0, 0, 0}, maxtime[4] = {0, 0, 0, 0},
		mintime[4] = {FLT_MAX,FLT_MAX,FLT_MAX,FLT_MAX};

static char	*label[4] = {"Copy:      ", "Scale:     ",
    "Add:       ", "Triad:     "};

static double	bytes[4] = {
    2 * sizeof(double) * N,
    2 * sizeof(double) * N,
    3 * sizeof(double) * N,
    3 * sizeof(double) * N
    };

//extern double mysecond();
extern double myclock();
extern void checkSTREAMresults();
#ifdef TUNED
extern void tuned_STREAM_Copy();
extern void tuned_STREAM_Scale(double scalar);
extern void tuned_STREAM_Add();
extern void tuned_STREAM_Triad(double scalar);
#endif
#ifdef _OPENMP
extern int omp_get_num_threads();
#endif

// Test Configuration
typedef int bool;
#define false 0
#define true  1

enum
{
    g_testAll   = 0,
    g_testCopy  = 1,
    g_testScale = 2,
    g_testAdd   = 3,
    g_testTriad = 4,
    g_testMax   = 5
};

typedef struct TestOpts_ {
   const char*   name;
}TestOpts;

void
streamv1_usage(int argc, char *argv[], char* usage)
{
  printf("Usage: %s %s", argv[0], usage);
  exit(-1);
}

//
// Convert a human readable string into an integer.
// Examples:
//   easy_size("1")    --> 1
//   easy_size("100")  --> 100
//   easy_size("1k")   --> 1024
//   easy_size("100k") --> 102400
//   easy_size("1m")   --> 1048576
//   easy_size("1g")   --> 1073741824
//
unsigned long long easy_size(const char *s)
{
  char *end;
  unsigned long long val;

  val = strtoll(s, &end, 0);

  if (end == s)
  {
    // No characters were accepted.
    printf("Bad size: %s\n", s);
    //streamv1_usage();
  }

  switch (*end)
  {
  case 'g':
  case 'G':
    val *= 1024;
    /* fall through to next case */

  case 'm':
  case 'M':
    val *= 1024;
    /* fall through to next case */

  case 'k':
  case 'K':
    val *= 1024;
    ++end;
    /* fall through to next case */

  case '\0':
    break;

  default:
    // Unrecognized character
    printf("Bad size: %s\n", s);
    //streamv1_usage();
  }

  if (*end != '\0')
  {
    // Extra characters at the end.
    printf("Bad size: %s\n", s);
    //streamv1_usage();
  }

  return val;
}

int
main(int ac, char **av)
{
#ifdef RT_SCHED
    struct sched_param sp;
    int ret;

    sp.sched_priority = sched_get_priority_max(SCHED_FIFO);
    ret=sched_setscheduler(0, SCHED_FIFO, &sp);
    if(ret==-1){
	perror("sched_setscheduler 1 failed");
        exit(1);
    }
#endif

    int			quantum, checktick();
    int			BytesPerWord;
    register int	i, j, k;
    double		scalar, t;
    double              **times;
    int                 op;

    /* --- SETUP --- determine precision and check timing --- */

    printf(HLINE);
    printf("STREAM version $Revision: 5.9 $\n");
    printf(HLINE);
    BytesPerWord = sizeof(double);
    printf("This system uses %d bytes per DOUBLE PRECISION word.\n",
	BytesPerWord);

    printf(HLINE);
    const char   *test = "all";
    bool isEnableAngelSignals = false;
    char   *usage = "[-T <copy/scale/add/triad/all>] \
                    [-R <repetitions>] \
                    [-A <enable(1)/disable angel signals(0)>] \n";
    TestOpts g_testopts[] =
    {
      { "all" },
      { "copy" },
      { "scale" },
      { "add" },
      { "triad" },
    };
    bool  algo[g_testMax] = {false};

    while (( op = getopt(ac, av, "T:R:A:")) != EOF)
    {
        switch(op)
        {
        case 'T':
                test = optarg;
                char *opt = (char *)&test[0];
                bool found = false;

                for(i = 0; i < sizeof(g_testopts) / sizeof(TestOpts); i++)
                {
                  if (strcmp(opt, g_testopts[i].name) == 0)
                  {
                    algo[i] = true;
                    printf( "Testing %s\n", g_testopts[i].name);
                    found = true;
                    break;
                  }
                }

                if (!found)
                {
                  printf( "Invalid option: -T=%s\n", optarg);
                  printf( "Valid tests: -T=%s", g_testopts[0].name);
                  for(i = 1; i < sizeof(g_testopts) / sizeof(TestOpts); i++)
                  {
                    printf( ",%s", g_testopts[i].name);
                  }
                  printf( "\n");
                  exit(0);
                }
                break;
        case 'A':
                isEnableAngelSignals = atoi(optarg);
                break;
        /*case 'M':
                N = easy_size(optarg);
                break;*/
        case 'R':
                NTIMES = atoi(optarg);
                break;
        default:
                streamv1_usage(ac, av, usage);
                break;
        }
    }

   //Allocate and initialize times array
   times = malloc(sizeof(double) * 4);
   for(j=0; j<4; j++)
   {
      times[j] = malloc( sizeof(times[j][0]) * NTIMES );
   }
   for (k=0; k<NTIMES; k++)
   {
       for (j=0; j<4; j++)
       {
         times[j][k] = 0;
       }
   }

#ifdef NO_LONG_LONG
    printf("Array size = %d, Offset = %d\n" , N, OFFSET);
#else
    printf("Array size = %llu, Offset = %d\n", (unsigned long long) N, OFFSET);
#endif

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

    printf(HLINE);
#pragma omp parallel
    {
    printf ("Printing one line per active thread....\n");
    }

    /* Get initial value for system clock. */
#pragma omp parallel for
    for (j=0; j<N; j++)
    {
	a[j] = 1.0;
	b[j] = 2.0;
	c[j] = 0.0;
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

    t = myclock();
#pragma omp parallel for
    for (j = 0; j < N; j++)
    {
	a[j] = 2.0E0 * a[j];
    }
    t = 1.0E6 * (myclock() - t);

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
          if( algo[g_testAll] || algo[g_testCopy] )
          {
#ifdef ANGEL
            if(isEnableAngelSignals && (k==0))
            {
              workload_ckpt_begin();
            }
#endif /* ANGEL  */
	    times[0][k] = myclock();
#ifdef TUNED
            tuned_STREAM_Copy();
#else
#pragma omp parallel for
	    for (j=0; j<N; j++)
	      c[j] = a[j];
#endif
	     times[0][k] = myclock() - times[0][k];
	     //printf("Copy time[%d][%d]=%8.8f and N=%d\n", 0, k, times[0][k], N);
          }

          if( algo[g_testAll] || algo[g_testScale] )
          {
#ifdef ANGEL
            if(isEnableAngelSignals && (k==0) && !algo[g_testAll])
            {
              workload_ckpt_begin();
            }
#endif /* ANGEL  */
	    times[1][k] = myclock();
#ifdef TUNED
            tuned_STREAM_Scale(scalar);
#else
#pragma omp parallel for
	    for (j=0; j<N; j++)
	      b[j] = scalar*c[j];
#endif
	    times[1][k] = myclock() - times[1][k];
          }

          if( algo[g_testAll] || algo[g_testAdd] )
          {
#ifdef ANGEL
            if(isEnableAngelSignals && (k==0) && !algo[g_testAll])
            {
              workload_ckpt_begin();
            }
#endif /* ANGEL  */
	    times[2][k] = myclock();
#ifdef TUNED
            tuned_STREAM_Add();
#else
#pragma omp parallel for
	   for (j=0; j<N; j++)
	     c[j] = a[j]+b[j];
#endif
	   times[2][k] = myclock() - times[2][k];
          }

	  if( algo[g_testAll] || algo[g_testTriad] )
          {
#ifdef ANGEL
            if(isEnableAngelSignals && (k==0) && !algo[g_testAll])
            {
              workload_ckpt_begin();
            }
#endif /* ANGEL  */
	    times[3][k] = myclock();
#ifdef TUNED
            tuned_STREAM_Triad(scalar);
#else
#pragma omp parallel for
	    for (j=0; j<N; j++)
	      a[j] = b[j]+scalar*c[j];
#endif
	    times[3][k] = myclock() - times[3][k];
          }
    } /*End of NTIMES loop*/

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

    printf("Function      Bytes        AvgClks        MinClks       MaxClks\n");
    for (j=0; j<4; j++)
    {
        avgtime[j] = avgtime[j]/(double)(NTIMES-1);

        printf("%s%11.4f  %11.11f  %11.11f  %11.11f\n", label[j],
               bytes[j],
               avgtime[j],
               mintime[j],
               maxtime[j]);
    }
    printf(HLINE);

#ifdef ANGEL
  if (isEnableAngelSignals)
  {
      workload_ckpt_end();
  }
#endif /* ANGEL  */

#ifndef ANGEL
    if (algo[g_testAll])
    {
      /* --- Check Results --- */
      checkSTREAMresults();
      printf(HLINE);
    }
#endif

    return 0;
}

# define	M	20

int checktick()
{
    int		i, minDelta, Delta;
    double	t1, t2, timesfound[M];

/*  Collect a sequence of M unique time values from the system. */

/* Making assumption of 500 MHz */

    for (i = 0; i < M; i++)
    {
	t1 = myclock();
	while( ((t2=myclock()) - t1) < 1.0E-6 );
	timesfound[i] = t1 = t2;
    }

/*
 * Determine the minimum difference between these M values.
 * This result will be our estimate (in microseconds) for the
 * clock granularity.
 */

    minDelta = 1000000;
    for (i = 1; i < M; i++)
    {
	Delta = (int)( 1.0E6 * (timesfound[i]-timesfound[i-1]));
	minDelta = MIN(minDelta, MAX(Delta,0));
    }

   return(minDelta);
}



/* A gettimeofday routine to give access to the wall
   clock timer on most UNIX-like systems.  */

//double mysecond()
//{
//        struct timeval tp;
//        struct timezone tzp;
//        int i;
//
//        i = gettimeofday(&tp,&tzp);
//        return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
//}

double myclock() {
#ifdef ANGEL
  angel_control_call(PRINT_CYCLES_ARG);
#endif
  return (double)clock() / CLOCKS_PER_SEC;
}

void checkSTREAMresults ()
{
	double aj,bj,cj,scalar;
	double asum,bsum,csum;
	double epsilon;
	int	j,k;

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
	aj = aj * (double) (N);
	bj = bj * (double) (N);
	cj = cj * (double) (N);

	asum = 0.0;
	bsum = 0.0;
	csum = 0.0;
	for (j=0; j<N; j++) {
		asum += a[j];
		bsum += b[j];
		csum += c[j];
	}
#ifdef VERBOSE
	printf ("Results Comparison: \n");
	printf ("        Expected  : %f %f %f \n",aj,bj,cj);
	printf ("        Observed  : %f %f %f \n",asum,bsum,csum);
#endif

#ifndef abs
#define abs(a) ((a) >= 0 ? (a) : -(a))
#endif
	epsilon = 1.e-8;

	if (abs(aj-asum)/asum > epsilon) {
		printf ("Failed Validation on array a[]\n");
		printf ("        Expected  : %f \n",aj);
		printf ("        Observed  : %f \n",asum);
	}
	else if (abs(bj-bsum)/bsum > epsilon) {
		printf ("Failed Validation on array b[]\n");
		printf ("        Expected  : %f \n",bj);
		printf ("        Observed  : %f \n",bsum);
	}
	else if (abs(cj-csum)/csum > epsilon) {
		printf ("Failed Validation on array c[]\n");
		printf ("        Expected  : %f \n",cj);
		printf ("        Observed  : %f \n",csum);
	}
	else {
		printf ("Solution Validates\n");
	}
}

void tuned_STREAM_Copy()
{
	int j;
#pragma omp parallel for
        for (j=0; j<N; j++)
            c[j] = a[j];
}

void tuned_STREAM_Scale(double scalar)
{
	int j;
#pragma omp parallel for
	for (j=0; j<N; j++)
	    b[j] = scalar*c[j];
}

void tuned_STREAM_Add()
{
	int j;
#pragma omp parallel for
	for (j=0; j<N; j++)
	    c[j] = a[j]+b[j];
}

void tuned_STREAM_Triad(double scalar)
{
	int j;
#pragma omp parallel for
	for (j=0; j<N; j++)
	    a[j] = b[j]+scalar*c[j];
}
