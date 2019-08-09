// *************************************************
// * Program:  STREAM2                             *
// * Revision: 0.1, 99.10.26                       *
// * Author:   John McCalpin                       *
// *           john@mccalpin.com                   *
// *************************************************
// * This program measures sustained bandwidth     *
// * using four computational kernels:             *
// *                                               *
// *       FILL:   a(i) = 0                        *
// *       COPY:   a(i) = b(i)                     *
// *       DAXPY:  a(i) = a(i) + q*b(i)            *
// *       SUM:    sum += a(i)                     *
// *                                               *
// * Results are presented in MB/s, assuming       *
// *   8 Bytes per iteration for FILL and SUM,     *
// *  16 Bytes per iteration for COPY, and         *
// *  24 Bytes per iteration for DAXPY             *
// *************************************************
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

#undef IA32
#undef AIX
#ifdef IA32
//IA32
extern double compute_tdelta ( int startlow, int starthigh, int endlow, int endhigh );
extern double compute_time ( int low, int high );
#include <asm/msr.h>
int	startlow,starthigh,endlow,endhigh,low,high;
#define START_TIMER rdtsc(startlow,starthigh);
#define STOP_TIMER rdtsc(endlow,endhigh);
#define COMPUTE_DELTA_TIME tdelta = compute_tdelta ( startlow, starthigh, endlow, endhigh);
#elif AIX
//AIX
#include <sys/time.h>
#include <sys/systemcfg.h>
timebasestruct_type start, finish;
#define START_TIMER read_real_time (&start, TIMEBASE_SZ);
#define STOP_TIMER read_real_time (&finish, TIMEBASE_SZ);
#define COMPUTE_DELTA_TIME {timebase_to_time(&start, TIMEBASE_SZ);
                            timebase_to_time(&finish, TIMEBASE_SZ);
#else
//x86_64/AARCH64
static double	avgtime[4] = {0, 0, 0, 0}, maxtime[4] = {0, 0, 0, 0},
		mintime[4] = {FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX};

#ifndef CLOCKS_PER_SEC
//#define CLOCKS_PER_SEC 1000000000UL
#define CLOCKS_PER_SEC 1UL
#endif
clock_t clock(void);      //forward declaration of sting function
int     checktick();
extern double myclock();
double	tstart, tstop;
#define START_TIMER tstart = myclock();
#define STOP_TIMER tstop = myclock();
#define COMPUTE_DELTA_TIME tdelta = tstop - tstart;

# define HLINE "-------------------------------------------------------------\n"

# ifndef MIN
# define MIN(x,y) ((x)<(y)?(x):(y))
# endif
# ifndef MAX
# define MAX(x,y) ((x)>(y)?(x):(y))
# endif
#endif

// 	program stream2
// 	integer NMIN, NMAX, NTIMES, NUMSIZES
// 	parameter (NMIN=30,NMAX=2 000 000)
// 	parameter (NTIMES=5,NUMSIZES=32)
// 	integer NPAD
// 	parameter (NPAD=5)
#undef COMPUTE_ITER_SIZE
#ifdef COMPUTE_ITER_SIZE
/*Default*/
#define NMIN		(30)
#define NMAX		(2000000)
#define NTIMES	(5)
#define NUMSIZES	(32)
#define NPAD		(5)
#else
/*QDT Modified*/
int       NTIMES    = 0;
#define   NUMSIZES	 (1)
#pragma message ("Size of array(bytes)=" STRINGIFY(NMAX))
#endif

static char	*label[4] = {"FILL1:     ", "COPY:      ",
    "DAXPY:     ", "SUM:       "};

//
// 	real*8 a(NMAX+NPAD),b(NMAX+NPAD)
double a[NMAX+NPAD], b[NMAX+NPAD], c[NMAX+NPAD];
/* Make these global, because most systems will not allow a
   sufficiently large stack by default. */

// Test Configuration
typedef int bool;
#define false 0
#define true  1

enum
{
  g_testAll      = 0,
  g_testFill     = 1,
  g_testCopy     = 2,
  g_testDaxpy    = 3,
  g_testSum      = 4,
  g_testFillZero = 5,
  g_testMax      = 6
};

typedef struct TestOpts_ {
  const char*   name;
}TestOpts;

void
streamv2_usage(int argc, char *argv[], char* usage)
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
    //streamv2_usage();
  }

  if (*end != '\0')
  {
    // Extra characters at the end.
    printf("Bad size: %s\n", s);
    //streamv2_usage();
  }

  return val;
}

int main( int ac, char **av )
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

  int                quantum;
  double                   t;
  double              **time;
  int                 op;
  double sum = 0.0, scalar;
  // 	real*8 sum0,sum1,sum2,sum3,sum4,sum5,start,finish
  double  sum0, sum1, sum2, sum3, sum4, sum5, start, finish;
  // 	real*8 rate(4),besttime(4),bytes(4)
  // 	real*8 exp,tdelta
  double expt, tdelta;
  // 	integer i,j,k,l,M
  int i, j, k, l, M, n, inner;
  // 	external second
  //
  // 	data bytes/8,16,24,8/
  static double	bytes[4] = {
      1 * sizeof(double) * NMAX,
      2 * sizeof(double) * NMAX,
      3 * sizeof(double) * NMAX,
      1 * sizeof(double) * NMAX
  };

  const char   *test = "all";
  bool isEnableAngelSignals = false;
  char   *usage = "[-T <fill/copy/daxpy/sum/fillzero/all>]  [-M <size(K/M)>]\
                  [-P <NPAD (offset. Match R)>] [-R <repetitions>] \
                  [-A <enable(1)/disable angel signals(0)>] \n";
  TestOpts g_testopts[] =
  {
    { "all" },
    { "fill" },
    { "copy" },
    { "daxpy" },
    { "sum" },
    { "fillzero" }
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
    case 'M':
        /*N = easy_size(optarg);*/
        break;
    case 'P':
        /*NPAD = atoi(optarg);*/
        break;
    case 'R':
        NTIMES = atoi(optarg);
        break;
    default:
        streamv2_usage(ac, av, usage);
        break;
    }
  }

  //Allocate and initialize time array
  time = malloc(sizeof(double) * 4);
  for(j=0; j<4; j++)
  {
    time[j] = malloc( sizeof(time[j][0]) * NTIMES );
  }
  for (k=0; k<NTIMES; k++)
  {
    for (j=0; j<4; j++)
    {
      time[j][k] = 0;
    }
  }

  printf("Array size = %llu, Padding = %d\n", (unsigned long long) NMAX, NPAD);

  printf("Total memory required = %.1f MB.\n",
	(3.0 * sizeof(double)) * ( (double) NMAX / 1048576.0));
  printf("Each test is run %d times, but only\n", NTIMES);
  printf("the *best* time for each is used.\n");

  printf(HLINE);
  //for (j=0; j<(NMAX+NPAD); j++)
  //{
    //a[j] = 1.0;
    //b[j] = 2.0;
    //c[j] = 0.0;
  //}

  if  ( (quantum = checktick()) >= 1)
  {
    printf("Your clock granularity/precision appears to be "
        "%d microseconds.\n", quantum);
  }
  else
  {
    printf("Your clock granularity appears to be "
        "less than one microsecond.\n");
    quantum = 1;
  }

  t = myclock();
  for (j = 0; j < NMAX+NPAD; j++)
  {
    a[j] = 2.0E0 * a[j];
  }
  t = 1.0E6 * (myclock() - t);

  printf("Each test below will take on the order"" of %d microseconds.\n", (int) t  );
  printf("   (= %d clock ticks)\n", (int) (t/quantum) );
  printf("Increase the size of the arrays if this shows that\n");
  printf("you are not getting at least 20 clock ticks per test.\n");

  printf(HLINE);

  printf("WARNING -- The above is only a rough guideline.\n");
  printf("For best results, please be sure you know the\n");
  printf("precision of your system timer.\n");
  printf(HLINE);

  if( algo[g_testFillZero] )
  {
    scalar = 0.0;
  }
  else
  {
    scalar = 1.0;  // Avoid FP overflows.
  }

  // * Loop over problem size
  // 	do j=1,NUMSIZES
  for( j = 0; j < NUMSIZES; ++j )
  {
      //printf("Value of NUMSIZES=%d\n", NUMSIZES);

  // 	    exp = log10(dble(NMIN)) + dble(j-1)/dble(NUMSIZES-1)*
  //      $          (log10(dble(NMAX))-log10(dble(NMIN)))
  // 	    M = NINT(10.**exp)
  //
#ifdef COMPUTE_ITER_SIZE
    expt = log( (double)(NMIN)) +
          ( (double) j / (double) (NUMSIZES-1) )
          * (log((double) NMAX)-log((double) NMIN));

    M = 0.5 + exp( expt );
#else
    //This is set to array size
    M = NMAX;
#endif

  //
  // * Initialize Arrays
  //
  // 	    do i=1,M
  // 	        a(i) = 0.0d0
  // 	        b(i) = 0.0d0
  // 	    end do

    for( i = 0; i < M+NPAD; ++i )  /* should this go up to M+NPAD ?? */
    {
      a[i] = 0.0;
      b[i] = 0.0;
      c[i] = 0.0;
    }

    //
    // 	    do k=1,NTIMES
    for( k = 0; k < NTIMES; ++k )
    {
      //inner = NMAX/M
      /* inner = NMAX / M; */
      inner = 1;

  //
  //      start = second()
  // 	    do l=1,inner
  // 	        scalar = a(M)+dble(k+l)
  // 	        do i=1,M
  // 	         a(i) = scalar
  // 	        end do
  // 	    end do
  // 	    finish = second()
  // 	    time(1,k) = (finish-start)/dble(inner)
      if( algo[g_testAll] || algo[g_testFill] || algo[g_testFillZero] )
      {
#ifdef ANGEL
        if(isEnableAngelSignals && (k==0))
        {
          workload_ckpt_begin();
        }
#endif /* ANGEL  */
        START_TIMER;
        for( l = 0; l < inner; ++l )
        {
          // scalar = a[M-1] + (double) (k+l+2);
          //printf("Value of inner=%d and M=%d and scalar=%f\n", inner, M, scalar);
          for( i = 0; i < M; ++i )
          {
            a[i] = scalar;
          }
        }
        STOP_TIMER;
        COMPUTE_DELTA_TIME;
        time[0][k] = tdelta/(double) (inner);
        //printf("Fill time[%d][%d]=%8.8f\n", 0 , k, time[0][k]);
      }

  //
  // 	    start = second()
  // 		do l=1,inner
  // 	        scalar = b(M)+dble(k+l)
  // 		    a(l) = scalar
  // 	        do i=1,M
  // 		        b(i) = a(i)
  // 	        end do
  // 	    end do
  // 	    finish = second()
  // 	    time(2,k) = (finish-start)/dble(inner)
      if( algo[g_testAll] || algo[g_testCopy] )
      {
#ifdef ANGEL
        if(isEnableAngelSignals && (k==0) && !algo[g_testAll])
        {
          workload_ckpt_begin();
        }
#endif /* ANGEL  */
        START_TIMER;
        for( l = 0; l < inner; ++l )
        {
          // scalar = b[M-1] + (double) (k+l+2);
          a[l] = scalar;

          for( i = 0; i < M; ++i )
          {
            b[i] = a[i];
          }
        }
        STOP_TIMER;
        COMPUTE_DELTA_TIME;
        time[1][k] = tdelta/(double) (inner);
        //printf("Copy time[%d][%d]=%8.8f\n", 1, k, time[1][k]);
      }


  //
  // 	    start = second()
  // 		do l=1,inner
  // 	        scalar = b(M)+dble(k+l)
  // 		    a(l) = scalar
  // 	        do i=1,M
  // 		        b(i) = b(i) + scalar*a(i)
  // 	        end do
  // 	    end do
  // 	    finish = second()
  // 	    time(3,k) = (finish-start)/dble(inner)
      if( algo[g_testAll] || algo[g_testDaxpy] )
      {
#ifdef ANGEL
        if(isEnableAngelSignals && (k==0) && !algo[g_testAll])
        {
          workload_ckpt_begin();
        }
#endif /* ANGEL  */
        START_TIMER;
        for( l = 0; l < inner; ++l )
        {
          // scalar = b[M-1] + (double) (k+l+2);
          a[l] = scalar;

          for( i = 0; i < M; ++i )
          {
            b[i] = b[i] + scalar * a[i];
          }
        }
        STOP_TIMER;
        COMPUTE_DELTA_TIME;
        time[2][k] = tdelta/(double) (inner);
      }


  //
  // 	    start = second()
  // 		do l=1,inner
  // 	        scalar = b(M)+dble(k+l)
  // 		    b(l) = scalar
  // 	        sum0 = 0.0d0
  // 	        sum1 = 0.0d0
  // 	        sum2 = 0.0d0
  // 	        sum3 = 0.0d0
  // 	        sum4 = 0.0d0
  // 	        sum5 = 0.0d0
  // 	        do i=1,M,6
  // 		        sum0 = sum0 + b(i+0)
  // 		        sum1 = sum1 + b(i+1)
  // 		        sum2 = sum2 + b(i+2)
  // 		        sum3 = sum3 + b(i+3)
  // 		        sum4 = sum4 + b(i+4)
  // 		        sum5 = sum5 + b(i+5)
  // 	        end do
  // 	    end do
  // 		sum = sum0 + sum1 + sum2 + sum3 + sum4 + sum5
  // 	    finish = second()
  // 	    time(4,k) = (finish-start);
      if( algo[g_testAll] || algo[g_testSum] )
      {
#ifdef ANGEL
        if(isEnableAngelSignals && (k==0) && !algo[g_testAll])
        {
          workload_ckpt_begin();
        }
#endif /* ANGEL  */
        START_TIMER;
        for( l = 0; l < inner; ++l )
        {
          // scalar = b[M-1] + (double) (k+l+2);
          b[l] = scalar;	/* WTF is this ?!? */
          sum0 = 0.0;
          sum1 = 0.0;
          sum2 = 0.0;
          sum3 = 0.0;
          sum4 = 0.0;
          sum5 = 0.0;

          // This code accesses beyond M elements.  Has the data been fixed
          // to cope ?
          for( i = 0; i < M-5; i += 6 )
          {
            sum0 += b[i+0];
            sum1 += b[i+1];
            sum2 += b[i+2];
            sum3 += b[i+3];
            sum4 += b[i+4];
            sum5 += b[i+5];
          }
          sum = sum0 + sum1 + sum2 + sum3 + sum4 + sum5;
        }

        STOP_TIMER;
        COMPUTE_DELTA_TIME;
        time[3][k] = tdelta/(double) (inner);
        //printf("Sum time[%d][%d]=%8.8f and sum=%8.8f\n", 3, k, time[3][k], sum);
      }
    }//End of NTIMES for...loop

    /*	--- SUMMARY --- */

    for (k=1; k<NTIMES; k++) /* note -- skip first iteration */
    {
      for (n=0; n<4; n++)
      {
        avgtime[n] = avgtime[n] + time[n][k];
        mintime[n] = MIN(mintime[n], time[n][k]);
        maxtime[n] = MAX(maxtime[n], time[n][k]);
      }
    }

    printf("Function      Bytes        AvgClks        MinClks       MaxClks\n");
    for (n=0; n<4; n++)
    {
      avgtime[n] = avgtime[n]/(double)(NTIMES-1);

      if ( algo[g_testFillZero] && (n == 0))
      {
        printf("%s%11.4f  %11.11f  %11.11f  %11.11f\n", "FILL0:     ",
               bytes[n],
               avgtime[n],
               mintime[n],
               maxtime[n]);
      }
      else
      {
        printf("%s%11.4f  %11.11f  %11.11f  %11.11f\n", label[n],
               bytes[n],
               avgtime[n],
               mintime[n],
               maxtime[n]);
      }
    }
    //Below print is to prevent compiler optimization for SUM
    printf("Total sum for SUM: %14.8E\n", sum);

  //
  // 	    open (unit=3,form='unformatted')
  // 	    write (3) sum
  // 	    close (unit=3)
  //
  // 	end do
  }//End of NUM_SIZES for..loop
//
//     1	format (1x,i8,2x,i4,1x,5(f8.1,2x))
// 	end
#ifdef ANGEL
  if (isEnableAngelSignals)
  {
    workload_ckpt_end();
  }
#endif /* ANGEL  */
  return 0;
}

/* A Fortran-callable gettimeofday routine to give access
   to the wall clock timer.

   This subroutine may need to be modified slightly to get
   it to link with Fortran on your computer.
*/

//#include <sys/time.h>
/* int gettimeofday(struct timeval *tp, struct timezone *tzp); */

//double second()
//{
/* struct timeval { long	tv_sec;
	    long	tv_usec;	};

struct timezone { int	tz_minuteswest;
	     int	tz_dsttime;	 };	*/

	//struct timeval tp;
	//struct timezone tzp;
	//int i;

	//i = gettimeofday(&tp,&tzp);
	//return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
//}

# define	T	20

int checktick()
{
    int		i, minDelta, Delta;
    double	t1, t2, timesfound[T];

/*  Collect a sequence of T unique time values from the system. */

/* Making assumption of 500 MHz */

    for (i = 0; i < T; i++)
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
    for (i = 1; i < T; i++)
    {
	Delta = (int)( 1.0E6 * (timesfound[i]-timesfound[i-1]));
	minDelta = MIN(minDelta, MAX(Delta,0));
    }

   return(minDelta);
}


double myclock()
{
#ifdef ANGEL
  angel_control_call(PRINT_CYCLES_ARG);
#endif
  return (double)clock() / CLOCKS_PER_SEC;
}

#ifdef IA32
double compute_tdelta( int startlow, int starthigh, int endlow, int endhigh )
{
	double time;

	time = ((double) (endhigh-starthigh) * 4294967296. + (double) (endlow-startlow)) / 1600000000.;
	return (time);
}
double compute_time( int low, int high )
{
	double time;

	time = ((double) high * 4294967296. + (double) low) / 1600000000.;
	return (time);
}
#elif AIX
double compute_tdelta( timebasestruct finish, timebasestruct start )
{
	int secs, nsecs;

	timebase_to_time(&start, TIMEBASE_SZ);
        timebase_to_time(&finish, TIMEBASE_SZ);

	secs = finish.tb_high - start.tb_high;
	nsecs = finish.tb_low - start.tb_low;

	if (nsecs < 0) {		/* check for carry and fix if needed */
		secs--;
		nsecs += 1000000000;
	}
	return ( (double) secs + (double) nsecs * 0.000000001);
}
double compute_time( timebasestruct start )
{
	double time;

	timebase_to_time(&start, TIMEBASE_SZ);

	time = ((double) high + (double) low * 0.000000001);
	return (time);
}
#endif

