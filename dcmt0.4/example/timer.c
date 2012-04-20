#include "timer.h"
#include <assert.h>
#include <stdlib.h>
#include <sys/time.h>

static inline double
get_seconds_gettimeofday ()
{
  struct timeval tv;
  assert (0 == gettimeofday (&tv, NULL));
  return (double) tv.tv_sec + 1.0e-6 * ((double) tv.tv_usec);
}

static double
gettimeofday_resolution ()
{
  double t0, t1;
  t0 = get_seconds_gettimeofday ();
  while ((t1 = get_seconds_gettimeofday ()) == t0)
    ;
  assert (t1 > t0);
  return t1 - t0;
}

#if defined(__ppc__) || defined(__powerpc__)
#  include <math.h>
#  include <pthread.h>
#  include <unistd.h>
#  include <stdio.h>

typedef unsigned long long ticks;
static pthread_mutex_t calibrate_timer_mutex = PTHREAD_MUTEX_INITIALIZER;
static double timer_factor = 1.0;

static inline ticks 
getticks ()
{
  unsigned int tbl, tbu0, tbu1;

   do {
     __asm__ __volatile__ ("mftbu %0" : "=r"(tbu0));
     __asm__ __volatile__ ("mftb %0" : "=r"(tbl));
     __asm__ __volatile__ ("mftbu %0" : "=r"(tbu1));
   } while (tbu0 != tbu1);

   return (((unsigned long long)tbu0) << 32) | tbl;
}

static void
calibrate_ticks ()
{
  ticks t0, t1;
  double s0, s1;
  t0 = getticks ();
  s0 = get_seconds_gettimeofday ();
  sleep (1); 
  t1 = getticks ();
  s1 = get_seconds_gettimeofday ();
  assert (s1 - s0 > 0.0);
  assert ((double)t1 - (double)t0 > 0.0);

  assert (0 == pthread_mutex_lock (&calibrate_timer_mutex));
  timer_factor = (s1 - s0) / ((double) t1 - (double) t0);
  assert (0 == pthread_mutex_unlock (&calibrate_timer_mutex));
}

double
get_seconds ()
{
  return (double) getticks () * timer_factor;
}

static double
timer_resolution ()
{
  double t0, t1;
  t0 = get_seconds ();
  while ((t1 = get_seconds ()) == t0)
    ;
  assert (t1 > t0);
  return t1 - t0;
}

static int
compare_double (const void* px, const void* py)
{
  const double x = *((double*) px);
  const double y = *((double*) py);
  if (x < y)
    return -1;
  else if (x > y)
    return +1;
  else
    return 0;
}

void
init_timer ()
{
  const int ntrials = 100;
  double resolutions[100];
  double mean = 0.0;
  double variance = 0.0;
  int i;

  calibrate_ticks ();

  for (i = 0; i < ntrials; i++)
    resolutions[i] = timer_resolution ();
  qsort (resolutions, ntrials, sizeof (double), &compare_double);
  for (i = 0; i < ntrials; i++)
    mean += resolutions[i] / (double) ntrials;
  for (i = 0; i < ntrials; i++)
    {
      const double diff = resolutions[i] - mean;
      variance += diff * diff / (double) ntrials;
    }
  printf ("Min PPC timer resolution: %g seconds\n", 
	  resolutions[0]);
  printf ("Max PPC timer resolution: %g seconds\n", 
	  resolutions[ntrials - 1]);
  printf ("PPC timer stddev: %g\n", sqrt (variance));

  for (i = 0; i < ntrials; i++)
    resolutions[i] = gettimeofday_resolution ();
  qsort (resolutions, ntrials, sizeof (double), &compare_double);
  for (i = 0; i < ntrials; i++)
    mean += resolutions[i] / (double) ntrials;
  for (i = 0; i < ntrials; i++)
    {
      const double diff = resolutions[i] - mean;
      variance += diff * diff / (double) ntrials;
    }
  printf ("Min gettimeofday() resolution: %g seconds\n", 
	  resolutions[0]);
  printf ("Max gettimeofday() resolution: %g seconds\n", 
	  resolutions[ntrials - 1]);
  printf ("gettimeofday() timer stddev: %g\n", sqrt (variance));
}


#else
#  warning "Defaulting to gettimeofday() timer"
double
get_seconds ()
{
  return get_seconds_gettimeofday ();
}

void
init_timer ()
{
  (void) get_seconds ();
}

static double
timer_resolution ()
{
  double t0, t1;
  t0 = get_seconds ();
  while ((t1 = get_seconds ()) == t0)
    ;
  assert (t1 > t0);
  return t1 - t0;
}

void
init_timer ()
{
  const int ntrials = 100;
  double resolutions[100];
  double mean = 0.0;
  double variance = 0.0;
  int i;

  for (i = 0; i < ntrials; i++)
    resolutions[i] = timer_resolution ();
  qsort (resolutions, ntrials, sizeof (double), &compare_double);
  for (i = 0; i < ntrials; i++)
    mean += resolutions[i] / (double) ntrials;
  for (i = 0; i < ntrials; i++)
    {
      const double diff = resolutions[i] - mean;
      variance += diff * diff / (double) ntrials;
    }
  printf ("Min gettimeofday() timer resolution: %g seconds\n", 
	  resolutions[0]);
  printf ("Max gettimeofday() timer resolution: %g seconds\n", 
	  resolutions[ntrials - 1]);
  printf ("gettimeofday() timer stddev: %g\n", sqrt (variance));
}


#endif

