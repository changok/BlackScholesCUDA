/*
 * This is a third example for the Dynamic Creator library, with
 * serialization extensions by Mark Hoemmen <mhoemmen AT cs DOT
 * berkeley DOT edu>.
 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "../include/dc.h"
#include "../include/timer.h"
#include "../include/serialize.h"
#include "log.h"


static int
mt_struct_equalp (mt_struct* mt1,
		  mt_struct* mt2)
{
  int i;
  int failedp = 0;

#define mt_struct_equalp_test( thing ) do { \
    if (mt1->thing != mt2->thing)	    \
      {					    \
	mfh_log (0, "*** mt1->%s != mt2->%s ***\n", #thing, #thing); \
	failedp = 1;                                                 \
      }						                     \
  } while(0)

  mt_struct_equalp_test (aaa);
  mt_struct_equalp_test (mm);
  mt_struct_equalp_test (nn);
  mt_struct_equalp_test (rr);
  mt_struct_equalp_test (ww);
  mt_struct_equalp_test (wmask);
  mt_struct_equalp_test (umask);
  mt_struct_equalp_test (shift0);
  mt_struct_equalp_test (shift1);
  mt_struct_equalp_test (shiftB);
  mt_struct_equalp_test (shiftC);
  mt_struct_equalp_test (maskB);
  mt_struct_equalp_test (maskC);
  mt_struct_equalp_test (i);
  if (failedp)
    return 0;

  for (i = 0; i < mt1->nn; i++)
    if (mt1->state[i] != mt2->state[i])
      {
	mfh_log (0, "*** mt1->state[%d] != mt2->state[%d] ***\n", i, i);
	failedp = 1;
      }
  if (failedp)
    return 0;
  else 
    return 1;
}

static void
test_random_numbers (mt_struct* mt1,
		     mt_struct* mt2,
		     const int howmany)
{
  int i;
  for (i = 0; i < howmany; i++)
    {
      if (genrand_mt (mt1) != genrand_mt (mt2))
	{
	  mfh_log (2, "*** genrand_mt(mt1) != genrand_mt(mt2) at i = %d ***\n", i);
	  exit (EXIT_FAILURE);
	}
    }
}

static void
test_distribution_timed (double* timing,
			 mt_struct* mt,
			 const int ntrials)
{
  double* trials = NULL;
  double t1 = 0.0;
  double mean = 0.0;
  double variance = 0.0;
  int i;

  mfh_log (2, "=== test_distribution_timed ===\n");

  if (ntrials < 1)
    return;
  trials = calloc (ntrials, sizeof (double));
  assert (trials != NULL);
  
  t1 = get_seconds ();
  for (i = 0; i < ntrials; i++)
    trials[i] = genrand_double (mt);
  *timing = get_seconds () - t1;

  /* Test if the numbers are in the range [0,1) */
  mfh_log (2, "Testing range of double-precision random numbers...\n");
  for (i = 0; i < ntrials; i++)
    {
      if (trials[i] < 0.0 || trials[i] >= 1.0)
	{
	  mfh_log (0, "*** test_distribution_timed: trials[%d]"
		   " = %g is out of range [0,1) ***\n", i, trials[i]);
	  exit (EXIT_FAILURE);
	}
    }
  mfh_log (2, "All numbers in range...\n");

  /* Compute the mean */
  for (i = 0; i < ntrials; i++)
    mean += trials[i] / (double) ntrials;
  mfh_log (2, "Arithmetic mean: %g\n", mean);
  
  /* Compute the standard deviation */
  for (i = 0; i < ntrials; i++)
    {
      const double diff = trials[i] - mean;
      variance += diff * diff / (double) ntrials;
    }
  mfh_log (2, "Standard deviation: %g\n", sqrt (variance));
  
  mfh_log (2, "=== Done with test_distribution_timed ===\n");
}

static int
to_int (const char* s)
{
  int val;

  errno = 0;
  val = strtol (s, NULL, 10);
  assert (errno == 0);
  
  return val;
}

static mt_struct*
get_and_init_mt_timed (double* get_timing,
		       double* init_timing,
		       const int id, 
		       const uint32_t seed)
{
  double t1, t2;
  mt_struct* mt = NULL;

  mfh_log (2, " Calling get_mt_parameter_id...");
  t1 = get_seconds ();
  mt = get_mt_parameter_id (32,521,id);
  assert (mt != NULL);
  t2 = get_seconds ();
  *get_timing = t2 - t1;
  mfh_log (2, "done.  ");

  mfh_log (2, "Calling sgenrand_mt...");
  t1 = get_seconds ();
  sgenrand_mt (seed, mt);
  t2 = get_seconds ();
  *init_timing = t2 - t1;
  mfh_log (2, "done.  ");

  return mt;
}

static void
create_twisters_timed (mt_struct* mts[], 
		       double get_timings[],
		       double init_timings[],
		       const int nstructs,
		       const uint32_t seed)
{
  /* Even if the seed is the same for each MT, the sequences should be
     independent, as long as their ids are different. */
  int i;
  for (i = 0; i < nstructs; i++)
    {
      mt_struct* mt = NULL;
      mfh_log (2, "\tcreate_twisters_timed: i = %d\n", i);
      mt = get_and_init_mt_timed (&get_timings[i],
				      &init_timings[i],
				      i, 
				      seed);
      mfh_log (2, "\tcreate_twisters_timed: i = %d done\n", i);
      mts[i] = mt;
    }
}
		       
static void
serialize_mt_struct_timed (double* timing, FILE* out, mt_struct* mt)
{
  double t1, t2;

  t1 = get_seconds ();
  serialize_mt_struct (out, mt);
  t2 = get_seconds ();
  *timing = t2 - t1;
}

static void
serialize_twisters_timed (double timings[],
			  FILE* out, 
			  mt_struct* mts[], 
			  const int nstructs)
{
  int i;
  for (i = 0; i < nstructs; i++)
    serialize_mt_struct_timed (&timings[i], out, mts[i]);
}

static mt_struct*
deserialize_mt_struct_timed (double* timing, FILE* in)
{
  double t1, t2;
  mt_struct* mt = NULL;

  t1 = get_seconds ();
  mt = deserialize_mt_struct (in);
  t2 = get_seconds ();
  *timing = t2 - t1;

  return mt;
}

static void
deserialize_twisters_timed (double timings[],
			    mt_struct* mts[], 
			    FILE* in,
			    const int nstructs)
{
  int i;
  for (i = 0; i < nstructs; i++)
    mts[i] = deserialize_mt_struct_timed (&timings[i], in);
}

static void
print_timings (FILE* out,
	       double get_timings[],
	       double init_timings[],
	       double serialize_timings[],
	       double deserialize_timings[],
	       const int nstructs)
{
  int i;
  fprintf (out, "Get, init, save, restore (s)\n");
  for (i = 0; i < nstructs; i++)
    {
      fprintf (out, "%g, %g, %g, %g\n", 
	       get_timings[i],
	       init_timings[i],
	       serialize_timings[i],
	       deserialize_timings[i]);
    }
}



int
main (int argc, char* argv[])
{
  int nstructs = 3;
  mt_struct** mts1 = NULL;
  mt_struct** mts2 = NULL;
  double* get_timings = NULL;
  double* init_timings = NULL;
  double* serialize_timings = NULL;
  double* deserialize_timings = NULL;
  double other_timing = 0.0;
  int ntrials = 1000000;
  FILE* in = NULL;
  FILE* out = NULL;
  int i;

  mfh_log_init (2, "benchmark.log", 1);
  
  mfh_log (2, "initializing timer...");
  init_timer ();
  mfh_log (2, "done.\n");
  mfh_log (2, "initializing dynamic creator...");
  init_dc (get_random_seed ());
  /* init_dc (4172); */
  mfh_log (2, "done.\n");

  mfh_log (2, "reading command-line args...");
  if (argc > 1)
    nstructs = to_int (argv[1]);
  assert (nstructs >= 1);
  mfh_log (2, "got nstructs = %d, done.\n", nstructs);

  mfh_log (2, "allocating space for timing data...");
  get_timings = calloc (nstructs, sizeof(double));
  init_timings = calloc (nstructs, sizeof(double));
  serialize_timings = calloc (nstructs, sizeof(double));
  deserialize_timings = calloc (nstructs, sizeof(double));
  assert (get_timings != NULL && init_timings != NULL && 
	  serialize_timings != NULL && deserialize_timings != NULL);
  mfh_log (2, "done.\n");

  mfh_log (2, "opening serialization output file...");
  /* Create an output file for serialization */
  out = fopen ("mt.out", "w");
  assert (out != NULL);
  mfh_log (2, "done.\n");

  /* Let's compute some Mersenne Twisters */
  mts1 = calloc (nstructs, sizeof(mt_struct));
  assert (mts1 != NULL);
  mfh_log (2, "calling create_twisters_timed...");
  create_twisters_timed (mts1, get_timings, init_timings, nstructs, 1234);
  mfh_log (2, "done.\n");

  /* Now let's serialize them... */
  mfh_log (2, "calling serialize_twisters_timed...");
  serialize_twisters_timed (serialize_timings,
			    out, 
			    mts1,
			    nstructs);
  mfh_log (2, "done.\n");
  mfh_log (2, "closing output file...");
  /* ...close the output file... */
  assert (0 == fclose (out));
  mfh_log (2, "done.\n");

  /* ...open input file for serialization... */
  mfh_log (2, "opening input file...");
  in = fopen ("mt.out", "r");
  assert (in != NULL);
  mfh_log (2, "done.\n");

  /* ...and read them back in. */
  mfh_log (2, "reading twisters back in...");
  mts2 = calloc (nstructs, sizeof(mt_struct));
  assert (mts2 != NULL);
  deserialize_twisters_timed (deserialize_timings,
			      mts2,
			      in,
			      nstructs);
  for (i = 0; i < nstructs; i++)
    assert (mts2[i] != NULL);
  mfh_log (2, "done.\n");

  /* Make sure that the deserialized mt_struct objects produce the
     same random numbers as the originals */
  for (i = 0; i < nstructs; i++)
    {
      mfh_log (2, "Testing deserialized mt_struct %d\n", i);
      if (! mt_struct_equalp (mts1[i], mts2[i]))
	{
	  mfh_log (2, "*** mts1[%d] != mts2[%d] ***\n", i, i);
	  exit (EXIT_FAILURE);
	}
      test_random_numbers (mts1[i], mts2[i], 20);
   }

  /* 
   * Benchmark the generation of some random numbers in [0,1).  Test
   * their distribution properties to make sure that it's uniform and
   * in the right range.
   */
  test_distribution_timed (&other_timing, mts1[0], ntrials);
  printf ("Generated %d doubles: timing is %g s, which is %g s per number\n\n", 
	  ntrials, other_timing, other_timing / (double) ntrials);

  /* Print the resulting timings */
  mfh_log (2, "printing timings...");
  print_timings (stdout,
		 get_timings, 
		 init_timings, 
		 serialize_timings, 
		 deserialize_timings, 
		 nstructs);
  mfh_log (2, "done.\n");

  /* Clean up */
  mfh_log (2, "cleaning up...");
  for (i = 0; i < nstructs; i++)
    {
      if (mts1[i] != NULL)
	free_mt_struct (mts1[i]);
      if (mts2[i] != NULL)
	free_mt_struct (mts2[i]);
    }
  if (mts1 != NULL)
    free (mts1);
  if (mts2 != NULL)
    free (mts2);
  mfh_log (2, "done.\n");
  return 0;
}


