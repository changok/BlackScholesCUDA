#include "black_scholes.h"
#include "parser.h"
#include "random.h"
#include "timer.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * Usage: ./hw1.x <filename> <nthreads>
 *
 * <filename> (don't include the angle brackets) is the name of 
 * a data file in the current directory containing the parameters
 * for the Black-Scholes simulation.  It has exactly six lines 
 * with no white space.  Put each parameter one to a line, with
 * an endline after it.  Here are the parameters:
 *
 * S
 * E
 * r
 * sigma
 * T
 * M
 *
 * <nthreads> (don't include the angle brackets) is the number of
 * worker threads to use at a time in the benchmark.  The sequential
 * code which we supply to you doesn't use this argument; your code
 * will.
 */
int
main (int argc, char* argv[])
{
  confidence_interval_t interval;
  double S, E, r, sigma, T;
  int M = 0;
  char* filename = NULL;
  int nthreads = 1;
  double t1, t2, prng_stream_spawn_time;
  
  if (argc < 3)
    {
      fprintf (stderr, 
	       "Usage: ./hw1.x <filename> <nthreads>\n\n");
      exit (EXIT_FAILURE);
    }
  filename = argv[1];
  nthreads = to_int (argv[2]);
  parse_parameters (&S, &E, &r, &sigma, &T, &M, filename);

  /* 
   * Make sure init_timer() is only called by one thread,
   * before all the other threads run!
   */
  init_timer ();

  /* Same goes for initializing the PRNG */
  init_prng (random_seed ());

  /*
   * Run the benchmark and time it.
   */
  t1 = get_seconds ();
  /* 
   * In the parallel case, you may want to set prng_stream_spawn_time to 
   * the max of all the prng_stream_spawn_times, or just take a representative
   * sample... 
   */
  black_scholes (&interval, &prng_stream_spawn_time, S, E, r, sigma, T, M);
  t2 = get_seconds ();

  /*
   * A fun fact about C string literals (i.e., strings enclosed in
   * double quotes) is that the C preprocessor automatically
   * concatenates them if they are separated only by whitespace.
   */
  printf ("Black-Scholes benchmark:\n"
	  "------------------------\n"
	  "S        %g\n"
	  "E        %g\n"
	  "r        %g\n"
	  "sigma    %g\n"
	  "T        %g\n"
	  "M        %d\n",
	  S, E, r, sigma, T, M);
  printf ("Confidence interval: (%g, %g)\n", interval.min, interval.max);
  printf ("Total simulation time: %g seconds\n", t2 - t1);
  printf ("PRNG stream spawn time: %g seconds\n", prng_stream_spawn_time);

  return 0;
}



