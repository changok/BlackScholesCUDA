#include "../include/dc.h"
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>

/**
 * Non-thread-safe implementation of get_random_seed().  The actual
 * get_random_seed() function protects the __get_random_seed() call 
 * with a mutex.  (The non-thread-safe part is when we read from 
 * /dev/random or /dev/urandom, or call time() (???).)
 */
static uint32_t
__get_random_seed ()
{
  const int nread = 4;
  FILE* r = NULL;
  uint32_t seed = (uint32_t) 0;
  char buffer[4];

  /* Make sure that a uint32_t is actually four bytes */
  assert ((size_t) nread == sizeof (uint32_t));

  /* 
   * Try to read some random bits from /dev/random.  The Darwin (MacOS
   * X) 10.4 documentation says that "[o]n Linux, /dev/urandom will
   * produce lower quality output if the entropy pool drains, while
   * /dev/random will prefer to block and wait for additional entropy
   * to be collected."  MacOS X allows (root) to write to /dev/random
   * in order to add what one thinks is entropy.
   */
  r = fopen ("/dev/random", "r");
  if (r == NULL)
    {
      /* Try to read some random bits from /dev/urandom */
      r = fopen ("/dev/urandom", "r");
      if (r == NULL)
	{
	  time_t current_time;
	  fprintf (stderr, "*** Warning: Failed to open /dev/urandom "
		   "or /dev/urandom, so falling back on system clock for"
		   " random seed.  This is undesirable because the time"
		   " at which a program is run isn't usually random ***");
	  current_time = time (NULL);
	  return (uint32_t) current_time;
	}
    }

  assert ((size_t) nread == fread (buffer, sizeof(char), (size_t) nread, r));
  seed = *((uint32_t*) buffer);
  fclose (r);
  return seed;
}

static pthread_mutex_t random_seed_mutex = PTHREAD_MUTEX_INITIALIZER;

uint32_t
get_random_seed ()
{
  uint32_t retval;
  assert (0 == pthread_mutex_lock (&random_seed_mutex));
  retval = __get_random_seed ();
  assert (0 == pthread_mutex_unlock (&random_seed_mutex));
  return retval;
}
