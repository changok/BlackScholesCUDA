#include "../include/serialize.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>


#ifdef EXPECTPOS( cmd )
#undef EXPECTPOS( cmd )
#endif

#define EXPECTPOS( cmd ) do { \
    assert ( (cmd) > 0 );     \
  } while(0)

#ifdef EXPECTN( num, cmd )
#undef EXPECTN( num, cmd )
#endif

#define EXPECTN( num, cmd ) do { \
    assert ( num == (cmd) );	 \
  } while(0)


void
serialize_mt_struct (FILE* out, const mt_struct* mt)
{
  int k;

  EXPECTPOS( fprintf (out, "%u\n", mt->aaa) );
  EXPECTPOS( fprintf (out, "%d %d %d %d\n", 
		      mt->mm, mt->nn, mt->rr, mt->ww) );
  EXPECTPOS( fprintf (out, "%u %u %u\n", 
		      mt->wmask, mt->umask, mt->lmask) );
  EXPECTPOS( fprintf (out, "%d %d %d %d\n", 
		      mt->shift0, mt->shift1, mt->shiftB, mt->shiftC) );
  EXPECTPOS( fprintf (out, "%u %u\n", mt->maskB, mt->maskC) );
  /* i: current position in the state vector */
  EXPECTPOS( fprintf (out, "%d\n", mt->i) );

  if (mt->nn > 0)
    EXPECTPOS( fprintf (out, "%u", mt->state[0]) );
  for (k = 1; k < mt->nn; k++)
    EXPECTPOS( fprintf (out, " %u", mt->state[k]) );
  EXPECTPOS( fprintf (out, "\n") );
}

/**
 * Reconstruct (load in) an mt_struct object from a stream.
 */ 
mt_struct*
deserialize_mt_struct (FILE* in)
{
  int k;
  mt_struct* mt = malloc (sizeof (mt_struct));
  assert (mt != NULL);

  EXPECTN( 1, fscanf(in, "%u\n", &(mt->aaa)) );
  EXPECTN( 4, fscanf(in, 
		     "%d %d %d %d\n", 
		     &(mt->mm), &(mt->nn), &(mt->rr), &(mt->ww)) );

  EXPECTN( 3, fscanf(in, 
		     "%u %u %u\n", 
		     &(mt->wmask), &(mt->umask), &(mt->lmask)) );
  EXPECTN( 4, fscanf(in, 
		     "%d %d %d %d\n", 
		     &(mt->shift0), 
		     &(mt->shift1), 
		     &(mt->shiftB), 
		     &(mt->shiftC)) );
  EXPECTN( 2, fscanf(in, "%u %u\n", &(mt->maskB), &(mt->maskC)) );
  EXPECTN( 1, fscanf(in, "%d\n", &(mt->i)) );

  assert (mt->nn >= 0);
  if (mt->nn > 0)
    {
      mt->state = calloc (mt->nn, sizeof (uint32_t));
      assert (mt->state != NULL);
    }
  else
    mt->state = NULL;

  if (mt->nn > 0)
    EXPECTN( 1, fscanf(in, "%u", &(mt->state[0])) );
  for (k = 1; k < mt->nn; k++)
    EXPECTN( 1, fscanf(in, " %u", &(mt->state[k])) );

  return mt;
}
