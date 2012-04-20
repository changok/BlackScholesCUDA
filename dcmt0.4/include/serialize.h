#ifndef _mfh_serialize_h
#define _mfh_serialize_h

#include "dc.h"
#include <stdio.h>

/**
 * Serialize the mt_struct object.  This may be somewhat
 * implementation-dependent, depending on sizeof(int).
 */
void
serialize_mt_struct (FILE* out, const mt_struct* mt);

/**
 * Reconstruct (load in) an mt_struct object from a stream.
 */ 
mt_struct*
deserialize_mt_struct (FILE* in);

#endif /* _mfh_serialize_h */


