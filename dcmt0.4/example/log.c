#include "log.h"
#include <assert.h>
#include <pthread.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


static pthread_mutex_t log_level_lock = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t log_write_lock = PTHREAD_MUTEX_INITIALIZER;
static int __mfh_log_level = 0;
static FILE* __mfh_log_stream = NULL;

static void
mfh_log_helper (const char* fmt, va_list ap)
{
  vfprintf (__mfh_log_stream, fmt, ap);
}



void
mfh_set_log_level (const int level)
{
  assert (0 == pthread_mutex_lock (&log_level_lock));
  __mfh_log_level = level;
  assert (0 == pthread_mutex_unlock (&log_level_lock));
}

void
mfh_log (const int level, const char* fmt, ...)
{
  int current_level = 0;

  assert (0 == pthread_mutex_lock (&log_level_lock));
  current_level = __mfh_log_level;
  assert (0 == pthread_mutex_unlock (&log_level_lock));

  if (level <= current_level)
    {
      va_list ap;

      va_start (ap, fmt);
      mfh_log_helper (fmt, ap);
      va_end (ap);
      assert (0 == pthread_mutex_unlock (&log_write_lock));
    }
}

static void
mfh_log_close ()
{
  assert (0 == pthread_mutex_lock (&log_write_lock));    
  if (__mfh_log_stream != NULL)
    assert (0 == fclose (__mfh_log_stream));
  assert (0 == pthread_mutex_unlock (&log_write_lock));
}

/**
 * DO NOT protect externally by locking log_write_lock, because it's
 * not a recursive lock.
 */
static void
mfh_log_open (const char* filename, const int appendp)
{
  time_t tod; 
    
  assert (0 == pthread_mutex_lock (&log_write_lock));  
  if (appendp)
    __mfh_log_stream = fopen (filename, "a");
  else
    __mfh_log_stream = fopen (filename, "w");
  assert (__mfh_log_stream != NULL);
  assert (0 == pthread_mutex_unlock (&log_write_lock));

  /* 
   * Log an initial message to the beginning of the file, with the
   * current time.
   */
  if (-1 == time (&tod))
    mfh_log (0, "\n\n### BEGIN LOGGING (can\'t get current time)\n\n");
  else
    {
      char ctime_result[50];
      ctime_r (&tod, ctime_result);
      mfh_log (0, "\n\n### BEGIN LOGGING %s\n\n", ctime_result);
    }
} 

void
mfh_log_init (const int level, const char* filename, const int appendp)
{
  mfh_log_open (filename, appendp);

  /* Register mfh_log_close() to run at exit() or return from main() */
  assert (0 == atexit (&mfh_log_close));

  mfh_set_log_level (level);
}
