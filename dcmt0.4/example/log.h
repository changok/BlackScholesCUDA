#ifndef _mfh_log_h
#define _mfh_log_h

void
mfh_set_log_level (const int level);

void
mfh_log (const int level, const char* fmt, ...);

void
mfh_log_init (const int level, const char* filename, const int appendp);

#endif /* _mfh_log_h */
