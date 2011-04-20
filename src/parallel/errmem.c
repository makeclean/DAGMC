#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

extern void fail(const char *fmt, ...);

void fail(const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
  exit(1);
}

