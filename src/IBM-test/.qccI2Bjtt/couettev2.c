#line 0 "couettev2-cpp.c"
#line 1 "/home/spencer/basilisk/CTFL/src/IBM-test/.qccI2Bjtt//"
#line 0 "<built-in>"
#line 0 "<command-line>"
#line 1 "/usr/include/stdc-predef.h"
#line 0 "<command-line>"
#line 1 "couettev2-cpp.c"
#if 700 < 700
  #undef 700
  #define 700 700
#endif
#if 1
#include <stdint.h>
#include <string.h>
#include <fenv.h>
#endif


#line 1 "/home/spencer/basilisk/src/common.h"
#line 1 "/home/spencer/basilisk/src/ast/std/stdlib.h"
#include <stdlib.h>
#line 2 "/home/spencer/basilisk/src/common.h"
#line 1 "/home/spencer/basilisk/src/ast/std/stdio.h"
#include <stdio.h>
#line 3 "/home/spencer/basilisk/src/common.h"
#line 1 "/home/spencer/basilisk/src/ast/std/stddef.h"
#include <stddef.h>
#line 4 "/home/spencer/basilisk/src/common.h"
#line 1 "/home/spencer/basilisk/src/ast/std/stdbool.h"
#include <stdbool.h>
#line 5 "/home/spencer/basilisk/src/common.h"
#line 1 "/home/spencer/basilisk/src/ast/std/stdarg.h"
#include <stdarg.h>
#line 6 "/home/spencer/basilisk/src/common.h"
#line 1 "/home/spencer/basilisk/src/ast/std/string.h"
#include <string.h>
#line 7 "/home/spencer/basilisk/src/common.h"
#line 1 "/home/spencer/basilisk/src/ast/std/float.h"
#include <float.h>
#line 8 "/home/spencer/basilisk/src/common.h"
#line 1 "/home/spencer/basilisk/src/ast/std/limits.h"
#include <limits.h>
#line 9 "/home/spencer/basilisk/src/common.h"
#line 1 "/home/spencer/basilisk/src/ast/std/math.h"
#include <math.h>
#line 10 "/home/spencer/basilisk/src/common.h"
#line 1 "/home/spencer/basilisk/src/ast/std/time.h"
#include <time.h>
#line 11 "/home/spencer/basilisk/src/common.h"
#line 1 "/home/spencer/basilisk/src/ast/std/sys/time.h"
#include <sys/time.h>
#line 12 "/home/spencer/basilisk/src/common.h"
#line 1 "/home/spencer/basilisk/src/ast/std/sys/resource.h"
#include <sys/resource.h>
#line 13 "/home/spencer/basilisk/src/common.h"

#if _OPENMP
# include <omp.h>
# define OMP(x) _Pragma(#x)
#elif _MPI

# define OMP(x)

# include <mpi.h>
static int mpi_rank, mpi_npe;
# define tid() mpi_rank
# define pid() mpi_rank
# define npe() mpi_npe

#else

# define OMP(x)

#endif
#line 49 "/home/spencer/basilisk/src/common.h"
#define _NVARMAX 65536
#define is_constant(v) ((v).i >= _NVARMAX)
#define constant(v) (is_constant(v) ? _constant[(v).i - _NVARMAX] : 1e30)

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))
#define sq(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))
#define sign(x) ((x) > 0 ? 1 : -1)
#define sign2(x) ((x) > 0 ? 1 : (x) < 0 ? -1 : 0)
#define noise() (1. - 2.*rand()/(double)RAND_MAX)
#define clamp(x,a,b) ((x) < (a) ? (a) : (x) > (b) ? (b) : (x))

#define unmap(x,y)

#define trash(x)


#define systderr stderr
#define systdout stdout

#if _MPI
FILE * qstderr (void);
FILE * qstdout (void);
FILE * ferr = NULL, * fout = NULL;
#define not_mpi_compatible()\
do {\
  if (npe() > 1) {\
    fprintf (ferr, "%s() is not compatible with MPI (yet)\n", __func__);\
    exit (1);\
  }\
} while(0)\

#line 81

# define system(command) (pid() == 0 ? system(command) : 0)
#else
# define qstderr() stderr
# define qstdout() stdout
# define ferr stderr
# define fout stdout
# define not_mpi_compatible()
#endif



static inline void qassert (const char * file, int line, const char * cond) {
  fprintf (ferr, "%s:%d: Assertion `%s' failed.\n", file, line, cond);
  abort();
}
#line 105 "/home/spencer/basilisk/src/common.h"
#define sysmalloc malloc
#define syscalloc calloc
#define sysrealloc realloc
#define sysfree free
#define systrdup strdup

#if MTRACE

struct {
  FILE * fp;
  size_t total, max;
  size_t overhead, maxoverhead;
  size_t nr;
  size_t startrss, maxrss;
  char * fname;
} pmtrace;

typedef struct {
  char * func, * file;
  size_t max, total;
  int line, id;
} pmfunc;

typedef struct {
  size_t id, size;
} pmdata;

static pmfunc * pmfuncs = NULL;
static int pmfuncn = 0;

static int pmfunc_index (const char * func, const char * file, int line)
{
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn; i++, p++)
    if (p->line == line && !strcmp(func, p->func) && !strcmp(file, p->file))
      return p->id;
  pmfuncn++;
  pmfunc