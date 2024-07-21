# 0 "cylinder-osc-cpp.c"
# 0 "<built-in>"
# 0 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 0 "<command-line>" 2
# 1 "cylinder-osc-cpp.c"
@if 700 < 700
  @undef 700
  @define 700 700
@endif
@if 1
@include <stdint.h>
@include <string.h>
@include <fenv.h>
@endif


# 1 "/home/spencer/basilisk/src/common.h" 1
# 1 "/home/spencer/basilisk/src/ast/std/stdlib.h" 1
@include <stdlib.h>
# 2 "/home/spencer/basilisk/src/common.h" 2
# 1 "/home/spencer/basilisk/src/ast/std/stdio.h" 1
@include <stdio.h>
# 3 "/home/spencer/basilisk/src/common.h" 2
# 1 "/home/spencer/basilisk/src/ast/std/stddef.h" 1
@include <stddef.h>
# 4 "/home/spencer/basilisk/src/common.h" 2
# 1 "/home/spencer/basilisk/src/ast/std/stdbool.h" 1
@include <stdbool.h>
# 5 "/home/spencer/basilisk/src/common.h" 2
# 1 "/home/spencer/basilisk/src/ast/std/stdarg.h" 1
@include <stdarg.h>
# 6 "/home/spencer/basilisk/src/common.h" 2
# 1 "/home/spencer/basilisk/src/ast/std/string.h" 1
@include <string.h>
# 7 "/home/spencer/basilisk/src/common.h" 2
# 1 "/home/spencer/basilisk/src/ast/std/float.h" 1
@include <float.h>
# 8 "/home/spencer/basilisk/src/common.h" 2
# 1 "/home/spencer/basilisk/src/ast/std/limits.h" 1
@include <limits.h>
# 9 "/home/spencer/basilisk/src/common.h" 2
# 1 "/home/spencer/basilisk/src/ast/std/math.h" 1
@include <math.h>
# 10 "/home/spencer/basilisk/src/common.h" 2
# 1 "/home/spencer/basilisk/src/ast/std/time.h" 1
@include <time.h>
# 11 "/home/spencer/basilisk/src/common.h" 2
# 1 "/home/spencer/basilisk/src/ast/std/sys/time.h" 1
@include <sys/time.h>
# 12 "/home/spencer/basilisk/src/common.h" 2
# 1 "/home/spencer/basilisk/src/ast/std/sys/resource.h" 1
@include <sys/resource.h>
# 13 "/home/spencer/basilisk/src/common.h" 2

@if _OPENMP
@ include <omp.h>
@ define OMP(x) Pragma(#x)
@elif _MPI

@ define OMP(x)

@ include <mpi.h>
static int mpi_rank, mpi_npe;
@ define tid() mpi_rank
@ define pid() mpi_rank
@ define npe() mpi_npe

@else

@ define OMP(x)

@endif
# 49 "/home/spencer/basilisk/src/common.h"
@define _NVARMAX 65536
@define is_constant(v) ((v).i >= _NVARMAX)
@define constant(v) (is_constant(v) ? _constant[(v).i - _NVARMAX] : 1e30)

@define max(a,b) ((a) > (b) ? (a) : (b))
@define min(a,b) ((a) < (b) ? (a) : (b))
@define sq(x) ((x)*(x))
@define cube(x) ((x)*(x)*(x))
@define sign(x) ((x) > 0 ? 1 : -1)
@define sign2(x) ((x) > 0 ? 1 : (x) < 0 ? -1 : 0)
@define noise() (1. - 2.*rand()/(double)RAND_MAX)
@define clamp(x,a,b) ((x) < (a) ? (a) : (x) > (b) ? (b) : (x))

@define unmap(x,y)

@define trash(x)


@define systderr stderr
@define systdout stdout

@if _MPI
FILE * qstderr (void);
FILE * qstdout (void);
FILE * ferr = NULL, * fout = NULL;
@ def not_mpi_compatible()
do {
  if (npe() > 1) {
    fprintf (ferr, "%s() is not compatible with MPI (yet)\n", __func__);
    exit (1);
  }
} while(0)
@
@ define system(command) (pid() == 0 ? system(command) : 0)
@else
@ define qstderr() stderr
@ define qstdout() stdout
@ define ferr stderr
@ define fout stdout
@ define not_mpi_compatible()
@endif



static inline void qassert (const char * file, int line, const char * cond) {
  fprintf (ferr, "%s:%d: Assertion `%s' failed.\n", file, line, cond);
  abort();
}
# 105 "/home/spencer/basilisk/src/common.h"
@define sysmalloc malloc
@define syscalloc calloc
@define sysrealloc realloc
@define sysfree free
@define systrdup strdup

@if MTRACE

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
  pmfuncs = (pmfunc *) sysrealloc (pmfuncs, pmfuncn*sizeof(pmfunc));
  p = &pmfuncs[pmfuncn - 1];
  memset (p, 0, sizeof(pmfunc));
  p->func = systrdup(func);
  p->file = systrdup(file);
  p->line = line;
  p->id = pmfuncn;
  if (pmtrace.fp)
    fprintf (pmtrace.fp, "@ %d %s %s %d\n", pmfuncn, func, file, line);
  return pmfuncn;
}

static void pmfunc_trace (pmfunc * f, char c)
{
  if (pmtrace.fp)
    fprintf (pmtrace.fp, "%c %d %ld %ld %ld",
      c, f->id, pmtrace.nr, pmtrace.total, f->total);
@if 1
  if (pmtrace.nr % 1 == 0) {
    struct rusage usage;
    getrusage (RUSAGE_SELF, &usage);
    if (pmtrace.fp)
      fprintf (pmtrace.fp, " %ld", usage.ru_maxrss*1024);
    if (!pmtrace.nr)
      pmtrace.startrss = usage.ru_maxrss;
    if (usage.ru_maxrss > pmtrace.maxrss)
      pmtrace.maxrss = usage.ru_maxrss;
  }
@endif
  if (pmtrace.fp)
    fputc ('\n', pmtrace.fp);
  pmtrace.nr++;
}

static void * pmfunc_alloc (pmdata * d, size_t size,
       const char * func, const char * file, int line,
       char c)
{
  if (!(d != NULL)) qassert ("/home/spencer/basilisk/src/common.h", 0, "d != NULL");
  OMP (omp critical)
  {
    d->id = pmfunc_index(func, file, line);
    d->size = size;
    pmfunc * f = &pmfuncs[d->id - 1];
    f->total += size;
    if (f->total > f->max)
      f->max = f->total;
    pmtrace.total += size;
    pmtrace.overhead += sizeof(pmdata);
    if (pmtrace.total > pmtrace.max) {
      pmtrace.max = pmtrace.total;
      pmtrace.maxoverhead = pmtrace.overhead;
    }
    pmfunc_trace (f, c);
  }
  return ((char *)d) + sizeof(pmdata);
}

static void * pmfunc_free (void * ptr, char c)
{
  if (!ptr)
    return ptr;
  pmdata * d = (pmdata *) (((char *)ptr) - sizeof(pmdata));
  if (d->id < 1 || d->id > pmfuncn) {
    fputs ("*** MTRACE: ERROR!: corrupted free()", ferr);
    if (d->size == 0)
      fputs (", possible double free()", ferr);
    else
      fputs (", not traced?", ferr);
    fputs (", aborting...\n", ferr);
    abort();
    return ptr;
  }
  else
  OMP (omp critical)
  {
    pmfunc * f = &pmfuncs[d->id - 1];
    if (f->total < d->size) {
      fprintf (ferr, "*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\n",
        f->total, d->size);
      abort();
    }
    else
      f->total -= d->size;
    if (pmtrace.total < d->size) {
      fprintf (ferr, "*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\n",
        pmtrace.total, d->size);
      abort();
    }
    else {
      pmtrace.total -= d->size;
      pmtrace.overhead -= sizeof(pmdata);
    }
    d->id = 0;
    d->size = 0;
    pmfunc_trace (f, c);
  }
  return d;
}

static void * pmalloc (size_t size,
         const char * func, const char * file, int line)
{
  return pmfunc_alloc ((pmdata *) sysmalloc (sizeof(pmdata) + size),
         size, func, file, line, '+');
}

static void * pcalloc (size_t nmemb, size_t size,
         const char * func, const char * file, int line)
{
  void * p = pmalloc (nmemb*size, func, file, line);
  return memset (p, 0, nmemb*size);
}

static void * prealloc (void * ptr, size_t size,
   const char * func, const char * file, int line)
{
  return pmfunc_alloc ((pmdata *) sysrealloc (pmfunc_free(ptr, '<'),
           sizeof(pmdata) + size),
         size, func, file, line, '>');
}

static void pfree (void * ptr,
     const char * func, const char * file, int line)
{
  sysfree (pmfunc_free (ptr, '-'));
}

static char * pstrdup (const char * s,
         const char * func, const char * file, int line)
{
  char * d = (char *) pmalloc (strlen(s) + 1, func, file, line);
  return strcpy (d, s);
}

@if MTRACE < 3
static int pmaxsort (const void * a, const void * b) {
  const pmfunc * p1 = a, * p2 = b;
  return p1->max < p2->max;
}
@endif

static int ptotalsort (const void * a, const void * b) {
  const pmfunc * p1 = (const pmfunc *) a, * p2 = (const pmfunc *) b;
  return p1->total < p2->total;
}

static void pmfuncs_free()
{
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn; i++, p++) {
    sysfree (p->func);
    sysfree (p->file);
  }
  sysfree (pmfuncs);
}

void pmuntrace (void)
{
@if MTRACE < 3
  fprintf (ferr,
    "*** MTRACE: max resident  set size: %10ld bytes\n"
    "*** MTRACE: max traced memory size: %10ld bytes"
    " (tracing overhead %.1g%%)\n"
    "%10s    %20s   %s\n",
    pmtrace.maxrss*1024,
    pmtrace.max, pmtrace.maxoverhead*100./pmtrace.max,
    "max bytes", "function", "file");
  qsort (pmfuncs, pmfuncn, sizeof(pmfunc), pmaxsort);
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn && p->max > 0; i++, p++)
    fprintf (ferr, "%10ld    %20s   %s:%d\n",
      p->max, p->func, p->file, p->line);

  if (pmtrace.fp) {
    char * fname = pmtrace.fname, * s;
    while ((s = strchr(fname,'/')))
      fname = s + 1;

    fputs ("load(\"`echo $BASILISK`/mtrace.plot\")\n", pmtrace.fp);
    fprintf (pmtrace.fp,
      "plot '%s' u 3:($6-%g) w l t 'ru_maxrss - %.3g',"
      "total(\"%s\") w l t 'total'",
      fname,
      pmtrace.startrss*1024.,
      pmtrace.startrss*1024.,
      fname);
    pmfunc * p = pmfuncs;
    for (int i = 0; i < pmfuncn && p->max > 0.01*pmtrace.max; i++, p++)
      fprintf (pmtrace.fp,
        ",func(\"%s\",%d) w l t '%s'",
        fname, p->id, p->func);
    fputc ('\n', pmtrace.fp);
    fprintf (ferr,
      "*** MTRACE: To get a graph use: tail -n 2 %s | gnuplot -persist\n",
      fname);
    fclose (pmtrace.fp);
    pmtrace.fp = NULL;
    sysfree (pmtrace.fname);
  }
@endif

  if (pmtrace.total > 0) {
    qsort (pmfuncs, pmfuncn, sizeof(pmfunc), ptotalsort);
    pmfunc * p = pmfuncs;
    for (int i = 0; i < pmfuncn && p->total > 0; i++, p++)
      fprintf (ferr, "%s:%d: error: %ld bytes leaked here\n",
        p->file, p->line, p->total);
    pmfuncs_free();
    exit(1);
  }
  else {
@if MTRACE < 3
    fputs ("*** MTRACE: No memory leaks\n", ferr);
@endif
    pmfuncs_free();
  }
}

@else
@ define pmalloc(s,func,file,line) malloc(s)
@ define pcalloc(n,s,func,file,line) calloc(n,s)
@ define prealloc(p,s,func,file,line) realloc(p,s)
@ define pfree(p,func,file,line) free(p)
@ define pstrdup(s,func,file,line) strdup(s)
@endif







typedef struct {
  void * p;
  long max, len;
} Array;

Array * array_new()
{
  Array * a = ((Array *) pmalloc ((1)*sizeof(Array),__func__,__FILE__,0));
  a->p = NULL;
  a->max = a->len = 0;
  return a;
}

void array_free (Array * a)
{
  pfree (a->p,__func__,__FILE__,0);
  pfree (a,__func__,__FILE__,0);
}

void array_append (Array * a, void * elem, size_t size)
{
  if (a->len + size >= a->max) {
    a->max += max (size, 4096);
    a->p = prealloc (a->p, a->max,__func__,__FILE__,0);
  }
  memcpy (((char *)a->p) + a->len, elem, size);
  a->len += size;
}

void * array_shrink (Array * a)
{
  void * p = prealloc (a->p, a->len,__func__,__FILE__,0);
  pfree (a,__func__,__FILE__,0);
  return p;
}



@if TRACE == 1
@include <extrae_user_events.h>

typedef struct {
  Array index, stack;
  extrae_type_t type;
} Trace;

Trace trace_func = {
  {NULL, 0, 0}, {NULL, 0, 0},
  60000010,
};

Trace trace_mpi_func = {
  {NULL, 0, 0}, {NULL, 0, 0},
  60000011,
};

static int lookup_func (Array * a, const char * func)
{
  for (int i = 0; i < a->len/sizeof(char *); i++) {
    char * s = ((char **)a->p)[i];
    if (!strcmp (func, s))
      return i + 1;
  }
  char * s = pstrdup (func,__func__,__FILE__,0);
  array_append (a, &s, sizeof(char *));
  return a->len;
}

static void trace_push (Trace * t, const char * func)
{
  int value = lookup_func (&t->index, func);
  Extrae_eventandcounters (t->type, value);
  array_append (&t->stack, &value, sizeof(int));
}

static void trace_pop (Trace * t, const char * func)
{
  if (!(t->stack.len > 0)) qassert ("/home/spencer/basilisk/src/common.h", 0, "t->stack.len > 0");
  t->stack.len -= sizeof(int);
  int value = t->stack.len > 0 ?
    ((int *)t->stack.p)[t->stack.len/sizeof(int) - 1] : 0;
  Extrae_eventandcounters (t->type, value);
}

static void trace_define (Trace * t, char * description)
{
  if (t->index.len > 0) {
    extrae_value_t values[t->index.len/sizeof(char *) + 1];
    char * names[t->index.len/sizeof(char *) + 1],
      ** func = (char **) t->index.p;
    names[0] = "OTHER";
    values[0] = 0;
    unsigned len = 1;
    for (int i = 0; i < t->index.len/sizeof(char *); i++, func++) {
      names[len] = *func;
      values[len++] = i + 1;
    }
    Extrae_define_event_type (&t->type, description, &len, values, names);
  }
}

static void trace_free (Trace * t)
{
  char ** func = (char **) t->index.p;
  for (int i = 0; i < t->index.len/sizeof(char *); i++, func++)
    pfree (*func,__func__,__FILE__,0);
  pfree (t->index.p,__func__,__FILE__,0);
  pfree (t->stack.p,__func__,__FILE__,0);
}

static void trace_off()
{
  trace_define (&trace_func, "Basilisk functions");
  trace_define (&trace_mpi_func, "Basilisk functions (MPI-related)");
  trace_free (&trace_func);
  trace_free (&trace_mpi_func);
}






@ define tracing(func, file, line) trace_push (&trace_func, func)
@ define end_tracing(func, file, line) trace_pop (&trace_func, func)

@elif TRACE

typedef struct {
  char * func, * file;
  int line, calls;
  double total, self;
@if _MPI
  double min, max;
@endif
} TraceIndex;

struct {
  Array stack, index;
  double t0;
} Trace = {
  {NULL, 0, 0}, {NULL, 0, 0},
  -1
};

static void trace_add (const char * func, const char * file, int line,
         double total, double self)
{
  TraceIndex * t = (TraceIndex *) Trace.index.p;
  int i, len = Trace.index.len/sizeof(TraceIndex);
  for (i = 0; i < len; i++, t++)
    if (t->line == line && !strcmp (func, t->func) && !strcmp (file, t->file))
      break;
  if (i == len) {
    TraceIndex t = {pstrdup(func,__func__,__FILE__,0), pstrdup(file,__func__,__FILE__,0), line, 1, total, self};
    array_append (&Trace.index, &t, sizeof(TraceIndex));
  }
  else
    t->calls++, t->total += total, t->self += self;
}

static void tracing (const char * func, const char * file, int line)
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  if (Trace.t0 < 0)
    Trace.t0 = tv.tv_sec + tv.tv_usec/1e6;
  double t[2] = {(tv.tv_sec - Trace.t0) + tv.tv_usec/1e6, 0.};
  array_append (&Trace.stack, t, 2*sizeof(double));




}

static void end_tracing (const char * func, const char * file, int line)
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  double te = (tv.tv_sec - Trace.t0) + tv.tv_usec/1e6;
  double * t = (double *) Trace.stack.p;
  if (!(Trace.stack.len >= 2*sizeof(double))) qassert ("/home/spencer/basilisk/src/common.h", 0, "Trace.stack.len >= 2*sizeof(double)");
  t += Trace.stack.len/sizeof(double) - 2;
  Trace.stack.len -= 2*sizeof(double);
  double dt = te - t[0];




  trace_add (func, file, line, dt, dt - t[1]);
  if (Trace.stack.len >= 2*sizeof(double)) {
    t -= 2;
    t[1] += dt;
  }
}

static int compar_self (const void * p1, const void * p2)
{
  const TraceIndex * t1 = p1, * t2 = p2;
  return t1->self < t2->self;
}

@if _MPI
static int compar_func (const void * p1, const void * p2)
{
  const TraceIndex * t1 = p1, * t2 = p2;
  if (t1->line != t2->line)
    return t1->line < t2->line;
  return strcmp (t1->file, t2->file);
}
@endif

void trace_print (FILE * fp, double threshold)
{
  int i, len = Trace.index.len/sizeof(TraceIndex);
  double total = 0.;
  TraceIndex * t;
  Array * index = array_new();
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    array_append (index, t, sizeof(TraceIndex)), total += t->self;
@if _MPI
  qsort (index->p, len, sizeof(TraceIndex), compar_func);
  double tot[len], self[len], min[len], max[len];
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    tot[i] = t->total, self[i] = t->self;
  MPI_Reduce (self, min, len, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce (self, max, len, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce (pid() ? self : MPI_IN_PLACE,
       self, len, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (pid() ? tot : MPI_IN_PLACE,
       tot, len, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  total = 0.;
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    t->total = tot[i]/npe(), t->self = self[i]/npe(),
      t->max = max[i], t->min = min[i], total += t->self;
@endif
  qsort (index->p, len, sizeof(TraceIndex), compar_self);
  fprintf (fp, "   calls    total     self   %% total   function\n");
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    if (t->self*100./total > threshold) {
      fprintf (fp, "%8d   %6.2f   %6.2f     %4.1f%%",
        t->calls, t->total, t->self, t->self*100./total);
@if _MPI
      fprintf (fp, " (%4.1f%% - %4.1f%%)", t->min*100./total, t->max*100./total);
@endif
      fprintf (fp, "   %s():%s:%d\n", t->func, t->file, t->line);
    }
  fflush (fp);
  array_free (index);
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    t->calls = t->total = t->self = 0.;
}

static void trace_off()
{
  trace_print (fout, 0.);

  int i, len = Trace.index.len/sizeof(TraceIndex);
  TraceIndex * t;
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    pfree (t->func,__func__,__FILE__,0), pfree (t->file,__func__,__FILE__,0);

  pfree (Trace.index.p,__func__,__FILE__,0);
  Trace.index.p = NULL;
  Trace.index.len = Trace.index.max = 0;

  pfree (Trace.stack.p,__func__,__FILE__,0);
  Trace.stack.p = NULL;
  Trace.stack.len = Trace.stack.max = 0;
}

@else
@ define tracing(...)
@ define end_tracing(...)
@endif



@if _OPENMP

@define tid() omp_get_thread_num()
@define pid() 0
@define npe() omp_get_num_threads()
@define mpi_all_reduce(v,type,op)
@define mpi_all_reduce_array(v,type,op,elem)

@elif _MPI

static bool in_prof = false;
static double prof_start, _prof;
@def prof_start(name)
  if (!(!in_prof)) qassert ("/home/spencer/basilisk/src/common.h", 0, "!in_prof"); in_prof = true;
  prof_start = MPI_Wtime();
@
@def prof_stop()
  if (!(in_prof)) qassert ("/home/spencer/basilisk/src/common.h", 0, "in_prof"); in_prof = false;
  _prof = MPI_Wtime();
  mpi_time += _prof - prof_start;
@

@if FAKE_MPI
@define mpi_all_reduce(v,type,op)
@define mpi_all_reduce_array(v,type,op,elem)
@else
     
int mpi_all_reduce0 (void *sendbuf, void *recvbuf, int count,
       MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{tracing("mpi_all_reduce0","/home/spencer/basilisk/src/common.h",0);
  { int _ret= MPI_Allreduce (sendbuf, recvbuf, count, datatype, op, comm);end_tracing("mpi_all_reduce0","/home/spencer/basilisk/src/common.h",0);return _ret;}
end_tracing("mpi_all_reduce0","/home/spencer/basilisk/src/common.h",0);}
@def mpi_all_reduce(v,type,op) {
  prof_start ("mpi_all_reduce");
  union { int a; float b; double c;} global;
  mpi_all_reduce0 (&(v), &global, 1, type, op, MPI_COMM_WORLD);
  memcpy (&(v), &global, sizeof (v));
  prof_stop();
}
@
@def mpi_all_reduce_array(v,type,op,elem) {
  prof_start ("mpi_all_reduce");
  type * global = malloc ((elem)*sizeof(type)), * tmp = malloc ((elem)*sizeof(type));
  for (int i = 0; i < elem; i++)
    tmp[i] = (v)[i];
  MPI_Datatype datatype;
  if (!strcmp(#type, "double")) datatype = MPI_DOUBLE;
  else if (!strcmp(#type, "int")) datatype = MPI_INT;
  else if (!strcmp(#type, "long")) datatype = MPI_LONG;
  else if (!strcmp(#type, "bool")) datatype = MPI_C_BOOL;
  else if (!strcmp(#type, "unsigned char")) datatype = MPI_UNSIGNED_CHAR;
  else {
    fprintf (stderr, "unknown reduction type '%s'\n", #type);
    fflush (stderr);
    abort();
  }
  mpi_all_reduce0 (tmp, global, elem, datatype, op, MPI_COMM_WORLD);
  for (int i = 0; i < elem; i++)
    (v)[i] = global[i];
  free (global), free (tmp);
  prof_stop();
}
@

@endif

@define QFILE FILE

FILE * qstderr (void)
{
  static QFILE * fp = NULL;
  if (!fp) {
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "log-%d", mpi_rank);
      fp = fopen (name, "w");
    }
    else
      fp = systderr;
  }
  return fp;
}

FILE * qstdout (void)
{
  static QFILE * fp = NULL;
  if (!fp) {
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "out-%d", mpi_rank);
      fp = fopen (name, "w");
    }
    else
      fp = systdout;
  }
  return fp;
}

static void finalize (void)
{
  MPI_Finalize();
}

void mpi_init()
{
  int initialized;
  MPI_Initialized (&initialized);
  if (!initialized) {
    MPI_Init (NULL, NULL);
    MPI_Comm_set_errhandler (MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
    atexit (finalize);
  }
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &mpi_npe);
  srand (mpi_rank + 1);
  if (ferr == NULL) {
    if (mpi_rank > 0) {
      ferr = fopen ("/dev/null", "w");
      fout = fopen ("/dev/null", "w");
    }
    else {
      ferr = systderr;
      fout = systdout;
    }
    char * etrace = getenv ("MALLOC_TRACE"), name[80];
    if (etrace && mpi_rank > 0) {
      sprintf (name, "%s-%d", etrace, mpi_rank);
      setenv ("MALLOC_TRACE", name, 1);
    }
@if MTRACE == 1
    etrace = getenv ("MTRACE");
    if (!etrace)
      etrace = "mtrace";
    if (mpi_rank > 0) {
      sprintf (name, "%s-%d", etrace, mpi_rank);
      pmtrace.fp = fopen (name, "w");
      pmtrace.fname = systrdup(name);
    }
    else {
      pmtrace.fp = fopen (etrace, "w");
      pmtrace.fname = systrdup(etrace);
    }
@endif
  }
}

@else

@define tid() 0
@define pid() 0
@define npe() 1
@define mpi_all_reduce(v,type,op)
@define mpi_all_reduce_array(v,type,op,elem)

@endif

@define OMP_PARALLEL() OMP(omp parallel)

@define NOT_UNUSED(x) (void)(x)

@define VARIABLES ;
@define _index(a,m) (a.i)
@define val(a,k,l,m) data(k,l,m)[_index(a,m)]

double _val_higher_dimension = 0.;
# 826 "/home/spencer/basilisk/src/common.h"
@if (1 || __APPLE__) && !_OPENMP && !_CADNA
double undefined;
@ if __APPLE__
@ include <stdint.h>
@ include "fp_osx.h"
@ endif
@ define enable_fpe(flags) feenableexcept (flags)
@ define disable_fpe(flags) fedisableexcept (flags)
static void set_fpe (void) {
  int64_t lnan = 0x7ff0000000000001;
  if (!(sizeof (int64_t) == sizeof (double))) qassert ("/home/spencer/basilisk/src/common.h", 0, "sizeof (int64_t) == sizeof (double)");
  memcpy (&undefined, &lnan, sizeof (double));
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}
@else
@ define undefined ((double) DBL_MAX)
@ define enable_fpe(flags)
@ define disable_fpe(flags)
static void set_fpe (void) {}
@endif


typedef struct {
  long n;
  long tn;
  int depth;
  int maxdepth;
} Grid;
Grid * grid = NULL;

double X0 = 0., Y0 = 0., Z0 = 0.;

double L0 = 1.;


int N = 64;




typedef struct { int i; } scalar;

typedef struct {
  scalar x;

  scalar y;




} vector;

typedef struct {
  scalar * x;

  scalar * y;




} vectorl;

typedef struct {
  vector x;

  vector y;




} tensor;

struct { int x, y, z; } Period = {false, false, false};

typedef struct {
  double x, y, z;
} coord;

OMP(omp declare reduction (+ : coord :
      omp_out.x += omp_in.x,
      omp_out.y += omp_in.y,
      omp_out.z += omp_in.z))
# 920 "/home/spencer/basilisk/src/common.h"
void normalize (coord * n)
{
  double norm = 0.;
  
    norm += sq(n->x);
    
#line 924
norm += sq(n->y);
  norm = sqrt(norm);
  
    n->x /= norm;
    
#line 927
n->y /= norm;
}

void origin (double x, double y, double z) {
  X0 = x; Y0 = y; Z0 = z;
}

void size (double L) {
  L0 = L;
}

double zero (double s0, double s1, double s2) { return 0.; }






  enum { right, left, top, bottom };



int nboundary = 2*2;



@define _dirichlet(expr, ...) (2.*(expr) - val(_s,0,0,0))
@define _dirichlet_homogeneous(...) (- val(_s,0,0,0))
@define _dirichlet_face(expr,...) (expr)
@define _dirichlet_face_homogeneous(...) (0.)
@define _neumann(expr,...) (Delta*(expr) + val(_s,0,0,0))
@define _neumann_homogeneous(...) (val(_s,0,0,0))

double * _constant = NULL;
size_t datasize = 0;
typedef struct _Point Point;

# 1 "/home/spencer/basilisk/src/grid/boundaries.h" 1


typedef struct _Boundary Boundary;

struct _Boundary {
  void (* destroy) (Boundary * b);
  void (* level) (const Boundary * b, scalar * list, int l);

  void (* restriction) (const Boundary * b, scalar * list, int l);
};

static Boundary ** boundaries = NULL;

void add_boundary (Boundary * b) {
  int len = 0;
  if (boundaries) {
    Boundary ** i = boundaries;
    while (*i++) len++;
  }
  boundaries = (Boundary * *) prealloc (boundaries, (len + 2)*sizeof(Boundary *),__func__,__FILE__,0);
  boundaries[len] = b;
  boundaries[len+1] = NULL;
}

void free_boundaries() {
  if (!boundaries)
    return;
  Boundary ** i = boundaries, * b;
  while ((b = *i++))
    if (b->destroy)
      b->destroy (b);
    else
      pfree (b,__func__,__FILE__,0);
  pfree (boundaries,__func__,__FILE__,0);
  boundaries = NULL;
}
# 47 "/home/spencer/basilisk/src/grid/boundaries.h"
typedef struct {
  Boundary parent;
  int d;
} BoxBoundary;
# 965 "/home/spencer/basilisk/src/common.h" 2



typedef struct {
  double (** boundary) (Point, Point, scalar, void *);
  double (** boundary_homogeneous) (Point, Point, scalar, void *);
  double (* gradient) (double, double, double);
  void (* delete) (scalar);
  char * name;
  struct {
    int x;

    int y;




  } d;
  vector v;
  int face;
  bool nodump, freed;
  int block;
  scalar * depends;

  
#line 19 "/home/spencer/basilisk/src/grid/stencils.h"
bool input, output;
  int width;
  int dirty;
  
#line 18 "/home/spencer/basilisk/src/grid/multigrid-common.h"
void (* prolongation) (Point, scalar);
  void (* restriction) (Point, scalar);
  
#line 9 "/home/spencer/basilisk/src/grid/tree-common.h"
void (* refine) (Point, scalar);
  
#line 97
void (* coarsen) (Point, scalar);
  
#line 81 "/home/spencer/basilisk/src/fractions.h"
vector n;

#line 988 "/home/spencer/basilisk/src/common.h"
} _Attributes;

static _Attributes * _attribute = NULL;






int list_len (scalar * list)
{
  if (!list) return 0;
  int ns = 0;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ ns++;}}
  return ns;
}

scalar * list_append (scalar * list, scalar s)
{
  int len = list_len (list);
  list = (scalar *) prealloc (list, (len + 2)*sizeof(scalar),__func__,__FILE__,0);
  list[len] = s;
  list[len + 1].i = -1;
  return list;
}

scalar * list_prepend (scalar * list, scalar s)
{
  int len = list_len (list);
  list = (scalar *) prealloc (list, (len + 2)*sizeof(scalar),__func__,__FILE__,0);
  for (int i = len; i >= 1; i--)
    list[i] = list[i-1];
  list[0] = s;
  list[len + 1].i = -1;
  return list;
}

scalar * list_add (scalar * list, scalar s)
{
  {scalar*_i=(scalar*)( list);if(_i)for(scalar t=*_i;(&t)->i>=0;t=*++_i){
    if (t.i == s.i)
      return list;}}
  return list_append (list, s);
}

int list_lookup (scalar * l, scalar s)
{
  if (l != NULL)
    {scalar*_i=(scalar*)( l);if(_i)for(scalar s1=*_i;(&s1)->i>=0;s1=*++_i){
      if (s1.i == s.i)
 return true;}}
  return false;
}

scalar * list_copy (scalar * l)
{
  scalar * list = NULL;
  if (l != NULL)
    {scalar*_i=(scalar*)( l);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      list = list_append (list, s);}}
  return list;
}

scalar * list_concat (scalar * l1, scalar * l2)
{
  scalar * l3 = list_copy (l1);
  {scalar*_i=(scalar*)( l2);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    l3 = list_append (l3, s);}}
  return l3;
}

void list_print (scalar * l, FILE * fp)
{
  int i = 0;
  {scalar*_i=(scalar*)( l);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    fprintf (fp, "%s%s", i++ == 0 ? "{" : ",", _attribute[s.i].name);}}
  fputs (i > 0 ? "}\n" : "{}\n", fp);
}

int vectors_len (vector * list)
{
  if (!list) return 0;
  int nv = 0;
  {vector*_i=(vector*)( list);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){ nv++;}}
  return nv;
}

vector * vectors_append (vector * list, vector v)
{
  int len = vectors_len (list);
  list = (vector *) prealloc (list, (len + 2)*sizeof(vector),__func__,__FILE__,0);
  list[len] = v;
  list[len + 1] = (vector){{-1}};
  return list;
}

vector * vectors_add (vector * list, vector v)
{
  {vector*_i=(vector*)( list);if(_i)for(vector w=*_i;(&w)->x.i>=0;w=*++_i){ {
    bool id = true;
    
      if (w.x.i != v.x.i)
 id = false;
      
#line 1089
if (w.y.i != v.y.i)
 id = false;
    if (id)
      return list;
  }}}
  return vectors_append (list, v);
}

vector * vectors_copy (vector * l)
{
  vector * list = NULL;
  if (l != NULL)
    {vector*_i=(vector*)( l);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
      list = vectors_append (list, v);}}
  return list;
}

vector * vectors_from_scalars (scalar * s)
{
  vector * list = NULL;
  while (s->i >= 0) {
    vector v;
     {
      if (!(s->i >= 0)) qassert ("/home/spencer/basilisk/src/common.h", 0, "s->i >= 0");
      v.x = *s++;
    } 
#line 1111
{
      if (!(s->i >= 0)) qassert ("/home/spencer/basilisk/src/common.h", 0, "s->i >= 0");
      v.y = *s++;
    }
    list = vectors_append (list, v);
  }
  return list;
}

int tensors_len (tensor * list)
{
  if (!list) return 0;
  int nt = 0;
  {tensor*_i=(tensor*)( list);if(_i)for(tensor t=*_i;(&t)->x.x.i>=0;t=*++_i){ nt++;}}
  return nt;
}

tensor * tensors_append (tensor * list, tensor t)
{
  int len = tensors_len (list);
  list = (tensor *) prealloc (list, (len + 2)*sizeof(tensor),__func__,__FILE__,0);
  list[len] = t;
  list[len + 1] = (tensor){{{-1}}};
  return list;
}

tensor * tensors_from_vectors (vector * v)
{
  tensor * list = NULL;
  while (v->x.i >= 0) {
    tensor t;
     {
      if (!(v->x.i >= 0)) qassert ("/home/spencer/basilisk/src/common.h", 0, "v->x.i >= 0");
      t.x = *v++;
    } 
#line 1142
{
      if (!(v->y.i >= 0)) qassert ("/home/spencer/basilisk/src/common.h", 0, "v->x.i >= 0");
      t.y = *v++;
    }
    list = tensors_append (list, t);
  }
  return list;
}

static inline bool is_vertex_scalar (scalar s)
{
  
    if (_attribute[s.i].d.x != -1)
      return false;
    
#line 1154
if (_attribute[s.i].d.y != -1)
      return false;
  return true;
}

scalar * all = NULL;
scalar * baseblock = NULL;



scalar (* init_scalar) (scalar, const char *);
scalar (* init_vertex_scalar) (scalar, const char *);
vector (* init_vector) (vector, const char *);
vector (* init_face_vector) (vector, const char *);
tensor (* init_tensor) (tensor, const char *);
void (* scalar_clone) (scalar, scalar);





typedef struct _Event Event;
typedef int (* Expr) (int *, double *, Event *);

struct _Event {
  int last, nexpr;
  int (* action) (const int, const double, Event *);
  Expr expr[3];
  int * arrayi;
  double * arrayt;
  char * file;
  int line;
  char * name;
  double t;
  int i, a;
  void * data;
  Event * next;
};

static Event * Events = NULL;

int iter = 0, inext = 0;
double t = 0, tnext = 0;
void init_events (void);
void event_register (Event event);
static void _init_solver (void);

void init_solver()
{
  Events = pmalloc (sizeof (Event),__func__,__FILE__,0);
  Events[0].last = 1;
  _attribute = pcalloc (datasize/sizeof(double), sizeof (_Attributes),__func__,__FILE__,0);
  int n = datasize/sizeof(double);
  all = (scalar *) pmalloc (sizeof (scalar)*(n + 1),__func__,__FILE__,0);
  baseblock = (scalar *) pmalloc (sizeof (scalar)*(n + 1),__func__,__FILE__,0);
  for (int i = 0; i < n; i++)
    baseblock[i].i = all[i].i = i;
  baseblock[n].i = all[n].i = -1;
@if _CADNA
  cadna_init (-1);
@endif
@if _MPI
  mpi_init();
@elif MTRACE == 1
  char * etrace = getenv ("MTRACE");
  pmtrace.fp = fopen (etrace ? etrace : "mtrace", "w");
  pmtrace.fname = systrdup (etrace ? etrace : "mtrace");
@endif
}



@if _MPI
static double mpi_time = 0.;
@endif

typedef struct {
  clock_t c;
  struct timeval tv;
  double tm;
} timer;

timer timer_start (void)
{
  timer t;
  t.c = clock();
  gettimeofday (&t.tv, NULL);
@if _MPI
  t.tm = mpi_time;
@endif
  return t;
}

double timer_elapsed (timer t)
{
  struct timeval tvend;
  gettimeofday (&tvend, NULL);
  return ((tvend.tv_sec - t.tv.tv_sec) +
   (tvend.tv_usec - t.tv.tv_usec)/1e6);
}



const vector zerof = {{_NVARMAX+0},{_NVARMAX+1}};
const vector unityf = {{_NVARMAX+2},{_NVARMAX+3}};
const scalar unity = {_NVARMAX+4};
const scalar zeroc = {_NVARMAX+5};



const vector unityf0 = {{_NVARMAX+6},{_NVARMAX+7}};
const scalar unity0 = {_NVARMAX+8};
        vector fm = {{_NVARMAX+6},{_NVARMAX+7}};
        scalar cm = {_NVARMAX+8};
# 1281 "/home/spencer/basilisk/src/common.h"
static FILE ** qpopen_pipes = NULL;

FILE * qpopen (const char * command, const char * type)
{
  if (pid() > 0)
    return fopen ("/dev/null", type);
  FILE * fp = popen (command, type);
  if (fp) {
    FILE ** i = qpopen_pipes;
    int n = 0;
    while (i && *i) { n++; i++; }
    qpopen_pipes = (FILE * *) prealloc (qpopen_pipes, (n + 2)*sizeof(FILE *),__func__,__FILE__,0);
    qpopen_pipes[n] = fp;
    qpopen_pipes[n+1] = NULL;
  }
  return fp;
}

int qpclose (FILE * fp)
{
  if (pid() > 0)
    return fclose (fp);
  FILE ** i = qpopen_pipes;
  while (i && *i) {
    if (*i == fp)
      *i = (FILE *) 1;
    i++;
  }
  return pclose (fp);
}

static void qpclose_all()
{
  FILE ** i = qpopen_pipes;
  while (i && *i) {
    if (*i != (FILE *) 1)
      pclose (*i);
    i++;
  }
  pfree (qpopen_pipes,__func__,__FILE__,0);
  qpopen_pipes = NULL;
}






FILE * lfopen (const char * name, const char * mode)
{
  char fname[80];
  sprintf (fname, "%s-%d", name, pid());
  return fopen (fname, mode);
}



void * matrix_new (int n, int p, size_t size)
{
  void ** m = ((void * *) pmalloc ((n)*sizeof(void *),__func__,__FILE__,0));
  char * a = ((char *) pmalloc ((n*p*size)*sizeof(char),__func__,__FILE__,0));
  for (int i = 0; i < n; i++)
    m[i] = a + i*p*size;
  return m;
}

double matrix_inverse (double ** m, int n, double pivmin)
{
  int indxc[n], indxr[n], ipiv[n];
  int i, icol = 0, irow = 0, j, k, l, ll;
  double big, dum, pivinv, minpiv = 1e30;

  for (j = 0; j < n; j++)
    ipiv[j] = -1;

  for (i = 0; i < n; i++) {
    big = 0.0;
    for (j = 0; j < n; j++)
      if (ipiv[j] != 0)
 for (k = 0; k < n; k++) {
   if (ipiv[k] == -1) {
     if (fabs (m[j][k]) >= big) {
       big = fabs (m[j][k]);
       irow = j;
       icol = k;
     }
   }
 }
    ipiv[icol]++;
    if (irow != icol)
      for (l = 0; l < n; l++)
 do { double __tmp = m[irow][l]; m[irow][l] = m[icol][l]; m[icol][l] = __tmp; } while(0);
    indxr[i] = irow;
    indxc[i] = icol;
    if (fabs (m[icol][icol]) <= pivmin)
      return 0.;
    if (fabs (m[icol][icol]) < minpiv)
      minpiv = fabs (m[icol][icol]);
    pivinv = 1.0/m[icol][icol];
    m[icol][icol] = 1.0;
    for (l = 0; l < n; l++) m[icol][l] *= pivinv;
    for (ll = 0; ll < n; ll++)
      if (ll != icol) {
 dum = m[ll][icol];
 m[ll][icol] = 0.0;
 for (l = 0; l < n; l++)
   m[ll][l] -= m[icol][l]*dum;
      }
  }
  for (l = n - 1; l >= 0; l--) {
    if (indxr[l] != indxc[l])
      for (k = 0; k < n; k++)
 do { double __tmp = m[k][indxr[l]]; m[k][indxr[l]] = m[k][indxc[l]]; m[k][indxc[l]] = __tmp; } while(0);
  }
  return minpiv;
}

void matrix_free (void * m)
{
  pfree (((void **) m)[0],__func__,__FILE__,0);
  pfree (m,__func__,__FILE__,0);
}



typedef void (* free_solver_func) (void);

static Array * free_solver_funcs = NULL;

void free_solver_func_add (free_solver_func func)
{
  if (!free_solver_funcs)
    free_solver_funcs = array_new();
  array_append (free_solver_funcs, &func, sizeof(free_solver_func));
}



static char * display_defaults = NULL;

static void free_display_defaults() {
  pfree (display_defaults,__func__,__FILE__,0);
}

void display (const char * commands, bool overwrite)
{
  if (display_defaults == NULL)
    free_solver_func_add (free_display_defaults);
  if (overwrite) {
    pfree (display_defaults,__func__,__FILE__,0);
    display_defaults = pmalloc (strlen(commands) + 2,__func__,__FILE__,0);
    strcpy (display_defaults, "@");
    strcat (display_defaults, commands);
  }
  else {
    if (!display_defaults)
      display_defaults = pstrdup ("@",__func__,__FILE__,0);
    display_defaults =
      prealloc (display_defaults,
        strlen(display_defaults) + strlen(commands) + 1,__func__,__FILE__,0);
    strcat (display_defaults, commands);
  }
}



typedef struct {
  double x;

  double y;




} _coord;





# 1 "/home/spencer/basilisk/src/grid/stencils.h" 1
# 17 "/home/spencer/basilisk/src/grid/stencils.h"










typedef struct {
  char * name;
  char * type;
  void * pointer;
  int * dimensions;
  int nd;
  char reduct;
  scalar data;
} NonLocal;

typedef struct {
  const char * fname;
  int line;
  int first;
  int face;
  bool vertex;
  int parallel;
  scalar * listc;
  vectorl listf;
  scalar * dirty;
  void * data;
} ForeachData;


@def foreach_stencil(...) {
  static ForeachData _loop = {
    .fname = S__FILE__, .line = S_LINENO, .first = 1
  };
  if (baseblock) for (scalar s = baseblock[0], * i = baseblock;
  s.i >= 0; i++, s = *i) {
    _attribute[s.i].input = _attribute[s.i].output = false;
    _attribute[s.i].width = 0;
  }
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
  Point point = {0}; NOT_UNUSED (point);
@

@def end_foreach_stencil()
  check_stencil (&_loop);
  boundary_stencil (&_loop);
  _loop.first = 0;
}
@

@define foreach_vertex_stencil(...) foreach_stencil(S__VA_ARGS__) _loop.vertex = true;
@define end_foreach_vertex_stencil() end_foreach_stencil()

@define foreach_face_stencil(...) foreach_stencil(S__VA_ARGS__)
@define end_foreach_face_stencil() end_foreach_stencil()

@define foreach_visible_stencil(...) foreach_stencil(S__VA_ARGS__)
@define end_foreach_visible_stencil() end_foreach_stencil()

@define foreach_point_stencil(...) foreach_stencil(S__VA_ARGS__)
@define end_foreach_point_stencil() end_foreach_stencil()

@define foreach_region_stencil(...) foreach_stencil(S__VA_ARGS__)
@define end_foreach_region_stencil() end_foreach_stencil()

@define _stencil_is_face_x() { _loop.face |= (1 << 0);
@define end__stencil_is_face_x() }
@define _stencil_is_face_y() { _loop.face |= (1 << 1);
@define end__stencil_is_face_y() }
@define _stencil_is_face_z() { _loop.face |= (1 << 2);
@define end__stencil_is_face_z() }

void stencil_val (Point p, scalar s, int i, int j, int k,
    const char * file, int line, bool overflow);
void stencil_val_a (Point p, scalar s, int i, int j, int k, bool input,
      const char * file, int line);

@def _stencil_val(a,_i,_j,_k)
  stencil_val (point, a, _i, _j, _k, S__FILE__, S_LINENO, false)
@
@def _stencil_val_o(a,_i,_j,_k)
  stencil_val (point, a, _i, _j, _k, S__FILE__, S_LINENO, true)
@
@def _stencil_val_a(a,_i,_j,_k)
  stencil_val_a (point, a, _i, _j, _k, false, S__FILE__, S_LINENO)
@
@def _stencil_val_r(a,_i,_j,_k)
  stencil_val_a (point, a, _i, _j, _k, true, S__FILE__, S_LINENO)
@

@define _stencil_fine(a,_i,_j,_k) _stencil_val(a,_i,_j,_k)
@define _stencil_fine_a(a,_i,_j,_k) _stencil_val_a(a,_i,_j,_k)
@define _stencil_fine_r(a,_i,_j,_k) _stencil_val_r(a,_i,_j,_k)

@define _stencil_coarse(a,_i,_j,_k) _stencil_val(a,_i,_j,_k)
@define _stencil_coarse_a(a,_i,_j,_k) _stencil_val_a(a,_i,_j,_k)
@define _stencil_coarse_r(a,_i,_j,_k) _stencil_val_r(a,_i,_j,_k)

@define r_assign(x)
@define _assign(x)

@define _stencil_neighbor(i,j,k)
@define _stencil_child(i,j,k)
@define _stencil_aparent(i,j,k)
@define _stencil_aparent_a(i,j,k)
@define _stencil_aparent_r(i,j,k)

@define _stencil_neighborp(i,j,k) neighborp(i,j,k)

int _stencil_nop;
@define _stencil_val_higher_dimension (_stencil_nop = 1)
@define _stencil__val_constant(a,_i,_j,_k) (_stencil_nop = 1)

typedef void _stencil_undefined;

@define o_stencil -2







static inline bool scalar_is_dirty (scalar s)
{
  if (_attribute[s.i].dirty)
    return true;
  scalar * depends = _attribute[s.i].depends;
  {scalar*_i=(scalar*)( depends);if(_i)for(scalar d=*_i;(&d)->i>=0;d=*++_i){
    if (_attribute[d.i].dirty)
      return true;}}
  return false;
}




static inline bool scalar_depends_from (scalar a, scalar b)
{
  scalar * depends = _attribute[a.i].depends;
  {scalar*_i=(scalar*)( depends);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (s.i == b.i)
      return true;}}
  return false;
}







void boundary_internal (scalar * list, const char * fname, int line);
void (* boundary_face) (vectorl);







void check_stencil (ForeachData * loop)
{
  loop->listf = (vectorl){NULL};




  {scalar*_i=(scalar*)( baseblock);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    bool write = _attribute[s.i].output, read = _attribute[s.i].input;




    {





      if (read && scalar_is_dirty (s)) {





 if (_attribute[s.i].face) {
   if (_attribute[s.i].width > 0)
     loop->listc = list_append (loop->listc, s);
   else if (!write) {
     scalar sn = _attribute[s.i].v.x.i >= 0 ? _attribute[s.i].v.x : s;
     
       if (_attribute[s.i].v.x.i == s.i) {




  if (_attribute[sn.i].boundary[left] || _attribute[sn.i].boundary[right])
    loop->listc = list_append (loop->listc, s);
  else if (_attribute[s.i].dirty != 2)
    loop->listf.x = list_append (loop->listf.x, s);
       }
       
#line 213
if (_attribute[s.i].v.y.i == s.i) {




  if (_attribute[sn.i].boundary[bottom] || _attribute[sn.i].boundary[top])
    loop->listc = list_append (loop->listc, s);
  else if (_attribute[s.i].dirty != 2)
    loop->listf.y = list_append (loop->listf.y, s);
       }
   }
 }





 else if (_attribute[s.i].width > 0)
   loop->listc = list_append (loop->listc, s);
      }





      if (write) {
 if (2 > 1 && !loop->vertex && loop->first) {
   bool vertex = true;
   
     if (_attribute[s.i].d.x != -1)
       vertex = false;
     
#line 242
if (_attribute[s.i].d.y != -1)
       vertex = false;
   if (vertex)
     fprintf (ferr,
       "%s:%d: warning: vertex scalar '%s' should be assigned with"
       " a foreach_vertex() loop\n",
       loop->fname, loop->line, _attribute[s.i].name);
 }
 if (_attribute[s.i].face) {
   if (loop->face == 0 && loop->first)
     fprintf (ferr,
       "%s:%d: warning: face vector '%s' should be assigned with"
       " a foreach_face() loop\n",
       loop->fname, loop->line, _attribute[s.i].name);
 }
 else if (loop->face) {
   if (_attribute[s.i].v.x.i < 0) {
     int d = 1, i = 0;
      {
       if (loop->face == d) {
  _attribute[s.i].face = 2, _attribute[s.i].v.x.i = s.i;
  _attribute[s.i].boundary[left] = _attribute[s.i].boundary[right] = NULL;





       }
       d *= 2, i++;
     } 
#line 260
{
       if (loop->face == d) {
  _attribute[s.i].face = 2, _attribute[s.i].v.y.i = s.i;
  _attribute[s.i].boundary[bottom] = _attribute[s.i].boundary[top] = NULL;





       }
       d *= 2, i++;
     }
     if (!_attribute[s.i].face && loop->first)
       fprintf (ferr,
         "%s:%d: warning: scalar '%s' should be assigned with "
         "a foreach_face(x|y|z) loop\n",
         loop->fname, loop->line, _attribute[s.i].name);
   }
   else {
     char * name = NULL;
     if (_attribute[s.i].name) {
       name = pstrdup (_attribute[s.i].name,__func__,__FILE__,0);
       char * s = name + strlen(name) - 1;
       while (s != name && *s != '.') s--;
       if (s != name) *s = '\0';
     }
     struct { int x, y, z; } input, output;
     vector v = _attribute[s.i].v;

     
       input.x = _attribute[v.x.i].input, output.x = _attribute[v.x.i].output;
       
#line 290
input.y = _attribute[v.y.i].input, output.y = _attribute[v.y.i].output;

     init_face_vector (v, name);


     
       _attribute[v.x.i].input = input.x, _attribute[v.x.i].output = output.x;
       
#line 296
_attribute[v.y.i].input = input.y, _attribute[v.y.i].output = output.y;





     pfree (name,__func__,__FILE__,0);
   }
 }
 else if (loop->vertex) {
   bool vertex = true;
   
     if (_attribute[s.i].d.x != -1)
       vertex = false;
     
#line 308
if (_attribute[s.i].d.y != -1)
       vertex = false;
   if (!vertex) {
     char * name = NULL;
     if (_attribute[s.i].name) name = pstrdup (_attribute[s.i].name,__func__,__FILE__,0);
     init_vertex_scalar (s, name);
     
       _attribute[s.i].v.x.i = -1;
       
#line 315
_attribute[s.i].v.y.i = -1;




     pfree (name,__func__,__FILE__,0);
   }
 }





 loop->dirty = list_append (loop->dirty, s);
 {scalar*_i=(scalar*)( baseblock);if(_i)for(scalar d=*_i;(&d)->i>=0;d=*++_i){
   if (scalar_depends_from (d, s))
     loop->dirty = list_append (loop->dirty, d);}}
      }
    }
  }}}
}




void boundary_stencil (ForeachData * loop)
{
  bool flux = false;
  
    if (loop->listf.x)
      flux = true;
    
#line 344
if (loop->listf.y)
      flux = true;
  if (flux) {
# 359 "/home/spencer/basilisk/src/grid/stencils.h"
    boundary_face (loop->listf);
    
      pfree (loop->listf.x,__func__,__FILE__,0), loop->listf.x = NULL;
      
#line 361
pfree (loop->listf.y,__func__,__FILE__,0), loop->listf.y = NULL;
  }




  if (loop->listc) {






    boundary_internal (loop->listc, loop->fname, loop->line);
    pfree (loop->listc,__func__,__FILE__,0), loop->listc = NULL;
  }





  if (loop->dirty) {






    {scalar*_i=(scalar*)( loop->dirty);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = true;}}
    pfree (loop->dirty,__func__,__FILE__,0), loop->dirty = NULL;
  }
}
# 1462 "/home/spencer/basilisk/src/common.h" 2
# 13 "cylinder-osc-cpp.c" 2
# 1 "grid/quadtree.h" 1
# 1 "/home/spencer/basilisk/src/grid/quadtree.h"


# 1 "grid/tree.h" 1
# 1 "/home/spencer/basilisk/src/grid/tree.h"
# 1 "grid/mempool.h" 1
# 1 "/home/spencer/basilisk/src/grid/mempool.h"





typedef struct _Pool Pool;

struct _Pool {
  Pool * next;
};

typedef struct {
  char * first, * lastb;
  size_t size;
  size_t poolsize;
  Pool * pool, * last;
} Mempool;

typedef struct {
  char * next;
} FreeBlock;

Mempool * mempool_new (size_t poolsize, size_t size)
{

  if (!(poolsize % 8 == 0)) qassert ("/home/spencer/basilisk/src/grid/mempool.h", 0, "poolsize % 8 == 0");
  if (!(size >= sizeof(FreeBlock))) qassert ("/home/spencer/basilisk/src/grid/mempool.h", 0, "size >= sizeof(FreeBlock)");


  poolsize = min(1 << 20, poolsize + sizeof(Pool));
  Mempool * m = ((Mempool *) pcalloc (1, sizeof(Mempool),__func__,__FILE__,0));
  m->poolsize = poolsize;
  m->size = size;
  return m;
}

void mempool_destroy (Mempool * m)
{
  Pool * p = m->pool;
  while (p) {
    Pool * next = p->next;
    pfree (p,__func__,__FILE__,0);
    p = next;
  }
  pfree (m,__func__,__FILE__,0);
}

void * mempool_alloc (Mempool * m)
{
  if (!m->first) {

    Pool * p = (Pool *) pmalloc (m->poolsize,__func__,__FILE__,0);
    p->next = NULL;
    if (m->last)
      m->last->next = p;
    else
      m->pool = p;
    m->last = p;
    m->first = m->lastb = ((char *)m->last) + sizeof(Pool);
    FreeBlock * b = (FreeBlock *) m->first;
    b->next = NULL;
  }
  void * ret = m->first;
  FreeBlock * b = (FreeBlock *) ret;
  char * next = b->next;
  if (!next) {
    m->lastb += m->size;
    next = m->lastb;
    if (next + m->size > ((char *) m->last) + m->poolsize)
      next = NULL;
    else {
      FreeBlock * b = (FreeBlock *) next;
      b->next = NULL;
    }
  }
  m->first = next;
@if TRASH
  double * v = (double *) ret;
  for (int i = 0; i < m->size/sizeof(double); i++)
    v[i] = undefined;
@endif
  return ret;
}

void * mempool_alloc0 (Mempool * m)
{
  void * ret = mempool_alloc (m);
  memset (ret, 0, m->size);
  return ret;
}

void mempool_free (Mempool * m, void * p)
{
@if TRASH
  double * v = (double *) p;
  for (int i = 0; i < m->size/sizeof(double); i++)
    v[i] = undefined;
@endif
  FreeBlock * b = (FreeBlock *) p;
  b->next = m->first;
  m->first = (char *) p;
}
# 2 "/home/spencer/basilisk/src/grid/tree.h" 2




# 1 "grid/memindex/range.h" 1
# 1 "/home/spencer/basilisk/src/grid/memindex/range.h"
# 15 "/home/spencer/basilisk/src/grid/memindex/range.h"
typedef struct {
  void ** p;
  int size;
} Memalloc;

typedef struct {
  int start, end;
} Memrange;
# 34 "/home/spencer/basilisk/src/grid/memindex/range.h"
void memrange_alloc (Memrange * r, Memalloc * mem, int i)
{
  if (r->start == r->end) {
    r->start = i;
    r->end = i + 1;
    for (Memalloc * m = mem; m->p; m++) {
      *m->p = pcalloc (1, m->size,__func__,__FILE__,0);
      *m->p = (char *)(*m->p) - i*m->size;
    }
  }
  else if (i >= r->end) {
    for (Memalloc * m = mem; m->p; m++) {
      *m->p = prealloc ((char *)(*m->p) + r->start*m->size,
         m->size*(i + 1 - r->start),__func__,__FILE__,0);
      *m->p = (char *)(*m->p) - r->start*m->size;
      memset ((char *)(*m->p) + r->end*m->size, 0, (i - r->end + 1)*m->size);
    }
    r->end = i + 1;
  }
  else if (i < r->start) {
    for (Memalloc * m = mem; m->p; m++) {
      *m->p = prealloc ((char *)(*m->p) + r->start*m->size, m->size*(r->end - i),__func__,__FILE__,0);
      memmove ((char *)(*m->p) + (r->start - i)*m->size, *m->p,
        m->size*(r->end - r->start));
      memset ((char *)(*m->p), 0, (r->start - i)*m->size);
      *m->p = (char *)(*m->p) - i*m->size;
    }
    r->start = i;
  }
}
# 73 "/home/spencer/basilisk/src/grid/memindex/range.h"
bool memrange_free (Memrange * r, Memalloc * mem, int i)
{
  if (i == r->start) {
    if (i == r->end - 1) {
      for (Memalloc * m = mem; m->p; m++) {
 pfree ((char *)(*m->p) + r->start*m->size,__func__,__FILE__,0);
 *m->p = NULL;
      }
      r->start = r->end = 0;
      return true;
    }
    else {
      for (i = i + 1; i < r->end &&
      !*(void **)((char *)(*mem->p) + i*mem->size); i++);
      for (Memalloc * m = mem; m->p; m++) {
 memmove ((char *)(*m->p) + r->start*m->size,
   (char *)(*m->p) + i*m->size, m->size*(r->end - i));
 *m->p = prealloc ((char *)(*m->p) + r->start*m->size,
    m->size*(r->end - i),__func__,__FILE__,0);
 *m->p = (char *)(*m->p) - i*m->size;
      }
      r->start = i;
    }
  }
  else if (i == r->end - 1) {
    for (i = i - 1; i >= r->start &&
    !*(void **)((char *)(*mem->p) + i*mem->size); i--);
    r->end = i + 1;
    for (Memalloc * m = mem; m->p; m++) {
      *m->p = prealloc ((char *)(*m->p) + r->start*m->size,
         m->size*(r->end - r->start),__func__,__FILE__,0);
      *m->p = (char *)(*m->p) - r->start*m->size;
    }
  }
  else {
    if (!(i > r->start && i < r->end)) qassert ("/home/spencer/basilisk/src/grid/memindex/range.h", 0, "i > r->start && i < r->end");
    for (Memalloc * m = mem; m->p; m++)
      memset ((char *)(*m->p) + i*m->size, 0, m->size);
  }
  return false;
}







struct _Memindex {
  Memrange r1;

  Memrange * r2;







  char *** b;



};
# 171 "/home/spencer/basilisk/src/grid/memindex/range.h"
struct _Memindex * mem_new (int len)
{
  struct _Memindex * m = pcalloc (1, sizeof (struct _Memindex),__func__,__FILE__,0);
  return m;
}





void mem_destroy (struct _Memindex * m, int len)
{

  for (int i = m->r1.start; i < m->r1.end; i++)
    if (m->b[i]) {






      pfree (m->b[i] + m->r2[i].start,__func__,__FILE__,0);
    }
  if (m->b) {
    pfree (m->r2 + m->r1.start,__func__,__FILE__,0);



  }

  if (m->b)
    pfree (m->b + m->r1.start,__func__,__FILE__,0);
  pfree (m,__func__,__FILE__,0);
}
# 218 "/home/spencer/basilisk/src/grid/memindex/range.h"
void mem_assign (struct _Memindex * m, int i, int j, int len, void * b)
{
  Memalloc mem[] = {{(void **)&m->b, sizeof(char **)},
      {(void **)&m->r2, sizeof(Memrange)},
      {NULL}};
  memrange_alloc (&m->r1, mem, i);
  Memalloc mem1[] = {{(void **)&m->b[i], sizeof(char *)},
       {NULL}};
  memrange_alloc (&m->r2[i], mem1, j);
  ((m)->b[i][j]) = b;
}
# 259 "/home/spencer/basilisk/src/grid/memindex/range.h"
void mem_free (struct _Memindex * m, int i, int j, int len)
{
  Memalloc mem[] = {{(void **)&m->b[i], sizeof(char *)},
      {NULL}};
  if (memrange_free (&m->r2[i], mem, j)) {
    Memalloc mem[] = {{(void **)&m->b, sizeof(char **)},
        {(void **)&m->r2, sizeof(Memrange)},
        {NULL}};
    memrange_free (&m->r1, mem, i);
  }
}
# 305 "/home/spencer/basilisk/src/grid/memindex/range.h"
@def foreach_mem(_m, _len, _i) {
  Point point = {0};
  for (point.i = max(Period.x*2, (_m)->r1.start);
       point.i < min(_len - Period.x*2, (_m)->r1.end);
       point.i += _i)
    if ((_m)->b[point.i])
      for (point.j = max(Period.y*2, (_m)->r2[point.i].start);
    point.j < min(_len - Period.y*2, (_m)->r2[point.i].end);
    point.j += _i)
 if ((_m)->b[point.i][point.j]) {
@
@define end_foreach_mem() }}
# 7 "/home/spencer/basilisk/src/grid/tree.h" 2





@ define BGHOSTS 1
# 24 "/home/spencer/basilisk/src/grid/tree.h"
typedef struct {
  unsigned short flags;

  unsigned short neighbors;
  int pid;
} Cell;

enum {
  active = 1 << 0,
  leaf = 1 << 1,
  border = 1 << 2,
  vertex = 1 << 3,
  user = 4,

  face_x = 1 << 0

  , face_y = 1 << 1




};

@define is_active(cell) ((cell).flags & active)
@define is_leaf(cell) ((cell).flags & leaf)
@define is_coarse() ((cell).neighbors > 0)
@define is_border(cell) ((cell).flags & border)
@define is_local(cell) ((cell).pid == pid())
@define is_vertex(cell) ((cell).flags & vertex)



typedef struct {
  int i;

  int j;




} IndexLevel;

typedef struct {
  IndexLevel * p;
  int n, nm;
} CacheLevel;

typedef struct {
  int i;

  int j;




  int level, flags;
} Index;

typedef struct {
  Index * p;
  int n, nm;
} Cache;



typedef struct {
  struct _Memindex * m;
  Mempool * pool;
  long nc;
  int len;
} Layer;

static size_t _size (size_t depth)
{
  return (1 << depth) + 2*2;
}

static size_t poolsize (size_t depth, size_t size)
{




  return sq(_size(depth))*size;



}

static Layer * new_layer (int depth)
{
  Layer * l = ((Layer *) pmalloc ((1)*sizeof(Layer),__func__,__FILE__,0));
  l->len = _size (depth);
  if (depth == 0)
    l->pool = NULL;
  else {
    size_t size = sizeof(Cell) + datasize;


    l->pool = mempool_new (poolsize (depth, size), (1 << 2)*size);
  }
  l->m = mem_new (l->len);
  l->nc = 0;
  return l;
}

static void destroy_layer (Layer * l)
{
  if (l->pool)
    mempool_destroy (l->pool);
  mem_destroy (l->m, l->len);
  pfree (l,__func__,__FILE__,0);
}



typedef struct {
  Grid g;
  Layer ** L;

  Cache leaves;
  Cache faces;
  Cache vertices;
  Cache refined;
  CacheLevel * active;
  CacheLevel * prolongation;
  CacheLevel * boundary;

  CacheLevel * restriction;

  bool dirty;
} Tree;



struct _Point {

  int i;

  int j;




  int level;
@ifdef foreach_block
  int l;
  @define _BLOCK_INDEX , point.l
@else
  @define _BLOCK_INDEX
@endif
};
static Point last_point;



static void cache_level_append (CacheLevel * c, Point p)
{
  if (c->n >= c->nm) {
    c->nm += 128;
    c->p = (IndexLevel *) prealloc (c->p, (c->nm)*sizeof(IndexLevel),__func__,__FILE__,0);
  }
  c->p[c->n].i = p.i;

  c->p[c->n].j = p.j;




  c->n++;
}

static void cache_level_shrink (CacheLevel * c)
{
  if (c->nm > (c->n/128 + 1)*128) {
    c->nm = (c->n/128 + 1)*128;
    if (!(c->nm > c->n)) qassert ("/home/spencer/basilisk/src/grid/tree.h", 0, "c->nm > c->n");
    c->p = (IndexLevel *) prealloc (c->p, sizeof (Index)*c->nm,__func__,__FILE__,0);
  }
}

static void cache_append (Cache * c, Point p, unsigned short flags)
{
  if (c->n >= c->nm) {
    c->nm += 128;
    c->p = (Index *) prealloc (c->p, (c->nm)*sizeof(Index),__func__,__FILE__,0);
  }
  c->p[c->n].i = p.i;

  c->p[c->n].j = p.j;




  c->p[c->n].level = p.level;
  c->p[c->n].flags = flags;
  c->n++;
}

void cache_shrink (Cache * c)
{
  cache_level_shrink ((CacheLevel *)c);
}
# 243 "/home/spencer/basilisk/src/grid/tree.h"
@def allocated(k,l,n) (((point.i+k) >= (((Tree *)grid)->L[point.level]->m)->r1.start && (point.i+k) < (((Tree *)grid)->L[point.level]->m->r1.end) && (((Tree *)grid)->L[point.level]->m)->b[point.i+k] && (point.j+l) >= (((Tree *)grid)->L[point.level]->m)->r2[point.i+k].start && (point.j+l) < (((Tree *)grid)->L[point.level]->m)->r2[point.i+k].end && (((Tree *)grid)->L[point.level]->m)->b[point.i+k][point.j+l])
                               )
@
@def NEIGHBOR(k,l,n) (((((Tree *)grid)->L[point.level]->m)->b[point.i+k][point.j+l])
                            )
@
@def PARENT(k,l,n) (((((Tree *)grid)->L[point.level-1]->m)->b[(point.i+2)/2+k][(point.j+2)/2+l])
                                                    )
@
@def allocated_child(k,l,n) (level < depth() &&
         ((2*point.i-2 +k) >= (((Tree *)grid)->L[point.level+1]->m)->r1.start && (2*point.i-2 +k) < (((Tree *)grid)->L[point.level+1]->m->r1.end) && (((Tree *)grid)->L[point.level+1]->m)->b[2*point.i-2 +k] && (2*point.j-2 +l) >= (((Tree *)grid)->L[point.level+1]->m)->r2[2*point.i-2 +k].start && (2*point.j-2 +l) < (((Tree *)grid)->L[point.level+1]->m)->r2[2*point.i-2 +k].end && (((Tree *)grid)->L[point.level+1]->m)->b[2*point.i-2 +k][2*point.j-2 +l])

                             )
@
@def CHILD(k,l,n) (((((Tree *)grid)->L[point.level+1]->m)->b[2*point.i-2 +k][2*point.j-2 +l])
                                                )
@
# 284 "/home/spencer/basilisk/src/grid/tree.h"
@define CELL(m) (*((Cell *)(m)))


@define depth() (grid->depth)
@define aparent(k,l,n) CELL(PARENT(k,l,n))
@define child(k,l,n) CELL(CHILD(k,l,n))


@define cell CELL(NEIGHBOR(0,0,0))
@define neighbor(k,l,n) CELL(NEIGHBOR(k,l,n))
@def neighborp(l,m,n) (Point) {
    point.i + l,

    point.j + m,




    point.level
    _BLOCK_INDEX
}
@


@define data(k,l,n) ((double *) (NEIGHBOR(k,l,n) + sizeof(Cell)))
@define fine(a,k,p,n) ((double *) (CHILD(k,p,n) + sizeof(Cell)))[_index(a,n)]
@define coarse(a,k,p,n) ((double *) (PARENT(k,p,n) + sizeof(Cell)))[_index(a,n)]

@def POINT_VARIABLES
  VARIABLES
  int level = point.level; NOT_UNUSED(level);



  struct { int x, y; } child = {
    2*((point.i+2)%2)-1, 2*((point.j+2)%2)-1
  };





  NOT_UNUSED(child);
  Point parent = point; NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + 2)/2;

  parent.j = (point.j + 2)/2;
# 341 "/home/spencer/basilisk/src/grid/tree.h"
@

# 1 "grid/foreach_cell.h" 1
# 1 "/home/spencer/basilisk/src/grid/foreach_cell.h"
# 66 "/home/spencer/basilisk/src/grid/foreach_cell.h"
@def foreach_cell_root(root)
  {
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
    Point point = {0};



    struct { int l, i, j, stage; } stack[20];




    int _s = -1;
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].stage = 0; };
    while (_s >= 0) {
      int stage;
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; stage = stack[_s].stage; _s--; };
      if (!allocated (0,0,0))
 continue;
      switch (stage) {
      case 0: {
 POINT_VARIABLES;

@
@def end_foreach_cell_root()
        if (point.level < grid->depth) {
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 1; };
          { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].stage = 0; };
        }
        break;
      }



      case 1: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 2; };
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].stage = 0; }; break;
      case 2: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 3; };
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].stage = 0; }; break;
      case 3: { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].stage = 0; }; break;
# 120 "/home/spencer/basilisk/src/grid/foreach_cell.h"
      }
    }
  }
@

@def foreach_cell() {



  Point root = {2,2,0};



  foreach_cell_root (root)
@
@define end_foreach_cell() end_foreach_cell_root() }

@def foreach_cell_all() {
  Point root = {0};
  for (root.i = 2*Period.x; root.i <= 2*(2 - Period.x); root.i++)

    for (root.j = 2*Period.y; root.j <= 2*(2 - Period.y); root.j++)




 foreach_cell_root (root)
@
@define end_foreach_cell_all() end_foreach_cell_root() }

@def foreach_cell_post_root(condition, root)
  {
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
    Point point = {0};



    struct { int l, i, j, stage; } stack[20];




    int _s = -1;
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].stage = 0; };
    while (_s >= 0) {
      int stage;
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; stage = stack[_s].stage; _s--; };
      if (!allocated (0,0,0))
 continue;
      switch (stage) {
      case 0: {
        POINT_VARIABLES;
 if (point.level == grid->depth) {
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 8; };
 }
 else {
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 1; };
   if (condition)
     { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].stage = 0; };
 }
 break;
      }







      case 1:
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 2; };
 if (condition)
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].stage = 0; };
 break;
      case 2:
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 3; };
 if (condition)
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].stage = 0; };
 break;
      case 3:
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 4; };
 if (condition)
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].stage = 0; };
 break;
# 241 "/home/spencer/basilisk/src/grid/foreach_cell.h"
      default: {
        POINT_VARIABLES;

@
@def end_foreach_cell_post_root()
      }
      }
    }
  }
@

@def foreach_cell_post(condition)
  {



    Point root = {2,2,0};



    foreach_cell_post_root(condition, root)
@
@define end_foreach_cell_post() end_foreach_cell_post_root() }

@def foreach_cell_post_all(condition) {
  Point root = {0};
  for (root.i = 0; root.i <= 2*2; root.i++)

    for (root.j = 0; root.j <= 2*2; root.j++)




 foreach_cell_post_root (condition, root)
@
@define end_foreach_cell_post_all() end_foreach_cell_post_root() }

@def foreach_leaf() foreach_cell()
  if (is_leaf (cell)) {
    if (is_active(cell) && is_local(cell)) {
@
@define end_foreach_leaf() } continue; } end_foreach_cell()
# 344 "/home/spencer/basilisk/src/grid/tree.h" 2
# 361 "/home/spencer/basilisk/src/grid/tree.h"
@def foreach_child() {
  int _i = 2*point.i - 2, _j = 2*point.j - 2;
  point.level++;
  for (int _k = 0; _k < 2; _k++) {
    point.i = _i + _k;
    for (int _l = 0; _l < 2; _l++) {
      point.j = _j + _l;
      POINT_VARIABLES;
@
@def end_foreach_child()
    }
  }
  point.i = (_i + 2)/2; point.j = (_j + 2)/2;
  point.level--;
}
@
@define foreach_child_break() _k = _l = 2
# 407 "/home/spencer/basilisk/src/grid/tree.h"
@def is_refined_check() ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) &&
    point.i > 0 && point.i < (1 << level) + 2*2 - 1

    && point.j > 0 && point.j < (1 << level) + 2*2 - 1




    )
@

@def foreach_cache(_cache) {
  OMP_PARALLEL() {
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
  Point point = {0};
  point.i = 2;

  point.j = 2;




  int _k; unsigned short _flags; NOT_UNUSED(_flags);
  OMP(omp for schedule(static))
  for (_k = 0; _k < _cache.n; _k++) {
    point.i = _cache.p[_k].i;

    point.j = _cache.p[_k].j;




    point.level = _cache.p[_k].level;
    _flags = _cache.p[_k].flags;
    POINT_VARIABLES;
@
@define end_foreach_cache() } } }

@def foreach_cache_level(_cache,_l) {
  OMP_PARALLEL() {
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
  Point point = {0};
  point.i = 2;

  point.j = 2;




  point.level = _l;
  int _k;
  OMP(omp for schedule(static))
  for (_k = 0; _k < _cache.n; _k++) {
    point.i = _cache.p[_k].i;

    point.j = _cache.p[_k].j;




    POINT_VARIABLES;
@
@define end_foreach_cache_level() } } }

@def foreach_boundary_level(_l) {
  if (_l <= depth()) {
    { if (((Tree *)grid)->dirty) update_cache_f(); };
    CacheLevel _boundary = ((Tree *)grid)->boundary[_l];
    foreach_cache_level (_boundary,_l)
@
@define end_foreach_boundary_level() end_foreach_cache_level(); }}



@def foreach_boundary(_b) {
  for (int _l = depth(); _l >= 0; _l--)
    foreach_boundary_level(_l) {
      if ((- cell.pid - 1) == _b)
 for (int _d = 0; _d < 2; _d++) {
   for (int _i = -1; _i <= 1; _i += 2) {
     if (_d == 0) ig = _i; else if (_d == 1) jg = _i; else kg = _i;
     if (allocated(-ig,-jg,-kg) &&
  is_leaf (neighbor(-ig,-jg,-kg)) &&
  !(neighbor(-ig,-jg,-kg).pid < 0) &&
  is_local(neighbor(-ig,-jg,-kg))) {
       point.i -= ig; x -= ig*Delta/2.;

       point.j -= jg; y -= jg*Delta/2.;




@
@def end_foreach_boundary()
       point.i += ig; x += ig*Delta/2.;

       point.j += jg; y += jg*Delta/2.;




            }
   }
   ig = jg = kg = 0;
 }
    } end_foreach_boundary_level(); }
@

@def foreach_halo(_name,_l) {
  if (_l <= depth()) {
    { if (((Tree *)grid)->dirty) update_cache_f(); };
    CacheLevel _cache = ((Tree *)grid)->_name[_l];
    foreach_cache_level (_cache, _l)
@
@define end_foreach_halo() end_foreach_cache_level(); }}

# 1 "grid/neighbors.h" 1
# 1 "/home/spencer/basilisk/src/grid/neighbors.h"
# 17 "/home/spencer/basilisk/src/grid/neighbors.h"
@def foreach_neighbor(_s) {
  int _nn = _s + 0 ? _s + 0 : 2;
  int _i = point.i, _j = point.j;
  for (int _k = - _nn; _k <= _nn; _k++) {
    point.i = _i + _k;
    for (int _l = - _nn; _l <= _nn; _l++) {
      point.j = _j + _l;
      POINT_VARIABLES;
@
@def end_foreach_neighbor()
    }
  }
  point.i = _i; point.j = _j;
}
@
@define foreach_neighbor_break() _k = _l = _nn + 1
# 524 "/home/spencer/basilisk/src/grid/tree.h" 2

static inline bool has_local_children (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  {foreach_child()
    if (is_local(cell))
      return true;end_foreach_child()}
  return false;
}

static inline void cache_append_face (Point point, unsigned short flags)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  Tree * q = ((Tree *)grid);
  cache_append (&q->faces, point, flags);

  if (!is_vertex(cell)) {
    cache_append (&q->vertices, point, 0);
    cell.flags |= vertex;
  }
  
    if ((flags & face_y) && !is_vertex(neighbor(1,0,0))) {
      cache_append (&q->vertices, neighborp(1,0,0), 0);
      neighbor(1,0,0).flags |= vertex;
    }
    
#line 543
if ((flags & face_x) && !is_vertex(neighbor(0,1,0))) {
      cache_append (&q->vertices, neighborp(0,1,0), 0);
      neighbor(0,1,0).flags |= vertex;
    }
# 557 "/home/spencer/basilisk/src/grid/tree.h"
}



static void update_cache_f (void)
{
  Tree * q = ((Tree *)grid);

  {foreach_cache (q->vertices)
    if (level <= depth() && allocated(0,0,0))
      cell.flags &= ~vertex;end_foreach_cache();}


  q->leaves.n = q->faces.n = q->vertices.n = 0;
  for (int l = 0; l <= depth(); l++)
    q->active[l].n = q->prolongation[l].n =
      q->boundary[l].n = q->restriction[l].n = 0;

  const unsigned short fboundary = 1 << user;
  {foreach_cell() {



    if (is_local(cell) && is_active(cell)) {


      cache_level_append (&q->active[level], point);
    }
# 601 "/home/spencer/basilisk/src/grid/tree.h"
    if (!(cell.pid < 0)) {

      {foreach_neighbor (BGHOSTS)
 if (allocated(0,0,0) && (cell.pid < 0) && !(cell.flags & fboundary)) {
   cache_level_append (&q->boundary[level], point);
   cell.flags |= fboundary;
 }end_foreach_neighbor()}
    }

    else if (level > 0 && is_local(aparent(0,0,0)))
      cache_level_append (&q->restriction[level], point);

    if (is_leaf (cell)) {
      if (is_local(cell)) {
 cache_append (&q->leaves, point, 0);

 unsigned short flags = 0;
 
   if ((neighbor(-1,0,0).pid < 0) || (!is_leaf(neighbor(-1,0,0)) && !neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0) ||
       is_leaf(neighbor(-1,0,0)))
     flags |= face_x;
   
#line 619
if ((neighbor(0,-1,0).pid < 0) || (!is_leaf(neighbor(0,-1,0)) && !neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0) ||
       is_leaf(neighbor(0,-1,0)))
     flags |= face_y;
 if (flags)
   cache_append (&q->faces, point, flags);
 
   if ((neighbor(1,0,0).pid < 0) || (!is_leaf(neighbor(1,0,0)) && !neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0) ||
       (!is_local(neighbor(1,0,0)) && is_leaf(neighbor(1,0,0))))
     cache_append (&q->faces, neighborp(1,0,0), face_x);
   
#line 625
if ((neighbor(0,1,0).pid < 0) || (!is_leaf(neighbor(0,1,0)) && !neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0) ||
       (!is_local(neighbor(0,1,0)) && is_leaf(neighbor(0,1,0))))
     cache_append (&q->faces, neighborp(0,1,0), face_y);

 for (int i = 0; i <= 1; i++)

   for (int j = 0; j <= 1; j++)




       if (!is_vertex(neighbor(i,j,k))) {
  cache_append (&q->vertices, neighborp(i,j,k), 0);
  neighbor(i,j,k).flags |= vertex;
       }

        if (cell.neighbors > 0)
   cache_level_append (&q->prolongation[level], point);
      }
      else if (!(cell.pid < 0) || is_local(aparent(0,0,0))) {

 unsigned short flags = 0;
 
   if (allocated(-1,0,0) &&
       is_local(neighbor(-1,0,0)) && (!is_leaf(neighbor(-1,0,0)) && !neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0))
     flags |= face_x;
   
#line 648
if (allocated(0,-1,0) &&
       is_local(neighbor(0,-1,0)) && (!is_leaf(neighbor(0,-1,0)) && !neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0))
     flags |= face_y;
 if (flags)
   cache_append_face (point, flags);
 
   if (allocated(1,0,0) && is_local(neighbor(1,0,0)) &&
       (!is_leaf(neighbor(1,0,0)) && !neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0))
     cache_append_face (neighborp(1,0,0), face_x);
   
#line 654
if (allocated(0,1,0) && is_local(neighbor(0,1,0)) &&
       (!is_leaf(neighbor(0,1,0)) && !neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0))
     cache_append_face (neighborp(0,1,0), face_y);
      }

      continue;

    }
  }end_foreach_cell();}


  cache_shrink (&q->leaves);
  cache_shrink (&q->faces);
  cache_shrink (&q->vertices);
  for (int l = 0; l <= depth(); l++) {
    cache_level_shrink (&q->active[l]);
    cache_level_shrink (&q->prolongation[l]);
    cache_level_shrink (&q->boundary[l]);
    cache_level_shrink (&q->restriction[l]);
}

  q->dirty = false;


  for (int l = depth(); l >= 0; l--)
    {foreach_boundary_level (l)
      cell.flags &= ~fboundary;end_foreach_boundary_level();}



  grid->n = q->leaves.n;

@if !_MPI
  grid->tn = grid->n;
  grid->maxdepth = grid->depth;
@endif
}

@define foreach() { if (((Tree *)grid)->dirty) update_cache_f(); }; foreach_cache(((Tree *)grid)->leaves)
@define end_foreach() end_foreach_cache()

@def foreach_face_generic()
  { if (((Tree *)grid)->dirty) update_cache_f(); };
  foreach_cache(((Tree *)grid)->faces) @
@define end_foreach_face_generic() end_foreach_cache()

@define is_face_x() { int ig = -1; VARIABLES; if (_flags & face_x) {
@define end_is_face_x() }}


@define is_face_y() { int jg = -1; VARIABLES; if (_flags & face_y) {
@define end_is_face_y() }}






@def foreach_vertex()
  { if (((Tree *)grid)->dirty) update_cache_f(); };
  foreach_cache(((Tree *)grid)->vertices) {
    x -= Delta/2.;

    y -= Delta/2.;




@
@define end_foreach_vertex() } end_foreach_cache()
# 734 "/home/spencer/basilisk/src/grid/tree.h"
@def foreach_level(l) {
  if (l <= depth()) {
    { if (((Tree *)grid)->dirty) update_cache_f(); };
    CacheLevel _active = ((Tree *)grid)->active[l];
    foreach_cache_level (_active,l)
@
@define end_foreach_level() end_foreach_cache_level(); }}

@define foreach_coarse_level(l) foreach_level(l) if (!is_leaf(cell)) {
@define end_foreach_coarse_level() } end_foreach_level()

@def foreach_level_or_leaf(l) {
  for (int _l1 = l; _l1 >= 0; _l1--)
    foreach_level(_l1)
      if (_l1 == l || is_leaf (cell)) {
@
@define end_foreach_level_or_leaf() } end_foreach_level(); }

@if TRASH
@ undef trash
@ define trash(list) reset(list, undefined)
@endif

void reset (void * alist, double val)
{
  scalar * list = (scalar *) alist;
  Tree * q = ((Tree *)grid);

  for (int l = 0; l <= depth(); l++) {
    Layer * L = q->L[l];
    {foreach_mem (L->m, L->len, 1) {
      point.level = l;
      {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
 if (!is_constant(s))
   for (int b = 0; b < _attribute[s.i].block; b++)
     data(0,0,0)[s.i + b] = val;
      }}}
    }end_foreach_mem();}
  }
}

static CacheLevel * cache_level_resize (CacheLevel * name, int a)
{
  for (int i = 0; i <= depth() - a; i++)
    pfree (name[i].p,__func__,__FILE__,0);
  pfree (name,__func__,__FILE__,0);
  return ((CacheLevel *) pcalloc (depth() + 1, sizeof(CacheLevel),__func__,__FILE__,0));
}

static void update_depth (int inc)
{
  Tree * q = ((Tree *)grid);
  grid->depth += inc;
  q->L = &(q->L[-1]);
  q->L = (Layer * *) prealloc (q->L, (grid->depth + 2)*sizeof(Layer *),__func__,__FILE__,0);
  q->L = &(q->L[1]);
  if (inc > 0)
    q->L[grid->depth] = new_layer (grid->depth);
  q->active = cache_level_resize (q->active, inc);
  q->prolongation = cache_level_resize (q->prolongation, inc);
  q->boundary = cache_level_resize (q->boundary, inc);
  q->restriction = cache_level_resize (q->restriction, inc);
}
# 823 "/home/spencer/basilisk/src/grid/tree.h"
typedef void (* PeriodicFunction) (struct _Memindex *, int, int, int, void *);

static void periodic_function (struct _Memindex * m, int i, int j, int len, void * b,
          PeriodicFunction f)
{
  f(m, i, j, len, b);
  if (Period.x) {
    int nl = len - 2*2;
    for (int l = - 1; l <= 1; l += 2)
      for (int n = i + l*nl; n >= 0 && n < len; n += l*nl)
 f(m, n, j, len, b);
    if (Period.y)
      for (int l = - 1; l <= 1; l += 2)
 for (int n = j + l*nl; n >= 0 && n < len; n += l*nl) {
   f(m, i, n, len, b);
   for (int o = - 1; o <= 1; o += 2)
     for (int p = i + o*nl; p >= 0 && p < len; p += o*nl)
       f(m, p, n, len, b);
 }
  }
  else if (Period.y) {
    int nl = len - 2*2;
    for (int l = - 1; l <= 1; l += 2)
      for (int n = j + l*nl; n >= 0 && n < len; n += l*nl)
 f(m, i, n, len, b);
  }
}

static void assign_periodic (struct _Memindex * m, int i, int j, int len, void * b)
{
  periodic_function (m, i, j, len, b, mem_assign);
}

static void free_periodic (struct _Memindex * m, int i, int j, int len)
{
  periodic_function (m, i, j, len, NULL, (PeriodicFunction) mem_free);
}
# 938 "/home/spencer/basilisk/src/grid/tree.h"
static void alloc_children (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (point.level == grid->depth)
    update_depth (+1);
  else if (allocated_child(0,0,0))
    return;


  Layer * L = ((Tree *)grid)->L[point.level + 1];
  L->nc++;
  size_t len = sizeof(Cell) + datasize;
  char * b = (char *) mempool_alloc0 (L->pool);
  int i = 2*point.i - 2;
  for (int k = 0; k < 2; k++, i++) {




    int j = 2*point.j - 2;
    for (int l = 0; l < 2; l++, j++) {
      assign_periodic (L->m, i, j, L->len, b);
      b += len;
    }
# 971 "/home/spencer/basilisk/src/grid/tree.h"
  }

  int pid = cell.pid;
  {foreach_child() {
    cell.pid = pid;
@if TRASH
    {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      val(s,0,0,0) = undefined;}}
@endif
  }end_foreach_child()}
}
# 1000 "/home/spencer/basilisk/src/grid/tree.h"
static void free_children (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

  Layer * L = ((Tree *)grid)->L[point.level + 1];
  int i = 2*point.i - 2, j = 2*point.j - 2;
  if (!(((L->m)->b[i][j]))) qassert ("/home/spencer/basilisk/src/grid/tree.h", 0, "mem_data (L->m,i,j)");
  mempool_free (L->pool, ((L->m)->b[i][j]));
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++)
      free_periodic (L->m, i + k, j + l, L->len);
  if (--L->nc == 0) {
    destroy_layer (L);
    if (!(point.level + 1 == grid->depth)) qassert ("/home/spencer/basilisk/src/grid/tree.h", 0, "point.level + 1 == grid->depth");
    update_depth (-1);
  }
}
# 1041 "/home/spencer/basilisk/src/grid/tree.h"
void increment_neighbors (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  ((Tree *)grid)->dirty = true;
  if (cell.neighbors++ == 0)
    alloc_children (point);
  {foreach_neighbor (2/2)
    if (cell.neighbors++ == 0)
      alloc_children (point);end_foreach_neighbor()}
  cell.neighbors--;
}

void decrement_neighbors (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  ((Tree *)grid)->dirty = true;
  {foreach_neighbor (2/2)
    if (allocated(0,0,0)) {
      cell.neighbors--;
      if (cell.neighbors == 0)
 free_children (point);
    }end_foreach_neighbor()}
  if (cell.neighbors) {
    int pid = cell.pid;
    {foreach_child() {
      cell.flags = 0;
      cell.pid = pid;
    }end_foreach_child()}
  }
}

void realloc_scalar (int size)
{

  Tree * q = ((Tree *)grid);
  size_t oldlen = sizeof(Cell) + datasize;
  size_t newlen = oldlen + size;
  datasize += size;

  Layer * L = q->L[0];
  {foreach_mem (L->m, L->len, 1) {




    char * p = (char *) prealloc (((L->m)->b[point.i][point.j]),
     newlen*sizeof(char),__func__,__FILE__,0);
    assign_periodic (L->m, point.i, point.j, L->len, p);





  }end_foreach_mem();}

  for (int l = 1; l <= depth(); l++) {
    Layer * L = q->L[l];
    Mempool * oldpool = L->pool;
    L->pool = mempool_new (poolsize (l, newlen), (1 << 2)*newlen);
    {foreach_mem (L->m, L->len, 2) {
      char * new = (char *) mempool_alloc (L->pool);







      for (int k = 0; k < 2; k++)
 for (int o = 0; o < 2; o++) {
   memcpy (new, ((L->m)->b[point.i + k][point.j + o]), oldlen);
   assign_periodic (L->m, point.i + k, point.j + o, L->len, new);
   new += newlen;
 }
# 1124 "/home/spencer/basilisk/src/grid/tree.h"
    }end_foreach_mem();}
    mempool_destroy (oldpool);
  }
}



@define VN v.x
@define VT v.y
@define VR v.z




@if _MPI
@ define disable_fpe_for_mpi() disable_fpe (FE_DIVBYZERO|FE_INVALID)
@ define enable_fpe_for_mpi() enable_fpe (FE_DIVBYZERO|FE_INVALID)
@else
@ define disable_fpe_for_mpi()
@ define enable_fpe_for_mpi()
@endif

static inline void no_restriction (Point point, scalar s);

static bool normal_neighbor (Point point, scalar * scalars, vector * vectors)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  for (int k = 1; k <= BGHOSTS; k++)
    {
      for (int i = -k; i <= k; i += 2*k)
 if ((allocated(i,0,0) && !(neighbor(i,0,0).pid < 0))) {
   Point neighbor = neighborp(i,0,0);
   int id = (- cell.pid - 1);
   {scalar*_i=(scalar*)( scalars);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    
       val(s,0,0,0) = _attribute[s.i].boundary[id](neighbor, point, s, NULL);}}
   {vector*_i=(vector*)( vectors);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
     {
       scalar vn = VN;
       val(v.x,0,0,0) = _attribute[vn.i].boundary[id](neighbor, point, v.x, NULL);

       scalar vt = VT;
       val(v.y,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.y, NULL);





     }}}
   return true;
 }
      
#line 1152
for (int i = -k; i <= k; i += 2*k)
 if ((allocated(0,i,0) && !(neighbor(0,i,0).pid < 0))) {
   Point neighbor = neighborp(0,i,0);
   int id = (- cell.pid - 1);
   {scalar*_i=(scalar*)( scalars);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    
       val(s,0,0,0) = _attribute[s.i].boundary[id](neighbor, point, s, NULL);}}
   {vector*_i=(vector*)( vectors);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
     {
       scalar vn = VN;
       val(v.y,0,0,0) = _attribute[vn.i].boundary[id](neighbor, point, v.y, NULL);

       scalar vt = VT;
       val(v.x,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.x, NULL);





     }}}
   return true;
 }}
  return false;
}

static bool diagonal_neighbor_2D (Point point,
      scalar * scalars, vector * vectors)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

  for (int k = 1; k <= BGHOSTS; k++)



      for (int i = -k; i <= k; i += 2*k)
 for (int j = -k; j <= k; j += 2*k)
   if (allocated(i,j,0) && (allocated(i,j,0) && !(neighbor(i,j,0).pid < 0)) &&
       allocated(i,0,0) && (neighbor(i,0,0).pid < 0) &&
       allocated(0,j,0) && (neighbor(0,j,0).pid < 0)) {
     Point n = neighborp(i,j,0),
       n1 = neighborp(i,0,0), n2 = neighborp(0,j,0);
     int id1 = (- neighbor(i,0,0).pid - 1), id2 = (- neighbor(0,j,0).pid - 1);
     {scalar*_i=(scalar*)( scalars);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      
  val(s,0,0,0) = (_attribute[s.i].boundary[id1](n,n1,s,NULL) +
         _attribute[s.i].boundary[id2](n,n2,s,NULL) -
         val(s,i,j,0));}}
     {vector*_i=(vector*)( vectors);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
       {
  scalar vt = VT, vn = VN;
  val(v.x,0,0,0) = (_attribute[vt.i].boundary[id1](n,n1,v.x,NULL) +
    _attribute[vn.i].boundary[id2](n,n2,v.x,NULL) -
    val(v.x,i,j,0));
  val(v.y,0,0,0) = (_attribute[vn.i].boundary[id1](n,n1,v.y,NULL) +
    _attribute[vt.i].boundary[id2](n,n2,v.y,NULL) -
    val(v.y,i,j,0));






       }}}
     return true;
   }

  return false;
}

static bool diagonal_neighbor_3D (Point point,
      scalar * scalars, vector * vectors)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
# 1266 "/home/spencer/basilisk/src/grid/tree.h"
  return false;
}



static Point tangential_neighbor_x (Point point, bool * zn)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  for (int k = 1; k <= BGHOSTS; k++)
    for (int j = -k; j <= k; j += 2*k) {
      if ((allocated(0,j,0) && !(neighbor(0,j,0).pid < 0)) || (allocated(-1,j,0) && !(neighbor(-1,j,0).pid < 0))) {
 *zn = false;
 return neighborp(0,j,0);
      }







    }
  return (Point){.level = -1};
}

#line 1271
static Point tangential_neighbor_y (Point point, bool * zn)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  for (int k = 1; k <= BGHOSTS; k++)
    for (int j = -k; j <= k; j += 2*k) {
      if ((allocated(j,0,0) && !(neighbor(j,0,0).pid < 0)) || (allocated(j,-1,0) && !(neighbor(j,-1,0).pid < 0))) {
 *zn = false;
 return neighborp(j,0,0);
      }







    }
  return (Point){.level = -1};
}


static inline bool is_boundary_point (Point point) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  return (cell.pid < 0);
}

static void box_boundary_level (const Boundary * b, scalar * list, int l)
{
  disable_fpe_for_mpi();
  scalar * scalars = NULL;
  vector * vectors = NULL, * faces = NULL;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant(s) && _attribute[s.i].refine != no_restriction) {
      if (_attribute[s.i].v.x.i == s.i) {
 if (_attribute[s.i].face)
   faces = vectors_add (faces, _attribute[s.i].v);
 else
   vectors = vectors_add (vectors, _attribute[s.i].v);
      }
      else if (_attribute[s.i].v.x.i < 0 && _attribute[s.i].boundary[0])
 scalars = list_add (scalars, s);
    }}}

  {foreach_boundary_level (l) {
    if (!normal_neighbor (point, scalars, vectors) &&
 !diagonal_neighbor_2D (point, scalars, vectors) &&
 !diagonal_neighbor_3D (point, scalars, vectors)) {

      {scalar*_i=(scalar*)( scalars);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){

   val(s,0,0,0) = undefined;}}
      {vector*_i=(vector*)( vectors);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){

   {
     val(v.x,0,0,0) = undefined;
     
#line 1323
val(v.y,0,0,0) = undefined;}}}
    }
    if (faces) {
      int id = (- cell.pid - 1);
      
 for (int i = -1; i <= 1; i += 2) {

   if ((allocated(i,0,0) && !(neighbor(i,0,0).pid < 0))) {
     Point neighbor = neighborp(i,0,0);
     {vector*_i=(vector*)( faces);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){ {
       scalar vn = VN;
       if (_attribute[vn.i].boundary[id])
 
    val(v.x,(i + 1)/2,0,0) = _attribute[vn.i].boundary[id](neighbor, point, v.x, NULL);
     }}}
   }

   else if (i == -1) {

     bool zn;
     Point neighbor = tangential_neighbor_x (point, &zn);
     if (neighbor.level >= 0) {
       int id = is_boundary_point (neighbor) ?
  (- neighbor(-1,0,0).pid - 1) : (- cell.pid - 1);
       {vector*_i=(vector*)( faces);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){ {

  scalar vt = VT;



 
    val(v.x,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.x, NULL);
       }}}
     }
     else

       {vector*_i=(vector*)( faces);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
 
    val(v.x,0,0,0) = 0.;}}
   }

 }
 
#line 1328
for (int i = -1; i <= 1; i += 2) {

   if ((allocated(0,i,0) && !(neighbor(0,i,0).pid < 0))) {
     Point neighbor = neighborp(0,i,0);
     {vector*_i=(vector*)( faces);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){ {
       scalar vn = VN;
       if (_attribute[vn.i].boundary[id])
 
    val(v.y,0,(i + 1)/2,0) = _attribute[vn.i].boundary[id](neighbor, point, v.y, NULL);
     }}}
   }

   else if (i == -1) {

     bool zn;
     Point neighbor = tangential_neighbor_y (point, &zn);
     if (neighbor.level >= 0) {
       int id = is_boundary_point (neighbor) ?
  (- neighbor(0,-1,0).pid - 1) : (- cell.pid - 1);
       {vector*_i=(vector*)( faces);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){ {

  scalar vt = VT;



 
    val(v.y,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.y, NULL);
       }}}
     }
     else

       {vector*_i=(vector*)( faces);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
 
    val(v.y,0,0,0) = 0.;}}
   }

 }
    }
  }end_foreach_boundary_level();}

  pfree (scalars,__func__,__FILE__,0);
  pfree (vectors,__func__,__FILE__,0);
  pfree (faces,__func__,__FILE__,0);
  enable_fpe_for_mpi();
}



@undef VN
@undef VT
@define VN _attribute[s.i].v.x
@define VT _attribute[s.i].v.y

static double masked_average (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  double sum = 0., n = 0.;
  {foreach_child()
    if (!(cell.pid < 0) && val(s,0,0,0) != 1e30)
      sum += val(s,0,0,0), n++;end_foreach_child()}
  return n ? sum/n : 1e30;
}


static double masked_average_x (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  double sum = 0., n = 0.;
  {foreach_child()
    if (child.x < 0 && (!(cell.pid < 0) || !(neighbor(1,0,0).pid < 0)) &&
 val(s,1,0,0) != 1e30)
      sum += val(s,1,0,0), n++;end_foreach_child()}
  return n ? sum/n : 1e30;
}

#line 1391
static double masked_average_y (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  double sum = 0., n = 0.;
  {foreach_child()
    if (child.y < 0 && (!(cell.pid < 0) || !(neighbor(0,1,0).pid < 0)) &&
 val(s,0,1,0) != 1e30)
      sum += val(s,0,1,0), n++;end_foreach_child()}
  return n ? sum/n : 1e30;
}

static void masked_boundary_restriction (const Boundary * b,
      scalar * list, int l)
{
  scalar * scalars = NULL;
  vector * faces = NULL;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant(s) && _attribute[s.i].refine != no_restriction) {
      if (_attribute[s.i].v.x.i == s.i && _attribute[s.i].face)
 faces = vectors_add (faces, _attribute[s.i].v);
      else
 scalars = list_add (scalars, s);
    }}}

  {foreach_halo (restriction, l) {
    {scalar*_i=(scalar*)( scalars);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      val(s,0,0,0) = masked_average (parent, s);}}
    {vector*_i=(vector*)( faces);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
      { {
 double average = masked_average_x (parent, v.x);
 if ((neighbor(-1,0,0).pid < 0))
   val(v.x,0,0,0) = average;
 if ((neighbor(1,0,0).pid < 0))
   val(v.x,1,0,0) = average;
      } 
#line 1418
{
 double average = masked_average_y (parent, v.y);
 if ((neighbor(0,-1,0).pid < 0))
   val(v.y,0,0,0) = average;
 if ((neighbor(0,1,0).pid < 0))
   val(v.y,0,1,0) = average;
      }}}}
  }end_foreach_halo();}

  pfree (scalars,__func__,__FILE__,0);
  pfree (faces,__func__,__FILE__,0);
}
# 1454 "/home/spencer/basilisk/src/grid/tree.h"
static void free_cache (CacheLevel * c)
{
  for (int l = 0; l <= depth(); l++)
    pfree (c[l].p,__func__,__FILE__,0);
  pfree (c,__func__,__FILE__,0);
}

void free_grid (void)
{
  if (!grid)
    return;
  free_boundaries();
  Tree * q = ((Tree *)grid);
  pfree (q->leaves.p,__func__,__FILE__,0);
  pfree (q->faces.p,__func__,__FILE__,0);
  pfree (q->vertices.p,__func__,__FILE__,0);
  pfree (q->refined.p,__func__,__FILE__,0);


  Layer * L = q->L[0];
  {foreach_mem (L->m, L->len, 1) {



    pfree (((L->m)->b[point.i][point.j]),__func__,__FILE__,0);



  }end_foreach_mem();}
  for (int l = 0; l <= depth(); l++)
    destroy_layer (q->L[l]);
  q->L = &(q->L[-1]);
  pfree (q->L,__func__,__FILE__,0);
  free_cache (q->active);
  free_cache (q->prolongation);
  free_cache (q->boundary);
  free_cache (q->restriction);
  pfree (q,__func__,__FILE__,0);
  grid = NULL;
}

static void refine_level (int depth);

     
void init_grid (int n)
{tracing("init_grid","/home/spencer/basilisk/src/grid/tree.h",0);

  if (!(sizeof(Cell) % 8 == 0)) qassert ("/home/spencer/basilisk/src/grid/tree.h", 0, "sizeof(Cell) % 8 == 0");

  free_grid();
  int depth = 0;
  while (n > 1) {
    if (n % 2) {
      fprintf (ferr, "tree: N must be a power-of-two\n");
      exit (1);
    }
    n /= 2;
    depth++;
  }
  Tree * q = ((Tree *) pcalloc (1, sizeof(Tree),__func__,__FILE__,0));
  grid = (Grid *) q;
  grid->depth = 0;


  q->L = ((Layer * *) pmalloc ((2)*sizeof(Layer *),__func__,__FILE__,0));

  q->L[0] = NULL; q->L = &(q->L[1]);

  Layer * L = new_layer (0);
  q->L[0] = L;
# 1537 "/home/spencer/basilisk/src/grid/tree.h"
  for (int i = Period.x*2; i < L->len - Period.x*2; i++)
    for (int j = Period.y*2; j < L->len - Period.y*2; j++)
      assign_periodic (L->m, i, j, L->len,
         (char *) pcalloc (1, sizeof(Cell) + datasize,__func__,__FILE__,0));
  CELL(((L->m)->b[2][2])).flags |= leaf;
  if (pid() == 0)
    CELL(((L->m)->b[2][2])).flags |= active;
  for (int k = - 2*(1 - Period.x); k <= 2*(1 - Period.x); k++)
    for (int l = -2*(1 - Period.y); l <= 2*(1 - Period.y); l++)
      CELL(((L->m)->b[2 +k][2 +l])).pid =
 (k < 0 ? -1 - left :
  k > 0 ? -1 - right :
  l > 0 ? -1 - top :
  l < 0 ? -1 - bottom :
  0);
  CELL(((L->m)->b[2][2])).pid = 0;
# 1575 "/home/spencer/basilisk/src/grid/tree.h"
  q->active = ((CacheLevel *) pcalloc (1, sizeof(CacheLevel),__func__,__FILE__,0));
  q->prolongation = ((CacheLevel *) pcalloc (1, sizeof(CacheLevel),__func__,__FILE__,0));
  q->boundary = ((CacheLevel *) pcalloc (1, sizeof(CacheLevel),__func__,__FILE__,0));
  q->restriction = ((CacheLevel *) pcalloc (1, sizeof(CacheLevel),__func__,__FILE__,0));
  q->dirty = true;
  N = 1 << depth;
@if _MPI
  void mpi_boundary_new();
  mpi_boundary_new();
@endif

  Boundary * b = ((Boundary *) pcalloc (1, sizeof(Boundary),__func__,__FILE__,0));
  b->level = box_boundary_level;
  b->restriction = masked_boundary_restriction;
  add_boundary (b);
  refine_level (depth);
  reset (all, 0.);
  { if (((Tree *)grid)->dirty) update_cache_f(); };
end_tracing("init_grid","/home/spencer/basilisk/src/grid/tree.h",0);}


void check_two_one (void)
{
  {foreach_leaf()
    if (level > 0)
      for (int k = -1; k <= 1; k++)
 for (int l = -1; l <= 1; l++) {

   int i = (point.i + 2)/2 + k;
   int j = (point.j + 2)/2 + l;
   double x = ((i - 2 + 0.5)*(1./(1 << point.level))*2. - 0.5);
   double y = ((j - 2 + 0.5)*(1./(1 << point.level))*2. - 0.5);
   if (x > -0.5 && x < 0.5 && y > -0.5 && y < 0.5 &&
       !(aparent(k,l,0).flags & active)) {
     FILE * fp = fopen("check_two_one_loc", "w");
     fprintf (fp,
       "# %d %d\n"
       "%g %g\n%g %g\n",
       k, l,
       (((point.i - 2) + 0.5)*(1./(1 << point.level)) - 0.5),
       (((point.j - 2) + 0.5)*(1./(1 << point.level)) - 0.5),
       x, y);
     fclose (fp);





     if (!(false)) qassert ("/home/spencer/basilisk/src/grid/tree.h", 0, "false");
   }
 }end_foreach_leaf();}
}


Point locate (double xp, double yp, double zp)
{
  for (int l = depth(); l >= 0; l--) {
    Point point = {0};int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
    point.level = l;
    int n = 1 << point.level;
    point.i = (xp - X0)/L0*n + 2;

    point.j = (yp - Y0)/L0*n + 2;




    if (point.i >= 0 && point.i < n + 2*2

 && point.j >= 0 && point.j < n + 2*2




 ) {
      if (allocated(0,0,0) && is_local(cell) && is_leaf(cell))
 return point;
    }
    else
      break;
  }
  Point point = {0};int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  point.level = -1;
  return point;
}



bool tree_is_full()
{
  { if (((Tree *)grid)->dirty) update_cache_f(); };
  return (grid->tn == 1L << grid->maxdepth*2);
}

# 1 "grid/tree-common.h" 1
# 1 "/home/spencer/basilisk/src/grid/tree-common.h"



# 1 "grid/multigrid-common.h" 1
# 1 "/home/spencer/basilisk/src/grid/multigrid-common.h"


# 1 "grid/cartesian-common.h" 1
# 1 "/home/spencer/basilisk/src/grid/cartesian-common.h"
# 1 "grid/events.h" 1
# 1 "/home/spencer/basilisk/src/grid/events.h"




static int END_EVENT = 1234567890;
static double TEND_EVENT = 1234567890;
static double TEPS = 1e-9;

static void event_error (Event * ev, const char * s)
{
  fprintf (ferr, "%s:%d: error: %s\n", ev->file, ev->line, s);
  exit (1);
}

static void init_event (Event * ev)
{
  if (ev->arrayi || ev->arrayt) {
    ev->i = -1; ev->t = - TEND_EVENT;
    if (ev->arrayi)
      ev->i = ev->arrayi[0];
    else
      ev->t = ev->arrayt[0];
    ev->a = 1;
    ev->expr[1] = NULL;
  }
  else {
    if (ev->nexpr > 0) {
      Expr init = NULL, cond = NULL, inc = NULL;
      for (int j = 0; j < ev->nexpr; j++) {
 int i = -123456; double t = - TEND_EVENT;
 (* ev->expr[j]) (&i, &t, ev);
 if (i == -123456 && t == - TEND_EVENT) {

   if (cond)
     event_error (ev, "events can only use a single condition");
   cond = ev->expr[j];
 }
 else {

   int i1 = i; double t1 = t;
   (* ev->expr[j]) (&i1, &t1, ev);
   if (i1 == i && t1 == t) {


     if (init)
       event_error (ev, "events can only use a single initialisation");
     init = ev->expr[j];
   }
   else {

     if (inc)
       event_error (ev, "events can only use a single increment");
     inc = ev->expr[j];
   }
 }
      }
      ev->expr[0] = init;
      ev->expr[1] = cond;
      ev->expr[2] = inc;
      ev->nexpr = 0;
    }
    ev->i = -1; ev->t = - TEND_EVENT;
    if (ev->expr[0]) {
      (* ev->expr[0]) (&ev->i, &ev->t, ev);
      if (ev->i == END_EVENT || ev->t == TEND_EVENT) {
 ev->i = END_EVENT; ev->t = - TEND_EVENT;
      }
    }
    else if (ev->expr[2]) {
      (* ev->expr[2]) (&ev->i, &ev->t, ev);
      if (ev->i != -1)
 ev->i = 0;
      if (ev->t != - TEND_EVENT)
 ev->t = 0;
    }
  }
}

enum { event_done, event_alive, event_stop };

static int event_finished (Event * ev)
{
  ev->i = -1; ev->t = - TEND_EVENT;
  return event_done;
}

void event_register (Event event) {
  if (!(Events)) qassert ("/home/spencer/basilisk/src/grid/events.h", 0, "Events");
  if (!(!event.last)) qassert ("/home/spencer/basilisk/src/grid/events.h", 0, "!event.last");
  int n = 0, parent = -1;
  for (Event * ev = Events; !ev->last; ev++) {
    if (!strcmp (event.name, ev->name)) {
      if (!(parent < 0)) qassert ("/home/spencer/basilisk/src/grid/events.h", 0, "parent < 0");
      parent = n;
    }
    n++;
  }
  if (parent < 0) {
    Events = (Event *) prealloc (Events, (n + 2)*sizeof(Event),__func__,__FILE__,0);
    Events[n] = event;
    Events[n].next = NULL;
    Events[n + 1].last = true;
    init_event (&Events[n]);
  }
  else {
    Event * ev = ((Event *) pcalloc (1, sizeof(Event),__func__,__FILE__,0));
    *ev = Events[parent];
    Events[parent] = event;
    Events[parent].next = ev;
    init_event (&Events[parent]);
  }
}

static int event_cond (Event * ev, int i, double t)
{
  if (!ev->expr[1])
    return true;
  return (* ev->expr[1]) (&i, &t, ev);
}
# 136 "/home/spencer/basilisk/src/grid/events.h"
static bool overload_event() { return true; }

static int event_do (Event * ev, bool action)
{
  if ((iter > ev->i && t > ev->t) || !event_cond (ev, iter, t))
    return event_finished (ev);
  if (!overload_event() || iter == ev->i || fabs (t - ev->t) <= TEPS*t) {
    if (action) {
      bool finished = false;
      for (Event * e = ev; e; e = e->next) {



 if ((* e->action) (iter, t, e))
   finished = true;
      }
      if (finished) {
 event_finished (ev);
 return event_stop;
      }
    }
    if (ev->arrayi) {
      ev->i = ev->arrayi[ev->a++];
      if (ev->i < 0)
 return event_finished (ev);
    }
    if (ev->arrayt) {
      ev->t = ev->arrayt[ev->a++];
      if (ev->t < 0)
 return event_finished (ev);
    }
    else if (ev->expr[2]) {
      int i0 = ev->i;
      (* ev->expr[2]) (&ev->i, &ev->t, ev);
      if (i0 == -1 && ev->i != i0)
 ev->i += iter + 1;
      if (!event_cond (ev, iter + 1, ev->t))
 return event_finished (ev);
    }
    else if (ev->expr[0] && !ev->expr[1])
      return event_finished (ev);
  }
  return event_alive;
}

static void end_event_do (bool action)
{




  for (Event * ev = Events; !ev->last; ev++)
    if (ev->i == END_EVENT && action)
      for (Event * e = ev; e; e = e->next) {



 e->action (iter, t, e);
      }
}

int events (bool action)
{





  if (iter == 0)
    for (Event * ev = Events; !ev->last; ev++)
      init_event (ev);

  int cond = 0, cond1 = 0;
  inext = END_EVENT; tnext = 1e30;
  for (Event * ev = Events; !ev->last && !cond; ev++)
    if (ev->i != END_EVENT &&
 (ev->expr[1] || (ev->expr[0] && !ev->expr[1] && !ev->expr[2]) || ev->arrayi || ev->arrayt))
      cond = 1;
  for (Event * ev = Events; !ev->last; ev++) {
    int status = event_do (ev, action);
    if (status == event_stop) {
      end_event_do (action);
      return 0;
    }
    if (status == event_alive && ev->i != END_EVENT &&
 (ev->expr[1] || (ev->expr[0] && !ev->expr[1] && !ev->expr[2]) || ev->arrayi || ev->arrayt))
      cond1 = 1;
    if (ev->t > t && ev->t < tnext)
      tnext = ev->t;
    if (ev->i > iter && ev->i < inext)
      inext = ev->i;
  }
  if (overload_event() && (!cond || cond1) && (tnext != 1e30 || inext != END_EVENT)) {
    inext = iter + 1;
    return 1;
  }
  end_event_do (action);
  return 0;
}

void event (const char * name)
{
  for (Event * ev = Events; !ev->last; ev++)
    if (!strcmp (ev->name, name))
      for (Event * e = ev; e; e = e->next) {



 (* e->action) (0, 0, e);
      }
}

double dtnext (double dt)
{
  if (tnext != 1e30 && tnext > t) {
    unsigned int n = (tnext - t)/dt;
    if (!(n < INT_MAX)) qassert ("/home/spencer/basilisk/src/grid/events.h", 0, "n < INT_MAX");
    if (n == 0)
      dt = tnext - t;
    else {
      double dt1 = (tnext - t)/n;
      if (dt1 > dt*(1. + TEPS))
 dt = (tnext - t)/(n + 1);
      else if (dt1 < dt)
 dt = dt1;
      tnext = t + dt;
    }
  }
  else
    tnext = t + dt;
  return dt;
}
# 2 "/home/spencer/basilisk/src/grid/cartesian-common.h" 2

void (* debug) (Point);

@define _val_constant(a,k,l,m) ((const double) _constant[a.i -_NVARMAX])
@define diagonalize(a)
@define val_diagonal(a,k,l,m) ((k) == 0 && (l) == 0 && (m) == 0)

@undef VARIABLES
@def VARIABLES
  double Delta = L0*(1./(1 << point.level));
  double Delta_x = Delta;

  double Delta_y = Delta;





  double x = (ig/2. + (point.i - 2) + 0.5)*Delta + X0; NOT_UNUSED(x);

  double y = (jg/2. + (point.j - 2) + 0.5)*Delta + Y0;



 NOT_UNUSED(y);



  double z = 0.;

  NOT_UNUSED(z);

  NOT_UNUSED(Delta);
  NOT_UNUSED(Delta_x);

  NOT_UNUSED(Delta_y);





  ;
@

# 1 "grid/fpe.h" 1
# 1 "/home/spencer/basilisk/src/grid/fpe.h"


@include <signal.h>
@include <unistd.h>

static int gdb()
{
  if (last_point.level >= 0) {
    debug (last_point);
    fputc ('\n', ferr);
    fflush (ferr);
  }
  char command[80];
  sprintf (command, "exec xterm -e 'gdb -p %d' & xterm -e 'gnuplot plot -'",
    getpid());
  return system (command);
}

static void caught_abort (int sig)
{
  fprintf (ferr, "Caught signal %d (Aborted)\n", sig);
  gdb();
}

static void caught_fpe (int sig)
{
  fprintf (ferr, "Caught signal %d (Floating Point Exception)\n", sig);
  gdb();
  exit (1);
}

static void caught_segfault (int sig)
{
  fprintf (ferr, "Caught signal %d (Segmentation Fault)\n", sig);
  gdb();
  exit (2);
}

void catch_fpe (void)
{
  struct sigaction act;
  act.sa_handler = caught_fpe;
  sigemptyset (&act.sa_mask);
  act.sa_flags = 0;
  last_point.level = -1;
  sigaction (8, &act, NULL);
  act.sa_handler = caught_segfault;
  sigaction (11, &act, NULL);
  act.sa_handler = caught_abort;
  act.sa_flags = SA_RESETHAND;
  sigaction (6, &act, NULL);
}
# 47 "/home/spencer/basilisk/src/grid/cartesian-common.h" 2

@define end_foreach_face()

@def foreach_point(...)
{
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
  coord _p = { S__VA_ARGS__ };
  Point point = locate (_p.x, _p.y, _p.z);
  if (point.level >= 0) {
    POINT_VARIABLES
@
@define end_foreach_point() }}

@def foreach_region(p, box, n)
  OMP_PARALLEL() { NOT_UNUSED (p);
    coord p = {0, 0, box[0].z};
    OMP(omp for schedule(static))
      for (int _i = 0; _i < (int) n.x; _i++) {
 p.x = box[0].x + (box[1].x - box[0].x)/n.x*(_i + 0.5);
 for (int _j = 0; _j < (int) n.y; _j++) {
   p.y = box[0].y + (box[1].y - box[0].y)/n.y*(_j + 0.5);
   Point point = locate (p.x, p.y, p.z);
   if (point.level >= 0) {
     int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
     POINT_VARIABLES
@
@define end_foreach_region() }}}}





static void init_block_scalar (scalar sb, const char * name, const char * ext,
          int n, int block)
{
  char bname[strlen(name) + strlen(ext) + 10];
  if (n == 0) {
    sprintf (bname, "%s%s", name, ext);
    _attribute[sb.i].block = block;
    baseblock = list_append (baseblock, sb);
  }
  else {
    sprintf (bname, "%s%d%s", name, n, ext);
    _attribute[sb.i].block = - n;
  }
  _attribute[sb.i].name = pstrdup (bname,__func__,__FILE__,0);
  all = list_append (all, sb);
}

@define interpreter_set_int(...)
@define interpreter_reset_scalar(...)

scalar alloc_block_scalar (const char * name, const char * ext, int block)
{
  interpreter_set_int (&block);
  int nvar = datasize/sizeof(double);

  scalar s = {0};
  while (s.i < nvar) {
    int n = 0;
    scalar sb = s;
    while (sb.i < nvar && n < block && _attribute[sb.i].freed)
      n++, sb.i++;
    if (n >= block) {
      memset (&_attribute[s.i], 0, block*sizeof (_Attributes));
      for (sb.i = s.i, n = 0; n < block; n++, sb.i++) {
 init_block_scalar (sb, name, ext, n, block);
 interpreter_reset_scalar (sb);
      }
      trash (((scalar []){s, {-1}}));
      return s;
    }
    s.i = sb.i + 1;
  }


  s = (scalar){nvar};
  if (!(nvar + block <= _NVARMAX)) qassert ("/home/spencer/basilisk/src/grid/cartesian-common.h", 0, "nvar + block <= _NVARMAX");

  if (_attribute == NULL)
    _attribute = (_Attributes *) pcalloc (nvar + block + 1, sizeof (_Attributes),__func__,__FILE__,0);
  else
    _attribute = (_Attributes *)
      prealloc (_attribute, (nvar + block + 1)*sizeof (_Attributes),__func__,__FILE__,0);
  memset (&_attribute[nvar], 0, block*sizeof (_Attributes));
  for (int n = 0; n < block; n++, nvar++) {
    scalar sb = (scalar){nvar};
    init_block_scalar (sb, name, ext, n, block);
  }

  realloc_scalar (block*sizeof(double));
  trash (((scalar []){s, {-1}}));
  return s;
}

scalar new_block_scalar (const char * name, const char * ext, int block)
{
  scalar s = alloc_block_scalar (name, ext, block), sb;
  int n = 0;
  for (sb.i = s.i, n = 0; n < block; n++, sb.i++)
    init_scalar (sb, NULL);
  return s;
}

scalar new_scalar (const char * name)
{
  return init_scalar (alloc_block_scalar (name, "", 1), NULL);
}

scalar new_vertex_scalar (const char * name)
{
  return init_vertex_scalar (alloc_block_scalar (name, "", 1), NULL);
}

static vector alloc_block_vector (const char * name, int block)
{
  vector v;
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
  
    v.x = alloc_block_scalar (name, ext.x, block);
    
#line 166
v.y = alloc_block_scalar (name, ext.y, block);
  return v;
}

vector new_vector (const char * name)
{
  vector v = alloc_block_vector (name, 1);
  init_vector (v, NULL);
  return v;
}

vector new_face_vector (const char * name)
{
  vector v = alloc_block_vector (name, 1);
  init_face_vector (v, NULL);
  return v;
}

vector new_block_vector (const char * name, int block)
{
  vector v = alloc_block_vector (name, block);
  for (int i = 0; i < block; i++) {
    vector vb;
    
      vb.x.i = v.x.i + i;
      
#line 190
vb.y.i = v.y.i + i;
    init_vector (vb, NULL);
    
      _attribute[vb.x.i].block = - i;
      
#line 193
_attribute[vb.y.i].block = - i;
  }
  
    _attribute[v.x.i].block = block;
    
#line 196
_attribute[v.y.i].block = block;
  return v;
}

vector new_block_face_vector (const char * name, int block)
{
  vector v = alloc_block_vector (name, block);
  for (int i = 0; i < block; i++) {
    vector vb;
    
      vb.x.i = v.x.i + i;
      
#line 206
vb.y.i = v.y.i + i;
    init_face_vector (vb, NULL);
    
      _attribute[vb.x.i].block = - i;
      
#line 209
_attribute[vb.y.i].block = - i;
  }
  
    _attribute[v.x.i].block = block;
    
#line 212
_attribute[v.y.i].block = block;
  return v;
}

tensor new_tensor (const char * name)
{
  char cname[strlen(name) + 3];
  struct { char * x, * y, * z; } ext = {"%s.x", "%s.y", "%s.z"};
  tensor t;
   {
    sprintf (cname, ext.x, name);
    t.x = alloc_block_vector (cname, 1);
  } 
#line 221
{
    sprintf (cname, ext.y, name);
    t.y = alloc_block_vector (cname, 1);
  }
  init_tensor (t, NULL);
  return t;
}

tensor new_symmetric_tensor (const char * name)
{
  struct { char * x, * y, * z; } ext = {".x.x", ".y.y", ".z.z"};
  tensor t;
  
    t.x.x = alloc_block_scalar (name, ext.x, 1);
    
#line 234
t.y.y = alloc_block_scalar (name, ext.y, 1);

    t.x.y = alloc_block_scalar (name, ".x.y", 1);
    t.y.x = t.x.y;
# 248 "/home/spencer/basilisk/src/grid/cartesian-common.h"
  init_tensor (t, NULL);
  return t;
}

static int nconst = 0;

void init_const_scalar (scalar s, const char * name, double val)
{
  if (s.i - _NVARMAX >= nconst) {
    nconst = s.i - _NVARMAX + 1;
    _constant = (double *) prealloc (_constant, (nconst)*sizeof(double),__func__,__FILE__,0);
  }
  _constant[s.i - _NVARMAX] = val;
}

scalar new_const_scalar (const char * name, int i, double val)
{
  scalar s = (scalar){i + _NVARMAX};
  init_const_scalar (s, name, val);
  return s;
}

void init_const_vector (vector v, const char * name, double * val)
{
  
    init_const_scalar (v.x, name, *val++);
    
#line 273
init_const_scalar (v.y, name, *val++);
}

vector new_const_vector (const char * name, int i, double * val)
{
  vector v;
  
    v.x.i = _NVARMAX + i++;
    
#line 280
v.y.i = _NVARMAX + i++;
  init_const_vector (v, name, val);
  return v;
}

static void cartesian_scalar_clone (scalar clone, scalar src)
{
  char * cname = _attribute[clone.i].name;
  double (** boundary) (Point, Point, scalar, void *) = _attribute[clone.i].boundary;
  double (** boundary_homogeneous) (Point, Point, scalar, void *) =
    _attribute[clone.i].boundary_homogeneous;
  if (!(_attribute[src.i].block > 0 && _attribute[clone.i].block == _attribute[src.i].block)) qassert ("/home/spencer/basilisk/src/grid/cartesian-common.h", 0, "src.block > 0 && clone.block == src.block");
  pfree (_attribute[clone.i].depends,__func__,__FILE__,0);
  _attribute[clone.i] = _attribute[src.i];
  _attribute[clone.i].name = cname;
  _attribute[clone.i].boundary = boundary;
  _attribute[clone.i].boundary_homogeneous = boundary_homogeneous;
  for (int i = 0; i < nboundary; i++) {
    _attribute[clone.i].boundary[i] = _attribute[src.i].boundary[i];
    _attribute[clone.i].boundary_homogeneous[i] = _attribute[src.i].boundary_homogeneous[i];
  }
  _attribute[clone.i].depends = list_copy (_attribute[src.i].depends);
}

scalar * list_clone (scalar * l)
{
  scalar * list = NULL;
  int nvar = datasize/sizeof(double), map[nvar];
  for (int i = 0; i < nvar; i++)
    map[i] = -1;
  {scalar*_i=(scalar*)( l);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    scalar c = _attribute[s.i].block > 1 ? new_block_scalar("c", "", _attribute[s.i].block) : new_scalar("c");
    scalar_clone (c, s);
    map[s.i] = c.i;
    list = list_append (list, c);
  }}}
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    {
      if (_attribute[s.i].v.x.i >= 0 && map[_attribute[s.i].v.x.i] >= 0)
 _attribute[s.i].v.x.i = map[_attribute[s.i].v.x.i];
      
#line 318
if (_attribute[s.i].v.y.i >= 0 && map[_attribute[s.i].v.y.i] >= 0)
 _attribute[s.i].v.y.i = map[_attribute[s.i].v.y.i];}}}
  return list;
}

void delete (scalar * list)
{
  if (all == NULL)
    return;

  {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){ {
    for (int i = 0; i < _attribute[f.i].block; i++) {
      scalar fb = {f.i + i};
      if (_attribute[f.i].delete)
 _attribute[f.i].delete (fb);
      pfree (_attribute[fb.i].name,__func__,__FILE__,0); _attribute[fb.i].name = NULL;
      pfree (_attribute[fb.i].boundary,__func__,__FILE__,0); _attribute[fb.i].boundary = NULL;
      pfree (_attribute[fb.i].boundary_homogeneous,__func__,__FILE__,0); _attribute[fb.i].boundary_homogeneous = NULL;
      pfree (_attribute[fb.i].depends,__func__,__FILE__,0); _attribute[fb.i].depends = NULL;
      _attribute[fb.i].freed = true;
    }
  }}}

  if (list == all) {
    all[0].i = -1;
    baseblock[0].i = -1;
    return;
  }

  trash (list);
  {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){ {
    if (_attribute[f.i].block > 0) {
      scalar * s;
      for (s = all; s->i >= 0 && s->i != f.i; s++);
      if (s->i == f.i) {
 for (; s[_attribute[f.i].block].i >= 0; s++)
   s[0] = s[_attribute[f.i].block];
 s->i = -1;
      }
      for (s = baseblock; s->i >= 0 && s->i != f.i; s++);
      if (s->i == f.i) {
 for (; s[1].i >= 0; s++)
   s[0] = s[1];
 s->i = -1;
      }
    }
  }}}
}

void free_solver()
{
  if (!(_val_higher_dimension == 0.)) qassert ("/home/spencer/basilisk/src/grid/cartesian-common.h", 0, "_val_higher_dimension == 0.");

  if (free_solver_funcs) {
    free_solver_func * a = (free_solver_func *) free_solver_funcs->p;
    for (int i = 0; i < free_solver_funcs->len/sizeof(free_solver_func); i++)
      a[i] ();
    array_free (free_solver_funcs);
  }

  delete (all);
  pfree (all,__func__,__FILE__,0); all = NULL;
  pfree (baseblock,__func__,__FILE__,0); baseblock = NULL;
  for (Event * ev = Events; !ev->last; ev++) {
    Event * e = ev->next;
    while (e) {
      Event * next = e->next;
      pfree (e,__func__,__FILE__,0);
      e = next;
    }
  }

  pfree (Events,__func__,__FILE__,0); Events = NULL;
  pfree (_attribute,__func__,__FILE__,0); _attribute = NULL;
  pfree (_constant,__func__,__FILE__,0); _constant = NULL;
  free_grid();
  qpclose_all();
@if TRACE
  trace_off();
@endif
@if MTRACE
  pmuntrace();
@endif
@if _CADNA
  cadna_end();
@endif
}



void (* boundary_level) (scalar *, int l);
void (* boundary_face) (vectorl);




void boundary_flux (vector * list) __attribute__ ((deprecated));

void boundary_flux (vector * list)
{
  vectorl list1 = {NULL};
  {vector*_i=(vector*)( list);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
    {
      list1.x = list_append (list1.x, v.x);
      
#line 421
list1.y = list_append (list1.y, v.y);}}}
  boundary_face (list1);
  
    pfree (list1.x,__func__,__FILE__,0);
    
#line 424
pfree (list1.y,__func__,__FILE__,0);
}

static scalar * list_add_depends (scalar * list, scalar s)
{
  {scalar*_i=(scalar*)( list);if(_i)for(scalar t=*_i;(&t)->i>=0;t=*++_i){
    if (t.i == s.i)
      return list;}}
  scalar * list1 = list;
  {scalar*_i=(scalar*)( _attribute[s.i].depends);if(_i)for(scalar d=*_i;(&d)->i>=0;d=*++_i){
    if (_attribute[d.i].dirty)
      list1 = list_add_depends (list1, d);}}
  return list_append (list1, s);
}

     
void boundary_internal (scalar * list, const char * fname, int line)
{tracing("boundary_internal","/home/spencer/basilisk/src/grid/cartesian-common.h",0);
  if (list == NULL)
    {end_tracing("boundary_internal","/home/spencer/basilisk/src/grid/cartesian-common.h",0);return;}
  scalar * listc = NULL;
  vectorl listf = {NULL};
  bool flux = false;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant(s) && _attribute[s.i].block > 0) {
      if (scalar_is_dirty (s)) {
 if (_attribute[s.i].face && _attribute[s.i].dirty != 2)
   {
     if (_attribute[s.i].v.x.i == s.i)
       listf.x = list_add (listf.x, s), flux = true;
     
#line 452
if (_attribute[s.i].v.y.i == s.i)
       listf.y = list_add (listf.y, s), flux = true;}
 if (!is_constant(cm) && _attribute[cm.i].dirty)
   listc = list_add_depends (listc, cm);
 if (_attribute[s.i].face != 2)
   listc = list_add_depends (listc, s);
      }




    }}}
  if (flux) {
    boundary_face (listf);
    
      pfree (listf.x,__func__,__FILE__,0);
      
#line 467
pfree (listf.y,__func__,__FILE__,0);
  }
  if (listc) {
    boundary_level (listc, -1);
    {scalar*_i=(scalar*)( listc);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = false;}}
    pfree (listc,__func__,__FILE__,0);
  }
end_tracing("boundary_internal","/home/spencer/basilisk/src/grid/cartesian-common.h",0);}

void cartesian_boundary_level (scalar * list, int l)
{
  { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list, l); };
}

void cartesian_boundary_face (vectorl list)
{
  
    {scalar*_i=(scalar*)( list.x);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = 2;}}
    
#line 485
{scalar*_i=(scalar*)( list.y);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = 2;}}
}

static double symmetry (Point point, Point neighbor, scalar s, void * data)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  return val(s,0,0,0);
}

static double antisymmetry (Point point, Point neighbor, scalar s, void * data)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  return -val(s,0,0,0);
}

double (* default_scalar_bc[]) (Point, Point, scalar, void *) = {
  symmetry, symmetry, symmetry, symmetry, symmetry, symmetry
};

scalar cartesian_init_scalar (scalar s, const char * name)
{

  char * pname;
  if (name) {
    pfree (_attribute[s.i].name,__func__,__FILE__,0);
    pname = pstrdup (name,__func__,__FILE__,0);
  }
  else
    pname = _attribute[s.i].name;
  int block = _attribute[s.i].block;
  double (** boundary) (Point, Point, scalar, void *) = _attribute[s.i].boundary;
  double (** boundary_homogeneous) (Point, Point, scalar, void *) =
    _attribute[s.i].boundary_homogeneous;
  _attribute[s.i].name = pname;
  if (block < 0)
    _attribute[s.i].block = block;
  else
    _attribute[s.i].block = block > 0 ? block : 1;

  _attribute[s.i].boundary = boundary ? boundary :
    (double (**)(Point, Point, scalar, void *))
    pmalloc (nboundary*sizeof (void (*)()),__func__,__FILE__,0);
  _attribute[s.i].boundary_homogeneous = boundary_homogeneous ? boundary_homogeneous :
    (double (**)(Point, Point, scalar, void *))
    pmalloc (nboundary*sizeof (void (*)()),__func__,__FILE__,0);
  for (int b = 0; b < nboundary; b++)
    _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b] =
      b < 2*2 ? default_scalar_bc[b] : symmetry;
  _attribute[s.i].gradient = NULL;
   {
    _attribute[s.i].d.x = 0;
    _attribute[s.i].v.x.i = -1;
  } 
#line 533
{
    _attribute[s.i].d.y = 0;
    _attribute[s.i].v.y.i = -1;
  }
  _attribute[s.i].face = false;
  return s;
}

scalar cartesian_init_vertex_scalar (scalar s, const char * name)
{
  cartesian_init_scalar (s, name);
  
    _attribute[s.i].d.x = -1;
    
#line 545
_attribute[s.i].d.y = -1;
  for (int d = 0; d < nboundary; d++)
    _attribute[s.i].boundary[d] = _attribute[s.i].boundary_homogeneous[d] = NULL;
  return s;
}

double (* default_vector_bc[]) (Point, Point, scalar, void *) = {
  antisymmetry, antisymmetry,
  antisymmetry, antisymmetry,
  antisymmetry, antisymmetry
};

vector cartesian_init_vector (vector v, const char * name)
{
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
   {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.x);
      cartesian_init_scalar (v.x, cname);
    }
    else
      cartesian_init_scalar (v.x, NULL);
    _attribute[v.x.i].v = v;
  } 
#line 560
{
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.y);
      cartesian_init_scalar (v.y, cname);
    }
    else
      cartesian_init_scalar (v.y, NULL);
    _attribute[v.y.i].v = v;
  }

  for (int d = 0; d < nboundary; d++)
    _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] =
      d < 2*2 ? default_vector_bc[d] : antisymmetry;
  return v;
}

vector cartesian_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_vector (v, name);
   {
    _attribute[v.x.i].d.x = -1;
    _attribute[v.x.i].face = true;
  } 
#line 580
{
    _attribute[v.y.i].d.y = -1;
    _attribute[v.y.i].face = true;
  }
  for (int d = 0; d < nboundary; d++)
    _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] = NULL;
  return v;
}

tensor cartesian_init_tensor (tensor t, const char * name)
{
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
   {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.x);
      cartesian_init_vector (t.x, cname);
    }
    else
      cartesian_init_vector (t.x, NULL);
  } 
#line 592
{
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.y);
      cartesian_init_vector (t.y, cname);
    }
    else
      cartesian_init_vector (t.y, NULL);
  }






    for (int b = 0; b < nboundary; b++) {
      _attribute[t.x.x.i].boundary[b] = _attribute[t.y.x.i].boundary[b] =
 _attribute[t.x.x.i].boundary_homogeneous[b] = _attribute[t.y.y.i].boundary_homogeneous[b] =
 b < 2*2 ? default_scalar_bc[b] : symmetry;
      _attribute[t.x.y.i].boundary[b] = _attribute[t.y.y.i].boundary[b] =
 _attribute[t.x.y.i].boundary_homogeneous[b] = _attribute[t.y.x.i].boundary_homogeneous[b] =
 b < 2*2 ? default_vector_bc[b] : antisymmetry;
    }



  return t;
}

void output_cells (FILE * fp, coord c, double size)
{
  {foreach() {
    bool inside = true;
    coord o = {x,y,z};
    
      if (inside && size > 0. &&
   (o.x > c.x + size || o.x < c.x - size))
 inside = false;
      
#line 627
if (inside && size > 0. &&
   (o.y > c.y + size || o.y < c.y - size))
 inside = false;
    if (inside) {
      Delta /= 2.;



      fprintf (fp, "%g %g\n%g %g\n%g %g\n%g %g\n%g %g\n\n",
        x - Delta, y - Delta,
        x - Delta, y + Delta,
        x + Delta, y + Delta,
        x + Delta, y - Delta,
        x - Delta, y - Delta);
# 655 "/home/spencer/basilisk/src/grid/cartesian-common.h"
    }
  }end_foreach();}
  fflush (fp);
}
# 667 "/home/spencer/basilisk/src/grid/cartesian-common.h"
static char * replace_ (const char * vname)
{
  char * name = pstrdup (vname,__func__,__FILE__,0), * c = name;
  while (*c != '\0') {
    if (*c == '.')
      *c = '_';
    c++;
  }
  return name;
}

static void debug_plot (FILE * fp, const char * name, const char * cells,
   const char * stencil)
{
  char * vname = replace_ (name);
  fprintf (fp,
    "  load 'debug.plot'\n"
    "  v=%s\n"




    "  plot '%s' w l lc 0, "
    "'%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 1 title columnhead(3+3*v)",





    vname, cells, stencil);
  pfree (vname,__func__,__FILE__,0);
}

void cartesian_debug (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  char name[80] = "cells";
  if (pid() > 0)
    sprintf (name, "cells-%d", pid());
  FILE * fp = fopen (name, "w");
  output_cells (fp, (coord){x,y,z}, 4.*Delta);
  fclose (fp);

  char stencil[80] = "stencil";
  if (pid() > 0)
    sprintf (stencil, "stencil-%d", pid());
  fp = fopen (stencil, "w");
  {scalar*_i=(scalar*)( all);if(_i)for(scalar v=*_i;(&v)->i>=0;v=*++_i){



    fprintf (fp, "x y %s ", _attribute[v.i].name);}}



  fputc ('\n', fp);
# 734 "/home/spencer/basilisk/src/grid/cartesian-common.h"
    for (int k = -2; k <= 2; k++)
      for (int l = -2; l <= 2; l++) {
 {scalar*_i=(scalar*)( all);if(_i)for(scalar v=*_i;(&v)->i>=0;v=*++_i){ {
   fprintf (fp, "%g %g ",
     x + k*Delta + _attribute[v.i].d.x*Delta/2.,
     y + l*Delta + _attribute[v.i].d.y*Delta/2.);
   if (allocated(k,l,0))
     fprintf (fp, "%g ", val(v,k,l,0));
   else
     fputs ("n/a ", fp);
 }}}
 fputc ('\n', fp);
      }
# 764 "/home/spencer/basilisk/src/grid/cartesian-common.h"
  fclose (fp);

  fp = fopen ("debug.plot", "w");
  fprintf (fp,
    "set term x11\n"
    "set size ratio -1\n"
    "set key outside\n");
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    char * name = replace_ (_attribute[s.i].name);
    fprintf (fp, "%s = %d\n", name, s.i);
    pfree (name,__func__,__FILE__,0);
  }}}
  fclose (fp);

  fprintf (ferr, "Last point stencils can be displayed using (in gnuplot)\n");
  debug_plot (ferr, _attribute[0].name, name, stencil);
  fflush (ferr);

  fp = fopen ("plot", "w");
  debug_plot (fp, _attribute[0].name, name, stencil);
  fclose (fp);
}

void cartesian_methods()
{
  init_scalar = cartesian_init_scalar;
  init_vertex_scalar = cartesian_init_vertex_scalar;
  init_vector = cartesian_init_vector;
  init_face_vector = cartesian_init_face_vector;
  init_tensor = cartesian_init_tensor;
  boundary_level = cartesian_boundary_level;
  boundary_face = cartesian_boundary_face;
  scalar_clone = cartesian_scalar_clone;
  debug = cartesian_debug;
}

tensor init_symmetric_tensor (tensor t, const char * name)
{
  return init_tensor (t, name);
}

static double interpolate_linear (Point point, scalar v,
      double xp, double yp, double zp)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;







  x = (xp - x)/Delta - _attribute[v.i].d.x/2.;
  y = (yp - y)/Delta - _attribute[v.i].d.y/2.;
  int i = sign(x), j = sign(y);
  x = fabs(x); y = fabs(y);

  return ((val(v,0,0,0)*(1. - x) + val(v,i,0,0)*x)*(1. - y) +
   (val(v,0,j,0)*(1. - x) + val(v,i,j,0)*x)*y);
# 834 "/home/spencer/basilisk/src/grid/cartesian-common.h"
}


#line 805
static void _stencil_interpolate_linear (Point point, scalar v,
_stencil_undefined * xp,_stencil_undefined * yp,_stencil_undefined * zp)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;      







        
        
  
       

_stencil_val(v,0,0,0);_stencil_val(v, o_stencil,0,0);
_stencil_val(v,0,o_stencil,0); _stencil_val(v,o_stencil,o_stencil,0);

  
#line 820
return         
    ;
# 834 "/home/spencer/basilisk/src/grid/cartesian-common.h"
}

     
double interpolate (scalar v, double xp, double yp, double zp,
      bool linear)
{tracing("interpolate","/home/spencer/basilisk/src/grid/cartesian-common.h",0);
  double val = 1e30;
  foreach_point_stencil (1,{(NonLocal[]){{"v","scalar",(void *)&v,NULL,0},{"linear","bool",(void *)&linear,NULL,0},{"val","double",(void *)&val,NULL,0,'m'},{"zp","double",(void *)&zp,NULL,0},{"yp","double",(void *)&yp,NULL,0},{"xp","double",(void *)&xp,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_("\n\n\n// #line 805 \"/home/spencer/basilisk/src/grid/cartesian-common.h\"\nstatic real interpolate_linear (Point point, scalar v,\n      real xp, real yp, real zp)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n\n\n\n\n\n\n\n  x = (xp - x)/Delta - v.d.x/2.;\n  y = (yp - y)/Delta - v.d.y/2.;\n  int i = sign(x), j = sign(y);\n  x = fabs(x); y = fabs(y);\n\n  return ((val(v,0,0,0)*(1. - x) + val(v,i,0,0)*x)*(1. - y) +\n   (val(v,0,j,0)*(1. - x) + val(v,i,j,0)*x)*y);\n// # 834 \"/home/spencer/basilisk/src/grid/cartesian-common.h\"\n}"),_("\n    \n// #line 842 \"/home/spencer/basilisk/src/grid/cartesian-common.h\"\nval = linear ? interpolate_linear (point, v, xp, yp, zp) : val(v,0,0,0);")})
    { _stencil_interpolate_linear (point, v, NULL, NULL, NULL); _stencil_val(v,0,0,0);    }end_foreach_point_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction (min:val)){
#line 841
foreach_point (xp, yp, zp)
    val = linear ? interpolate_linear (point, v, xp, yp, zp) : val(v,0,0,0);end_foreach_point();mpi_all_reduce_array(&val,double,MPI_MIN,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
  
#line 843
{end_tracing("interpolate","/home/spencer/basilisk/src/grid/cartesian-common.h",0);return val;}
end_tracing("interpolate","/home/spencer/basilisk/src/grid/cartesian-common.h",0);}

     
void interpolate_array (scalar * list, coord * a, int n, double * v,
   bool linear)
{tracing("interpolate_array","/home/spencer/basilisk/src/grid/cartesian-common.h",0);
  int len = 0;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    len++;}}
  for (int i = 0; i < n; i++) {
    double * w = v;
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      *(w++) = 1e30;}}
    foreach_point_stencil (1,{(NonLocal[]){{"linear","bool",(void *)&linear,NULL,0},{"v","double",(void *)v,NULL,1},{"list","scalar",(void *)list,NULL,1},{"len","int",(void *)&len,NULL,0},{"i","int",(void *)&i,NULL,0},{"a","coord",(void *)a,NULL,1},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_("\n\n\n// #line 805 \"/home/spencer/basilisk/src/grid/cartesian-common.h\"\nstatic real interpolate_linear (Point point, scalar v,\n      real xp, real yp, real zp)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n\n\n\n\n\n\n\n  x = (xp - x)/Delta - v.d.x/2.;\n  y = (yp - y)/Delta - v.d.y/2.;\n  int i = sign(x), j = sign(y);\n  x = fabs(x); y = fabs(y);\n\n  return ((val(v,0,0,0)*(1. - x) + val(v,i,0,0)*x)*(1. - y) +\n   (val(v,0,j,0)*(1. - x) + val(v,i,j,0)*x)*y);\n// # 834 \"/home/spencer/basilisk/src/grid/cartesian-common.h\"\n}"),_(" \n// #line 857 \"/home/spencer/basilisk/src/grid/cartesian-common.h\"\n{\n      int j = 0;\n      {forin (scalar, s , list)\n v[j++] = !linear ? val(s,0,0,0) : interpolate_linear (point, s, a[i].x, a[i].y, a[i].z); endforin()}\n    }")}) {   
      
      {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 { _stencil_val(s,0,0,0); _stencil_interpolate_linear (point, s, NULL, NULL, NULL);    }}}
    }end_foreach_point_stencil();
    
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(min:v[:len])){
#line 857
foreach_point (a[i].x, a[i].y, a[i].z) {
      int j = 0;
      {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 v[j++] = !linear ? val(s,0,0,0) : interpolate_linear (point, s, a[i].x, a[i].y, a[i].z);}}
    }end_foreach_point();mpi_all_reduce_array(v,double,MPI_MIN,len);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
    
#line 862
v = w;
  }
end_tracing("interpolate_array","/home/spencer/basilisk/src/grid/cartesian-common.h",0);}



typedef int bid;

bid new_bid()
{
  int b = nboundary++;
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    _attribute[s.i].boundary = (double (**)(Point, Point, scalar, void *))
      prealloc (_attribute[s.i].boundary, nboundary*sizeof (void (*)()),__func__,__FILE__,0);
    _attribute[s.i].boundary_homogeneous = (double (**)(Point, Point, scalar, void *))
      prealloc (_attribute[s.i].boundary_homogeneous, nboundary*sizeof (void (*)()),__func__,__FILE__,0);
  }}}
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    if (_attribute[s.i].v.x.i < 0)
      _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b] = symmetry;
    else if (_attribute[s.i].v.x.i == s.i) {
      vector v = _attribute[s.i].v;
      
 _attribute[v.y.i].boundary[b] = _attribute[v.y.i].boundary_homogeneous[b] = symmetry;
 
#line 885
_attribute[v.x.i].boundary[b] = _attribute[v.x.i].boundary_homogeneous[b] = symmetry;
      _attribute[v.x.i].boundary[b] = _attribute[v.x.i].boundary_homogeneous[b] =
 _attribute[v.x.i].face ? NULL : antisymmetry;
    }
  }}}
  return b;
}



static double periodic_bc (Point point, Point neighbor, scalar s, void * data)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  return 1e30;
}

static void periodic_boundary (int d)
{

  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (is_vertex_scalar (s))
      _attribute[s.i].boundary[d] = _attribute[s.i].boundary_homogeneous[d] = NULL;
    else
      _attribute[s.i].boundary[d] = _attribute[s.i].boundary_homogeneous[d] = periodic_bc;}}

  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (_attribute[s.i].face) {
      vector v = _attribute[s.i].v;
      _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] = NULL;
    }}}

  default_scalar_bc[d] = periodic_bc;
  default_vector_bc[d] = periodic_bc;
}

void periodic (int dir)
{



    if (!(dir <= bottom)) qassert ("/home/spencer/basilisk/src/grid/cartesian-common.h", 0, "dir <= bottom");




  int c = dir/2;
  periodic_boundary (2*c);
  periodic_boundary (2*c + 1);
  (&Period.x)[c] = true;
}


double getvalue (Point point, scalar s, int i, int j, int k)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  return val(s,i,j,k);
}

void default_stencil (Point p, scalar * list)
{
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    _attribute[s.i].input = true, _attribute[s.i].width = 2;}}
}




static void write_stencil_index (int * index)
{
  fprintf (qstderr(), "[%d", index[0]);
  for (int d = 1; d < 2; d++)
    fprintf (qstderr(), ",%d", index[d]);
  fputs ("]", qstderr());
}

void stencil_val (Point p, scalar s, int i, int j, int k,
    const char * file, int line, bool overflow)
{
  if (is_constant(s) || s.i < 0)
    return;
  int index[] = {i, j, k};
  for (int d = 0; d < 2; d++)
    index[d] += (&p.i)[d];
  bool central = true;
  for (int d = 0; d < 2; d++) {
    if (!overflow && (index[d] > 2 || index[d] < - 2)) {
      fprintf (qstderr(), "%s:%d: error: stencil overflow: %s",
        file, line, _attribute[s.i].name);
      write_stencil_index (index);
      fprintf (qstderr(), "\n");
      fflush (qstderr());
      abort();
    }
    if (index[d] != 0)
      central = false;
  }
  if (central) {
    if (!_attribute[s.i].output)
      _attribute[s.i].input = true;
  }
  else {
    _attribute[s.i].input = true;
    int d = 0;
     {
      if ((!_attribute[s.i].face || _attribute[s.i].v.x.i != s.i) && abs(index[d]) > _attribute[s.i].width)
 _attribute[s.i].width = abs(index[d]);
      d++;
    } 
#line 986
{
      if ((!_attribute[s.i].face || _attribute[s.i].v.y.i != s.i) && abs(index[d]) > _attribute[s.i].width)
 _attribute[s.i].width = abs(index[d]);
      d++;
    }
  }
}

void stencil_val_a (Point p, scalar s, int i, int j, int k, bool input,
      const char * file, int line)
{
  if (is_constant(s) || s.i < 0)
    abort();
  int index[] = {i, j, k};
  for (int d = 0; d < 2; d++)
    index[d] += (&p.i)[d];
  for (int d = 0; d < 2; d++)
    if (index[d] != 0) {
      fprintf (qstderr(), "%s:%d: error: illegal write: %s",
        file, line, _attribute[s.i].name);
      write_stencil_index (index);
      fprintf (qstderr(), "\n");
      fflush (qstderr());
      abort();
    }
  if (input && !_attribute[s.i].output)
    _attribute[s.i].input = true;
  _attribute[s.i].output = true;
}




@define dimensional(...)

@define show_dimension_internal(...)
@define display_value(...)
@define interpreter_verbosity(...)
# 4 "/home/spencer/basilisk/src/grid/multigrid-common.h" 2

@ifndef foreach_level_or_leaf
@ define foreach_level_or_leaf foreach_level
@ define end_foreach_level_or_leaf end_foreach_level
@endif

@ifndef foreach_coarse_level
@ define foreach_coarse_level foreach_level
@ define end_foreach_coarse_level end_foreach_level
@endif










void (* restriction) (scalar *);

static inline void restriction_average (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  double sum = 0.;
  {foreach_child()
    sum += val(s,0,0,0);end_foreach_child()}
  val(s,0,0,0) = sum/(1 << 2);
}

static inline void restriction_volume_average (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

#line 35
if(!is_constant(cm)){{
  double sum = 0.;
  {foreach_child()
    sum += val(cm,0,0,0)*val(s,0,0,0);end_foreach_child()}
  val(s,0,0,0) = sum/(1 << 2)/(val(cm,0,0,0) + 1e-30);
}}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);

#line 35
{
  double sum = 0.;
  {foreach_child()
    sum += _const_cm*val(s,0,0,0);end_foreach_child()}
  val(s,0,0,0) = sum/(1 << 2)/(_const_cm + 1e-30);
}}

#line 40
}

static inline void face_average (Point point, vector v)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
   {




      val(v.x,0,0,0) = (fine(v.x,0,0,0) + fine(v.x,0,1,0))/2.;
      val(v.x,1,0,0) = (fine(v.x,2,0,0) + fine(v.x,2,1,0))/2.;






  } 
#line 44
{




      val(v.y,0,0,0) = (fine(v.y,0,0,0) + fine(v.y,1,0,0))/2.;
      val(v.y,0,1,0) = (fine(v.y,0,2,0) + fine(v.y,1,2,0))/2.;






  }
}

static inline void restriction_face (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  face_average (point, _attribute[s.i].v);
}

static inline void restriction_vertex (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  for (int i = 0; i <= 1; i++) {
    val(s,i,0,0) = fine(s,2*i,0,0);

    val(s,i,1,0) = fine(s,2*i,2,0);





  }
}

static inline void no_restriction (Point point, scalar s) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;}

static inline void no_data (Point point, scalar s) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  {foreach_child()
    val(s,0,0,0) = 1e30;end_foreach_child()}
}

void wavelet (scalar s, scalar w)
{
  restriction (((scalar[]){s,{-1}}));
  for (int l = grid->maxdepth - 1; l >= 0; l--) {
    {foreach_coarse_level (l) {
      {foreach_child()
        val(w,0,0,0) = val(s,0,0,0);end_foreach_child()}
      _attribute[s.i].prolongation (point, s);
      {foreach_child() {
        double sp = val(s,0,0,0);
        val(s,0,0,0) = val(w,0,0,0);

        val(w,0,0,0) -= sp;
      }end_foreach_child()}
    }end_foreach_coarse_level();}
    boundary_level (((scalar[]){w,{-1}}), l + 1);
  }

  {foreach_level(0)
    val(w,0,0,0) = val(s,0,0,0);end_foreach_level();}
  boundary_level (((scalar[]){w,{-1}}), 0);
}

void inverse_wavelet (scalar s, scalar w)
{
  {foreach_level(0)
    val(s,0,0,0) = val(w,0,0,0);end_foreach_level();}
  boundary_level (((scalar[]){s,{-1}}), 0);
  for (int l = 0; l <= grid->maxdepth - 1; l++) {
    {foreach_coarse_level (l) {
      _attribute[s.i].prolongation (point, s);
      {foreach_child()
        val(s,0,0,0) += val(w,0,0,0);end_foreach_child()}
    }end_foreach_coarse_level();}
    boundary_level (((scalar[]){s,{-1}}), l + 1);
  }
}

static inline double bilinear (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;



    return (9.*coarse(s,0,0,0) +
     3.*(coarse(s,child.x,0,0) + coarse(s,0,child.y,0)) +
     coarse(s,child.x,child.y,0))/16.;
# 140 "/home/spencer/basilisk/src/grid/multigrid-common.h"
}

static inline void refine_bilinear (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  {foreach_child()
    val(s,0,0,0) = bilinear (point, s);end_foreach_child()}
}

static inline double quadratic (double a, double b, double c)
{
  return (30.*a + 5.*b - 3.*c)/32.;
}

static inline double biquadratic (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;



  return
    quadratic (quadratic (coarse(s,0,0,0),
     coarse(s,child.x,0,0),
     coarse(s,-child.x,0,0)),
        quadratic (coarse(s,0,child.y,0),
     coarse(s,child.x,child.y,0),
     coarse(s,-child.x,child.y,0)),
        quadratic (coarse(s,0,-child.y,0),
     coarse(s,child.x,-child.y,0),
     coarse(s,-child.x,-child.y,0)));




}

static inline double biquadratic_vertex (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;



  return (36.*val(s,0,0,0) + 18.*(val(s,-1,0,0) + val(s,0,-1,0)) - 6.*(val(s,1,0,0) + val(s,0,1,0)) +
   9.*val(s,-1,-1,0) - 3.*(val(s,1,-1,0) + val(s,-1,1,0)) + val(s,1,1,0))/64.;




}

static inline void refine_biquadratic (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  {foreach_child()
    val(s,0,0,0) = biquadratic (point, s);end_foreach_child()}
}

static inline void refine_linear (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

#line 194
if(!is_constant(cm)){{
  coord g;
  if (_attribute[s.i].gradient)
    {
      g.x = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0));
      
#line 198
g.y = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0));}
  else
    {
      g.x = (val(s,1,0,0) - val(s,-1,0,0))/2.;
      
#line 201
g.y = (val(s,0,1,0) - val(s,0,-1,0))/2.;}

  double sc = val(s,0,0,0), cmc = 4.*val(cm,0,0,0), sum = val(cm,0,0,0)*(1 << 2);
  {foreach_child() {
    val(s,0,0,0) = sc;
    
      val(s,0,0,0) += child.x*g.x*val(cm,-child.x,0,0)/cmc;
      
#line 207
val(s,0,0,0) += child.y*g.y*val(cm,0,-child.y,0)/cmc;
    sum -= val(cm,0,0,0);
  }end_foreach_child()}
  if (!(fabs(sum) < 1e-10)) qassert ("/home/spencer/basilisk/src/grid/multigrid-common.h", 0, "fabs(sum) < 1e-10");
}}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);

#line 194
{
  coord g;
  if (_attribute[s.i].gradient)
    {
      g.x = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0));
      
#line 198
g.y = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0));}
  else
    {
      g.x = (val(s,1,0,0) - val(s,-1,0,0))/2.;
      
#line 201
g.y = (val(s,0,1,0) - val(s,0,-1,0))/2.;}

  double sc = val(s,0,0,0), cmc = 4.*_const_cm, sum = _const_cm*(1 << 2);
  {foreach_child() {
    val(s,0,0,0) = sc;
    
      val(s,0,0,0) += child.x*g.x*_const_cm/cmc;
      
#line 207
val(s,0,0,0) += child.y*g.y*_const_cm/cmc;
    sum -= _const_cm;
  }end_foreach_child()}
  if (!(fabs(sum) < 1e-10)) qassert ("/home/spencer/basilisk/src/grid/multigrid-common.h", 0, "fabs(sum) < 1e-10");
}}

#line 211
}

static inline void refine_reset (Point point, scalar v)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  {foreach_child()
    val(v,0,0,0) = 0.;end_foreach_child()}
}

static inline void refine_injection (Point point, scalar v)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  double val = val(v,0,0,0);
  {foreach_child()
    val(v,0,0,0) = val;end_foreach_child()}
}

static scalar multigrid_init_scalar (scalar s, const char * name)
{
  s = cartesian_init_scalar (s, name);
  _attribute[s.i].prolongation = refine_bilinear;
  _attribute[s.i].restriction = restriction_average;
  return s;
}

static scalar multigrid_init_vertex_scalar (scalar s, const char * name)
{
  s = cartesian_init_vertex_scalar (s, name);
  _attribute[s.i].restriction = restriction_vertex;
  return s;
}

static void multigrid_setup_vector (vector v)
{
   {
    _attribute[v.x.i].prolongation = refine_bilinear;
    _attribute[v.x.i].restriction = restriction_average;
  } 
#line 243
{
    _attribute[v.y.i].prolongation = refine_bilinear;
    _attribute[v.y.i].restriction = restriction_average;
  }
}

static vector multigrid_init_vector (vector v, const char * name)
{
  v = cartesian_init_vector (v, name);
  multigrid_setup_vector (v);
  return v;
}

static vector multigrid_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_face_vector (v, name);
  
    _attribute[v.y.i].restriction = no_restriction;
    
#line 260
_attribute[v.x.i].restriction = no_restriction;
  _attribute[v.x.i].restriction = restriction_face;
  return v;
}

static tensor multigrid_init_tensor (tensor t, const char * name)
{
  t = cartesian_init_tensor (t, name);
  
    multigrid_setup_vector (t.x);
    
#line 269
multigrid_setup_vector (t.y);
  return t;
}

void multigrid_debug (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  cartesian_debug (point);

  FILE * plot = fopen ("plot", "a");
  if (point.level > 0) {
    char name[80] = "coarse";
    if (pid() > 0)
      sprintf (name, "coarse-%d", pid());
    FILE * fp = fopen (name, "w");
# 294 "/home/spencer/basilisk/src/grid/multigrid-common.h"
      double xc = x - child.x*Delta/2., yc = y - child.y*Delta/2.;
      for (int k = 0; k <= 1; k++)
 for (int l = 0; l <= 1; l++) {
   {scalar*_i=(scalar*)( all);if(_i)for(scalar v=*_i;(&v)->i>=0;v=*++_i){
     fprintf (fp, "%g %g %g ",
       xc + k*child.x*Delta*2. + _attribute[v.i].d.x*Delta,
       yc + l*child.y*Delta*2. + _attribute[v.i].d.y*Delta,
       coarse(v,k*child.x,l*child.y,0));}}
   fputc ('\n', fp);
 }
      fprintf (ferr, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 3 t ''", name);
      fprintf (plot, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 3 t ''", name);
# 325 "/home/spencer/basilisk/src/grid/multigrid-common.h"
    fclose (fp);
  }

  if (is_coarse()) {
    char name[80] = "fine";
    if (pid() > 0)
      sprintf (name, "fine-%d", pid());
    FILE * fp = fopen (name, "w");
# 347 "/home/spencer/basilisk/src/grid/multigrid-common.h"
      double xf = x - Delta/4., yf = y - Delta/4.;
      for (int k = -2; k <= 3; k++)
 for (int l = -2; l <= 3; l++) {
   {scalar*_i=(scalar*)( all);if(_i)for(scalar v=*_i;(&v)->i>=0;v=*++_i){ {
     fprintf (fp, "%g %g ",
       xf + k*Delta/2. + _attribute[v.i].d.x*Delta/4.,
       yf + l*Delta/2. + _attribute[v.i].d.y*Delta/4.);
     if (allocated_child(k,l,0))
       fprintf (fp, "%g ", fine(v,k,l,0));
     else
       fputs ("n/a ", fp);
   }}}
   fputc ('\n', fp);
 }
      fprintf (ferr, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 2 t ''", name);
      fprintf (plot, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 2 t ''", name);
# 385 "/home/spencer/basilisk/src/grid/multigrid-common.h"
    fclose (fp);
  }
  fflush (ferr);
  fclose (plot);
}

static void multigrid_restriction (scalar * list)
{
  scalar * listdef = NULL, * listc = NULL, * list2 = NULL;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant (s) && _attribute[s.i].block > 0) {
      if (_attribute[s.i].restriction == restriction_average) {
 listdef = list_add (listdef, s);
 list2 = list_add (list2, s);
      }
      else if (_attribute[s.i].restriction != no_restriction) {
 listc = list_add (listc, s);
 if (_attribute[s.i].face)
   {
     list2 = list_add (list2, _attribute[s.i].v.x);
     
#line 404
list2 = list_add (list2, _attribute[s.i].v.y);}
 else
   list2 = list_add (list2, s);
      }
    }}}

  if (listdef || listc) {
    for (int l = depth() - 1; l >= 0; l--) {
      {foreach_coarse_level(l) {
 {scalar*_i=(scalar*)( listdef);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
  
     restriction_average (point, s);}}
 {scalar*_i=(scalar*)( listc);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
  
     _attribute[s.i].restriction (point, s);
 }}}
      }end_foreach_coarse_level();}
      { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list2, l); };
    }
    pfree (listdef,__func__,__FILE__,0);
    pfree (listc,__func__,__FILE__,0);
    pfree (list2,__func__,__FILE__,0);
  }
}

void multigrid_methods()
{
  cartesian_methods();
  init_scalar = multigrid_init_scalar;
  init_vertex_scalar = multigrid_init_vertex_scalar;
  init_vector = multigrid_init_vector;
  init_face_vector = multigrid_init_face_vector;
  init_tensor = multigrid_init_tensor;
  restriction = multigrid_restriction;
  debug = multigrid_debug;
}







void subtree_size (scalar size, bool leaves)
{




  foreach_stencil(1,{(NonLocal[]){{"size","scalar",(void *)&size,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n    \n// #line 454 \"/home/spencer/basilisk/src/grid/multigrid-common.h\"\nval_out_(size,0,0,0) = 1;")})
    {_stencil_val_a(size,0,0,0);  }end_foreach_stencil();




  {
#line 453
foreach()
    val(size,0,0,0) = 1;end_foreach();}





  { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b,((scalar[]) {size,{-1}}), depth()); };
  for (int l = depth() - 1; l >= 0; l--) {
    {foreach_coarse_level(l) {
      double sum = !leaves;
      {foreach_child()
 sum += val(size,0,0,0);end_foreach_child()}
      val(size,0,0,0) = sum;
    }end_foreach_coarse_level();}
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b,((scalar[]) {size,{-1}}), l); };
  }
}
# 5 "/home/spencer/basilisk/src/grid/tree-common.h" 2




# 21 "/home/spencer/basilisk/src/grid/tree-common.h"
int refine_cell (Point point, scalar * list, int flag, Cache * refined)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  int nr = 0;


  if (level > 0)
    for (int k = 0; k != 2*child.x; k += child.x)

      for (int l = 0; l != 2*child.y; l += child.y)




   if (aparent(k,l,m).pid >= 0 && is_leaf(aparent(k,l,m))) {
     Point p = point;


     p.level = point.level - 1;
     p.i = (point.i + 2)/2 + k;
     do { if (p.i < 2) p.i += 1 << p.level; else if (p.i >= 2 + (1 << p.level)) p.i -= 1 << p.level; } while(0);

       p.j = (point.j + 2)/2 + l;
       do { if (p.j < 2) p.j += 1 << p.level; else if (p.j >= 2 + (1 << p.level)) p.j -= 1 << p.level; } while(0);





     nr += refine_cell (p, list, flag, refined);
     aparent(k,l,m).flags |= flag;
   }



  increment_neighbors (point);

  int cflag = is_active(cell) ? (active|leaf) : leaf;
  {foreach_child()
    cell.flags |= cflag;end_foreach_child()}


  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (is_local(cell) || _attribute[s.i].face)
      _attribute[s.i].refine (point, s);}}


  cell.flags &= ~leaf;

@if _MPI
  if (is_border(cell)) {
    {foreach_child() {
      bool bord = false;
      {foreach_neighbor() {
 if (!is_local(cell) || (level > 0 && !is_local(aparent(0,0,0)))) {
   bord = true; foreach_neighbor_break();
 }
 if (is_refined_check())
   {foreach_child()
     if (!is_local(cell)) {
       bord = true; foreach_child_break();
     }end_foreach_child()}
 if (bord)
   foreach_neighbor_break();
      }end_foreach_neighbor()}
      if (bord)
 cell.flags |= border;
    }end_foreach_child()}
    if (refined)
      cache_append (refined, point, cell.flags);
    nr++;
  }
@endif
  return nr;
}





bool coarsen_cell (Point point, scalar * list)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;



  int pid = cell.pid;
  {foreach_child()
    if (cell.neighbors || (cell.pid < 0 && cell.pid != pid))
      return false;end_foreach_child()}



  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    _attribute[s.i].restriction (point, s);
    if (_attribute[s.i].coarsen)
      _attribute[s.i].coarsen (point, s);
  }}}


  cell.flags |= leaf;


  decrement_neighbors (point);

@if _MPI
  if (!is_local(cell)) {
    cell.flags &= ~(active|border);
    {foreach_neighbor(1)
      if (cell.neighbors)
 {foreach_child()
   if (allocated(0,0,0) && is_local(cell) && !is_border(cell))
     cell.flags |= border;end_foreach_child()}end_foreach_neighbor()}
  }
@endif

  return true;
}

void coarsen_cell_recursive (Point point, scalar * list)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;


  {foreach_child()
    if (cell.neighbors)
      {foreach_neighbor(1)
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
   coarsen_cell_recursive (point, list);end_foreach_neighbor()}end_foreach_child()}

  if (!(coarsen_cell (point, list))) qassert ("/home/spencer/basilisk/src/grid/tree-common.h", 0, "coarsen_cell (point, list)");
}

void mpi_boundary_refine (scalar *);
void mpi_boundary_coarsen (int, int);
void mpi_boundary_update (scalar *);

static
scalar * list_add_depend (scalar * list, scalar s)
{
  if (is_constant(s) || _attribute[s.i].restriction == no_restriction)
    return list;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar t=*_i;(&t)->i>=0;t=*++_i){
    if (t.i == s.i)
      return list;}}
  {scalar*_i=(scalar*)( _attribute[s.i].depends);if(_i)for(scalar d=*_i;(&d)->i>=0;d=*++_i){
    list = list_add_depend (list, d);}}
  return list_append (list, s);
}

typedef struct {
  int nc, nf;
} astats;

     
astats adapt_wavelet (scalar * slist,
        double * max,
        int maxlevel,
        int minlevel,
        scalar * list)
{tracing("adapt_wavelet","/home/spencer/basilisk/src/grid/tree-common.h",0);
  scalar * ilist = list;

  if (is_constant(cm)) {
    if (list == NULL || list == all)
      list = list_copy (all);
    boundary_internal ((scalar *)list, "/home/spencer/basilisk/src/grid/tree-common.h", 0);
    restriction (slist);
  }
  else {
    if (list == NULL || list == all) {
      list = list_copy (((scalar[]){cm, fm.x, fm.y,{-1}}));
      {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 list = list_add (list, s);}}
    }
    boundary_internal ((scalar *)list, "/home/spencer/basilisk/src/grid/tree-common.h", 0);
    scalar * listr = list_concat (slist,((scalar[]) {cm,{-1}}));
    restriction (listr);
    pfree (listr,__func__,__FILE__,0);
  }

  astats st = {0, 0};
  scalar * listc = NULL;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    listc = list_add_depend (listc, s);}}


  if (minlevel < 1)
    minlevel = 1;
  ((Tree *)grid)->refined.n = 0;
  static const int refined = 1 << user, too_fine = 1 << (user + 1);
  {foreach_cell() {
    if (is_active(cell)) {
      static const int too_coarse = 1 << (user + 2);
      if (is_leaf (cell)) {
 if (cell.flags & too_coarse) {
   cell.flags &= ~too_coarse;
   refine_cell (point, listc, refined, &((Tree *)grid)->refined);
   st.nf++;
 }
 continue;
      }
      else {
 if (cell.flags & refined) {

   cell.flags &= ~too_coarse;
   continue;
 }

 bool local = is_local(cell);
 if (!local)
   {foreach_child()
     if (is_local(cell)) {
       local = true; foreach_child_break();
     }end_foreach_child()}
 if (local) {
   int i = 0;
   static const int just_fine = 1 << (user + 3);
   {scalar*_i=(scalar*)( slist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
     double emax = max[i++], sc[1 << 2];
     int c = 0;
     {foreach_child()
       sc[c++] = val(s,0,0,0);end_foreach_child()}
     _attribute[s.i].prolongation (point, s);
     c = 0;
     {foreach_child() {
       double e = fabs(sc[c] - val(s,0,0,0));
       if (e > emax && level < maxlevel) {
  cell.flags &= ~too_fine;
  cell.flags |= too_coarse;
       }
       else if ((e <= emax/1.5 || level > maxlevel) &&
         !(cell.flags & (too_coarse|just_fine))) {
  if (level >= minlevel)
    cell.flags |= too_fine;
       }
       else if (!(cell.flags & too_coarse)) {
  cell.flags &= ~too_fine;
  cell.flags |= just_fine;
       }
       val(s,0,0,0) = sc[c++];
     }end_foreach_child()}
   }}}
   {foreach_child() {
     cell.flags &= ~just_fine;
     if (!is_leaf(cell)) {
       cell.flags &= ~too_coarse;
       if (level >= maxlevel)
  cell.flags |= too_fine;
     }
     else if (!is_active(cell))
       cell.flags &= ~too_coarse;
   }end_foreach_child()}
 }
      }
    }
    else
      continue;
  }end_foreach_cell();}
  mpi_boundary_refine (listc);



  for (int l = depth(); l >= 0; l--) {
    {foreach_cell()
      if (!(cell.pid < 0)) {
 if (level == l) {
   if (!is_leaf(cell)) {
     if (cell.flags & refined)

       cell.flags &= ~(refined|too_fine);
     else if (cell.flags & too_fine) {
       if (is_local(cell) && coarsen_cell (point, listc))
  st.nc++;
       cell.flags &= ~too_fine;
     }
   }
   if (cell.flags & too_fine)
     cell.flags &= ~too_fine;
   else if (level > 0 && (aparent(0,0,0).flags & too_fine))
     aparent(0,0,0).flags &= ~too_fine;
   continue;
 }
 else if (is_leaf(cell))
   continue;
      }end_foreach_cell();}
    mpi_boundary_coarsen (l, too_fine);
  }
  pfree (listc,__func__,__FILE__,0);

  mpi_all_reduce (st.nf, MPI_INT, MPI_SUM);
  mpi_all_reduce (st.nc, MPI_INT, MPI_SUM);
  if (st.nc || st.nf)
    mpi_boundary_update (list);

  if (list != ilist)
    pfree (list,__func__,__FILE__,0);

  {end_tracing("adapt_wavelet","/home/spencer/basilisk/src/grid/tree-common.h",0);return st;}
end_tracing("adapt_wavelet","/home/spencer/basilisk/src/grid/tree-common.h",0);}
# 339 "/home/spencer/basilisk/src/grid/tree-common.h"
static void refine_level (int depth)
{
  int refined;
  do {
    refined = 0;
    ((Tree *)grid)->refined.n = 0;
    {foreach_leaf()
      if (level < depth) {
 refine_cell (point, NULL, 0, &((Tree *)grid)->refined);
 refined++;
 continue;
      }end_foreach_leaf();}
    mpi_all_reduce (refined, MPI_INT, MPI_SUM);
    if (refined) {
      mpi_boundary_refine (NULL);
      mpi_boundary_update (NULL);
    }
  } while (refined);
}
# 384 "/home/spencer/basilisk/src/grid/tree-common.h"
     
static void halo_face (vectorl vl)
{tracing("halo_face","/home/spencer/basilisk/src/grid/tree-common.h",0);
  
    {scalar*_i=(scalar*)( vl.x);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = 2;}}
    
#line 388
{scalar*_i=(scalar*)( vl.y);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = 2;}}

  for (int l = depth() - 1; l >= 0; l--)
    {foreach_halo (prolongation, l)
      {
        if (vl.x) {
# 403 "/home/spencer/basilisk/src/grid/tree-common.h"
   if ((!is_leaf (neighbor(-1,0,0)) && neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0))
     {scalar*_i=(scalar*)( vl.x);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
       val(s,0,0,0) = (fine(s,0,0,0) + fine(s,0,1,0))/2.;}}
   if ((!is_leaf (neighbor(1,0,0)) && neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0))
     {scalar*_i=(scalar*)( vl.x);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
       val(s,1,0,0) = (fine(s,2,0,0) + fine(s,2,1,0))/2.;}}
# 419 "/home/spencer/basilisk/src/grid/tree-common.h"
 }
        
#line 394
if (vl.y) {
# 403 "/home/spencer/basilisk/src/grid/tree-common.h"
   if ((!is_leaf (neighbor(0,-1,0)) && neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0))
     {scalar*_i=(scalar*)( vl.y);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
       val(s,0,0,0) = (fine(s,0,0,0) + fine(s,1,0,0))/2.;}}
   if ((!is_leaf (neighbor(0,1,0)) && neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0))
     {scalar*_i=(scalar*)( vl.y);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
       val(s,0,1,0) = (fine(s,0,2,0) + fine(s,1,2,0))/2.;}}
# 419 "/home/spencer/basilisk/src/grid/tree-common.h"
 }}end_foreach_halo();}
end_tracing("halo_face","/home/spencer/basilisk/src/grid/tree-common.h",0);}



static scalar tree_init_scalar (scalar s, const char * name)
{
  s = multigrid_init_scalar (s, name);
  _attribute[s.i].refine = _attribute[s.i].prolongation;
  return s;
}

static void prolongation_vertex (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

  fine(s,1,1,0) = (val(s,0,0,0) + val(s,1,0,0) + val(s,0,1,0) + val(s,1,1,0))/4.;





  for (int i = 0; i <= 1; i++) {
    for (int j = 0; j <= 1; j++)





      if (allocated_child(2*i,2*j,0))
 fine(s,2*i,2*j,0) = val(s,i,j,0);


    
      if (neighbor(i,0,0).neighbors) {

 fine(s,2*i,1,0) = (val(s,i,0,0) + val(s,i,1,0))/2.;
# 464 "/home/spencer/basilisk/src/grid/tree-common.h"
      }
      
#line 452
if (neighbor(0,i,0).neighbors) {

 fine(s,1,2*i,0) = (val(s,0,i,0) + val(s,1,i,0))/2.;
# 464 "/home/spencer/basilisk/src/grid/tree-common.h"
      }
  }
}

static scalar tree_init_vertex_scalar (scalar s, const char * name)
{
  s = multigrid_init_vertex_scalar (s, name);
  _attribute[s.i].refine = _attribute[s.i].prolongation = prolongation_vertex;
  return s;
}

static void tree_setup_vector (vector v)
{
  
    _attribute[v.x.i].refine = _attribute[v.x.i].prolongation;
    
#line 478
_attribute[v.y.i].refine = _attribute[v.y.i].prolongation;
}

static vector tree_init_vector (vector v, const char * name)
{
  v = multigrid_init_vector (v, name);
  tree_setup_vector (v);
  return v;
}


static void refine_face_x (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  vector v = _attribute[s.i].v;

  if (!(!is_leaf (neighbor(-1,0,0)) && neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0) &&
      (is_local(cell) || is_local(neighbor(-1,0,0)))) {
    double g1 = (val(v.x,0,+1,0) - val(v.x,0,-1,0))/8.;
    for (int j = 0; j <= 1; j++)
      fine(v.x,0,j,0) = val(v.x,0,0,0) + (2*j - 1)*g1;
  }
  if (!(!is_leaf (neighbor(1,0,0)) && neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0) && neighbor(1,0,0).neighbors &&
      (is_local(cell) || is_local(neighbor(1,0,0)))) {
    double g1 = (val(v.x,1,+1,0) - val(v.x,1,-1,0))/8.;
    for (int j = 0; j <= 1; j++)
      fine(v.x,2,j,0) = val(v.x,1,0,0) + (2*j - 1)*g1;
  }
  if (is_local(cell)) {
    double g1 = (val(v.x,0,+1,0) - val(v.x,0,-1,0) + val(v.x,1,+1,0) - val(v.x,1,-1,0))/16.;
    for (int j = 0; j <= 1; j++)
      fine(v.x,1,j,0) = (val(v.x,0,0,0) + val(v.x,1,0,0))/2. + (2*j - 1)*g1;
  }
# 535 "/home/spencer/basilisk/src/grid/tree-common.h"
}

#line 489
static void refine_face_y (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  vector v = _attribute[s.i].v;

  if (!(!is_leaf (neighbor(0,-1,0)) && neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0) &&
      (is_local(cell) || is_local(neighbor(0,-1,0)))) {
    double g1 = (val(v.y,+1,0,0) - val(v.y,-1,0,0))/8.;
    for (int j = 0; j <= 1; j++)
      fine(v.y,j,0,0) = val(v.y,0,0,0) + (2*j - 1)*g1;
  }
  if (!(!is_leaf (neighbor(0,1,0)) && neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0) && neighbor(0,1,0).neighbors &&
      (is_local(cell) || is_local(neighbor(0,1,0)))) {
    double g1 = (val(v.y,+1,1,0) - val(v.y,-1,1,0))/8.;
    for (int j = 0; j <= 1; j++)
      fine(v.y,j,2,0) = val(v.y,0,1,0) + (2*j - 1)*g1;
  }
  if (is_local(cell)) {
    double g1 = (val(v.y,+1,0,0) - val(v.y,-1,0,0) + val(v.y,+1,1,0) - val(v.y,-1,1,0))/16.;
    for (int j = 0; j <= 1; j++)
      fine(v.y,j,1,0) = (val(v.y,0,0,0) + val(v.y,0,1,0))/2. + (2*j - 1)*g1;
  }
# 535 "/home/spencer/basilisk/src/grid/tree-common.h"
}

void refine_face (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  vector v = _attribute[s.i].v;
  
    _attribute[v.x.i].prolongation (point, v.x);
    
#line 541
_attribute[v.y.i].prolongation (point, v.y);
}

void refine_face_solenoidal (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  refine_face (point, s);

  if (is_local(cell)) {

    vector v = _attribute[s.i].v;
    double d[1 << 2], p[1 << 2];
    int i = 0;
    {foreach_child() {
      d[i] = 0.;
      
 d[i] += val(v.x,1,0,0) - val(v.x,0,0,0);
 
#line 556
d[i] += val(v.y,0,1,0) - val(v.y,0,0,0);
      i++;
    }end_foreach_child()}

    p[0] = 0.;
    p[1] = (3.*d[3] + d[0])/4. + d[2]/2.;
    p[2] = (d[3] + 3.*d[0])/4. + d[2]/2.;
    p[3] = (d[3] + d[0])/2. + d[2];
    fine(v.x,1,1,0) += p[1] - p[0];
    fine(v.x,1,0,0) += p[3] - p[2];
    fine(v.y,0,1,0) += p[0] - p[2];
    fine(v.y,1,1,0) += p[1] - p[3];
# 595 "/home/spencer/basilisk/src/grid/tree-common.h"
  }

}

static vector tree_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_face_vector (v, name);
  
    _attribute[v.x.i].restriction = _attribute[v.x.i].refine = no_restriction;
    
#line 603
_attribute[v.y.i].restriction = _attribute[v.y.i].refine = no_restriction;
  _attribute[v.x.i].restriction = restriction_face;
  _attribute[v.x.i].refine = refine_face;
  
    _attribute[v.x.i].prolongation = refine_face_x;
    
#line 607
_attribute[v.y.i].prolongation = refine_face_y;
  return v;
}

static tensor tree_init_tensor (tensor t, const char * name)
{
  t = multigrid_init_tensor (t, name);
  
    tree_setup_vector (t.x);
    
#line 615
tree_setup_vector (t.y);
  return t;
}

     
static void tree_boundary_level (scalar * list, int l)
{tracing("tree_boundary_level","/home/spencer/basilisk/src/grid/tree-common.h",0);
  int depth = l < 0 ? depth() : l;

  if (tree_is_full()) {
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list, depth); };
    {end_tracing("tree_boundary_level","/home/spencer/basilisk/src/grid/tree-common.h",0);return;}
  }

  scalar * listdef = NULL, * listc = NULL, * list2 = NULL, * vlist = NULL;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant (s)) {
      if (_attribute[s.i].restriction == restriction_average) {
 listdef = list_add (listdef, s);
 list2 = list_add (list2, s);
      }
      else if (_attribute[s.i].restriction != no_restriction) {
 listc = list_add (listc, s);
 if (_attribute[s.i].face)
   {
     list2 = list_add (list2, _attribute[s.i].v.x);
     
#line 640
list2 = list_add (list2, _attribute[s.i].v.y);}
 else {
   list2 = list_add (list2, s);
   if (_attribute[s.i].restriction == restriction_vertex)
     vlist = list_add (vlist, s);
 }
      }
    }}}

  if (vlist)






    {foreach_vertex () {
      if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) || (!is_leaf (neighbor(-1,0,0)) && neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0) ||
   (!is_leaf (neighbor(0,-1,0)) && neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0) || (!is_leaf (neighbor(-1,-1,0)) && neighbor(-1,-1,0).neighbors && neighbor(-1,-1,0).pid >= 0)) {

 {scalar*_i=(scalar*)( vlist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
   val(s,0,0,0) = is_vertex (child(0,0,0)) ? fine(s,0,0,0) : 1e30;}}
      }
      else
 {
   if (child.y == 1 &&
       ((!is_leaf(cell) && !cell.neighbors && cell.pid >= 0) || (!is_leaf(neighbor(-1,0,0)) && !neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0))) {

     {scalar*_i=(scalar*)( vlist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
       val(s,0,0,0) = is_vertex(neighbor(0,-1,0)) && is_vertex(neighbor(0,1,0)) ?
  (val(s,0,-1,0) + val(s,0,1,0))/2. : 1e30;}}
   }
   
#line 665
if (child.x == 1 &&
       ((!is_leaf(cell) && !cell.neighbors && cell.pid >= 0) || (!is_leaf(neighbor(0,-1,0)) && !neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0))) {

     {scalar*_i=(scalar*)( vlist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
       val(s,0,0,0) = is_vertex(neighbor(-1,0,0)) && is_vertex(neighbor(1,0,0)) ?
  (val(s,-1,0,0) + val(s,1,0,0))/2. : 1e30;}}
   }}
    }end_foreach_vertex();}
# 705 "/home/spencer/basilisk/src/grid/tree-common.h"
  pfree (vlist,__func__,__FILE__,0);

  if (listdef || listc) {
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, list2, depth); };
    for (int l = depth - 1; l >= 0; l--) {
      {foreach_coarse_level(l) {
 {scalar*_i=(scalar*)( listdef);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
   restriction_average (point, s);}}
 {scalar*_i=(scalar*)( listc);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
   _attribute[s.i].restriction (point, s);}}
      }end_foreach_coarse_level();}
      { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, list2, l); };
    }
    pfree (listdef,__func__,__FILE__,0);
    pfree (listc,__func__,__FILE__,0);
    pfree (list2,__func__,__FILE__,0);
  }

  scalar * listr = NULL;
  vector * listf = NULL;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant (s) && _attribute[s.i].refine != no_restriction) {
      if (_attribute[s.i].face)
 listf = vectors_add (listf, _attribute[s.i].v);
      else
 listr = list_add (listr, s);
    }}}

  if (listr || listf) {
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list, 0); };
    for (int i = 0; i < depth; i++) {
      {foreach_halo (prolongation, i) {
 {scalar*_i=(scalar*)( listr);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
          _attribute[s.i].prolongation (point, s);}}
 {vector*_i=(vector*)( listf);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
   {
     _attribute[v.x.i].prolongation (point, v.x);
     
#line 741
_attribute[v.y.i].prolongation (point, v.y);}}}
      }end_foreach_halo();}
      { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list, i + 1); };
    }
    pfree (listr,__func__,__FILE__,0);
    pfree (listf,__func__,__FILE__,0);
  }
end_tracing("tree_boundary_level","/home/spencer/basilisk/src/grid/tree-common.h",0);}

double treex (Point point) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (level == 0)
    return 0;

  double i = 2*child.x - child.y;
  if (i <= 1 && i >= -1) i = -i;




  return treex(parent) + i/(1 << 2*(level - 1));
}

double treey (Point point) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (level == 0)
    return 0;
  return treey(parent) + 4./(1 << 2*(level - 1));
}

void output_tree (FILE * fp)
{
  {foreach_cell()
    if (cell.neighbors)
      {foreach_child()
 if (is_local(cell))
   fprintf (fp, "%g %g\n%g %g\n\n",
     treex(parent), treey(parent), treex(point), treey(point));end_foreach_child()}end_foreach_cell();}
}

     
void tree_check()
{tracing("tree_check","/home/spencer/basilisk/src/grid/tree-common.h",0);


  long nleaves = 0, nactive = 0;
  {foreach_cell_all() {
    if (is_leaf(cell)) {
      if (!(cell.pid >= 0)) qassert ("/home/spencer/basilisk/src/grid/tree-common.h", 0, "cell.pid >= 0");
      nleaves++;
    }
    if (is_local(cell))
      if (!(is_active(cell) || (!is_leaf(cell) && !cell.neighbors && cell.pid >= 0))) qassert ("/home/spencer/basilisk/src/grid/tree-common.h", 0, "is_active(cell) || is_prolongation(cell)");
    if (is_active(cell))
      nactive++;

    int neighbors = 0;
    {foreach_neighbor(1)
      if (allocated(0,0,0) && (!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
 neighbors++;end_foreach_neighbor()}
    if (!(cell.neighbors == neighbors)) qassert ("/home/spencer/basilisk/src/grid/tree-common.h", 0, "cell.neighbors == neighbors");


    if (!cell.neighbors)
      if (!(!allocated_child(0,0,0))) qassert ("/home/spencer/basilisk/src/grid/tree-common.h", 0, "!allocated_child(0)");
  }end_foreach_cell_all();}


  long reachable = 0;
  {foreach_cell() {
    if (is_active(cell))
      reachable++;
    else
      continue;
  }end_foreach_cell();}
  if (!(nactive == reachable)) qassert ("/home/spencer/basilisk/src/grid/tree-common.h", 0, "nactive == reachable");


  reachable = 0;
  {foreach_cell()
    if (is_leaf(cell)) {
      reachable++;
      continue;
    }end_foreach_cell();}
  if (!(nleaves == reachable)) qassert ("/home/spencer/basilisk/src/grid/tree-common.h", 0, "nleaves == reachable");
end_tracing("tree_check","/home/spencer/basilisk/src/grid/tree-common.h",0);}

     
static void tree_restriction (scalar * list) {tracing("tree_restriction","/home/spencer/basilisk/src/grid/tree-common.h",0);
  boundary_internal ((scalar *)list, "/home/spencer/basilisk/src/grid/tree-common.h", 0);
  if (tree_is_full())
    multigrid_restriction (list);
end_tracing("tree_restriction","/home/spencer/basilisk/src/grid/tree-common.h",0);}

void tree_methods()
{
  multigrid_methods();
  init_scalar = tree_init_scalar;
  init_vertex_scalar = tree_init_vertex_scalar;
  init_vector = tree_init_vector;
  init_face_vector = tree_init_face_vector;
  init_tensor = tree_init_tensor;
  boundary_level = tree_boundary_level;
  boundary_face = halo_face;
  restriction = tree_restriction;
}
# 1670 "/home/spencer/basilisk/src/grid/tree.h" 2


void tree_periodic (int dir)
{
  int depth = grid ? depth() : -1;
  if (grid)
    free_grid();
  periodic (dir);
  if (depth >= 0)
    init_grid (1 << depth);
}


@if _MPI
# 1 "grid/tree-mpi.h" 1
# 1 "/home/spencer/basilisk/src/grid/tree-mpi.h"

int debug_iteration = -1;

void debug_mpi (FILE * fp1);

typedef struct {
  CacheLevel * halo;
  void * buf;
  MPI_Request r;
  int depth;
  int pid;
  int maxdepth;
} Rcv;

typedef struct {
  Rcv * rcv;
  char * name;
  int npid;
} RcvPid;

typedef struct {
  RcvPid * rcv, * snd;
} SndRcv;

typedef struct {
  Boundary parent;

  SndRcv mpi_level, mpi_level_root, restriction;
  Array * send, * receive;
} MpiBoundary;

static void cache_level_init (CacheLevel * c)
{
  c->p = NULL;
  c->n = c->nm = 0;
}

static void rcv_append (Point point, Rcv * rcv)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (level > rcv->depth) {
    rcv->halo = (CacheLevel *) prealloc (rcv->halo, (level + 1)*sizeof(CacheLevel),__func__,__FILE__,0);
    for (int j = rcv->depth + 1; j <= level; j++)
      cache_level_init (&rcv->halo[j]);
    rcv->depth = level;
  }
  cache_level_append (&rcv->halo[level], point);
  if (level > rcv->maxdepth)
    rcv->maxdepth = level;
}

void rcv_print (Rcv * rcv, FILE * fp, const char * prefix)
{
  for (int l = 0; l <= rcv->depth; l++)
    if (rcv->halo[l].n > 0)
      {foreach_cache_level(rcv->halo[l], l)
 fprintf (fp, "%s%g %g %g %d %d\n", prefix, x, y, z, rcv->pid, level);end_foreach_cache_level();}
}

static void rcv_free_buf (Rcv * rcv)
{
  if (rcv->buf) {
    prof_start ("rcv_pid_receive");
    MPI_Wait (&rcv->r, MPI_STATUS_IGNORE);
    pfree (rcv->buf,__func__,__FILE__,0);
    rcv->buf = NULL;
    prof_stop();
  }
}

static void rcv_destroy (Rcv * rcv)
{
  rcv_free_buf (rcv);
  for (int i = 0; i <= rcv->depth; i++)
    if (rcv->halo[i].n > 0)
      pfree (rcv->halo[i].p,__func__,__FILE__,0);
  pfree (rcv->halo,__func__,__FILE__,0);
}

static RcvPid * rcv_pid_new (const char * name)
{
  RcvPid * r = ((RcvPid *) pcalloc (1, sizeof(RcvPid),__func__,__FILE__,0));
  r->name = pstrdup (name,__func__,__FILE__,0);
  return r;
}

static Rcv * rcv_pid_pointer (RcvPid * p, int pid)
{
  if (!(pid >= 0 && pid < npe())) qassert ("/home/spencer/basilisk/src/grid/tree-mpi.h", 0, "pid >= 0 && pid < npe()");

  int i;
  for (i = 0; i < p->npid; i++)
    if (pid == p->rcv[i].pid)
      break;

  if (i == p->npid) {
    p->rcv = (Rcv *) prealloc (p->rcv, (++p->npid)*sizeof(Rcv),__func__,__FILE__,0);
    Rcv * rcv = &p->rcv[p->npid-1];
    rcv->pid = pid;
    rcv->depth = rcv->maxdepth = 0;
    rcv->halo = ((CacheLevel *) pmalloc ((1)*sizeof(CacheLevel),__func__,__FILE__,0));
    rcv->buf = NULL;
    cache_level_init (&rcv->halo[0]);
  }
  return &p->rcv[i];
}

static void rcv_pid_append (RcvPid * p, int pid, Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  rcv_append (point, rcv_pid_pointer (p, pid));
}

static void rcv_pid_append_pids (RcvPid * p, Array * pids)
{

  for (int i = 0; i < p->npid; i++) {
    int pid = p->rcv[i].pid, j, * a;
    for (j = 0, a = pids->p; j < pids->len/sizeof(int); j++,a++)
      if (*a == pid)
 break;
    if (j == pids->len/sizeof(int))
      array_append (pids, &pid, sizeof(int));
  }
}

void rcv_pid_write (RcvPid * p, const char * name)
{
  for (int i = 0; i < p->npid; i++) {
    Rcv * rcv = &p->rcv[i];
    char fname[80];
    sprintf (fname, "%s-%d-%d", name, pid(), rcv->pid);
    FILE * fp = fopen (fname, "w");
    rcv_print (rcv, fp, "");
    fclose (fp);
  }
}

static void rcv_pid_print (RcvPid * p, FILE * fp, const char * prefix)
{
  for (int i = 0; i < p->npid; i++)
    rcv_print (&p->rcv[i], fp, prefix);
}

static void rcv_pid_destroy (RcvPid * p)
{
  for (int i = 0; i < p->npid; i++)
    rcv_destroy (&p->rcv[i]);
  pfree (p->rcv,__func__,__FILE__,0);
  pfree (p->name,__func__,__FILE__,0);
  pfree (p,__func__,__FILE__,0);
}

static Boundary * mpi_boundary = NULL;






void debug_mpi (FILE * fp1);

static void apply_bc (Rcv * rcv, scalar * list, scalar * listv,
        vector * listf, int l, MPI_Status s)
{
  double * b = rcv->buf;
  {foreach_cache_level(rcv->halo[l], l) {
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
      memcpy (&val(s,0,0,0), b, sizeof(double)*_attribute[s.i].block);
      b += _attribute[s.i].block;
    }}}
    {vector*_i=(vector*)( listf);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
      { {
 memcpy (&val(v.x,0,0,0), b, sizeof(double)*_attribute[v.x.i].block);
 b += _attribute[v.x.i].block;
 if (*b != 1e30 && allocated(1,0,0))
   memcpy (&val(v.x,1,0,0), b, sizeof(double)*_attribute[v.x.i].block);
 b += _attribute[v.x.i].block;
      } 
#line 171
{
 memcpy (&val(v.y,0,0,0), b, sizeof(double)*_attribute[v.y.i].block);
 b += _attribute[v.y.i].block;
 if (*b != 1e30 && allocated(0,1,0))
   memcpy (&val(v.y,0,1,0), b, sizeof(double)*_attribute[v.y.i].block);
 b += _attribute[v.y.i].block;
      }}}}
    {scalar*_i=(scalar*)( listv);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
      for (int i = 0; i <= 1; i++)
 for (int j = 0; j <= 1; j++)







          {
     if (*b != 1e30 && allocated(i,j,0))
       memcpy (&val(s,i,j,0), b, sizeof(double)*_attribute[s.i].block);
     b += _attribute[s.i].block;
          }

    }}}
  }end_foreach_cache_level();}
  size_t size = b - (double *) rcv->buf;
  pfree (rcv->buf,__func__,__FILE__,0);
  rcv->buf = NULL;

  int rlen;
  MPI_Get_count (&s, MPI_DOUBLE, &rlen);
  if (rlen != size) {
    fprintf (ferr,
      "rlen (%d) != size (%ld), %d receiving from %d at level %d\n"
      "Calling debug_mpi(NULL)...\n"
      "Aborting...\n",
      rlen, size, pid(), rcv->pid, l);
    fflush (ferr);
    debug_mpi (NULL);
    MPI_Abort (MPI_COMM_WORLD, -2);
  }
}
# 234 "/home/spencer/basilisk/src/grid/tree-mpi.h"
static void mpi_recv_check (void * buf, int count, MPI_Datatype datatype,
       int source, int tag,
       MPI_Comm comm, MPI_Status * status,
       const char * name)
{
# 269 "/home/spencer/basilisk/src/grid/tree-mpi.h"
  int errorcode = MPI_Recv (buf, count, datatype, source, tag, comm, status);
  if (errorcode != MPI_SUCCESS) {
    char string[MPI_MAX_ERROR_STRING];
    int resultlen;
    MPI_Error_string (errorcode, string, &resultlen);
    fprintf (ferr,
      "ERROR MPI_Recv \"%s\" (count = %d, source = %d, tag = %d):\n%s\n"
      "Calling debug_mpi(NULL)...\n"
      "Aborting...\n",
      name, count, source, tag, string);
    fflush (ferr);
    debug_mpi (NULL);
    MPI_Abort (MPI_COMM_WORLD, -1);
  }





}

     
static int mpi_waitany (int count, MPI_Request array_of_requests[], int *indx,
   MPI_Status *status)
{tracing("mpi_waitany","/home/spencer/basilisk/src/grid/tree-mpi.h",0);
  { int _ret= MPI_Waitany (count, array_of_requests, indx, status);end_tracing("mpi_waitany","/home/spencer/basilisk/src/grid/tree-mpi.h",0);return _ret;}
end_tracing("mpi_waitany","/home/spencer/basilisk/src/grid/tree-mpi.h",0);}

static int list_lenb (scalar * list) {
  int len = 0;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    len += _attribute[s.i].block;}}
  return len;
}

static int vectors_lenb (vector * list) {
  int len = 0;
  {vector*_i=(vector*)( list);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
    len += _attribute[v.x.i].block;}}
  return len;
}

static void rcv_pid_receive (RcvPid * m, scalar * list, scalar * listv,
        vector * listf, int l)
{
  if (m->npid == 0)
    return;

  prof_start ("rcv_pid_receive");

  int len = list_lenb (list) + 2*2*vectors_lenb (listf) +
    (1 << 2)*list_lenb (listv);

  MPI_Request r[m->npid];
  Rcv * rrcv[m->npid];
  int nr = 0;
  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
    if (l <= rcv->depth && rcv->halo[l].n > 0) {
      if (!(!rcv->buf)) qassert ("/home/spencer/basilisk/src/grid/tree-mpi.h", 0, "!rcv->buf");
      rcv->buf = pmalloc (sizeof (double)*rcv->halo[l].n*len,__func__,__FILE__,0);






      MPI_Irecv (rcv->buf, rcv->halo[l].n*len, MPI_DOUBLE, rcv->pid,
   (l), MPI_COMM_WORLD, &r[nr]);
      rrcv[nr++] = rcv;






    }
  }


  if (nr > 0) {
    int i;
    MPI_Status s;
    mpi_waitany (nr, r, &i, &s);
    while (i != MPI_UNDEFINED) {
      Rcv * rcv = rrcv[i];
      if (!(l <= rcv->depth && rcv->halo[l].n > 0)) qassert ("/home/spencer/basilisk/src/grid/tree-mpi.h", 0, "l <= rcv->depth && rcv->halo[l].n > 0");
      if (!(rcv->buf)) qassert ("/home/spencer/basilisk/src/grid/tree-mpi.h", 0, "rcv->buf");
      apply_bc (rcv, list, listv, listf, l, s);
      mpi_waitany (nr, r, &i, &s);
    }
  }

  prof_stop();
}

     
static void rcv_pid_wait (RcvPid * m)
{tracing("rcv_pid_wait","/home/spencer/basilisk/src/grid/tree-mpi.h",0);

  for (int i = 0; i < m->npid; i++)
    rcv_free_buf (&m->rcv[i]);
end_tracing("rcv_pid_wait","/home/spencer/basilisk/src/grid/tree-mpi.h",0);}

static void rcv_pid_send (RcvPid * m, scalar * list, scalar * listv,
     vector * listf, int l)
{
  if (m->npid == 0)
    return;

  prof_start ("rcv_pid_send");

  int len = list_lenb (list) + 2*2*vectors_lenb (listf) +
    (1 << 2)*list_lenb (listv);


  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
    if (l <= rcv->depth && rcv->halo[l].n > 0) {
      if (!(!rcv->buf)) qassert ("/home/spencer/basilisk/src/grid/tree-mpi.h", 0, "!rcv->buf");
      rcv->buf = pmalloc (sizeof (double)*rcv->halo[l].n*len,__func__,__FILE__,0);
      double * b = rcv->buf;
      {foreach_cache_level(rcv->halo[l], l) {
 {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
   memcpy (b, &val(s,0,0,0), sizeof(double)*_attribute[s.i].block);
   b += _attribute[s.i].block;
 }}}
 {vector*_i=(vector*)( listf);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
   { {
     memcpy (b, &val(v.x,0,0,0), sizeof(double)*_attribute[v.x.i].block);
     b += _attribute[v.x.i].block;
     if (allocated(1,0,0))
       memcpy (b, &val(v.x,1,0,0), sizeof(double)*_attribute[v.x.i].block);
     else
       *b = 1e30;
     b += _attribute[v.x.i].block;
   } 
#line 397
{
     memcpy (b, &val(v.y,0,0,0), sizeof(double)*_attribute[v.y.i].block);
     b += _attribute[v.y.i].block;
     if (allocated(0,1,0))
       memcpy (b, &val(v.y,0,1,0), sizeof(double)*_attribute[v.y.i].block);
     else
       *b = 1e30;
     b += _attribute[v.y.i].block;
   }}}}
 {scalar*_i=(scalar*)( listv);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
   for (int i = 0; i <= 1; i++)
     for (int j = 0; j <= 1; j++)
# 418 "/home/spencer/basilisk/src/grid/tree-mpi.h"
       {
  if (allocated(i,j,0))
    memcpy (b, &val(s,i,j,0), sizeof(double)*_attribute[s.i].block);
  else
    *b = 1e30;
  b += _attribute[s.i].block;
       }

 }}}
      }end_foreach_cache_level();}





      MPI_Isend (rcv->buf, (b - (double *) rcv->buf),
   MPI_DOUBLE, rcv->pid, (l), MPI_COMM_WORLD,
   &rcv->r);
    }
  }

  prof_stop();
}

static void rcv_pid_sync (SndRcv * m, scalar * list, int l)
{
  scalar * listr = NULL, * listv = NULL;
  vector * listf = NULL;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant(s) && _attribute[s.i].block > 0) {
      if (_attribute[s.i].face)
 listf = vectors_add (listf, _attribute[s.i].v);
      else if (_attribute[s.i].restriction == restriction_vertex)
 listv = list_add (listv, s);
      else
 listr = list_add (listr, s);
    }}}
  rcv_pid_send (m->snd, listr, listv, listf, l);
  rcv_pid_receive (m->rcv, listr, listv, listf, l);
  rcv_pid_wait (m->snd);
  pfree (listr,__func__,__FILE__,0);
  pfree (listf,__func__,__FILE__,0);
  pfree (listv,__func__,__FILE__,0);
}

static void snd_rcv_destroy (SndRcv * m)
{
  rcv_pid_destroy (m->rcv);
  rcv_pid_destroy (m->snd);
}

static void snd_rcv_init (SndRcv * m, const char * name)
{
  char s[strlen(name) + 5];
  strcpy (s, name);
  strcat (s, ".rcv");
  m->rcv = rcv_pid_new (s);
  strcpy (s, name);
  strcat (s, ".snd");
  m->snd = rcv_pid_new (s);
}

static void mpi_boundary_destroy (Boundary * b)
{
  MpiBoundary * m = (MpiBoundary *) b;
  snd_rcv_destroy (&m->mpi_level);
  snd_rcv_destroy (&m->mpi_level_root);
  snd_rcv_destroy (&m->restriction);
  array_free (m->send);
  array_free (m->receive);
  pfree (m,__func__,__FILE__,0);
}

     
static void mpi_boundary_level (const Boundary * b, scalar * list, int l)
{tracing("mpi_boundary_level","/home/spencer/basilisk/src/grid/tree-mpi.h",0);
  MpiBoundary * m = (MpiBoundary *) b;
  rcv_pid_sync (&m->mpi_level, list, l);
  rcv_pid_sync (&m->mpi_level_root, list, l);
end_tracing("mpi_boundary_level","/home/spencer/basilisk/src/grid/tree-mpi.h",0);}

     
static void mpi_boundary_restriction (const Boundary * b, scalar * list, int l)
{tracing("mpi_boundary_restriction","/home/spencer/basilisk/src/grid/tree-mpi.h",0);
  MpiBoundary * m = (MpiBoundary *) b;
  rcv_pid_sync (&m->restriction, list, l);
end_tracing("mpi_boundary_restriction","/home/spencer/basilisk/src/grid/tree-mpi.h",0);}

void mpi_boundary_new()
{
  mpi_boundary = (Boundary *) ((MpiBoundary *) pcalloc (1, sizeof(MpiBoundary),__func__,__FILE__,0));
  mpi_boundary->destroy = mpi_boundary_destroy;
  mpi_boundary->level = mpi_boundary_level;
  mpi_boundary->restriction = mpi_boundary_restriction;
  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;
  snd_rcv_init (&mpi->mpi_level, "mpi_level");
  snd_rcv_init (&mpi->mpi_level_root, "mpi_level_root");
  snd_rcv_init (&mpi->restriction, "restriction");
  mpi->send = array_new();
  mpi->receive = array_new();
  add_boundary (mpi_boundary);
}

static FILE * fopen_prefix (FILE * fp, const char * name, char * prefix)
{
  if (fp) {
    sprintf (prefix, "%s-%d ", name, pid());
    return fp;
  }
  else {
    strcpy (prefix, "");
    char fname[80];
    if (debug_iteration >= 0)
      sprintf (fname, "%s-%d-%d", name, debug_iteration, pid());
    else
      sprintf (fname, "%s-%d", name, pid());
    return fopen (fname, "w");
  }
}

void debug_mpi (FILE * fp1)
{
  void output_cells_internal (FILE * fp);

  char prefix[80];
  FILE * fp;


  if (fp1 == NULL) {
    char name[80];
    sprintf (name, "halo-%d", pid()); remove (name);
    sprintf (name, "cells-%d", pid()); remove (name);
    sprintf (name, "faces-%d", pid()); remove (name);
    sprintf (name, "vertices-%d", pid()); remove (name);
    sprintf (name, "neighbors-%d", pid()); remove (name);
    sprintf (name, "mpi-level-rcv-%d", pid()); remove (name);
    sprintf (name, "mpi-level-snd-%d", pid()); remove (name);
    sprintf (name, "mpi-level-root-rcv-%d", pid()); remove (name);
    sprintf (name, "mpi-level-root-snd-%d", pid()); remove (name);
    sprintf (name, "mpi-restriction-rcv-%d", pid()); remove (name);
    sprintf (name, "mpi-restriction-snd-%d", pid()); remove (name);
    sprintf (name, "mpi-border-%d", pid()); remove (name);
    sprintf (name, "exterior-%d", pid()); remove (name);
    sprintf (name, "depth-%d", pid()); remove (name);
    sprintf (name, "refined-%d", pid()); remove (name);
  }


  fp = fopen_prefix (fp1, "halo", prefix);
  for (int l = 0; l < depth(); l++)
    {foreach_halo (prolongation, l)
      {foreach_child()
        fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level);end_foreach_child()}end_foreach_halo();}
  if (!fp1)
    fclose (fp);

  if (!fp1) {
    fp = fopen_prefix (fp1, "cells", prefix);
    output_cells_internal (fp);
    fclose (fp);
  }

  fp = fopen_prefix (fp1, "faces", prefix);
  {foreach_face_generic(){is_face_x(){
    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level);}end_is_face_x()
#line 581
is_face_y(){
    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level);}end_is_face_y()}end_foreach_face_generic();}
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "vertices", prefix);
  {foreach_vertex()
    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level);end_foreach_vertex();}
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "neighbors", prefix);
  {foreach() {
    int n = 0;
    {foreach_neighbor(1)
      if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
 n++;end_foreach_neighbor()}
    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, cell.neighbors);
    if (!(cell.neighbors == n)) qassert ("/home/spencer/basilisk/src/grid/tree-mpi.h", 0, "cell.neighbors == n");
  }end_foreach();}
  if (!fp1)
    fclose (fp);

  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;

  fp = fopen_prefix (fp1, "mpi-level-rcv", prefix);
  rcv_pid_print (mpi->mpi_level.rcv, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-level-root-rcv", prefix);
  rcv_pid_print (mpi->mpi_level_root.rcv, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-restriction-rcv", prefix);
  rcv_pid_print (mpi->restriction.rcv, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-level-snd", prefix);
  rcv_pid_print (mpi->mpi_level.snd, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-level-root-snd", prefix);
  rcv_pid_print (mpi->mpi_level_root.snd, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-restriction-snd", prefix);
  rcv_pid_print (mpi->restriction.snd, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-border", prefix);
  {foreach_cell() {
    if (is_border(cell))
      fprintf (fp, "%s%g %g %g %d %d %d\n",
        prefix, x, y, z, level, cell.neighbors, cell.pid);
    else
      continue;
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "exterior", prefix);
  {foreach_cell() {
    if (!is_local(cell))
      fprintf (fp, "%s%g %g %g %d %d %d %d\n",
        prefix, x, y, z, level, cell.neighbors,
        cell.pid, cell.flags & leaf);






  }end_foreach_cell();}
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "depth", prefix);
  fprintf (fp, "depth: %d %d\n", pid(), depth());
  fprintf (fp, "======= mpi_level.snd ======\n");
  RcvPid * snd = mpi->mpi_level.snd;
  for (int i = 0; i < snd->npid; i++)
    fprintf (fp, "%d %d %d\n", pid(), snd->rcv[i].pid, snd->rcv[i].maxdepth);
  fprintf (fp, "======= mpi_level.rcv ======\n");
  snd = mpi->mpi_level.rcv;
  for (int i = 0; i < snd->npid; i++)
    fprintf (fp, "%d %d %d\n", pid(), snd->rcv[i].pid, snd->rcv[i].maxdepth);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "refined", prefix);
  {foreach_cache (((Tree *)grid)->refined)
    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level);end_foreach_cache();}
  if (!fp1)
    fclose (fp);
}

static void snd_rcv_free (SndRcv * p)
{
  char name[strlen(p->rcv->name) + 1];
  strcpy (name, p->rcv->name);
  rcv_pid_destroy (p->rcv);
  p->rcv = rcv_pid_new (name);
  strcpy (name, p->snd->name);
  rcv_pid_destroy (p->snd);
  p->snd = rcv_pid_new (name);
}

static bool is_root (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
    {foreach_child()
      if (is_local(cell))
 return true;end_foreach_child()}
  return false;
}


static bool is_local_prolongation (Point point, Point p)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

  struct { int x, y; } dp = {p.i - point.i, p.j - point.j};



   {
    if (dp.x == 0 && ((!is_leaf (neighbor(-1,0,0)) && neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0) || (!is_leaf (neighbor(1,0,0)) && neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0)))
      return true;
    if ((!is_leaf (neighbor(dp.x,0,0)) && neighbor(dp.x,0,0).neighbors && neighbor(dp.x,0,0).pid >= 0))
      return true;
  } 
#line 713
{
    if (dp.y == 0 && ((!is_leaf (neighbor(0,-1,0)) && neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0) || (!is_leaf (neighbor(0,1,0)) && neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0)))
      return true;
    if ((!is_leaf (neighbor(0,dp.y,0)) && neighbor(0,dp.y,0).neighbors && neighbor(0,dp.y,0).pid >= 0))
      return true;
  }
  return false;
}



static void append_pid (Array * pids, int pid)
{
  for (int i = 0, * p = (int *) pids->p; i < pids->len/sizeof(int); i++, p++)
    if (*p == pid)
      return;
  array_append (pids, &pid, sizeof(int));
}

static int locals_pids (Point point, Array * pids)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (is_leaf(cell)) {
    if (is_local(cell)) {
      Point p = point;
      {foreach_neighbor(1) {
 if ((cell.pid >= 0 && cell.pid != pid()) &&
     ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) || is_local_prolongation (point, p)))
   append_pid (pids, cell.pid);
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
   {foreach_child()
     if ((cell.pid >= 0 && cell.pid != pid()))
       append_pid (pids, cell.pid);end_foreach_child()}
      }end_foreach_neighbor()}
    }
  }
  else
    {foreach_neighbor(1) {
      if ((cell.pid >= 0 && cell.pid != pid()))
 append_pid (pids, cell.pid);
      if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
 {foreach_child()
   if ((cell.pid >= 0 && cell.pid != pid()))
     append_pid (pids, cell.pid);end_foreach_child()}
    }end_foreach_neighbor()}
  return pids->len/sizeof(int);
}

static int root_pids (Point point, Array * pids)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  {foreach_child()
    if ((cell.pid >= 0 && cell.pid != pid()))
      append_pid (pids, cell.pid);end_foreach_child()}
  return pids->len/sizeof(int);
}







static void rcv_pid_row (RcvPid * m, int l, int * row)
{
  for (int i = 0; i < npe(); i++)
    row[i] = 0;
  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
    if (l <= rcv->depth && rcv->halo[l].n > 0)
      row[rcv->pid] = rcv->halo[l].n;
  }
}

void check_snd_rcv_matrix (SndRcv * sndrcv, const char * name)
{
  int maxlevel = depth();
  mpi_all_reduce (maxlevel, MPI_INT, MPI_MAX);
  int * row = ((int *) pmalloc ((npe())*sizeof(int),__func__,__FILE__,0));
  for (int l = 0; l <= maxlevel; l++) {
    int status = 0;
    if (pid() == 0) {


      int ** send = matrix_new (npe(), npe(), sizeof(int));
      int ** receive = matrix_new (npe(), npe(), sizeof(int));
      rcv_pid_row (sndrcv->snd, l, row);
      MPI_Gather (row, npe(), MPI_INT, &send[0][0], npe(), MPI_INT, 0,
    MPI_COMM_WORLD);
      rcv_pid_row (sndrcv->rcv, l, row);
      MPI_Gather (row, npe(), MPI_INT, &receive[0][0], npe(), MPI_INT, 0,
    MPI_COMM_WORLD);

      int * astatus = ((int *) pmalloc ((npe())*sizeof(int),__func__,__FILE__,0));
      for (int i = 0; i < npe(); i++)
 astatus[i] = 0;
      for (int i = 0; i < npe(); i++)
 for (int j = 0; j < npe(); j++)
   if (send[i][j] != receive[j][i]) {
     fprintf (ferr, "%s: %d sends    %d to   %d at level %d\n",
       name, i, send[i][j], j, l);
     fprintf (ferr, "%s: %d receives %d from %d at level %d\n",
       name, j, receive[j][i], i, l);
     fflush (ferr);
     for (int k = i - 2; k <= i + 2; k++)
       if (k >= 0 && k < npe())
  astatus[k] = 1;
     for (int k = j - 2; k <= j + 2; k++)
       if (k >= 0 && k < npe())
  astatus[k] = 1;
   }
      MPI_Scatter (astatus, 1, MPI_INT, &status, 1, MPI_INT, 0, MPI_COMM_WORLD);
      pfree (astatus,__func__,__FILE__,0);

      matrix_free (send);
      matrix_free (receive);
    }
    else {
      rcv_pid_row (sndrcv->snd, l, row);
      MPI_Gather (row, npe(), MPI_INT, NULL, npe(), MPI_INT, 0, MPI_COMM_WORLD);
      rcv_pid_row (sndrcv->rcv, l, row);
      MPI_Gather (row, npe(), MPI_INT, NULL, npe(), MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Scatter (NULL, 1, MPI_INT, &status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    if (status) {
      fprintf (ferr,
        "check_snd_rcv_matrix \"%s\" failed\n"
        "Calling debug_mpi(NULL)...\n"
        "Aborting...\n",
        name);
      fflush (ferr);
      debug_mpi (NULL);
      MPI_Abort (MPI_COMM_WORLD, -3);
    }
  }
  pfree (row,__func__,__FILE__,0);
}

static bool has_local_child (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  {foreach_child()
    if (is_local(cell))
      return true;end_foreach_child()}
  return false;
}

     
void mpi_boundary_update_buffers()
{tracing("mpi_boundary_update_buffers","/home/spencer/basilisk/src/grid/tree-mpi.h",0);
  if (npe() == 1)
    {end_tracing("mpi_boundary_update_buffers","/home/spencer/basilisk/src/grid/tree-mpi.h",0);return;}

  prof_start ("mpi_boundary_update_buffers");

  MpiBoundary * m = (MpiBoundary *) mpi_boundary;
  SndRcv * mpi_level = &m->mpi_level;
  SndRcv * mpi_level_root = &m->mpi_level_root;
  SndRcv * restriction = &m->restriction;

  snd_rcv_free (mpi_level);
  snd_rcv_free (mpi_level_root);
  snd_rcv_free (restriction);

  static const unsigned short used = 1 << user;
  {foreach_cell() {
    if (is_active(cell) && !is_border(cell))



      continue;

    if (cell.neighbors) {

      Array pids = {NULL, 0, 0};
      int n = locals_pids (point, &pids);
      if (n) {
 {foreach_child()
   if (is_local(cell))
     for (int i = 0, * p = (int *) pids.p; i < n; i++, p++)
       rcv_pid_append (mpi_level->snd, *p, point);end_foreach_child()}
 pfree (pids.p,__func__,__FILE__,0);
      }

      bool locals = false;
      if (is_leaf(cell)) {
 if ((cell.pid >= 0 && cell.pid != pid())) {
   Point p = point;
   {foreach_neighbor(1)
     if ((is_local(cell) &&
   ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) || is_local_prolongation (point, p))) ||
  is_root(point)) {
       locals = true; foreach_neighbor_break();
     }end_foreach_neighbor()}
 }
      }
      else
 {foreach_neighbor(1)
   if (is_local(cell) || is_root(point)) {
     locals = true; foreach_neighbor_break();
   }end_foreach_neighbor()}
      if (locals)
 {foreach_child()
   if ((cell.pid >= 0 && cell.pid != pid()))
            rcv_pid_append (mpi_level->rcv, cell.pid, point),
       cell.flags |= used;end_foreach_child()}


      if (!is_leaf(cell)) {

 if (is_local(cell)) {
   Array pids = {NULL, 0, 0};

   int n = root_pids (point, &pids);
   if (n) {
     {foreach_neighbor()
       for (int i = 0, * p = (int *) pids.p; i < n; i++, p++)
  if (cell.pid >= 0 && cell.pid != *p)
    rcv_pid_append (mpi_level_root->snd, *p, point);end_foreach_neighbor()}

     for (int i = 0, * p = (int *) pids.p; i < n; i++, p++)
       rcv_pid_append (restriction->snd, *p, point);
     pfree (pids.p,__func__,__FILE__,0);
   }
 }

 else if ((cell.pid >= 0 && cell.pid != pid())) {
   bool root = false;
   {foreach_child()
     if (is_local(cell)) {
       root = true; foreach_child_break();
     }end_foreach_child()}
   if (root) {
     int pid = cell.pid;
     {foreach_neighbor()
       if ((cell.pid >= 0 && cell.pid != pid()))
  rcv_pid_append (mpi_level_root->rcv, pid, point),
    cell.flags |= used;end_foreach_neighbor()}

     rcv_pid_append (restriction->rcv, pid, point);
   }
 }
      }
    }


    if (level > 0) {
      if (is_local(cell)) {

 Array pids = {NULL, 0, 0};
 if ((aparent(0,0,0).pid >= 0 && aparent(0,0,0).pid != pid()))
   append_pid (&pids, aparent(0,0,0).pid);
 int n = root_pids (parent, &pids);
 if (n) {
   for (int i = 0, * p = (int *) pids.p; i < n; i++, p++)
     rcv_pid_append (restriction->snd, *p, point);
   pfree (pids.p,__func__,__FILE__,0);
 }
      }
      else if ((cell.pid >= 0 && cell.pid != pid())) {

 if (is_local(aparent(0,0,0)) || has_local_child (parent))
   rcv_pid_append (restriction->rcv, cell.pid, point);
      }
    }
  }end_foreach_cell();}





  static const unsigned short keep = 1 << (user + 1);
  for (int l = depth(); l >= 0; l--)
    {foreach_cell()
      if (level == l) {
 if (level > 0 && (cell.pid < 0 || is_local(cell) || (cell.flags & used)))
   aparent(0,0,0).flags |= keep;
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) && !(cell.flags & keep))
   coarsen_cell (point, NULL);
 cell.flags &= ~(used|keep);
 continue;
      }end_foreach_cell();}


  m->send->len = m->receive->len = 0;
  rcv_pid_append_pids (mpi_level->snd, m->send);
  rcv_pid_append_pids (mpi_level_root->snd, m->send);
  rcv_pid_append_pids (mpi_level->rcv, m->receive);
  rcv_pid_append_pids (mpi_level_root->rcv, m->receive);

  prof_stop();
# 1015 "/home/spencer/basilisk/src/grid/tree-mpi.h"
end_tracing("mpi_boundary_update_buffers","/home/spencer/basilisk/src/grid/tree-mpi.h",0);}

     
void mpi_boundary_refine (scalar * list)
{tracing("mpi_boundary_refine","/home/spencer/basilisk/src/grid/tree-mpi.h",0);
  prof_start ("mpi_boundary_refine");

  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;


  Array * snd = mpi->send;
  MPI_Request r[2*snd->len/sizeof(int)];
  int nr = 0;
  for (int i = 0, * dest = snd->p; i < snd->len/sizeof(int); i++,dest++) {
    int len = ((Tree *)grid)->refined.n;
    MPI_Isend (&((Tree *)grid)->refined.n, 1, MPI_INT, *dest,
        (128), MPI_COMM_WORLD, &r[nr++]);
    if (len > 0)
      MPI_Isend (((Tree *)grid)->refined.p, sizeof(Index)/sizeof(int)*len,
   MPI_INT, *dest, (128), MPI_COMM_WORLD, &r[nr++]);
  }



  Array * rcv = mpi->receive;
  Cache rerefined = {NULL, 0, 0};
  for (int i = 0, * source = rcv->p; i < rcv->len/sizeof(int); i++,source++) {
    int len;
    mpi_recv_check (&len, 1, MPI_INT, *source, (128),
      MPI_COMM_WORLD, MPI_STATUS_IGNORE,
      "mpi_boundary_refine (len)");
    if (len > 0) {
      Index p[len];
      mpi_recv_check (p, sizeof(Index)/sizeof(int)*len,
        MPI_INT, *source, (128),
        MPI_COMM_WORLD, MPI_STATUS_IGNORE,
        "mpi_boundary_refine (p)");
      Cache refined = {p, len, len};
      {foreach_cache (refined)
 if (level <= depth() && allocated(0,0,0)) {
   if (is_leaf(cell)) {
     bool neighbors = false;
     {foreach_neighbor()
       if (allocated(0,0,0) && (is_active(cell) || is_local(aparent(0,0,0)))) {
  neighbors = true; foreach_neighbor_break();
       }end_foreach_neighbor()}

     if (neighbors)
       refine_cell (point, list, 0, &rerefined);
   }
 }end_foreach_cache();}
    }
  }


  if (nr)
    MPI_Waitall (nr, r, MPI_STATUSES_IGNORE);


  pfree (((Tree *)grid)->refined.p,__func__,__FILE__,0);
  ((Tree *)grid)->refined = rerefined;

  prof_stop();



  mpi_all_reduce (rerefined.n, MPI_INT, MPI_SUM);
  if (rerefined.n)
    mpi_boundary_refine (list);
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    _attribute[s.i].dirty = true;}}
end_tracing("mpi_boundary_refine","/home/spencer/basilisk/src/grid/tree-mpi.h",0);}

static void check_depth()
{
# 1121 "/home/spencer/basilisk/src/grid/tree-mpi.h"
}

typedef struct {
  int refined, leaf;
} Remote;



     
void mpi_boundary_coarsen (int l, int too_fine)
{tracing("mpi_boundary_coarsen","/home/spencer/basilisk/src/grid/tree-mpi.h",0);
  if (npe() == 1)
    {end_tracing("mpi_boundary_coarsen","/home/spencer/basilisk/src/grid/tree-mpi.h",0);return;}

  check_depth();

  if (!(sizeof(Remote) == sizeof(double))) qassert ("/home/spencer/basilisk/src/grid/tree-mpi.h", 0, "sizeof(Remote) == sizeof(double)");

  scalar  remote=new_scalar("remote");
  {foreach_cell() {
    if (level == l) {
      if (is_local(cell)) {
 ((Remote *)&val(remote,0,0,0))->refined = (!is_leaf (cell) && cell.neighbors && cell.pid >= 0);
 ((Remote *)&val(remote,0,0,0))->leaf = is_leaf(cell);
      }
      else {
 ((Remote *)&val(remote,0,0,0))->refined = true;
 ((Remote *)&val(remote,0,0,0))->leaf = false;
      }
      continue;
    }
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}
  mpi_boundary_level (mpi_boundary,((scalar[]) {remote,{-1}}), l);

  {foreach_cell() {
    if (level == l) {
      if (!is_local(cell)) {
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) && !((Remote *)&val(remote,0,0,0))->refined)
   coarsen_cell_recursive (point, NULL);
 else if (is_leaf(cell) && cell.neighbors && ((Remote *)&val(remote,0,0,0))->leaf) {
   int pid = cell.pid;
   {foreach_child()
     cell.pid = pid;end_foreach_child()}
 }
      }
      continue;
    }
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}

  check_depth();

  if (l > 0) {
    {foreach_cell() {
      if (level == l) {
 val(remote,0,0,0) = is_local(cell) ? cell.neighbors : 0;
 continue;
      }
      if (is_leaf(cell))
 continue;
    }end_foreach_cell();}
    mpi_boundary_level (mpi_boundary,((scalar[]) {remote,{-1}}), l);
    {foreach_cell() {
      if (level == l)
 if (!is_local(cell) && is_local(aparent(0,0,0)) && val(remote,0,0,0)) {
   aparent(0,0,0).flags &= ~too_fine;
   continue;
 }
      if (is_leaf(cell))
 continue;
    }end_foreach_cell();}
  }delete((scalar*)((scalar[]){remote,{-1}}));
end_tracing("mpi_boundary_coarsen","/home/spencer/basilisk/src/grid/tree-mpi.h",0);}

static void flag_border_cells()
{
  {foreach_cell() {
    if (is_active(cell)) {
      short flags = cell.flags & ~border;
      {foreach_neighbor() {
 if (!is_local(cell) || (level > 0 && !is_local(aparent(0,0,0)))) {
   flags |= border; foreach_neighbor_break();
 }

 if (is_refined_check())
   {foreach_child()
     if (!is_local(cell)) {
       flags |= border; foreach_child_break();
     }end_foreach_child()}
 if (flags & border)
   foreach_neighbor_break();
      }end_foreach_neighbor()}
      cell.flags = flags;
    }
    else {
      cell.flags &= ~border;

    }
    if (is_leaf(cell)) {
      if (cell.neighbors) {
 {foreach_child()
   cell.flags &= ~border;end_foreach_child()}
 if (is_border(cell)) {
   bool remote = false;
   {foreach_neighbor (2/2)
     if (!is_local(cell)) {
       remote = true; foreach_neighbor_break();
     }end_foreach_neighbor()}
   if (remote)
     {foreach_child()
       cell.flags |= border;end_foreach_child()}
 }
      }
      continue;
    }
  }end_foreach_cell();}
}

static int balanced_pid (long index, long nt, int nproc)
{
  long ne = max(1, nt/nproc), nr = nt % nproc;
  int pid = index < nr*(ne + 1) ?
    index/(ne + 1) :
    nr + (index - nr*(ne + 1))/ne;
  return min(nproc - 1, pid);
}


     
void mpi_partitioning()
{tracing("mpi_partitioning","/home/spencer/basilisk/src/grid/tree-mpi.h",0);
  prof_start ("mpi_partitioning");

  long nt = 0;
  
#if _OPENMP
  #undef OMP
  #define OMP(x)
#endif
{
#line 1258
foreach ()
    nt++;end_foreach();}
#if _OPENMP
  #undef OMP
  #define OMP(x) _Pragma(#x)
#endif



  
#line 1262
long i = 0;
  ((Tree *)grid)->dirty = true;
  {foreach_cell_post (is_active (cell))
    if (is_active (cell)) {
      if (is_leaf (cell)) {
 cell.pid = balanced_pid (i++, nt, npe());
 if (cell.neighbors > 0) {
   int pid = cell.pid;
   {foreach_child()
     cell.pid = pid;end_foreach_child()}
 }
 if (!is_local(cell))
   cell.flags &= ~active;
      }
      else {
 cell.pid = child(0,0,0).pid;
 bool inactive = true;
 {foreach_child()
   if (is_active(cell)) {
     inactive = false; foreach_child_break();
   }end_foreach_child()}
 if (inactive)
   cell.flags &= ~active;
      }
    }end_foreach_cell_post();}

  flag_border_cells();

  prof_stop();

  mpi_boundary_update_buffers();
end_tracing("mpi_partitioning","/home/spencer/basilisk/src/grid/tree-mpi.h",0);}

void restore_mpi (FILE * fp, scalar * list1)
{
  long index = 0, nt = 0, start = ftell (fp);
  scalar  size=new_scalar("size"), * list = list_concat (((scalar[]){size,{-1}}), list1);;
  long offset = sizeof(double)*list_len(list);


  static const unsigned short set = 1 << user;
  scalar * listm = is_constant(cm) ? NULL : (scalar *)((vector[]){fm,{{-1},{-1}}});
  {foreach_cell()
    if (balanced_pid (index, nt, npe()) <= pid()) {
      unsigned flags;
      if (fread (&flags, sizeof(unsigned), 1, fp) != 1) {
 fprintf (ferr, "restore(): error: expecting 'flags'\n");
 exit (1);
      }
      {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
 double val;
 if (fread (&val, sizeof(double), 1, fp) != 1) {
   fprintf (ferr, "restore(): error: expecting scalar\n");
   exit (1);
 }
 if (s.i != INT_MAX)
   val(s,0,0,0) = val;
      }}}
      if (level == 0)
 nt = val(size,0,0,0);
      cell.pid = balanced_pid (index, nt, npe());
      cell.flags |= set;
      if (!(flags & leaf) && is_leaf(cell)) {
 if (balanced_pid (index + val(size,0,0,0) - 1, nt, npe()) < pid()) {
   fseek (fp, (sizeof(unsigned) + offset)*(val(size,0,0,0) - 1), SEEK_CUR);
   index += val(size,0,0,0);
   continue;
 }
 refine_cell (point, listm, 0, NULL);
      }
      index++;
      if (is_leaf(cell))
 continue;
    }end_foreach_cell();}


  fseek (fp, start, SEEK_SET);
  index = 0;
  {foreach_cell() {
    unsigned flags;
    if (fread (&flags, sizeof(unsigned), 1, fp) != 1) {
      fprintf (ferr, "restore(): error: expecting 'flags'\n");
      exit (1);
    }
    if (cell.flags & set)
      fseek (fp, offset, SEEK_CUR);
    else {
      {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
 double val;
 if (fread (&val, sizeof(double), 1, fp) != 1) {
   fprintf (ferr, "restore(): error: expecting a scalar\n");
   exit (1);
 }
 if (s.i != INT_MAX)
   val(s,0,0,0) = val;
      }}}
      cell.pid = balanced_pid (index, nt, npe());
      if (is_leaf(cell) && cell.neighbors) {
 int pid = cell.pid;
 {foreach_child()
   cell.pid = pid;end_foreach_child()}
      }
    }
    if (!(flags & leaf) && is_leaf(cell)) {
      bool locals = false;
      {foreach_neighbor(1)
 if ((cell.flags & set) && (is_local(cell) || is_root(point))) {
   locals = true; foreach_neighbor_break();
 }end_foreach_neighbor()}
      if (locals)
 refine_cell (point, listm, 0, NULL);
      else {
 fseek (fp, (sizeof(unsigned) + offset)*(val(size,0,0,0) - 1), SEEK_CUR);
 index += val(size,0,0,0);
 continue;
      }
    }
    index++;
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}


  {foreach_cell_post (is_active (cell)) {
    cell.flags &= ~set;
    if (is_active (cell)) {
      if (is_leaf (cell)) {
 if (cell.neighbors > 0) {
   int pid = cell.pid;
   {foreach_child()
     cell.pid = pid;end_foreach_child()}
 }
 if (!is_local(cell))
   cell.flags &= ~active;
      }
      else if (!is_local(cell)) {
 bool inactive = true;
 {foreach_child()
   if (is_active(cell)) {
     inactive = false; foreach_child_break();
   }end_foreach_child()}
 if (inactive)
   cell.flags &= ~active;
      }
    }
  }end_foreach_cell_post();}

  flag_border_cells();

  mpi_boundary_update (list);
  pfree (list,__func__,__FILE__,0);delete((scalar*)((scalar[]){size,{-1}}));
}
# 1435 "/home/spencer/basilisk/src/grid/tree-mpi.h"
     
double z_indexing (scalar index, bool leaves)
{tracing("z_indexing","/home/spencer/basilisk/src/grid/tree-mpi.h",0);



  scalar  size=new_scalar("size");
  subtree_size (size, leaves);






  double maxi = -1.;
  if (pid() == 0)
    {foreach_level(0)
      maxi = val(size,0,0,0) - 1.;end_foreach_level();}




  {foreach_level(0)
    val(index,0,0,0) = 0;end_foreach_level();}
  for (int l = 0; l < depth(); l++) {
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b,((scalar[]) {index,{-1}}), l); };
    {foreach_cell() {
      if (level == l) {
 if (is_leaf(cell)) {
   if (is_local(cell) && cell.neighbors) {
     int i = val(index,0,0,0);
     {foreach_child()
       val(index,0,0,0) = i;end_foreach_child()}
   }
 }
 else {
   bool loc = is_local(cell);
   if (!loc)
     {foreach_child()
       if (is_local(cell)) {
  loc = true; foreach_child_break();
       }end_foreach_child()}
   if (loc) {
     int i = val(index,0,0,0) + !leaves;
     {foreach_child() {
       val(index,0,0,0) = i;
       i += val(size,0,0,0);
     }end_foreach_child()}
   }
 }
 continue;
      }
      if (is_leaf(cell))
 continue;
    }end_foreach_cell();}
  }
  { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b,((scalar[]) {index,{-1}}), depth()); };

  {delete((scalar*)((scalar[]){size,{-1}}));{end_tracing("z_indexing","/home/spencer/basilisk/src/grid/tree-mpi.h",0);return maxi;}}delete((scalar*)((scalar[]){size,{-1}}));
end_tracing("z_indexing","/home/spencer/basilisk/src/grid/tree-mpi.h",0);}
# 1685 "/home/spencer/basilisk/src/grid/tree.h" 2
# 1 "grid/balance.h" 1
# 1 "/home/spencer/basilisk/src/grid/balance.h"


typedef struct {
  short leaf, prolongation;
  int pid;
} NewPid;



@if TRASH
@ define is_newpid() (!isnan(val(newpid,0,0,0)) && ((NewPid *)&val(newpid,0,0,0))->pid > 0)
@else
@ define is_newpid() (((NewPid *)&val(newpid,0,0,0))->pid > 0)
@endif

Array * linear_tree (size_t size, scalar newpid)
{
  const unsigned short sent = 1 << user, next = 1 << (user + 1);
  Array * a = array_new();

  {foreach_cell_post_all (true)
    if (level > 0 && (cell.flags & (sent|next)))
      aparent(0,0,0).flags |= next;end_foreach_cell_post_all();}

  bool empty = true;
  {foreach_cell_all() {
    if (cell.flags & sent) {
      array_append (a, &cell, size);
      cell.flags &= ~sent;
      empty = false;
    }
    else {
      if (cell.pid >= 0 && ((NewPid *)&val(newpid,0,0,0))->leaf)
 if (!(is_leaf(cell))) qassert ("/home/spencer/basilisk/src/grid/balance.h", 0, "is_leaf(cell)");
      if (is_refined_check()) {


 bool prolo = false;
 {foreach_child()
   if (((NewPid *)&val(newpid,0,0,0))->prolongation)
     prolo = true;end_foreach_child()}
 if (prolo) {

   cell.flags |= leaf;
   array_append (a, &cell, sizeof(Cell));
   cell.flags &= ~leaf;
 }
 else
   array_append (a, &cell, sizeof(Cell));
      }
      else
 array_append (a, &cell, sizeof(Cell));
    }
    if (cell.flags & next)
      cell.flags &= ~next;
    else
      continue;
  }end_foreach_cell_all();}

  if (empty)
    a->len = 0;
  return a;
}

@def foreach_tree(t, size, list)
{
  const unsigned short _sent = 1 << user, _next = 1 << (user + 1);
  scalar * _list = list;
  char * _i = (char *) (t)->p;
  foreach_cell_all() {
    Cell * c = (Cell *) _i;
    if (c->flags & _sent) {
      _i += size;
@

@def end_foreach_tree()
    }
    else
      _i += sizeof(Cell);
    if (c->flags & _next) {
      if (!(c->neighbors)) qassert ("/home/spencer/basilisk/src/grid/balance.h", 0, "c->neighbors");
      if (!(c->flags & leaf) && is_leaf(cell) &&
   (!is_newpid() || !((NewPid *)&val(newpid,0,0,0))->leaf))

 refine_cell (point, _list, 0, NULL);
      else if (!cell.neighbors)

 alloc_children (point);
    }
    else
      continue;
  } end_foreach_cell_all();
}
@

Array * neighborhood (scalar newpid, int nextpid, FILE * fp)
{
  const unsigned short sent = 1 << user;
  {foreach_cell() {

    bool root = false;
    if ((!is_local(cell) || ((NewPid *)&val(newpid,0,0,0))->pid - 1 != nextpid) && (!is_leaf (cell) && cell.neighbors && cell.pid >= 0)) {
      {foreach_child()
 if (is_local(cell) && ((NewPid *)&val(newpid,0,0,0))->pid - 1 == nextpid) {
   root = true; foreach_child_break();
 }end_foreach_child()}
      if (root && cell.pid != nextpid) {
 {foreach_neighbor()
   if (cell.pid != nextpid && is_newpid()) {
     if (fp)
       fprintf (fp, "%g %g %g %d %d root\n",
         x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid);
     cell.flags |= sent;
   }end_foreach_neighbor()}
      }
    }

    if ((is_local(cell) && ((NewPid *)&val(newpid,0,0,0))->pid - 1 == nextpid) || root) {
      {foreach_neighbor(1)
 if (cell.neighbors && cell.pid != nextpid)
   {foreach_child()
     if (cell.pid != nextpid && is_newpid()) {
       if (fp)
  fprintf (fp, "%g %g %g %d %d nextpid\n",
    x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid);
       cell.flags |= sent;
     }end_foreach_child()}end_foreach_neighbor()}
    }
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}

  return linear_tree (sizeof(Cell) + datasize, newpid);
}

static void send_tree (Array * a, int to, MPI_Request * r)
{
  MPI_Isend (&a->len, 1, MPI_LONG, to, (256), MPI_COMM_WORLD, &r[0]);
  if (a->len > 0) {
    MPI_Isend (a->p, a->len, MPI_BYTE, to, (256), MPI_COMM_WORLD, &r[1]);
    ((Tree *)grid)->dirty = true;
  }
}

static void receive_tree (int from, scalar newpid, FILE * fp)
{
  Array a;
  mpi_recv_check (&a.len, 1, MPI_LONG, from, (256),
    MPI_COMM_WORLD, MPI_STATUS_IGNORE, "receive_tree (len)");
  if (a.len > 0) {
    a.p = pmalloc (a.len,__func__,__FILE__,0);
    if (fp)
      fprintf (fp, "receiving %ld from %d\n", a.len, from);
    mpi_recv_check (a.p, a.len, MPI_BYTE, from, (256),
      MPI_COMM_WORLD, MPI_STATUS_IGNORE, "receive_tree (p)");

    {foreach_tree (&a, sizeof(Cell) + datasize, NULL) {
      memcpy (((char *)&cell) + sizeof(Cell), ((char *)c) + sizeof(Cell),
       datasize);
      if (!(((NewPid *)&val(newpid,0,0,0))->pid > 0)) qassert ("/home/spencer/basilisk/src/grid/balance.h", 0, "NEWPID()->pid > 0");
      if (fp)
 fprintf (fp, "%g %g %g %d %d %d %d %d %d recv\n",
   x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid,
   c->flags & leaf,
   cell.flags & leaf, from, ((NewPid *)&val(newpid,0,0,0))->leaf);
    }end_foreach_tree();}
    pfree (a.p,__func__,__FILE__,0);
    ((Tree *)grid)->dirty = true;
  }
}

static void wait_tree (Array * a, MPI_Request * r)
{
  MPI_Wait (&r[0], MPI_STATUS_IGNORE);
  if (a->len > 0)
    MPI_Wait (&r[1], MPI_STATUS_IGNORE);
}

static void check_flags()
{







}

struct {
  int min;
  bool leaves;

  int npe;
} mpi = {
  1,
  true
};

     
bool balance()
{tracing("balance","/home/spencer/basilisk/src/grid/balance.h",0);
  if (npe() == 1)
    {end_tracing("balance","/home/spencer/basilisk/src/grid/balance.h",0);return false;}

  if (!(sizeof(NewPid) == sizeof(double))) qassert ("/home/spencer/basilisk/src/grid/balance.h", 0, "sizeof(NewPid) == sizeof(double)");

  check_flags();

  long nl = 0, nt = 0;
  {foreach_cell() {
    if (is_local(cell)) {
      nt++;
      if (is_leaf(cell))
 nl++;
    }
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}

  grid->n = grid->tn = nl;
  grid->maxdepth = depth();
  long nmin = nl, nmax = nl;

  mpi_all_reduce (nmax, MPI_LONG, MPI_MAX);
  mpi_all_reduce (nmin, MPI_LONG, MPI_MIN);
  mpi_all_reduce (grid->tn, MPI_LONG, MPI_SUM);
  mpi_all_reduce (grid->maxdepth, MPI_INT, MPI_MAX);
  if (mpi.leaves)
    nt = grid->tn;
  else
    mpi_all_reduce (nt, MPI_LONG, MPI_SUM);

  long ne = max(1, nt/npe());

  if (ne < mpi.min) {
    mpi.npe = max(1, nt/mpi.min);
    ne = max(1, nt/mpi.npe);
  }
  else
    mpi.npe = npe();

  if (nmax - nmin <= 1)
    {end_tracing("balance","/home/spencer/basilisk/src/grid/balance.h",0);return false;}

  scalar  newpid=new_scalar("newpid");
  double zn = z_indexing (newpid, mpi.leaves);
  if (pid() == 0)
    if (!(zn + 1 == nt)) qassert ("/home/spencer/basilisk/src/grid/balance.h", 0, "zn + 1 == nt");

  FILE * fp = NULL;
# 261 "/home/spencer/basilisk/src/grid/balance.h"
  bool next = false, prev = false;
  {foreach_cell_all() {
    if (is_local(cell)) {
      int pid = balanced_pid (val(newpid,0,0,0), nt, mpi.npe);
      pid = clamp (pid, cell.pid - 1, cell.pid + 1);
      if (pid == pid() + 1)
 next = true;
      else if (pid == pid() - 1)
 prev = true;
      ((NewPid *)&val(newpid,0,0,0))->pid = pid + 1;
      ((NewPid *)&val(newpid,0,0,0))->leaf = is_leaf(cell);
      ((NewPid *)&val(newpid,0,0,0))->prolongation = (!is_leaf(cell) && !cell.neighbors && cell.pid >= 0);
      if (fp)
 fprintf (fp, "%g %g %d %d newpid\n", x, y, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid);
    }
    else
      val(newpid,0,0,0) = 0;
  }end_foreach_cell_all();}
  for (int l = 0; l <= depth(); l++)
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b,((scalar[]) {newpid,{-1}}), l); };
# 305 "/home/spencer/basilisk/src/grid/balance.h"
  Array * anext = next ? neighborhood (newpid, pid() + 1, fp) : array_new();
  Array * aprev = prev ? neighborhood (newpid, pid() - 1, fp) : array_new();

  if (fp)
    fflush (fp);

  check_flags();


  MPI_Request rprev[2], rnext[2];
  if (pid() > 0)
    send_tree (aprev, pid() - 1, rprev);
  if (pid() < npe() - 1)
    send_tree (anext, pid() + 1, rnext);


  if (pid() < npe() - 1)
    receive_tree (pid() + 1, newpid, fp);
  if (pid() > 0)
    receive_tree (pid() - 1, newpid, fp);


  if (pid() > 0)
    wait_tree (aprev, rprev);
  array_free (aprev);
  if (pid() < npe() - 1)
    wait_tree (anext, rnext);
  array_free (anext);

  if (fp)
    fflush (fp);


  int pid_changed = false;
  {foreach_cell_all() {
    if (cell.pid >= 0) {
      if (is_newpid()) {
 if (fp)
   fprintf (fp, "%g %g %g %d %d %d %d %d new\n",
     x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid,
     is_leaf(cell), cell.neighbors, ((NewPid *)&val(newpid,0,0,0))->leaf);
 if (cell.pid != ((NewPid *)&val(newpid,0,0,0))->pid - 1) {
   cell.pid = ((NewPid *)&val(newpid,0,0,0))->pid - 1;
   cell.flags &= ~(active|border);
   if (is_local(cell))
     cell.flags |= active;
   pid_changed = true;
 }
 if (((NewPid *)&val(newpid,0,0,0))->leaf && !is_leaf(cell) && cell.neighbors)
   coarsen_cell_recursive (point, NULL);
      }
      else if (level > 0 && ((NewPid *)&coarse(newpid,0,0,0))->leaf)
 cell.pid = aparent(0,0,0).pid;
    }

    if (!cell.neighbors && allocated_child(0,0,0)) {
      if (fp)
 fprintf (fp, "%g %g %g %d %d freechildren\n",
   x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid);
      free_children (point);
    }
  }end_foreach_cell_all();}

  if (((Tree *)grid)->dirty || pid_changed) {


    {foreach_cell_post (!is_leaf (cell))
      if (!is_leaf(cell) && !is_local(cell)) {
 unsigned short flags = cell.flags & ~active;
 {foreach_child()
   if (is_active(cell)) {
     flags |= active; foreach_child_break();
   }end_foreach_child()}
 cell.flags = flags;
      }end_foreach_cell_post();}

    flag_border_cells();
    pid_changed = true;
  }

  if (fp)
    fclose (fp);

  mpi_all_reduce (pid_changed, MPI_INT, MPI_MAX);
  if (pid_changed)
    mpi_boundary_update_buffers();

  {delete((scalar*)((scalar[]){newpid,{-1}}));{end_tracing("balance","/home/spencer/basilisk/src/grid/balance.h",0);return pid_changed;}}delete((scalar*)((scalar[]){newpid,{-1}}));
end_tracing("balance","/home/spencer/basilisk/src/grid/balance.h",0);}

void mpi_boundary_update (scalar * list)
{
  mpi_boundary_update_buffers();
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    _attribute[s.i].dirty = true;}}
  grid->tn = 0;
  boundary_internal ((scalar *)list, "/home/spencer/basilisk/src/grid/balance.h", 0);
  while (balance());
}
# 1686 "/home/spencer/basilisk/src/grid/tree.h" 2
@else
void mpi_boundary_refine (scalar * list){}
void mpi_boundary_coarsen (int a, int b){}
void mpi_boundary_update (scalar * list) {
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    _attribute[s.i].dirty = true;}}
  boundary_internal ((scalar *)list, "/home/spencer/basilisk/src/grid/tree.h", 0);
}
@endif
# 4 "/home/spencer/basilisk/src/grid/quadtree.h" 2

void quadtree_methods() {
  tree_methods();
}
# 14 "cylinder-osc-cpp.c" 2
# 1 "cylinder-osc.c"
# 1 "navier-stokes/centered.h" 1
# 1 "/home/spencer/basilisk/src/navier-stokes/centered.h"
# 27 "/home/spencer/basilisk/src/navier-stokes/centered.h"
# 1 "./run.h" 1
# 1 "/home/spencer/basilisk/src/run.h"
# 9 "/home/spencer/basilisk/src/run.h"
double dt = 1.;

# 1 "./utils.h" 1
# 1 "/home/spencer/basilisk/src/utils.h"







double DT = 1e30, CFL = 0.5;




struct {

  long nc;

  long tnc;

  double t;

  double speed;

  timer gt;
} perf = {0};





void update_perf() {
  perf.nc += grid->n;
  perf.tnc += grid->tn;
  perf.t = timer_elapsed (perf.gt);
  perf.speed = perf.tnc/perf.t;
}






typedef struct {
  double cpu;
  double real;
  double speed;
  double min;
  double avg;
  double max;
  size_t tnc;
  long mem;
} timing;






timing timer_timing (timer t, int i, size_t tnc, double * mpi)
{
  timing s;
@if _MPI
  s.avg = mpi_time - t.tm;
@endif
  clock_t end = clock();
  s.cpu = ((double) (end - t.c))/CLOCKS_PER_SEC;
  s.real = timer_elapsed (t);
  if (tnc == 0) {
    double n = 0;
    
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(+:n)){
#line 69
foreach() n++;end_foreach();mpi_all_reduce_array(&n,double,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
    
#line 70
s.tnc = n;
    tnc = n*i;
  }
  else
    s.tnc = tnc;
@if 1
  struct rusage usage;
  getrusage (RUSAGE_SELF, &usage);
  s.mem = usage.ru_maxrss;
@else
  s.mem = 0;
@endif
@if _MPI
  if (mpi)
    MPI_Allgather (&s.avg, 1, MPI_DOUBLE, mpi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
  s.max = s.min = s.avg;
  mpi_all_reduce (s.max, MPI_DOUBLE, MPI_MAX);
  mpi_all_reduce (s.min, MPI_DOUBLE, MPI_MIN);
  mpi_all_reduce (s.avg, MPI_DOUBLE, MPI_SUM);
  mpi_all_reduce (s.real, MPI_DOUBLE, MPI_SUM);
  mpi_all_reduce (s.mem, MPI_LONG, MPI_SUM);
  s.real /= npe();
  s.avg /= npe();
  s.mem /= npe();
@else
  s.min = s.max = s.avg = 0.;
@endif
  s.speed = s.real > 0. ? tnc/s.real : -1.;
  return s;
}




void timer_print (timer t, int i, size_t tnc)
{
  timing s = timer_timing (t, i, tnc, NULL);
  fprintf (fout,
    "\n# " "Quadtree"
    ", %d steps, %g CPU, %.4g real, %.3g points.step/s, %d var\n",
    i, s.cpu, s.real, s.speed, (int) (datasize/sizeof(double)));
@if _MPI
  fprintf (fout,
    "# %d procs, MPI: min %.2g (%.2g%%) "
    "avg %.2g (%.2g%%) max %.2g (%.2g%%)\n",
    npe(),
    s.min, 100.*s.min/s.real,
    s.avg, 100.*s.avg/s.real,
    s.max, 100.*s.max/s.real);
@endif
}







typedef struct {
  double avg, rms, max, volume;
} norm;

norm normf (scalar f)
{
  double avg = 0., rms = 0., max = 0., volume = 0.;
  if(!is_constant(cm)){
  
#line 135
foreach_stencil(1,{(NonLocal[]){{"rms","double",(void *)&rms,NULL,0,'+'},{"avg","double",(void *)&avg,NULL,0,'+'},{"volume","double",(void *)&volume,NULL,0,'+'},{"max","double",(void *)&max,NULL,0,'M'},{"cm","scalar",(void *)&cm,NULL,0},{"f","scalar",(void *)&f,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n    \n// #line 137 \"/home/spencer/basilisk/src/utils.h\"\nif (val(f,0,0,0) != 1e30 && (sq(Delta)*val(cm,0,0,0)) > 0.) {\n      real v = fabs(val(f,0,0,0));\n      if (v > max) max = v;\n      volume += (sq(Delta)*val(cm,0,0,0));\n      avg += (sq(Delta)*val(cm,0,0,0))*v;\n      rms += (sq(Delta)*val(cm,0,0,0))*sq(v);\n    }")}
)
    {_stencil_val(f,0,0,0);_stencil_val(cm,0,0,0); {   
      _stencil_val(f,0,0,0);   
         
_stencil_val(cm,0,0,0); 
       _stencil_val(cm,0,0,0); 
       _stencil_val(cm,0,0,0); 
       
    
#line 143
}       }end_foreach_stencil();
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(+:volume)
   reduction(+:rms) reduction(+:avg)reduction(max:max)){
#line 135
foreach(
)
    if (val(f,0,0,0) != 1e30 && (sq(Delta)*val(cm,0,0,0)) > 0.) {
      double v = fabs(val(f,0,0,0));
      if (v > max) max = v;
      volume += (sq(Delta)*val(cm,0,0,0));
      avg += (sq(Delta)*val(cm,0,0,0))*v;
      rms += (sq(Delta)*val(cm,0,0,0))*sq(v);
    }end_foreach();mpi_all_reduce_array(&volume,double,MPI_SUM,1);mpi_all_reduce_array(&rms,double,MPI_SUM,1);mpi_all_reduce_array(&avg,double,MPI_SUM,1);mpi_all_reduce_array(&max,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 143
}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#line 135
foreach_stencil(1,{(NonLocal[]){{"rms","double",(void *)&rms,NULL,0,'+'},{"avg","double",(void *)&avg,NULL,0,'+'},{"volume","double",(void *)&volume,NULL,0,'+'},{"max","double",(void *)&max,NULL,0,'M'},{"_const_cm","double",(void *)&_const_cm,NULL,0},{"f","scalar",(void *)&f,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n    \n// #line 137 \"/home/spencer/basilisk/src/utils.h\"\nif (val(f,0,0,0) != 1e30 && (sq(Delta)*_const_cm) > 0.) {\n      real v = fabs(val(f,0,0,0));\n      if (v > max) max = v;\n      volume += (sq(Delta)*_const_cm);\n      avg += (sq(Delta)*_const_cm)*v;\n      rms += (sq(Delta)*_const_cm)*sq(v);\n    }")}
)
    {_stencil_val(f,0,0,0);; {   
      _stencil_val(f,0,0,0);

;
;
; 
       
    
#line 143
}       }end_foreach_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(+:volume)
   reduction(+:rms) reduction(+:avg)reduction(max:max)){
#line 135
foreach(
)
    if (val(f,0,0,0) != 1e30 && (sq(Delta)*_const_cm) > 0.) {
      double v = fabs(val(f,0,0,0));
      if (v > max) max = v;
      volume += (sq(Delta)*_const_cm);
      avg += (sq(Delta)*_const_cm)*v;
      rms += (sq(Delta)*_const_cm)*sq(v);
    }end_foreach();mpi_all_reduce_array(&volume,double,MPI_SUM,1);mpi_all_reduce_array(&rms,double,MPI_SUM,1);mpi_all_reduce_array(&avg,double,MPI_SUM,1);mpi_all_reduce_array(&max,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 143
}
  norm n;
  n.avg = volume ? avg/volume : 0.;
  n.rms = volume ? sqrt(rms/volume) : 0.;
  n.max = max;
  n.volume = volume;
  return n;
}





typedef struct {
  double min, max, sum, stddev, volume;
} stats;

stats statsf (scalar f)
{
  double min = 1e100, max = -1e100, sum = 0., sum2 = 0., volume = 0.;
  if(!is_constant(cm)){
  
#line 163
foreach_stencil(1,{(NonLocal[]){{"min","double",(void *)&min,NULL,0,'m'},{"max","double",(void *)&max,NULL,0,'M'},{"sum2","double",(void *)&sum2,NULL,0,'+'},{"sum","double",(void *)&sum,NULL,0,'+'},{"volume","double",(void *)&volume,NULL,0,'+'},{"f","scalar",(void *)&f,NULL,0},{"cm","scalar",(void *)&cm,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n    \n// #line 165 \"/home/spencer/basilisk/src/utils.h\"\nif ((sq(Delta)*val(cm,0,0,0)) > 0. && val(f,0,0,0) != 1e30) {\n      volume += (sq(Delta)*val(cm,0,0,0));\n      sum += (sq(Delta)*val(cm,0,0,0))*val(f,0,0,0);\n      sum2 += (sq(Delta)*val(cm,0,0,0))*sq(val(f,0,0,0));\n      if (val(f,0,0,0) > max) max = val(f,0,0,0);\n      if (val(f,0,0,0) < min) min = val(f,0,0,0);\n    }")}
)
    {_stencil_val(cm,0,0,0); _stencil_val(f,0,0,0); {
_stencil_val(cm,0,0,0); 
       _stencil_val(cm,0,0,0);_stencil_val(f,0,0,0); 
       _stencil_val(cm,0,0,0);_stencil_val(f,0,0,0); 
       _stencil_val(f,0,0,0); { _stencil_val(f,0,0,0); }
_stencil_val(f,0,0,0); { _stencil_val(f,0,0,0); }
         
         
    
#line 171
}      }end_foreach_stencil();
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(min:min)
   reduction(max:max) reduction(+:volume) reduction(+:sum2)reduction(+:sum)){
#line 163
foreach(
)
    if ((sq(Delta)*val(cm,0,0,0)) > 0. && val(f,0,0,0) != 1e30) {
      volume += (sq(Delta)*val(cm,0,0,0));
      sum += (sq(Delta)*val(cm,0,0,0))*val(f,0,0,0);
      sum2 += (sq(Delta)*val(cm,0,0,0))*sq(val(f,0,0,0));
      if (val(f,0,0,0) > max) max = val(f,0,0,0);
      if (val(f,0,0,0) < min) min = val(f,0,0,0);
    }end_foreach();mpi_all_reduce_array(&min,double,MPI_MIN,1);mpi_all_reduce_array(&max,double,MPI_MAX,1);mpi_all_reduce_array(&volume,double,MPI_SUM,1);mpi_all_reduce_array(&sum2,double,MPI_SUM,1);mpi_all_reduce_array(&sum,double,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 171
}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#line 163
foreach_stencil(1,{(NonLocal[]){{"min","double",(void *)&min,NULL,0,'m'},{"max","double",(void *)&max,NULL,0,'M'},{"sum2","double",(void *)&sum2,NULL,0,'+'},{"sum","double",(void *)&sum,NULL,0,'+'},{"volume","double",(void *)&volume,NULL,0,'+'},{"f","scalar",(void *)&f,NULL,0},{"_const_cm","double",(void *)&_const_cm,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n    \n// #line 165 \"/home/spencer/basilisk/src/utils.h\"\nif ((sq(Delta)*_const_cm) > 0. && val(f,0,0,0) != 1e30) {\n      volume += (sq(Delta)*_const_cm);\n      sum += (sq(Delta)*_const_cm)*val(f,0,0,0);\n      sum2 += (sq(Delta)*_const_cm)*sq(val(f,0,0,0));\n      if (val(f,0,0,0) > max) max = val(f,0,0,0);\n      if (val(f,0,0,0) < min) min = val(f,0,0,0);\n    }")}
)
    {; _stencil_val(f,0,0,0); {
;
;_stencil_val(f,0,0,0);
;_stencil_val(f,0,0,0); 
       _stencil_val(f,0,0,0); { _stencil_val(f,0,0,0); }
_stencil_val(f,0,0,0); { _stencil_val(f,0,0,0); }
         
         
    
#line 171
}      }end_foreach_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(min:min)
   reduction(max:max) reduction(+:volume) reduction(+:sum2)reduction(+:sum)){
#line 163
foreach(
)
    if ((sq(Delta)*_const_cm) > 0. && val(f,0,0,0) != 1e30) {
      volume += (sq(Delta)*_const_cm);
      sum += (sq(Delta)*_const_cm)*val(f,0,0,0);
      sum2 += (sq(Delta)*_const_cm)*sq(val(f,0,0,0));
      if (val(f,0,0,0) > max) max = val(f,0,0,0);
      if (val(f,0,0,0) < min) min = val(f,0,0,0);
    }end_foreach();mpi_all_reduce_array(&min,double,MPI_MIN,1);mpi_all_reduce_array(&max,double,MPI_MAX,1);mpi_all_reduce_array(&volume,double,MPI_SUM,1);mpi_all_reduce_array(&sum2,double,MPI_SUM,1);mpi_all_reduce_array(&sum,double,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 171
}
  stats s;
  s.min = min, s.max = max, s.sum = sum, s.volume = volume;
  if (volume > 0.)
    sum2 -= sum*sum/volume;
  s.stddev = sum2 > 0. ? sqrt(sum2/volume) : 0.;
  return s;
}
# 187 "/home/spencer/basilisk/src/utils.h"
static double generic_limiter (double r, double beta)
{
  double v1 = min (r, beta), v2 = min (beta*r, 1.);
  v1 = max (0., v1);
  return max (v1, v2);
}

double minmod (double s0, double s1, double s2) {
  return s1 == s0 ? 0. : generic_limiter ((s2 - s1)/(s1 - s0), 1.)*(s1 - s0);
}

double superbee (double s0, double s1, double s2) {
  return s1 == s0 ? 0. : generic_limiter ((s2 - s1)/(s1 - s0), 2.)*(s1 - s0);
}

double sweby (double s0, double s1, double s2) {
  return s1 == s0 ? 0. : generic_limiter ((s2 - s1)/(s1 - s0), 1.5)*(s1 - s0);
}
# 213 "/home/spencer/basilisk/src/utils.h"
double theta = 1.3;

double minmod2 (double s0, double s1, double s2)
{
  if (s0 < s1 && s1 < s2) {
    double d1 = theta*(s1 - s0), d2 = (s2 - s0)/2., d3 = theta*(s2 - s1);
    if (d2 < d1) d1 = d2;
    return min(d1, d3);
  }
  if (s0 > s1 && s1 > s2) {
    double d1 = theta*(s1 - s0), d2 = (s2 - s0)/2., d3 = theta*(s2 - s1);
    if (d2 > d1) d1 = d2;
    return max(d1, d3);
  }
  return 0.;
}
# 237 "/home/spencer/basilisk/src/utils.h"
void gradients (scalar * f, vector * g)
{
  if (!(list_len(f) == vectors_len(g))) qassert ("/home/spencer/basilisk/src/utils.h", 0, "list_len(f) == vectors_len(g)");
  foreach_stencil(1,{(NonLocal[]){{"g","vector",(void *)g,NULL,1},{"f","scalar",(void *)f,NULL,1},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_(" \n// #line 240 \"/home/spencer/basilisk/src/utils.h\"\n{\n    scalar s; vector v;\n    {forin2 (s,v , f,g) {\n      if (s.gradient)\n { {\n\n\n\n\n\n     val_out_(v.x,0,0,0) = s.gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0))/Delta;\n } \n// #line 244\n{\n\n\n\n\n\n     val_out_(v.y,0,0,0) = s.gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0))/Delta;\n }}\n      else\n { {\n\n\n\n\n\n     val_out_(v.x,0,0,0) = (val(s,1,0,0) - val(s,-1,0,0))/(2.*Delta);\n } \n// #line 253\n{\n\n\n\n\n\n     val_out_(v.y,0,0,0) = (val(s,0,1,0) - val(s,0,-1,0))/(2.*Delta);\n }}\n    } endforin2()}\n  }")}) {
    scalar s; vector v;
    {vector*_i0=g;scalar*_i1= f;if(_i0)for(v=*_i0,s=*_i1;_i0->x.i>= 0;v=*++_i0,s=*++_i1){ {
      if (_attribute[s.i].gradient)
 { {





_stencil_val(s,-1,0,0); _stencil_val(s,0,0,0); _stencil_val(s,1,0,0);





     
#line 250
_stencil_val_a(v.x,0,0,0);   
 } 
#line 244
{





_stencil_val(s,0,-1,0); _stencil_val(s,0,0,0); _stencil_val(s,0,1,0);





     
#line 250
_stencil_val_a(v.y,0,0,0);   
 }}
      else
 { {





_stencil_val(s,1,0,0); _stencil_val(s,-1,0,0);





     
#line 259
_stencil_val_a(v.x,0,0,0);   
 } 
#line 253
{





_stencil_val(s,0,1,0); _stencil_val(s,0,-1,0);





     
#line 259
_stencil_val_a(v.y,0,0,0);   
 }}
    }}}
  }end_foreach_stencil();
  {
#line 240
foreach() {
    scalar s; vector v;
    {vector*_i0=g;scalar*_i1= f;if(_i0)for(v=*_i0,s=*_i1;_i0->x.i>= 0;v=*++_i0,s=*++_i1){ {
      if (_attribute[s.i].gradient)
 { {





     val(v.x,0,0,0) = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0))/Delta;
 } 
#line 244
{





     val(v.y,0,0,0) = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0))/Delta;
 }}
      else
 { {





     val(v.x,0,0,0) = (val(s,1,0,0) - val(s,-1,0,0))/(2.*Delta);
 } 
#line 253
{





     val(v.y,0,0,0) = (val(s,0,1,0) - val(s,0,-1,0))/(2.*Delta);
 }}
    }}}
  }end_foreach();}
}
# 280 "/home/spencer/basilisk/src/utils.h"
void vorticity (const vector u, scalar omega)
{
  if(!is_constant(fm.x) && !is_constant(cm)){
  
#line 282
foreach_stencil(1,{(NonLocal[]){{"cm","scalar",(void *)&cm,NULL,0},{"u","vector",(void *)&u,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"omega","scalar",(void *)&omega,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n    \n// #line 283 \"/home/spencer/basilisk/src/utils.h\"\nval_out_(omega,0,0,0) = ((val(fm.x,1,0,0) - val(fm.x,0,0,0))*val(u.y,0,0,0) +\n        val(fm.x,1,0,0)*val(u.y,1,0,0) - val(fm.x,0,0,0)*val(u.y,-1,0,0) -\n        (val(fm.y,0,1,0) - val(fm.y,0,0,0))*val(u.x,0,0,0) +\n        val(fm.y,0,0,0)*val(u.x,0,-1,0) - val(fm.y,0,1,0)*val(u.x,0,1,0))/(2.*(val(cm,0,0,0) + 0.)*Delta);")})
    {_stencil_val(fm.x,1,0,0); _stencil_val(fm.x,0,0,0);_stencil_val(u.y,0,0,0);
        _stencil_val(fm.x,1,0,0);_stencil_val(u.y,1,0,0); _stencil_val(fm.x,0,0,0);_stencil_val(u.y,-1,0,0);
_stencil_val(fm.y,0,1,0); _stencil_val(fm.y,0,0,0);_stencil_val(u.x,0,0,0);
        _stencil_val(fm.y,0,0,0);_stencil_val(u.x,0,-1,0); _stencil_val(fm.y,0,1,0);_stencil_val(u.x,0,1,0);_stencil_val(cm,0,0,0);
#line 283
_stencil_val_a(omega,0,0,0);      
             

}end_foreach_stencil();{
#line 282
foreach()
    val(omega,0,0,0) = ((val(fm.x,1,0,0) - val(fm.x,0,0,0))*val(u.y,0,0,0) +
        val(fm.x,1,0,0)*val(u.y,1,0,0) - val(fm.x,0,0,0)*val(u.y,-1,0,0) -
        (val(fm.y,0,1,0) - val(fm.y,0,0,0))*val(u.x,0,0,0) +
        val(fm.y,0,0,0)*val(u.x,0,-1,0) - val(fm.y,0,1,0)*val(u.x,0,1,0))/(2.*(val(cm,0,0,0) + 0.)*Delta);end_foreach();}}else if(is_constant(fm.x) && !is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 282
foreach_stencil(1,{(NonLocal[]){{"cm","scalar",(void *)&cm,NULL,0},{"u","vector",(void *)&u,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"omega","scalar",(void *)&omega,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n    \n// #line 283 \"/home/spencer/basilisk/src/utils.h\"\nval_out_(omega,0,0,0) = ((_const_fm.x - _const_fm.x)*val(u.y,0,0,0) +\n        _const_fm.x*val(u.y,1,0,0) - _const_fm.x*val(u.y,-1,0,0) -\n        (_const_fm.y - _const_fm.y)*val(u.x,0,0,0) +\n        _const_fm.y*val(u.x,0,-1,0) - _const_fm.y*val(u.x,0,1,0))/(2.*(val(cm,0,0,0) + 0.)*Delta);")})
    {;;_stencil_val(u.y,0,0,0);
;_stencil_val(u.y,1,0,0);;_stencil_val(u.y,-1,0,0);
;;_stencil_val(u.x,0,0,0);
;_stencil_val(u.x,0,-1,0);;_stencil_val(u.x,0,1,0);_stencil_val(cm,0,0,0);
#line 283
_stencil_val_a(omega,0,0,0);      
             

}end_foreach_stencil();
  {
#line 282
foreach()
    val(omega,0,0,0) = ((_const_fm.x - _const_fm.x)*val(u.y,0,0,0) +
        _const_fm.x*val(u.y,1,0,0) - _const_fm.x*val(u.y,-1,0,0) -
        (_const_fm.y - _const_fm.y)*val(u.x,0,0,0) +
        _const_fm.y*val(u.x,0,-1,0) - _const_fm.y*val(u.x,0,1,0))/(2.*(val(cm,0,0,0) + 0.)*Delta);end_foreach();}}else if(!is_constant(fm.x) && is_constant(cm)){double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#line 282
foreach_stencil(1,{(NonLocal[]){{"_const_cm","double",(void *)&_const_cm,NULL,0},{"u","vector",(void *)&u,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"omega","scalar",(void *)&omega,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n    \n// #line 283 \"/home/spencer/basilisk/src/utils.h\"\nval_out_(omega,0,0,0) = ((val(fm.x,1,0,0) - val(fm.x,0,0,0))*val(u.y,0,0,0) +\n        val(fm.x,1,0,0)*val(u.y,1,0,0) - val(fm.x,0,0,0)*val(u.y,-1,0,0) -\n        (val(fm.y,0,1,0) - val(fm.y,0,0,0))*val(u.x,0,0,0) +\n        val(fm.y,0,0,0)*val(u.x,0,-1,0) - val(fm.y,0,1,0)*val(u.x,0,1,0))/(2.*(_const_cm + 0.)*Delta);")})
    {_stencil_val(fm.x,1,0,0); _stencil_val(fm.x,0,0,0);_stencil_val(u.y,0,0,0);
        _stencil_val(fm.x,1,0,0);_stencil_val(u.y,1,0,0); _stencil_val(fm.x,0,0,0);_stencil_val(u.y,-1,0,0);
_stencil_val(fm.y,0,1,0); _stencil_val(fm.y,0,0,0);_stencil_val(u.x,0,0,0);
        _stencil_val(fm.y,0,0,0);_stencil_val(u.x,0,-1,0); _stencil_val(fm.y,0,1,0);_stencil_val(u.x,0,1,0);;
#line 283
_stencil_val_a(omega,0,0,0);      
             

}end_foreach_stencil();
  {
#line 282
foreach()
    val(omega,0,0,0) = ((val(fm.x,1,0,0) - val(fm.x,0,0,0))*val(u.y,0,0,0) +
        val(fm.x,1,0,0)*val(u.y,1,0,0) - val(fm.x,0,0,0)*val(u.y,-1,0,0) -
        (val(fm.y,0,1,0) - val(fm.y,0,0,0))*val(u.x,0,0,0) +
        val(fm.y,0,0,0)*val(u.x,0,-1,0) - val(fm.y,0,1,0)*val(u.x,0,1,0))/(2.*(_const_cm + 0.)*Delta);end_foreach();}}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#line 282
foreach_stencil(1,{(NonLocal[]){{"_const_cm","double",(void *)&_const_cm,NULL,0},{"u","vector",(void *)&u,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"omega","scalar",(void *)&omega,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n    \n// #line 283 \"/home/spencer/basilisk/src/utils.h\"\nval_out_(omega,0,0,0) = ((_const_fm.x - _const_fm.x)*val(u.y,0,0,0) +\n        _const_fm.x*val(u.y,1,0,0) - _const_fm.x*val(u.y,-1,0,0) -\n        (_const_fm.y - _const_fm.y)*val(u.x,0,0,0) +\n        _const_fm.y*val(u.x,0,-1,0) - _const_fm.y*val(u.x,0,1,0))/(2.*(_const_cm + 0.)*Delta);")})
    {;;_stencil_val(u.y,0,0,0);
;_stencil_val(u.y,1,0,0);;_stencil_val(u.y,-1,0,0);
;;_stencil_val(u.x,0,0,0);
;_stencil_val(u.x,0,-1,0);;_stencil_val(u.x,0,1,0);;
#line 283
_stencil_val_a(omega,0,0,0);      
             

}end_foreach_stencil();
  {
#line 282
foreach()
    val(omega,0,0,0) = ((_const_fm.x - _const_fm.x)*val(u.y,0,0,0) +
        _const_fm.x*val(u.y,1,0,0) - _const_fm.x*val(u.y,-1,0,0) -
        (_const_fm.y - _const_fm.y)*val(u.x,0,0,0) +
        _const_fm.y*val(u.x,0,-1,0) - _const_fm.y*val(u.x,0,1,0))/(2.*(_const_cm + 0.)*Delta);end_foreach();}}
}





double change (scalar s, scalar sn)
{
  double max = 0.;
  if(!is_constant(cm)){
  
#line 296
foreach_stencil(1,{(NonLocal[]){{"max","double",(void *)&max,NULL,0,'M'},{"sn","scalar",(void *)&sn,NULL,0},{"s","scalar",(void *)&s,NULL,0},{"cm","scalar",(void *)&cm,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_(" \n// #line 296 \"/home/spencer/basilisk/src/utils.h\"\n{\n    if ((sq(Delta)*val(cm,0,0,0)) > 0.) {\n      real ds = fabs (val(s,0,0,0) - val(sn,0,0,0));\n      if (ds > max)\n max = ds;\n    }\n    val_out_(sn,0,0,0) = val(s,0,0,0);\n  }")}) {
_stencil_val(cm,0,0,0); {     
       _stencil_val(sn,0,0,0);_stencil_val(s,0,0,0);   
       
  
    } 
_stencil_val(s,0,0,0);
       
    
#line 302
_stencil_val_a(sn,0,0,0); 
  }end_foreach_stencil();
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:max)){
#line 296
foreach() {
    if ((sq(Delta)*val(cm,0,0,0)) > 0.) {
      double ds = fabs (val(s,0,0,0) - val(sn,0,0,0));
      if (ds > max)
 max = ds;
    }
    val(sn,0,0,0) = val(s,0,0,0);
  }end_foreach();mpi_all_reduce_array(&max,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 303
}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#line 296
foreach_stencil(1,{(NonLocal[]){{"max","double",(void *)&max,NULL,0,'M'},{"sn","scalar",(void *)&sn,NULL,0},{"s","scalar",(void *)&s,NULL,0},{"_const_cm","double",(void *)&_const_cm,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_(" \n// #line 296 \"/home/spencer/basilisk/src/utils.h\"\n{\n    if ((sq(Delta)*_const_cm) > 0.) {\n      real ds = fabs (val(s,0,0,0) - val(sn,0,0,0));\n      if (ds > max)\n max = ds;\n    }\n    val_out_(sn,0,0,0) = val(s,0,0,0);\n  }")}) {
; {     
       _stencil_val(sn,0,0,0);_stencil_val(s,0,0,0);   
       
  
    } 
_stencil_val(s,0,0,0);
       
    
#line 302
_stencil_val_a(sn,0,0,0); 
  }end_foreach_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:max)){
#line 296
foreach() {
    if ((sq(Delta)*_const_cm) > 0.) {
      double ds = fabs (val(s,0,0,0) - val(sn,0,0,0));
      if (ds > max)
 max = ds;
    }
    val(sn,0,0,0) = val(s,0,0,0);
  }end_foreach();mpi_all_reduce_array(&max,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 303
}
  return max;
}





scalar lookup_field (const char * name)
{
  if (name)
    {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      if (!strcmp (_attribute[s.i].name, name))
 return s;}}
  return (scalar){-1};
}

vector lookup_vector (const char * name)
{
  if (name) {
    char component[strlen(name) + 3];
    strcpy (component, name);
    strcat (component, ".x");
    {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      if (!strcmp (_attribute[s.i].name, component))
 return _attribute[s.i].v;}}
  }
  return (vector){{-1}};
}
# 340 "/home/spencer/basilisk/src/utils.h"
@def foreach_segment(_S,_p) {
  coord t = {(_S)[1].x - (_S)[0].x, (_S)[1].y - (_S)[0].y};
  double norm = sqrt(sq(t.x) + sq(t.y));
  if (!(norm > 0.)) qassert ("/home/spencer/basilisk/src/utils.h", 0, "norm > 0.");
  t.x = t.x/norm + 1e-6, t.y = t.y/norm - 1.5e-6;
  double alpha = ((_S)[0].x*((_S)[1].y - (_S)[0].y) -
    (_S)[0].y*((_S)[1].x - (_S)[0].x))/norm;
  foreach()
    if (fabs(t.y*x - t.x*y - alpha) < 0.708*Delta) {
      coord _o = {x,y}, _p[2];
      int _n = 0;
 if (t.x)
   for (int _i = -1; _i <= 1 && _n < 2; _i += 2) {
     _p[_n].x = _o.x + _i*Delta/2.;
     double a = (_p[_n].x - (_S)[0].x)/t.x;
     _p[_n].y = (_S)[0].y + a*t.y;
     if (fabs(_p[_n].y - _o.y) <= Delta/2.) {
       a = clamp (a, 0., norm);
       _p[_n].x = (_S)[0].x + a*t.x, _p[_n].y = (_S)[0].y + a*t.y;
       if (fabs(_p[_n].x - _o.x) <= Delta/2. &&
    fabs(_p[_n].y - _o.y) <= Delta/2.)
  _n++;
     }
   }

 if (t.y)
   for (int _i = -1; _i <= 1 && _n < 2; _i += 2) {
     _p[_n].y = _o.y + _i*Delta/2.;
     double a = (_p[_n].y - (_S)[0].y)/t.y;
     _p[_n].x = (_S)[0].x + a*t.x;
     if (fabs(_p[_n].x - _o.x) <= Delta/2.) {
       a = clamp (a, 0., norm);
       _p[_n].y = (_S)[0].y + a*t.y, _p[_n].x = (_S)[0].x + a*t.x;
       if (fabs(_p[_n].y - _o.y) <= Delta/2. &&
    fabs(_p[_n].x - _o.x) <= Delta/2.)
  _n++;
     }
   }

      if (_n == 2) {
@
# 410 "/home/spencer/basilisk/src/utils.h"
@define end_foreach_segment() } } end_foreach(); }




void fields_stats()
{
  fprintf (ferr, "# t = %g, fields = {", t);
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    fprintf (ferr, " %s", _attribute[s.i].name);}}
  fputs (" }\n", ferr);
  fprintf (ferr, "# %12s: %12s %12s %12s %12s\n",
    "name", "min", "avg", "stddev", "max");
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    stats ss = statsf (s);
    fprintf (ferr, "# %12s: %12g %12g %12g %12g\n",
      _attribute[s.i].name, ss.min, ss.sum/ss.volume, ss.stddev, ss.max);
  }}}
}

# 1 "./output.h" 1
# 1 "/home/spencer/basilisk/src/output.h"
# 37 "/home/spencer/basilisk/src/output.h"
     
void output_field (scalar * list,
     FILE * fp,
     int n,
     bool linear,
     coord box[2])
{tracing("output_field","/home/spencer/basilisk/src/output.h",0);
  n++;
  int len = list_len (list);
  double Delta = 0.999999*(box[1].x - box[0].x)/(n - 1);
  int ny = (box[1].y - box[0].y)/Delta + 1;
  double ** field = (double **) matrix_new (n, ny, len*sizeof(double)), * v = field[0];
  for (int i = 0; i < n*ny*len; i++, v++)
    *v = 1e30;
  coord box1[2] = {{box[0].x - Delta/2., box[0].y - Delta/2.},
     {box[0].x + (n - 0.5)*Delta, box[0].y + (ny - 0.5)*Delta}};
  coord cn = {n, ny}, p;




  foreach_region_stencil (1,{(NonLocal[]){{"linear","bool",(void *)&linear,NULL,0},{"len","int",(void *)&len,NULL,0},{"list","scalar",(void *)list,NULL,1},{"field","double",(void *)field,NULL,2},{"cn","coord",(void *)&cn,NULL,0},{"box1","coord",(void *)box1,(int[]){2,0},0},{"p","coord",(void *)&p,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_("\n\n\n// #line 805 \"/home/spencer/basilisk/src/grid/cartesian-common.h\"\nstatic real interpolate_linear (Point point, scalar v,\n      real xp, real yp, real zp)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n\n\n\n\n\n\n\n  x = (xp - x)/Delta - v.d.x/2.;\n  y = (yp - y)/Delta - v.d.y/2.;\n  int i = sign(x), j = sign(y);\n  x = fabs(x); y = fabs(y);\n\n  return ((val(v,0,0,0)*(1. - x) + val(v,i,0,0)*x)*(1. - y) +\n   (val(v,0,j,0)*(1. - x) + val(v,i,j,0)*x)*y);\n// # 834 \"/home/spencer/basilisk/src/grid/cartesian-common.h\"\n}"),_("\n\n  \n// #line 60 \"/home/spencer/basilisk/src/output.h\"\n{\n    real ** alias = field;\n    int i = (p.x - box1[0].x)/(box1[1].x - box1[0].x)*cn.x;\n    int j = (p.y - box1[0].y)/(box1[1].y - box1[0].y)*cn.y;\n    int k = 0;\n    {forin (scalar, s , list)\n      alias[i][len*j + k++] = linear ? interpolate_linear (point, s, p.x, p.y, p.z) : val(s,0,0,0); endforin()}\n  }")})

  {                     
    
    
    
    
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      { _stencil_interpolate_linear (point, s, NULL, NULL, NULL); _stencil_val(s,0,0,0);      }}}
  }end_foreach_region_stencil();




  {
#line 58
foreach_region (p, box1, cn)

  {
    double ** alias = field;
    int i = (p.x - box1[0].x)/(box1[1].x - box1[0].x)*cn.x;
    int j = (p.y - box1[0].y)/(box1[1].y - box1[0].y)*cn.y;
    int k = 0;
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      alias[i][len*j + k++] = linear ? interpolate_linear (point, s, p.x, p.y, p.z) : val(s,0,0,0);}}
  }end_foreach_region();}

  if (pid() == 0) {
    fprintf (fp, "# 1:x 2:y");
    int i = 3;
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      fprintf (fp, " %d:%s", i++, _attribute[s.i].name);}}
    fputc('\n', fp);
    for (int i = 0; i < n; i++) {
      double x = Delta*i + box[0].x;
      for (int j = 0; j < ny; j++) {
 double y = Delta*j + box[0].y;

 fprintf (fp, "%g %g", x, y);
 int k = 0;
 {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
   fprintf (fp, " %g", field[i][len*j + k++]);}}
 fputc ('\n', fp);
      }
      fputc ('\n', fp);
    }
    fflush (fp);
  }

  matrix_free (field);
end_tracing("output_field","/home/spencer/basilisk/src/output.h",0);}
# 120 "/home/spencer/basilisk/src/output.h"
     
void output_matrix (scalar f, FILE * fp, int n, bool linear)
{tracing("output_matrix","/home/spencer/basilisk/src/output.h",0);
  float fn = n;
  float Delta = (float) L0/fn;
  fwrite (&fn, sizeof(float), 1, fp);
  for (int j = 0; j < n; j++) {
    float yp = (float) (Delta*j + X0 + Delta/2.);
    fwrite (&yp, sizeof(float), 1, fp);
  }
  for (int i = 0; i < n; i++) {
    float xp = (float) (Delta*i + X0 + Delta/2.);
    fwrite (&xp, sizeof(float), 1, fp);
    for (int j = 0; j < n; j++) {
      float yp = (float)(Delta*j + Y0 + Delta/2.), v;
      v = interpolate (f, xp, yp
#line 837 "/home/spencer/basilisk/src/grid/cartesian-common.h"
, 0.
#line 135 "/home/spencer/basilisk/src/output.h"
, linear);
      fwrite (&v, sizeof(float), 1, fp);
    }
  }
  fflush (fp);
end_tracing("output_matrix","/home/spencer/basilisk/src/output.h",0);}
# 149 "/home/spencer/basilisk/src/output.h"
typedef void (* Colormap) (double cmap[127][3]);

void jet (double cmap[127][3])
{
  for (int i = 0; i < 127; i++) {
    cmap[i][0] =
      i <= 46 ? 0. :
      i >= 111 ? -0.03125*(i - 111) + 1. :
      i >= 78 ? 1. :
      0.03125*(i - 46);
    cmap[i][1] =
      i <= 14 || i >= 111 ? 0. :
      i >= 79 ? -0.03125*(i - 111) :
      i <= 46 ? 0.03125*(i - 14) :
      1.;
    cmap[i][2] =
      i >= 79 ? 0. :
      i >= 47 ? -0.03125*(i - 79) :
      i <= 14 ? 0.03125*(i - 14) + 1.:
      1.;
  }
}

void cool_warm (double cmap[127][3])
{






  static double basemap[33][3] = {
    {0.2298057, 0.298717966, 0.753683153},
    {0.26623388, 0.353094838, 0.801466763},
    {0.30386891, 0.406535296, 0.84495867},
    {0.342804478, 0.458757618, 0.883725899},
    {0.38301334, 0.50941904, 0.917387822},
    {0.424369608, 0.558148092, 0.945619588},
    {0.46666708, 0.604562568, 0.968154911},
    {0.509635204, 0.648280772, 0.98478814},
    {0.552953156, 0.688929332, 0.995375608},
    {0.596262162, 0.726149107, 0.999836203},
    {0.639176211, 0.759599947, 0.998151185},
    {0.681291281, 0.788964712, 0.990363227},
    {0.722193294, 0.813952739, 0.976574709},
    {0.761464949, 0.834302879, 0.956945269},
    {0.798691636, 0.849786142, 0.931688648},
    {0.833466556, 0.860207984, 0.901068838},
    {0.865395197, 0.86541021, 0.865395561},
    {0.897787179, 0.848937047, 0.820880546},
    {0.924127593, 0.827384882, 0.774508472},
    {0.944468518, 0.800927443, 0.726736146},
    {0.958852946, 0.769767752, 0.678007945},
    {0.96732803, 0.734132809, 0.628751763},
    {0.969954137, 0.694266682, 0.579375448},
    {0.966811177, 0.650421156, 0.530263762},
    {0.958003065, 0.602842431, 0.481775914},
    {0.943660866, 0.551750968, 0.434243684},
    {0.923944917, 0.49730856, 0.387970225},
    {0.89904617, 0.439559467, 0.343229596},
    {0.869186849, 0.378313092, 0.300267182},
    {0.834620542, 0.312874446, 0.259301199},
    {0.795631745, 0.24128379, 0.220525627},
    {0.752534934, 0.157246067, 0.184115123},
    {0.705673158, 0.01555616, 0.150232812}
  };

  for (int i = 0; i < 127; i++) {
    double x = i*(32 - 1e-10)/(127 - 1);
    int j = x; x -= j;
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (1. - x)*basemap[j][k] + x*basemap[j+1][k];
  }
}

void gray (double cmap[127][3])
{
  for (int i = 0; i < 127; i++)
    for (int k = 0; k < 3; k++)
      cmap[i][k] = i/(127 - 1.);
}

void randomap (double cmap[127][3])
{
  srand(0);
  for (int i = 0; i < 127; i++)
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (noise() + 1.)/2.;
}

void blue_white_red (double cmap[127][3])
{
  for (int i = 0; i < (127 + 1)/2; i++) {
    cmap[i][0] = i/((127 - 1)/2.);
    cmap[i][1] = i/((127 - 1)/2.);
    cmap[i][2] = 1.;
  }
  for (int i = 0; i < (127 - 1)/2; i++) {
    cmap[i + (127 + 1)/2][0] = 1.;
    cmap[i + (127 + 1)/2][1] = cmap[(127 - 3)/2 - i][1];
    cmap[i + (127 + 1)/2][2] = cmap[(127 - 3)/2 - i][1];
  }
}





typedef struct {
  unsigned char r, g, b;
} Color;

Color colormap_color (double cmap[127][3],
        double val, double min, double max)
{
  Color c;
  if (val == 1e30) {
    c.r = c.g = c.b = 0;
    return c;
  }
  int i;
  double coef;
  if (max != min)
    val = (val - min)/(max - min);
  else
    val = 0.;
  if (val <= 0.) i = 0, coef = 0.;
  else if (val >= 1.) i = 127 - 2, coef = 1.;
  else {
    i = val*(127 - 1);
    coef = val*(127 - 1) - i;
  }
  if (!(i >= 0 && i < 127 - 1)) qassert ("/home/spencer/basilisk/src/output.h", 0, "i >= 0 && i < NCMAP - 1");
  unsigned char * c1 = (unsigned char *) &c;
  for (int j = 0; j < 3; j++)
    c1[j] = 255*(cmap[i][j]*(1. - coef) + cmap[i + 1][j]*coef);
  return c;
}
# 300 "/home/spencer/basilisk/src/output.h"
static const char * extension (const char * file, const char * ext) {
  int len = strlen(file);
  return len > 4 && !strcmp (file + len - 4, ext) ? file + len - 4 : NULL;
}

static const char * is_animation (const char * file) {
  const char * ext;
  if ((ext = extension (file, ".mp4")) ||
      (ext = extension (file, ".ogv")) ||
      (ext = extension (file, ".gif")))
    return ext;
  return NULL;
}

static struct {
  FILE ** fp;
  char ** names;
  int n;
} open_image_data = {NULL, NULL, 0};

static void open_image_cleanup()
{
  for (int i = 0; i < open_image_data.n; i++) {
    qpclose (open_image_data.fp[i]);
    pfree (open_image_data.names[i],__func__,__FILE__,0);
  }
  pfree (open_image_data.fp,__func__,__FILE__,0);
  pfree (open_image_data.names,__func__,__FILE__,0);
  open_image_data.fp = NULL;
  open_image_data.names = NULL;
  open_image_data.n = 0;
}

static FILE * open_image_lookup (const char * file)
{
  for (int i = 0; i < open_image_data.n; i++)
    if (!strcmp (file, open_image_data.names[i]))
      return open_image_data.fp[i];
  return NULL;
}

static bool which (const char * command)
{
  char * s = getenv ("PATH");
  if (!s)
    return false;
  char path[strlen(s) + 1];
  strcpy (path, s);
  s = strtok (path, ":");
  while (s) {
    char f[strlen(s) + strlen(command) + 2];
    strcpy (f, s);
    strcat (f, "/");
    strcat (f, command);
    FILE * fp = fopen (f, "r");
    if (fp) {
      fclose (fp);
      return true;
    }
    s = strtok (NULL, ":");
  }
  return false;
}

static FILE * ppm_fallback (const char * file, const char * mode)
{
  char filename[strlen(file) + 5];
  strcpy (filename, file);
  strcat (filename, ".ppm");
  FILE * fp = fopen (filename, mode);
  if (!fp) {
    perror (file);



    exit (1);
  }
  return fp;
}

FILE * open_image (const char * file, const char * options)
{
  if (!(pid() == 0)) qassert ("/home/spencer/basilisk/src/output.h", 0, "pid() == 0");
  const char * ext;
  if ((ext = is_animation (file))) {
    FILE * fp = open_image_lookup (file);
    if (fp)
      return fp;

    int len = strlen ("ppm2???    ") + strlen (file) +
      (options ? strlen (options) : 0);
    char command[len];
    strcpy (command, "ppm2"); strcat (command, ext + 1);

    static int has_ffmpeg = -1;
    if (has_ffmpeg < 0) {
      if (which (command) && (which ("ffmpeg") || which ("avconv")))
 has_ffmpeg = true;
      else {
 fprintf (ferr,
   "open_image(): cannot find '%s' or 'ffmpeg'/'avconv'\n"
   "  falling back to raw PPM outputs\n", command);
 has_ffmpeg = false;
      }
    }
    if (!has_ffmpeg)
      return ppm_fallback (file, "a");

    static bool added = false;
    if (!added) {
      free_solver_func_add (open_image_cleanup);
      added = true;
    }
    open_image_data.n++;
    open_image_data.names = (char * *) prealloc (open_image_data.names, (open_image_data.n)*sizeof(char *),__func__,__FILE__,0);
    open_image_data.names[open_image_data.n - 1] = pstrdup (file,__func__,__FILE__,0);

    if (options) {
      strcat (command, " ");
      strcat (command, options);
    }
    strcat (command, !strcmp (ext, ".mp4") ? " " : " > ");
    strcat (command, file);
    open_image_data.fp = (FILE * *) prealloc (open_image_data.fp, (open_image_data.n)*sizeof(FILE *),__func__,__FILE__,0);
    return open_image_data.fp[open_image_data.n - 1] = qpopen (command, "w");
  }
  else {
    static int has_convert = -1;
    if (has_convert < 0) {
      if (which ("convert"))
 has_convert = true;
      else {
 fprintf (ferr,
   "open_image(): cannot find 'convert'\n"
   "  falling back to raw PPM outputs\n");
 has_convert = false;
      }
    }
    if (!has_convert)
      return ppm_fallback (file, "w");

    int len = strlen ("convert ppm:-   ") + strlen (file) +
      (options ? strlen (options) : 0);
    char command[len];
    strcpy (command, "convert ppm:- ");
    if (options) {
      strcat (command, options);
      strcat (command, " ");
    }
    strcat (command, file);
    return qpopen (command, "w");
  }
}

void close_image (const char * file, FILE * fp)
{
  if (!(pid() == 0)) qassert ("/home/spencer/basilisk/src/output.h", 0, "pid() == 0");
  if (is_animation (file)) {
    if (!open_image_lookup (file))
      fclose (fp);
  }
  else if (which ("convert"))
    qpclose (fp);
  else
    fclose (fp);
}
# 534 "/home/spencer/basilisk/src/output.h"
     
void output_ppm (scalar f,
   FILE * fp,
   int n,
   char * file,
   double min, double max, double spread,
   double z,
   bool linear,
   coord box[2],
   scalar mask,
   Colormap map,
   char * opt)
{tracing("output_ppm","/home/spencer/basilisk/src/output.h",0);

  if (!min && !max) {
    stats s = statsf (f);
    if (spread < 0.)
      min = s.min, max = s.max;
    else {
      double avg = s.sum/s.volume;
      min = avg - spread*s.stddev; max = avg + spread*s.stddev;
    }
  }
  box[0].z = z, box[1].z = z;

  coord cn = {n}, p;
  double delta = (box[1].x - box[0].x)/n;
  cn.y = (int)((box[1].y - box[0].y)/delta);
  if (((int)cn.y) % 2) cn.y++;

  Color ** ppm = (Color **) matrix_new (cn.y, cn.x, sizeof(Color));
  unsigned char * ppm0 = &ppm[0][0].r;
  int len = 3*cn.x*cn.y;
  memset (ppm0, 0, len*sizeof (unsigned char));
  double cmap[127][3];
  (* map) (cmap);




  foreach_region_stencil (1,{(NonLocal[]){{"max","double",(void *)&max,NULL,0},{"min","double",(void *)&min,NULL,0},{"cmap","double",(void *)cmap,(int[]){1273,0},0},{"ferr","not implemented yet",(void *)ferr,NULL,1},{"ppm","Color",(void *)ppm,NULL,2},{"f","scalar",(void *)&f,NULL,0},{"linear","bool",(void *)&linear,NULL,0},{"mask","scalar",(void *)&mask,NULL,0},{"cn","coord",(void *)&cn,NULL,0},{"box","coord",(void *)box,(int[]){2,0},0},{"p","coord",(void *)&p,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_("\n@ def not_mpi_compatible()\ndo {\n  if (npe() > 1) {\n    fprintf (ferr, \"%s() is not compatible with MPI (yet)\\n\", __func__);\n    exit (1);\n  }\n} while(0)\n@\n@ define system(command) (pid() == 0 ? system(command) : 0)\n@else\n@ define qstderr() stderr\n@ define qstdout() stdout\n@ define ferr stderr\n@ define fout stdout\n@ define not_mpi_compatible()\n@endif\n\n\n\n\n// #line 93 \"/home/spencer/basilisk/src/common.h\"\nstatic inline void qassert (const char * file, int line, const char * cond) {\n  fprintf (ferr, \"%s:%d: Assertion `%s' failed.\\n\", file, line, cond);\n  abort();\n}\n\n\n// #line 261 \"/home/spencer/basilisk/src/output.h\"\nColor colormap_color (real cmap[127][3],\n        real val, real min, real max)\n{\n  Color c;\n  if (val == 1e30) {\n    c.r = c.g = c.b = 0;\n    return c;\n  }\n  int i;\n  real coef;\n  if (max != min)\n    val = (val - min)/(max - min);\n  else\n    val = 0.;\n  if (val <= 0.) i = 0, coef = 0.;\n  else if (val >= 1.) i = 127 - 2, coef = 1.;\n  else {\n    i = val*(127 - 1);\n    coef = val*(127 - 1) - i;\n  }\n  if (!(i >= 0 && i < 127 - 1)) qassert (\"/home/spencer/basilisk/src/output.h\", 0, \"i >= 0 && i < NCMAP - 1\");\n  unsigned char * c1 = (unsigned char *) &c;\n  for (int j = 0; j < 3; j++)\n    c1[j] = 255*(cmap[i][j]*(1. - coef) + cmap[i + 1][j]*coef);\n  return c;\n}\n\n\n// #line 805 \"/home/spencer/basilisk/src/grid/cartesian-common.h\"\nstatic real interpolate_linear (Point point, scalar v,\n      real xp, real yp, real zp)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n\n\n\n\n\n\n\n  x = (xp - x)/Delta - v.d.x/2.;\n  y = (yp - y)/Delta - v.d.y/2.;\n  int i = sign(x), j = sign(y);\n  x = fabs(x); y = fabs(y);\n\n  return ((val(v,0,0,0)*(1. - x) + val(v,i,0,0)*x)*(1. - y) +\n   (val(v,0,j,0)*(1. - x) + val(v,i,j,0)*x)*y);\n// # 834 \"/home/spencer/basilisk/src/grid/cartesian-common.h\"\n}"),_("\n\n  \n// #line 576 \"/home/spencer/basilisk/src/output.h\"\n{\n    real v;\n    if (mask.i >= 0) {\n      if (linear) {\n real m = interpolate_linear (point, mask, p.x, p.y, p.z);\n if (m < 0.)\n   v = 1e30;\n else\n   v = interpolate_linear (point, f, p.x, p.y, p.z);\n      }\n      else {\n if (val(mask,0,0,0) < 0.)\n   v = 1e30;\n else\n   v = val(f,0,0,0);\n      }\n    }\n    else if (linear)\n      v = interpolate_linear (point, f, p.x, p.y, p.z);\n    else\n      v = val(f,0,0,0);\n    int i = (p.x - box[0].x)/(box[1].x - box[0].x)*cn.x;\n    int j = (p.y - box[0].y)/(box[1].y - box[0].y)*cn.y;\n    Color ** alias = ppm;\n    alias[(int)cn.y - 1 - j][i] = colormap_color (cmap, v, min, max);\n  }")})

  { 
    
    if (mask.i >= 0) {
      if (linear) {  
  _stencil_interpolate_linear (point, mask, NULL, NULL, NULL);
{ 
    
   
{ _stencil_interpolate_linear (point, f, NULL, NULL, NULL); }}
    
 
      
#line 585
}
      else {
_stencil_val(mask,0,0,0);{
     
   
{ _stencil_val(f,0,0,0); }}
    
 
      
#line 591
}
    }
    else if (linear)
      { _stencil_interpolate_linear (point, f, NULL, NULL, NULL); }
    else
      { _stencil_val(f,0,0,0); }                  
    
    
         
         
  }end_foreach_region_stencil();




  {
#line 574
foreach_region (p, box, cn)

  {
    double v;
    if (mask.i >= 0) {
      if (linear) {
 double m = interpolate_linear (point, mask, p.x, p.y, p.z);
 if (m < 0.)
   v = 1e30;
 else
   v = interpolate_linear (point, f, p.x, p.y, p.z);
      }
      else {
 if (val(mask,0,0,0) < 0.)
   v = 1e30;
 else
   v = val(f,0,0,0);
      }
    }
    else if (linear)
      v = interpolate_linear (point, f, p.x, p.y, p.z);
    else
      v = val(f,0,0,0);
    int i = (p.x - box[0].x)/(box[1].x - box[0].x)*cn.x;
    int j = (p.y - box[0].y)/(box[1].y - box[0].y)*cn.y;
    Color ** alias = ppm;
    alias[(int)cn.y - 1 - j][i] = colormap_color (cmap, v, min, max);
  }end_foreach_region();}

  if (pid() == 0) {
    if (file)
      fp = open_image (file, opt);

    fprintf (fp, "P6\n%g %g 255\n", cn.x, cn.y);
    fwrite (ppm0, sizeof(unsigned char), 3*cn.x*cn.y, fp);

    if (file)
      close_image (file, fp);
    else
      fflush (fp);
  }

  matrix_free (ppm);
end_tracing("output_ppm","/home/spencer/basilisk/src/output.h",0);}
# 649 "/home/spencer/basilisk/src/output.h"
     
void output_grd (scalar f,
   FILE * fp,
   double Delta,
   bool linear,
   double box[2][2],
   scalar mask)
{tracing("output_grd","/home/spencer/basilisk/src/output.h",0);
  int nx = (box[1][0] - box[0][0])/Delta;
  int ny = (box[1][1] - box[0][1])/Delta;


  fprintf (fp, "ncols          %d\n", nx);
  fprintf (fp, "nrows          %d\n", ny);
  fprintf (fp, "xllcorner      %g\n", box[0][0]);
  fprintf (fp, "yllcorner      %g\n", box[0][1]);
  fprintf (fp, "cellsize       %g\n", Delta);
  fprintf (fp, "nodata_value   -9999\n");


  for (int j = ny-1; j >= 0; j--) {
    double yp = Delta*j + box[0][1] + Delta/2.;
    for (int i = 0; i < nx; i++) {
      double xp = Delta*i + box[0][0] + Delta/2., v;
      if (mask.i >= 0) {
 double m = interpolate (mask, xp, yp
#line 837 "/home/spencer/basilisk/src/grid/cartesian-common.h"
, 0.
#line 674 "/home/spencer/basilisk/src/output.h"
, linear);
 if (m < 0.)
   v = 1e30;
 else
   v = interpolate (f, xp, yp
#line 837 "/home/spencer/basilisk/src/grid/cartesian-common.h"
, 0.
#line 678 "/home/spencer/basilisk/src/output.h"
, linear);
      }
      else
 v = interpolate (f, xp, yp
#line 837 "/home/spencer/basilisk/src/grid/cartesian-common.h"
, 0.
#line 681 "/home/spencer/basilisk/src/output.h"
, linear);
      if (v == 1e30)
 fprintf (fp, "-9999 ");
      else
 fprintf (fp, "%f ", v);
    }
    fprintf (fp, "\n");
  }

  fflush (fp);
end_tracing("output_grd","/home/spencer/basilisk/src/output.h",0);}
# 718 "/home/spencer/basilisk/src/output.h"
static char * replace (const char * input, int target, int with,
         bool translate)
{
  if (translate) {
    if (!strcmp (input, "u.x"))
      return pstrdup ("U",__func__,__FILE__,0);
    if (!strcmp (input, "u.y"))
      return pstrdup ("V",__func__,__FILE__,0);
    if (!strcmp (input, "u.z"))
      return pstrdup ("W",__func__,__FILE__,0);
  }
  char * name = pstrdup (input,__func__,__FILE__,0), * i = name;
  while (*i != '\0') {
    if (*i == target)
      *i = with;
    i++;
  }
  return name;
}

     
void output_gfs (FILE * fp,
   scalar * list,
   char * file,
   bool translate)
{tracing("output_gfs","/home/spencer/basilisk/src/output.h",0);
  char * fname = file;

@if _MPI



  FILE * sfp = fp;
  if (file == NULL) {
    long pid = getpid();
    MPI_Bcast (&pid, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    fname = ((char *) pmalloc ((80)*sizeof(char),__func__,__FILE__,0));
    snprintf (fname, 80, ".output-%ld", pid);
    fp = NULL;
  }
@endif

  bool opened = false;
  if (fp == NULL) {
    if (fname == NULL)
      fp = fout;
    else if (!(fp = fopen (fname, "w"))) {
      perror (fname);
      exit (1);
    }
    else
      opened = true;
  }

  scalar * slist = list ? list : list_copy (all);

  restriction (slist);
  fprintf (fp,
    "1 0 GfsSimulation GfsBox GfsGEdge { binary = 1"
    " x = %g y = %g ",
    0.5 + X0/L0, 0.5 + Y0/L0);




  if (slist != NULL && slist[0].i != -1) {
    scalar s = slist[0];
    char * name = replace (_attribute[s.i].name, '.', '_', translate);
    fprintf (fp, "variables = %s", name);
    pfree (name,__func__,__FILE__,0);
    for (int i = 1; i < list_len(slist); i++) {
      scalar s = slist[i];
      if (_attribute[s.i].name) {
 char * name = replace (_attribute[s.i].name, '.', '_', translate);
 fprintf (fp, ",%s", name);
 pfree (name,__func__,__FILE__,0);
      }
    }
    fprintf (fp, " ");
  }
  fprintf (fp, "} {\n");
  fprintf (fp, "  Time { t = %g }\n", t);
  if (L0 != 1.)
    fprintf (fp, "  PhysicalParams { L = %g }\n", L0);
  fprintf (fp, "  VariableTracerVOF f\n");
  fprintf (fp, "}\nGfsBox { x = 0 y = 0 z = 0 } {\n");

@if _MPI
  long header;
  if ((header = ftell (fp)) < 0) {
    perror ("output_gfs(): error in header");
    exit (1);
  }
  int cell_size = sizeof(unsigned) + sizeof(double);
  {scalar*_i=(scalar*)( slist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (_attribute[s.i].name)
      cell_size += sizeof(double);}}
  scalar index = new_scalar("index");
  size_t total_size = header + (z_indexing (index, false) + 1)*cell_size;
@endif



  {foreach_cell() {
@if _MPI
    if (is_local(cell))
@endif
    {
@if _MPI
      if (fseek (fp, header + val(index,0,0,0)*cell_size, SEEK_SET) < 0) {
 perror ("output_gfs(): error while seeking");
 exit (1);
      }
@endif
      unsigned flags =
 level == 0 ? 0 :



      child.x == -1 && child.y == -1 ? 0 :
 child.x == -1 && child.y == 1 ? 1 :
 child.x == 1 && child.y == -1 ? 2 :
 3;
# 851 "/home/spencer/basilisk/src/output.h"
      if (is_leaf(cell))
 flags |= (1 << 4);
      fwrite (&flags, sizeof (unsigned), 1, fp);
      double a = -1;
      fwrite (&a, sizeof (double), 1, fp);
      {scalar*_i=(scalar*)( slist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 if (_attribute[s.i].name) {
   if (_attribute[s.i].v.x.i >= 0) {




     if (_attribute[s.i].v.x.i == s.i) {
       s = _attribute[s.i].v.y;
       a = is_local(cell) && val(s,0,0,0) != 1e30 ? val(s,0,0,0) : (double) DBL_MAX;
     }
     else if (_attribute[s.i].v.y.i == s.i) {
       s = _attribute[s.i].v.x;
       a = is_local(cell) && val(s,0,0,0) != 1e30 ? - val(s,0,0,0) : (double) DBL_MAX;
     }





   }
   else
     a = is_local(cell) && val(s,0,0,0) != 1e30 ? val(s,0,0,0) : (double) DBL_MAX;
   fwrite (&a, sizeof (double), 1, fp);
 }}}
    }
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}

@if _MPI
  delete (((scalar[]){index,{-1}}));
  if (!pid() && fseek (fp, total_size, SEEK_SET) < 0) {
    perror ("output_gfs(): error while finishing");
    exit (1);
  }
  if (!pid())
@endif
    fputs ("}\n", fp);
  fflush (fp);

  if (!list)
    pfree (slist,__func__,__FILE__,0);
  if (opened)
    fclose (fp);

@if _MPI
  if (file == NULL) {
    MPI_Barrier (MPI_COMM_WORLD);
    if (pid() == 0) {
      if (sfp == NULL)
 sfp = fout;
      fp = fopen (fname, "r");
      size_t l;
      unsigned char buffer[8192];
      while ((l = fread (buffer, 1, 8192, fp)) > 0)
 fwrite (buffer, 1, l, sfp);
      fflush (sfp);
      remove (fname);
    }
    pfree (fname,__func__,__FILE__,0);
  }
@endif
end_tracing("output_gfs","/home/spencer/basilisk/src/output.h",0);}
# 943 "/home/spencer/basilisk/src/output.h"
struct DumpHeader {
  double t;
  long len;
  int i, depth, npe, version;
  coord n;
};

static const int dump_version =

  170901;

static scalar * dump_list (scalar * lista)
{
  scalar * list = is_constant(cm) ? NULL : list_concat (((scalar[]){cm,{-1}}), NULL);
  {scalar*_i=(scalar*)( lista);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!_attribute[s.i].face && !_attribute[s.i].nodump && s.i != cm.i)
      list = list_add (list, s);}}
  return list;
}

static void dump_header (FILE * fp, struct DumpHeader * header, scalar * list)
{
  if (fwrite (header, sizeof(struct DumpHeader), 1, fp) < 1) {
    perror ("dump(): error while writing header");
    exit (1);
  }
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    unsigned len = strlen(_attribute[s.i].name);
    if (fwrite (&len, sizeof(unsigned), 1, fp) < 1) {
      perror ("dump(): error while writing len");
      exit (1);
    }
    if (fwrite (_attribute[s.i].name, sizeof(char), len, fp) < len) {
      perror ("dump(): error while writing s.name");
      exit (1);
    }
  }}}
  double o[4] = {X0,Y0,Z0,L0};
  if (fwrite (o, sizeof(double), 4, fp) < 4) {
    perror ("dump(): error while writing coordinates");
    exit (1);
  }
}

@if !_MPI
     
void dump (const char * file,
    scalar * list,
    FILE * fp,
    bool unbuffered)
{tracing("dump","/home/spencer/basilisk/src/output.h",0);
  char * name = NULL;
  if (!fp) {
    name = (char *) pmalloc (strlen(file) + 2,__func__,__FILE__,0);
    strcpy (name, file);
    if (!unbuffered)
      strcat (name, "~");
    if ((fp = fopen (name, "w")) == NULL) {
      perror (name);
      exit (1);
    }
  }
  if (!(fp)) qassert ("/home/spencer/basilisk/src/output.h", 0, "fp");

  scalar * dlist = dump_list (list);
  scalar  size=new_scalar("size");
  scalar * slist = list_concat (((scalar[]){size,{-1}}), dlist); pfree (dlist,__func__,__FILE__,0);
  struct DumpHeader header = { t, list_len(slist), iter, depth(), npe(),
          dump_version };
  dump_header (fp, &header, slist);

  subtree_size (size, false);

  {foreach_cell() {
    unsigned flags = is_leaf(cell) ? leaf : 0;
    if (fwrite (&flags, sizeof(unsigned), 1, fp) < 1) {
      perror ("dump(): error while writing flags");
      exit (1);
    }
    {scalar*_i=(scalar*)( slist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      if (fwrite (&val(s,0,0,0), sizeof(double), 1, fp) < 1) {
 perror ("dump(): error while writing scalars");
 exit (1);
      }}}
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}

  pfree (slist,__func__,__FILE__,0);
  if (file) {
    fclose (fp);
    if (!unbuffered)
      rename (name, file);
    pfree (name,__func__,__FILE__,0);
  }delete((scalar*)((scalar[]){size,{-1}}));
end_tracing("dump","/home/spencer/basilisk/src/output.h",0);}
@else
     
void dump (const char * file,
    scalar * list,
    FILE * fp,
    bool unbuffered)
{tracing("dump","/home/spencer/basilisk/src/output.h",0);
  if (fp != NULL || file == NULL) {
    fprintf (ferr, "dump(): must specify a file name when using MPI\n");
    exit(1);
  }

  char name[strlen(file) + 2];
  strcpy (name, file);
  if (!unbuffered)
    strcat (name, "~");
  FILE * fh = fopen (name, "w");
  if (fh == NULL) {
    perror (name);
    exit (1);
  }

  scalar * dlist = dump_list (list);
  scalar  size=new_scalar("size");
  scalar * slist = list_concat (((scalar[]){size,{-1}}), dlist); pfree (dlist,__func__,__FILE__,0);
  struct DumpHeader header = { t, list_len(slist), iter, depth(), npe(),
          dump_version };







  if (pid() == 0)
    dump_header (fh, &header, slist);

  scalar index = {-1};

  index = new_scalar("index");
  z_indexing (index, false);
  int cell_size = sizeof(unsigned) + header.len*sizeof(double);
  int sizeofheader = sizeof(header) + 4*sizeof(double);
  {scalar*_i=(scalar*)( slist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    sizeofheader += sizeof(unsigned) + sizeof(char)*strlen(_attribute[s.i].name);}}
  long pos = pid() ? 0 : sizeofheader;

  subtree_size (size, false);

  {foreach_cell() {

    if (is_local(cell)) {
      long offset = sizeofheader + val(index,0,0,0)*cell_size;
      if (pos != offset) {
 fseek (fh, offset, SEEK_SET);
 pos = offset;
      }
      unsigned flags = is_leaf(cell) ? leaf : 0;
      fwrite (&flags, 1, sizeof(unsigned), fh);
      {scalar*_i=(scalar*)( slist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 fwrite (&val(s,0,0,0), 1, sizeof(double), fh);}}
      pos += cell_size;
    }
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}

  delete (((scalar[]){index,{-1}}));

  pfree (slist,__func__,__FILE__,0);
  fclose (fh);
  if (!unbuffered && pid() == 0)
    rename (name, file);delete((scalar*)((scalar[]){size,{-1}}));
end_tracing("dump","/home/spencer/basilisk/src/output.h",0);}
@endif

     
bool restore (const char * file,
       scalar * list,
       FILE * fp)
{tracing("restore","/home/spencer/basilisk/src/output.h",0);
  if (!fp && (fp = fopen (file, "r")) == NULL)
    {end_tracing("restore","/home/spencer/basilisk/src/output.h",0);return false;}
  if (!(fp)) qassert ("/home/spencer/basilisk/src/output.h", 0, "fp");

  struct DumpHeader header = {0};
  if (fread (&header, sizeof(header), 1, fp) < 1) {
    fprintf (ferr, "restore(): error: expecting header\n");
    exit (1);
  }


  init_grid (1);
  {foreach_cell() {
    cell.pid = pid();
    cell.flags |= active;
  }end_foreach_cell();}
  ((Tree *)grid)->dirty = true;
# 1157 "/home/spencer/basilisk/src/output.h"
  bool restore_all = (list == all);
  scalar * slist = dump_list (list ? list : all);
  if (header.version == 161020) {
    if (header.len - 1 != list_len (slist)) {
      fprintf (ferr,
        "restore(): error: the list lengths don't match: "
        "%ld (file) != %d (code)\n",
        header.len - 1, list_len (slist));
      exit (1);
    }
  }
  else {
    if (header.version != dump_version) {
      fprintf (ferr,
        "restore(): error: file version mismatch: "
        "%d (file) != %d (code)\n",
        header.version, dump_version);
      exit (1);
    }

    scalar * input = NULL;
    for (int i = 0; i < header.len; i++) {
      unsigned len;
      if (fread (&len, sizeof(unsigned), 1, fp) < 1) {
 fprintf (ferr, "restore(): error: expecting len\n");
 exit (1);
      }
      char name[len + 1];
      if (fread (name, sizeof(char), len, fp) < 1) {
 fprintf (ferr, "restore(): error: expecting s.name\n");
 exit (1);
      }
      name[len] = '\0';

      if (i > 0) {
 bool found = false;
 {scalar*_i=(scalar*)( slist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
   if (!strcmp (_attribute[s.i].name, name)) {
     input = list_append (input, s);
     found = true; break;
   }}}
 if (!found) {
   if (restore_all) {
     scalar s = new_scalar("s");
     pfree (_attribute[s.i].name,__func__,__FILE__,0);
     _attribute[s.i].name = pstrdup (name,__func__,__FILE__,0);
     input = list_append (input, s);
   }
   else
     input = list_append (input, (scalar){INT_MAX});
 }
      }
    }
    pfree (slist,__func__,__FILE__,0);
    slist = input;

    double o[4];
    if (fread (o, sizeof(double), 4, fp) < 4) {
      fprintf (ferr, "restore(): error: expecting coordinates\n");
      exit (1);
    }
    origin (o[0], o[1], o[2]);
    size (o[3]);
  }
# 1232 "/home/spencer/basilisk/src/output.h"
  scalar * listm = is_constant(cm) ? NULL : (scalar *)((vector[]){fm,{{-1},{-1}}});



  {foreach_cell() {
    unsigned flags;
    if (fread (&flags, sizeof(unsigned), 1, fp) != 1) {
      fprintf (ferr, "restore(): error: expecting 'flags'\n");
      exit (1);
    }

    fseek (fp, sizeof(double), SEEK_CUR);
    {scalar*_i=(scalar*)( slist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
      double val;
      if (fread (&val, sizeof(double), 1, fp) != 1) {
 fprintf (ferr, "restore(): error: expecting a scalar\n");
 exit (1);
      }
      if (s.i != INT_MAX)
 val(s,0,0,0) = val;
    }}}
    if (!(flags & leaf) && is_leaf(cell))
      refine_cell (point, listm, 0, NULL);
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    _attribute[s.i].dirty = true;}}


  scalar * other = NULL;
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!list_lookup (slist, s) && !list_lookup (listm, s))
      other = list_append (other, s);}}
  reset (other, 0.);
  pfree (other,__func__,__FILE__,0);

  pfree (slist,__func__,__FILE__,0);
  if (file)
    fclose (fp);


  while (iter < header.i && events (false))
    iter = inext;
  events (false);
  while (t < header.t && events (false))
    t = tnext;
  t = header.t;
  events (false);

  {end_tracing("restore","/home/spencer/basilisk/src/output.h",0);return true;}
end_tracing("restore","/home/spencer/basilisk/src/output.h",0);}
# 431 "/home/spencer/basilisk/src/utils.h" 2


# 1 "./display.h" 1
# 1 "/home/spencer/basilisk/src/display.h"
# 86 "/home/spencer/basilisk/src/display.h"
# 1 "/home/spencer/basilisk/src/ast/std/netdb.h" 1
@include <netdb.h>
# 87 "/home/spencer/basilisk/src/display.h" 2
# 1 "/home/spencer/basilisk/src/wsServer/include/ws.h" 1
# 28 "/home/spencer/basilisk/src/wsServer/include/ws.h"
# 1 "/home/spencer/basilisk/src/ast/std/stdbool.h" 1
@include <stdbool.h>
# 29 "/home/spencer/basilisk/src/wsServer/include/ws.h" 2
# 1 "/home/spencer/basilisk/src/ast/std/stdint.h" 1
@include <stdint.h>
# 30 "/home/spencer/basilisk/src/wsServer/include/ws.h" 2
# 227 "/home/spencer/basilisk/src/wsServer/include/ws.h"
        extern int ws_socket_open (uint16_t port);
 extern int get_handshake_accept(char *wsKey, unsigned char **dest);
 extern int get_handshake_response(char *hsrequest, char **hsresponse);
 extern char *ws_getaddress(int fd);
        extern int ws_sendframe_init(int fd, ssize_t size, bool broadcast, int type);
        extern int ws_sendframe(
  int fd, const char *msg, ssize_t size, bool broadcast, int type);
 extern int ws_sendframe_txt(int fd, const char *msg, bool broadcast);
 extern int ws_sendframe_bin(int fd, const char *msg, size_t size, bool broadcast);
        extern ssize_t ws_send(int sockfd, const void *buf, size_t len);
 extern int ws_get_state(int fd);
 extern int ws_close_client(int fd);

        struct ws_message {
  int fd;
  char * msg;
  size_t size;
  int type;
 };

        extern struct ws_message * ws_socket_poll (int sock, int timeout);
# 88 "/home/spencer/basilisk/src/display.h" 2
#pragma autolink -L$BASILISK/wsServer -lws

# 1 "./view.h" 1
# 1 "/home/spencer/basilisk/src/view.h"
# 56 "/home/spencer/basilisk/src/view.h"
# 1 "/home/spencer/basilisk/src/gl/framebuffer.h" 1
typedef struct _framebuffer framebuffer;

framebuffer * framebuffer_new (unsigned width, unsigned height);
void framebuffer_destroy (framebuffer * p);
unsigned char * framebuffer_image (framebuffer * p);
float * framebuffer_depth (framebuffer * p);
# 57 "/home/spencer/basilisk/src/view.h" 2
# 1 "/home/spencer/basilisk/src/gl/trackball.h" 1
# 50 "/home/spencer/basilisk/src/gl/trackball.h"
void
gl_trackball(float q[4], float p1x, float p1y, float p2x, float p2y);
# 61 "/home/spencer/basilisk/src/gl/trackball.h"
void
gl_add_quats(float q1[4], float q2[4], float dest[4]);





void
gl_build_rotmatrix(float m[4][4], float q[4]);






void
gl_axis_to_quat(float a[3], float phi, float q[4]);
# 58 "/home/spencer/basilisk/src/view.h" 2
# 1 "/home/spencer/basilisk/src/gl/utils.h" 1
# 1 "./gl/tinygl.h" 1
# 1 "/home/spencer/basilisk/src/gl/tinygl.h"
# 122 "/home/spencer/basilisk/src/gl/tinygl.h"
typedef unsigned int GLenum;
typedef int GLint;
typedef unsigned int GLuint;
typedef float GLfloat;
typedef double GLdouble;
typedef int GLsizei;
typedef unsigned int GLbitfield;
typedef unsigned char GLubyte;

void glBegin (GLenum mode);
void glEnd (void);
void glClear (GLbitfield mask);
void glClearColor (float red, float green, float blue, float alpha);
void glBindTexture (GLenum target, GLuint texture);
void glColor3f (GLfloat red, GLfloat green, GLfloat blue);
void glColorMaterial (GLenum face, GLenum mode);
void glDisable (GLenum cap);
void glEnable (GLenum cap);
void glFinish (void);
void glGetDoublev (GLenum pname, GLdouble * params);
void glGetIntegerv (GLenum pname, GLint * params);
void glHint (GLenum target, GLenum mode);
void glLightfv (GLenum light, GLenum pname, const GLfloat *params);
void glGetLightfv (GLenum light, GLenum pname, GLfloat * params);
void glLightModeli (GLenum pname, GLint param);
void glLineWidth (GLfloat width);
void glPointSize (GLfloat size);
void glNormal3d (GLdouble nx, GLdouble ny, GLdouble nz);
void glOrtho (GLdouble left, GLdouble right, GLdouble bottom, GLdouble top,
       GLdouble nearVal, GLdouble farVal);
void glShadeModel (GLenum mode);
void glTexCoord1d (GLdouble s);
void glTexCoord2f (GLfloat s, GLfloat t);
void glTexImage1D (GLenum target, GLint level, GLint internalFormat,
     GLsizei width, GLint border, GLenum format, GLenum type,
     const void * data);
void glTexParameteri (GLenum target, GLenum pname, GLint param);
void glVertex3d (GLdouble x, GLdouble y, GLdouble z);
void glVertex3f (GLfloat x, GLfloat y, GLfloat z);

GLenum glGetError (void);
void glGetFloatv (GLenum pname, GLfloat * params);
void glMultMatrixf (const GLfloat * m);
void glLoadIdentity (void);
void glScalef (GLfloat x, GLfloat y, GLfloat z);
void glTranslatef (GLfloat x, GLfloat y, GLfloat z);
void glRotatef (GLfloat angle, GLfloat x, GLfloat y, GLfloat z);
void glMatrixMode (GLenum mode);
void glPopMatrix (void);
void glPushMatrix (void);
void glLoadMatrixd (const GLdouble * m);
# 2 "/home/spencer/basilisk/src/gl/utils.h" 2
# 1 "/home/spencer/basilisk/src/ast/std/stdio.h" 1
@include <stdio.h>
# 3 "/home/spencer/basilisk/src/gl/utils.h" 2
# 1 "/home/spencer/basilisk/src/ast/std/stdbool.h" 1
@include <stdbool.h>
# 4 "/home/spencer/basilisk/src/gl/utils.h" 2

void gl_write_image (FILE * fp, const GLubyte * buffer,
       unsigned width, unsigned height, unsigned samples);
void gl_write_image_png (FILE * fp, const GLubyte * buffer,
    unsigned width, unsigned height, unsigned samples);
void init_gl();

void matrix_multiply (float * m, const float * n);
void vector_multiply (float * v, const float * m);

typedef struct {
  float m[16], p[16];
  float n[6][3];
  float d[6];
  unsigned width;
} Frustum;

void gl_get_frustum (Frustum * f);
int sphere_in_frustum (double x, double y, double z, double r, Frustum * f);
float sphere_diameter (double x, double y, double z, double r, Frustum * f);
void gl_check_error();

int polygonize (const double val[8], double isolevel, double triangles[5][3][3]);
void gl_perspective (double fovy, double aspect, double zNear, double zFar);
int gl_project (float objx, float objy, float objz,
  const float modelMatrix[16],
  const float projMatrix[16],
  const GLint viewport[4],
  float *winx, float *winy, float *winz);

# 1 "/home/spencer/basilisk/src/gl/parser.h" 1
typedef struct _Node Node;

struct _Node {
  char type;
  union {
    char * id;
    double (* func) (double);
    double value;
  } d;
  int s;
  Node * e[3];
};

Node * parse_node (char * code);
void free_node (Node * n);
void print_node (Node * n, FILE * fp);
void reset_node_type (Node * n, char type);
# 35 "/home/spencer/basilisk/src/gl/utils.h" 2
# 59 "/home/spencer/basilisk/src/view.h" 2
#pragma autolink -L$BASILISK/gl -lglutils $OPENGLIBS

# 1 "./utils.h" 1
# 62 "/home/spencer/basilisk/src/view.h" 2
# 1 "./input.h" 1
# 1 "/home/spencer/basilisk/src/input.h"
# 16 "/home/spencer/basilisk/src/input.h"
void input_pgm (scalar s, FILE * fp,
  double ox, double oy, double width)
{
  char line[81];
  if (!fgets (line, 81, fp)) {
    fprintf (ferr, "input_pgm: could not read magic number\n");
    exit (1);
  }
  if (strcmp (line, "P2\n") && strcmp (line, "P5\n")) {
    fprintf (ferr, "input_pgm: magic number '%s' does not match PGM\n",
      line);
    exit (1);
  }
  int binary = !strcmp (line, "P5\n");
  if (!fgets (line, 81, fp)) {
    fprintf (ferr, "input_pgm: could not read width and height\n");
    exit (1);
  }
  int W, H;
  while (line[0] == '#' && fgets (line, 81, fp));
  if (line[0] == '#' || sscanf (line, "%d %d", &W, &H) != 2) {
    fprintf (ferr, "input_pgm: could not read width and height\n");
    exit (1);
  }
  if (!fgets (line, 81, fp)) {
    fprintf (ferr, "input_pgm: could not read maxval\n");
    exit (1);
  }
  int maxval;
  if (sscanf (line, "%d", &maxval) != 1) {
    fprintf (ferr, "input_pgm: could not read maxval\n");
    exit (1);
  }
  if (maxval < 256) {
    unsigned char * a = ((unsigned char *) pmalloc ((W*H)*sizeof(unsigned char),__func__,__FILE__,0));
    size_t n = 0;
    if (binary)
      n = fread (a, 1, W*H, fp);
    else {
      int v;
      while (n < W*H && fscanf (fp, "%d ", &v) == 1)
 a[n++] = v;
    }
    if (n != W*H) {
      fprintf (ferr, "input_pgm: read only %ld values\n", n);
      exit (1);
    }
    foreach_stencil(1,{(NonLocal[]){{"maxval","int",(void *)&maxval,NULL,0},{"a","not implemented yet",(void *)a,NULL,1},{"s","scalar",(void *)&s,NULL,0},{"H","int",(void *)&H,NULL,0},{"oy","double",(void *)&oy,NULL,0},{"width","double",(void *)&width,NULL,0},{"W","int",(void *)&W,NULL,0},{"ox","double",(void *)&ox,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_(" \n// #line 63 \"/home/spencer/basilisk/src/input.h\"\n{\n      int i = (x - ox)*W/width, j = (y - oy)*W/width;\n      if (i >= 0 && i < W && j >= 0 && j < H)\n val_out_(s,0,0,0) = 1. - a[(H - 1 - j)*W + i]/(real)maxval;\n      else\n val_out_(s,0,0,0) = 0.;\n    }")}) {          
      
{
 {_stencil_val_a(s,0,0,0);          }
 
{_stencil_val_a(s,0,0,0);  }}
                     
      
    
#line 69
}end_foreach_stencil();
    {
#line 63
foreach() {
      int i = (x - ox)*W/width, j = (y - oy)*W/width;
      if (i >= 0 && i < W && j >= 0 && j < H)
 val(s,0,0,0) = 1. - a[(H - 1 - j)*W + i]/(double)maxval;
      else
 val(s,0,0,0) = 0.;
    }end_foreach();}
    pfree (a,__func__,__FILE__,0);
  }
  else {
    unsigned short * a = ((unsigned short *) pmalloc ((W*H)*sizeof(unsigned short),__func__,__FILE__,0));
    size_t n = 0;
    if (binary)
      n = fread (a, 2, W*H, fp);
    else {
      int v;
      while (n < W*H && fscanf (fp, "%d ", &v) == 1)
 a[n++] = v;
    }
    if (n != W*H) {
      fprintf (ferr, "input_pgm: read only %ld values\n", n);
      exit (1);
    }
    foreach_stencil(1,{(NonLocal[]){{"maxval","int",(void *)&maxval,NULL,0},{"a","not implemented yet",(void *)a,NULL,1},{"s","scalar",(void *)&s,NULL,0},{"H","int",(void *)&H,NULL,0},{"oy","double",(void *)&oy,NULL,0},{"width","double",(void *)&width,NULL,0},{"W","int",(void *)&W,NULL,0},{"ox","double",(void *)&ox,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_(" \n// #line 86 \"/home/spencer/basilisk/src/input.h\"\n{\n      int i = (x - ox)*W/width, j = (y - oy)*W/width;\n      if (i >= 0 && i < W && j >= 0 && j < H)\n val_out_(s,0,0,0) = 1. - a[(H - 1 - j)*W + i]/(real)maxval;\n      else\n val_out_(s,0,0,0) = 0.;\n    }")}) {          
      
{
 {_stencil_val_a(s,0,0,0);          }
 
{_stencil_val_a(s,0,0,0);  }}
                     
      
    
#line 92
}end_foreach_stencil();
    {
#line 86
foreach() {
      int i = (x - ox)*W/width, j = (y - oy)*W/width;
      if (i >= 0 && i < W && j >= 0 && j < H)
 val(s,0,0,0) = 1. - a[(H - 1 - j)*W + i]/(double)maxval;
      else
 val(s,0,0,0) = 0.;
    }end_foreach();}
    pfree (a,__func__,__FILE__,0);
  }
}

static void next_char (FILE * fp, int target)
{
  int c = fgetc(fp), para = 0;
  while (c != EOF && (c != target || para > 0)) {
    if (c == '{') para++;
    if (c == '}') para--;
    c = fgetc(fp);
  }
  if (c != target) {
    fprintf (ferr, "input_gfs(): error: expecting '%c'\n", target);
    exit (1);
  }
}

static int next_string (FILE * fp, const char * target)
{
  int slen = strlen (target), para = 0;
  char s[slen + 1];
  s[slen] = '\0';
  int len = 0, c = fgetc (fp);
  while (c != EOF && len < slen) {
    if (c == '{') para++;
    if (c == '}') para--;
    s[len++] = c;
    c = fgetc (fp);
  }
  while (c != EOF && para >= 0) {
    if (!strcmp (s, target) && para == 0)
      break;
    if (c == '{') para++;
    if (c == '}') para--;
    for (int i = 0; i < slen - 1; i++)
      s[i] = s[i+1];
    s[slen - 1] = c;
    c = fgetc (fp);
  }
  if (strcmp (s, target))
    c = -1;
  return c;
}
# 156 "/home/spencer/basilisk/src/input.h"
     
void input_gfs (FILE * fp,
  scalar * list,
  char * file)
{tracing("input_gfs","/home/spencer/basilisk/src/input.h",0);
  not_mpi_compatible();

  if (file && !(fp = fopen (file, "r"))) {
    perror (file);
    exit (1);
  }

  bool input_all = (list == all);
  if (!list) list = all;


  init_grid (1);


  next_char (fp, '{');

  char * s = ((char *) pmalloc ((1)*sizeof(char),__func__,__FILE__,0));
  int len = 0;
  int c = fgetc(fp);
  while (c != EOF && c != '}') {
    s[len++] = c;
    s = (char *) prealloc (s, (len + 1)*sizeof(char),__func__,__FILE__,0);
    s[len] = '\0';
    c = fgetc(fp);
  }
  if (c != '}') {
    fprintf (ferr, "input_gfs(): error: expecting '}'\n");
    exit (1);
  }

  char * s1 = strstr (s, "variables");
  if (!s1) {
    fprintf (ferr, "input_gfs(): error: expecting 'variables'\n");
    exit (1);
  }

  s1 = strstr (s1, "=");
  if (!s1) {
    fprintf (ferr, "input_gfs(): error: expecting '='\n");
    exit (1);
  }
  s1++;

  while (strchr (" \t", *s1))
    s1++;

  scalar * input = NULL;
  s1 = strtok (s1, ", \t");
  while (s1) {
    char * name = replace (s1, '_', '.', false);
    bool found = false;
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      if (!is_constant(s) && _attribute[s.i].name && !strcmp (_attribute[s.i].name, name)) {
 input = list_append (input, s);
 found = true; break;
      }}}
    if (!found) {
      if (input_all) {
 scalar s = new_scalar("s");
 pfree (_attribute[s.i].name,__func__,__FILE__,0);
 _attribute[s.i].name = pstrdup (name,__func__,__FILE__,0);
 input = list_append (input, s);
      }
      else
 input = list_append (input, (scalar){INT_MAX});
    }
    pfree (name,__func__,__FILE__,0);
    s1 = strtok (NULL, ", \t");
  }
  pfree (s,__func__,__FILE__,0);

  next_char (fp, '{');
  double t1 = 0.;
  if (next_string (fp, "Time") >= 0) {
    next_char (fp, '{');
    next_char (fp, 't');
    next_char (fp, '=');
    if (fscanf (fp, "%lf", &t1) != 1) {
      fprintf (ferr, "input_gfs(): error: expecting 't'\n");
      exit (1);
    }
    next_char (fp, '}');
    next_char (fp, '}');
  }

  if (next_string (fp, "Box") < 0) {
    fprintf (ferr, "input_gfs(): error: expecting 'GfsBox'\n");
    exit (1);
  }

  next_char (fp, '{');
  next_char (fp, '{');
  next_char (fp, '\n');

  scalar * listm =((scalar[]) {cm,fm.x,fm.y,{-1}});
  scalar * listr = !is_constant(cm) ? listm : NULL;
  NOT_UNUSED (listr);

  {foreach_cell() {
    unsigned flags;
    if (fread (&flags, sizeof (unsigned), 1, fp) != 1) {
      fprintf (ferr, "input_gfs(): error: expecting 'flags'\n");
      exit (1);
    }
    if (!(flags & (1 << 4)) && is_leaf(cell))
      refine_cell (point, listr, 0, NULL);
    double a;
    if (fread (&a, sizeof (double), 1, fp) != 1 || a != -1) {
      fprintf (ferr, "input_gfs(): error: expecting '-1'\n");
      exit (1);
    }
    {scalar*_i=(scalar*)( input);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
      if (fread (&a, sizeof (double), 1, fp) != 1) {
 fprintf (ferr, "input_gfs(): error: expecting a scalar\n");
 exit (1);
      }
      if (s.i != INT_MAX) {
 if (_attribute[s.i].v.x.i >= 0) {



   if (_attribute[s.i].v.x.i == s.i) {
     s = _attribute[s.i].v.y;
     val(s,0,0,0) = a;
   }
   else if (_attribute[s.i].v.y.i == s.i) {
     s = _attribute[s.i].v.x;
     val(s,0,0,0) = - a;
   }





 }
 else
   val(s,0,0,0) = a;
      }
    }}}
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}
  {scalar*_i=(scalar*)( listm);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant(s))
      _attribute[s.i].dirty = true;}}
  {scalar*_i=(scalar*)( input);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant(s))
      _attribute[s.i].dirty = true;}}

  pfree (input,__func__,__FILE__,0);
  if (file)
    fclose (fp);


  while (t < t1 && events (false))
    t = tnext;
  events (false);
end_tracing("input_gfs","/home/spencer/basilisk/src/input.h",0);}
# 357 "/home/spencer/basilisk/src/input.h"
void input_grd (scalar s,
  FILE * fp, const char * file,
  double nodatavalue,
  bool linear, bool periodic, bool zero,
  int smooth)
{
  scalar input = s;

  if (file && !(fp = fopen (file, "r"))) {
    perror (file);
    exit (1);
  }


  double DeltaGRD;
  int nx, ny;
  double XG0, YG0, ndv;


  char waste[100];
  if (fscanf (fp, "%s %d", waste, &nx) != 2) {
    fprintf (ferr, "input_grd(): error reading 'nx'\n");
    if (file) fclose (fp);
    return;
  }
  if (fscanf (fp, "%s %d", waste, &ny) != 2) {
    fprintf (ferr, "input_grd(): error reading 'ny'\n");
    if (file) fclose (fp);
    return;
  }
  if (fscanf (fp, "%s %lf", waste, &XG0) != 2) {
    fprintf (ferr, "input_grd(): error reading 'XG0'\n");
    if (file) fclose (fp);
    return;
  }
  if (fscanf (fp, "%s %lf", waste, &YG0) != 2) {
    fprintf (ferr, "input_grd(): error reading 'YG0'\n");
    if (file) fclose (fp);
    return;
  }
  if (fscanf (fp, "%s %lf", waste, &DeltaGRD) != 2) {
    fprintf (ferr, "input_grd(): error reading 'DeltaGRD'\n");
    if (file) fclose (fp);
    return;
  }
  if (fscanf (fp, "%s %lf", waste, &ndv) != 2) {
    fprintf (ferr, "input_grd(): error reading 'ndv'\n");
    if (file) fclose (fp);
    return;
  }


  if (!nodatavalue)
    nodatavalue = ndv;


  double * value = ((double *) pmalloc ((nx*ny)*sizeof(double),__func__,__FILE__,0));
  for (int i = ny - 1; i >= 0; i--)
    for (int j = 0 ; j < nx; j++) {
      if (fscanf (fp, "%lf ", &value[j + i*nx]) != 1) {
 fprintf (ferr, "input_grd(): error reading value %d,%d\n", i, j);
 if (file) fclose (fp);
 pfree (value,__func__,__FILE__,0);
 return;
      }
      if (zero && value[j + i*nx] == ndv)
 value[j + i*nx] = 0.;
    }


  if (smooth > 0) {
    double * smoothed = ((double *) pmalloc ((nx*ny)*sizeof(double),__func__,__FILE__,0));
    for (int s = 0; s < smooth; s++) {
      for (int i = 0; i < ny; i++)
 for (int j = 0 ; j < nx; j++) {
   int n = 0;
   smoothed[j + i*nx] = 0.;
   for (int k = -1; k <= 1; k++)
     for (int l = -1; l <= 1; l++)
       if ((l != 0 || k != 0) &&
    i + k >= 0 && i + k < ny &&
    j + l >= 0 && j + l < nx &&
    value[j + l + (i + k)*nx] != ndv)
  smoothed[j + i*nx] += value[j + l + (i + k)*nx], n++;
   if (n == 0)
     smoothed[j + i*nx] = zero ? 0. : ndv;
   else
     smoothed[j + i*nx] /= n;
 }
      do { double * __tmp = value; value = smoothed; smoothed = __tmp; } while(0);
    }
    pfree (smoothed,__func__,__FILE__,0);
  }

  bool warning = false;
  foreach_stencil (0,{(NonLocal[]){{"warning","bool",(void *)&warning,NULL,0},{"ndv","double",(void *)&ndv,NULL,0},{"value","double",(void *)value,NULL,1},{"linear","bool",(void *)&linear,NULL,0},{"ny","int",(void *)&ny,NULL,0},{"YG0","double",(void *)&YG0,NULL,0},{"DeltaGRD","double",(void *)&DeltaGRD,NULL,0},{"nx","int",(void *)&nx,NULL,0},{"XG0","double",(void *)&XG0,NULL,0},{"right","ast_int",(void *)&right,NULL,0},{"input","scalar",(void *)&input,NULL,0},{"periodic","bool",(void *)&periodic,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},NULL}) {                   
     
  
           
           
  
     

    
    
{ {              

 
  
          
                 
           
                  
                   
         
                 
      

      
       





{
 {_stencil_val_a(input,0,0,0);  }
 
{_stencil_val_a(input,0,0,0);  }}
                               
                  
              
      
     
         
      
    
#line 486
} 
{
      _stencil_val_a(input,0,0,0);    
      
    }}
                   
    
  
#line 491
}end_foreach_stencil();
  
#if _OPENMP
  #undef OMP
  #define OMP(x)
#endif
{
#line 452
foreach () {
    if (periodic || _attribute[input.i].boundary[right] == periodic_bc) {
      if (x > XG0 + nx*DeltaGRD)
 x -= nx*DeltaGRD;
      else if (x < XG0)
 x += nx*DeltaGRD;
    }

    int j = (x - XG0 + DeltaGRD/2.)/DeltaGRD;
    int i = (y - YG0 + DeltaGRD/2.)/DeltaGRD;
    if (i >= 0 && i < ny && j >= 0 && j < nx) {
      double val;

      int j1 = (x - XG0)/DeltaGRD;
      int i1 = (y - YG0)/DeltaGRD;
      if (linear && i1 >= 0 && j1 >= 0 && i1 < ny - 1 && j1 < nx - 1 &&
   value[j1 + i1*nx] != ndv && value[j1 + 1 + i1*nx] != ndv &&
   value[j1 + (i1 + 1)*nx] != ndv && value[j1 + 1 + (i1 + 1)*nx] != ndv) {

 double dx = x - (j1*DeltaGRD + XG0);
 double dy = y - (i1*DeltaGRD + YG0);
 val = (value[j1 + i1*nx] +
        dx*(value[j1 + 1 + i1*nx] - value[j1 + i1*nx])/DeltaGRD +
        dy*(value[j1 + (i1 + 1)*nx] - value[j1 + i1*nx])/DeltaGRD +
        dx*dy*(value[j1 + i1*nx] + value[j1 + 1 + (i1 + 1)*nx] -
        value[j1 + (i1 + 1)*nx] - value[j1 + 1 + i1*nx])
        /sq(DeltaGRD));
      }
      else
 val = value[j + i*nx];
      if (val == ndv)
 val(input,0,0,0) = 1e30;
      else
 val(input,0,0,0) = val;
    }
    else {
      val(input,0,0,0) = 1e30;
      warning = true;
    }
  }end_foreach();}
#if _OPENMP
  #undef OMP
  #define OMP(x) _Pragma(#x)
#endif

  
#line 492
pfree (value,__func__,__FILE__,0);

  if (warning)
    fprintf (ferr,
      "input_grd(): Warning: Raster data is not covering all"
      " the simulation area\n");

  if (file)
    fclose (fp);
}
# 63 "/home/spencer/basilisk/src/view.h" 2







typedef struct {
  char * expr;
  scalar s;
} cexpr;

static scalar get_cexpr (cexpr * cache, const char * expr)
{
  cexpr * c = cache;
  while (c->expr) {
    if (!strcmp (c->expr, expr)) {


      cexpr tmp = *c;
      while ((c + 1)->expr)
 *c = *(c + 1), c++;
      *c = tmp;
      return c->s;
    }
    c++;
  }
  return (scalar){-1};
}

static cexpr * add_cexpr (cexpr * cache, int maxlen,
     const char * expr, scalar s)
{
  cexpr * c = cache;
  while (c->expr) c++;
  int len = c - cache;
  if (len < maxlen) {
    cache = prealloc (cache, sizeof(cexpr)*(len + 2),__func__,__FILE__,0);
    c = &cache[len];
  }
  else {

    c = cache;
    pfree (c->expr,__func__,__FILE__,0);
    scalar s = c->s;
    delete (((scalar[]){s,{-1}}));

    while ((c + 1)->expr)
      *c = *(c + 1), c++;
  }
  c->expr = pstrdup (expr,__func__,__FILE__,0);
  c->s = s;
  (c + 1)->expr = NULL;
  return cache;
}

static void free_cexpr (cexpr * cache)
{
  cexpr * c = cache;
  while (c->expr) {
    pfree (c->expr,__func__,__FILE__,0);
    scalar s = c->s;
    delete (((scalar[]){s,{-1}}));
    c++;
  }
  pfree (cache,__func__,__FILE__,0);
}






typedef void (* MapFunc) (coord *);

struct _bview {
  float tx, ty, sx, sy, sz;
  float quat[4];
  float fov;
  float tz, near, far;

  bool gfsview;
  bool reversed;

  float bg[3];
  float lc;
  float res;

  unsigned width, height, samples;

  framebuffer * fb;
  Frustum frustum;

  MapFunc map;

  int ni;

  bool active;

  cexpr * cache;
  int maxlen;
};

typedef struct _bview bview;




bview * bview_new()
{
  bview * p = ((bview *) pcalloc (1, sizeof(bview),__func__,__FILE__,0));

  p->tx = p->ty = 0;
  p->sx = p->sy = p->sz = 1.;
  p->quat[0] = p->quat[1] = p->quat[2] = 0; p->quat[3] = 1;
  p->fov = 24.;
  gl_trackball (p->quat, 0.0, 0.0, 0.0, 0.0);


  p->bg[0] = 1; p->bg[1] = 1; p->bg[2] = 1;



  p->res = 1.;
  p->lc = 0.004;

  p->samples = 4;
  p->width = 600*p->samples, p->height = 600*p->samples;


  disable_fpe (FE_DIVBYZERO|FE_INVALID);

  p->fb = framebuffer_new (p->width, p->height);

  init_gl();
  p->active = false;

  enable_fpe (FE_DIVBYZERO|FE_INVALID);

  return p;
}




void bview_destroy (bview * p)
{
  framebuffer_destroy (p->fb);
  if (p->cache)
    free_cexpr (p->cache);
  pfree (p,__func__,__FILE__,0);
}




static bview * _view = NULL;






static void destroy_view()
{
  if (!(_view)) qassert ("/home/spencer/basilisk/src/view.h", 0, "_view");
  bview_destroy (_view);
}

bview * get_view() {
  if (!_view) {
    _view = bview_new();
    free_solver_func_add (destroy_view);
  }
  return _view;
}




static void redraw() {
  bview * view = get_view();


  disable_fpe (FE_DIVBYZERO|FE_INVALID);

  glMatrixMode (0x1701);
  glLoadIdentity ();

  if (view->far <= view->near) {
    double max = 2.;
    gl_perspective (view->fov, view->width/(float)view->height, 1., 1. + 2.*max);

    glMatrixMode (0x1700);
    glLoadIdentity ();
    glTranslatef (view->tx, view->ty, - (1. + max));
  }
  else {
    gl_perspective (view->fov, view->width/(float)view->height,
      view->near, view->far);

    glMatrixMode (0x1700);
    glLoadIdentity ();
    glTranslatef (view->tx, view->ty, view->tz);
  }

  GLfloat m[4][4];
  gl_build_rotmatrix (m, view->quat);
  glMultMatrixf (&m[0][0]);

  if (view->gfsview) {
    m[0][0] = 0., m[0][1] = 0., m[0][2] = -1.;
    m[1][0] = 0., m[1][1] = -1., m[1][2] = 0.;
    m[2][0] = 1., m[2][1] = 0., m[2][2] = 0.;
    glMultMatrixf (&m[0][0]);
  }

  glScalef (view->sx/L0, view->sy/L0, view->sz/L0);

  glClearColor (view->bg[0], view->bg[1], view->bg[2], 0.);
  glClear (0x00004000|0x00000100);

  gl_get_frustum (&view->frustum);

  view->active = true;
  view->ni = 0;
}




bview * draw() {
  bview * view = get_view();
  if (!view->active)
    redraw();
  else


    disable_fpe (FE_DIVBYZERO|FE_INVALID);
  glMatrixMode (0x1701);
  glTranslatef (0, 0, - 1e-4);
  return view;
}







typedef void * pointer;


     
static pointer compose_image (bview * view) {tracing("compose_image","/home/spencer/basilisk/src/view.h",0);
  { pointer _ret= framebuffer_image((view)->fb);end_tracing("compose_image","/home/spencer/basilisk/src/view.h",0);return _ret;}
end_tracing("compose_image","/home/spencer/basilisk/src/view.h",0);}
# 414 "/home/spencer/basilisk/src/view.h"
# 1 "./vertexbuffer.h" 1
# 1 "/home/spencer/basilisk/src/vertexbuffer.h"
# 14 "/home/spencer/basilisk/src/vertexbuffer.h"
struct {

  Array * position, * normal, * color, * index;
  float modelview[16];
  int type;
  int dim;
  int vertex, nvertex;
  bool visible;


  int line_loop, lines, line_strip ;
  int quads, polygon, fan;
  int state;
} VertexBuffer = {
  .visible = false,
  .modelview = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 }
};

static void vertex_buffer_push_index (unsigned int i)
{
  i -= VertexBuffer.vertex;
  array_append (VertexBuffer.index, &i, sizeof(unsigned int));
}

void vertex_buffer_setup()
{
  VertexBuffer.nvertex = 0;
  VertexBuffer.type = -1;
  VertexBuffer.dim = -1;
  VertexBuffer.position = array_new();
  VertexBuffer.normal = array_new();
  VertexBuffer.color = array_new();
  VertexBuffer.index = array_new();
}

void vertex_buffer_free()
{
  array_free (VertexBuffer.position);
  VertexBuffer.position = NULL;
  array_free (VertexBuffer.normal);
  VertexBuffer.normal = NULL;
  array_free (VertexBuffer.color);
  VertexBuffer.color = NULL;
  array_free (VertexBuffer.index);
  VertexBuffer.index = NULL;
}

static void vertex_buffer_glBegin (unsigned int state)
{
  if (VertexBuffer.index) {

    glGetFloatv (0x0BA6, VertexBuffer.modelview);

    bview * view = get_view();

    float q[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,
      - view->tx, - view->ty, 3, 1 };
    matrix_multiply (q, VertexBuffer.modelview);
    for (int i = 0; i < 16; i++)
      VertexBuffer.modelview[i] = q[i];

    gl_build_rotmatrix ((float (*)[4])q, view->quat);
    do { float __tmp = q[1]; q[1] = q[4]; q[4] = __tmp; } while(0);
    do { float __tmp = q[2]; q[2] = q[8]; q[8] = __tmp; } while(0);
    do { float __tmp = q[6]; q[6] = q[9]; q[9] = __tmp; } while(0);
    matrix_multiply (q, VertexBuffer.modelview);
    for (int i = 0; i < 16; i++)
      VertexBuffer.modelview[i] = q[i];

    VertexBuffer.state = state;
    switch (state) {
    case 0x0002:
      VertexBuffer.line_loop = VertexBuffer.nvertex;
      break;
    case 0x0001:
      VertexBuffer.lines = VertexBuffer.nvertex;
      break;
    case 0x0003:
      VertexBuffer.line_strip = VertexBuffer.nvertex;
      break;
    case 0x0007:
      VertexBuffer.quads = VertexBuffer.nvertex;
      break;
    case 0x0009:
      VertexBuffer.polygon = VertexBuffer.nvertex;
      break;
    case 0x0006:
      VertexBuffer.fan = VertexBuffer.nvertex;
      break;
    default:
      fprintf (ferr, "glBegin (%d) not implemented yet\n", state);
      break;
    }
  }
  else
    glBegin (state);
}

static void vertex_buffer_glEnd()
{
  if (VertexBuffer.index) {
    int type = -1;
    switch (VertexBuffer.state) {

    case 0x0002:
      for (int i = VertexBuffer.line_loop; i < VertexBuffer.nvertex - 1; i++) {
 vertex_buffer_push_index (i);
 vertex_buffer_push_index (i + 1);
      }
      vertex_buffer_push_index (VertexBuffer.nvertex - 1);
      vertex_buffer_push_index (VertexBuffer.line_loop);
      type = 0;
      break;

    case 0x0001:
      for (int i = VertexBuffer.lines; i < VertexBuffer.nvertex; i += 2) {
 vertex_buffer_push_index (i);
 vertex_buffer_push_index (i + 1);
      }
      type = 0;
      break;

    case 0x0003:
      for (int i = VertexBuffer.line_strip; i < VertexBuffer.nvertex - 1; i++) {
 vertex_buffer_push_index (i);
 vertex_buffer_push_index (i + 1);
      }
      type = 0;
      break;

    case 0x0007:
      for (int i = VertexBuffer.quads; i < VertexBuffer.nvertex; i += 4)
 for (int j = 1; j <= 2; j++) {
   vertex_buffer_push_index (i);
   vertex_buffer_push_index (i + j);
   vertex_buffer_push_index (i + j + 1);
 }
      type = 1;
      break;

    case 0x0009:
      for (int j = 1; j <= VertexBuffer.nvertex - VertexBuffer.polygon - 2;
    j++) {
 vertex_buffer_push_index (VertexBuffer.polygon);
 vertex_buffer_push_index (VertexBuffer.polygon + j);
 vertex_buffer_push_index (VertexBuffer.polygon + j + 1);
      }
      type = 1;
      break;

    case 0x0006:
      for (int i = VertexBuffer.fan + 1; i < VertexBuffer.nvertex - 1; i++) {
 vertex_buffer_push_index (VertexBuffer.fan);
 vertex_buffer_push_index (i);
 vertex_buffer_push_index (i + 1);
      }
      type = 1;
      break;

    default:
      break;
    }
    VertexBuffer.state = 0;
    if (VertexBuffer.type >= 0 && type >= 0) {

      if (!(VertexBuffer.type == type)) qassert ("/home/spencer/basilisk/src/vertexbuffer.h", 0, "VertexBuffer.type == type");
    }
    else
      VertexBuffer.type = type;
  }
  else
    glEnd();
}

static void vertex_buffer_glColor3f (float r, float g, float b)
{
  if (VertexBuffer.color) {
    struct { float x, y, z; } color = {r, g, b};
    array_append (VertexBuffer.color, &color, 3*sizeof(float));
  }
  else
    glColor3f (r, g, b);
}

static void vertex_buffer_glNormal3d (double nx, double ny, double nz)
{
  if (VertexBuffer.normal) {
    struct { float x, y, z; } normal = {nx, ny, nz};
    array_append (VertexBuffer.normal, &normal, 3*sizeof(float));
  }
  else
    glNormal3d (nx, ny, nz);
}

static void vertex_buffer_glVertex3d (double x, double y, double z)
{
  if (VertexBuffer.position) {
    if (VertexBuffer.dim < 3)
      VertexBuffer.dim = 3;
    float v[4] = {x, y, z, 1.};
    vector_multiply (v, VertexBuffer.modelview);
    array_append (VertexBuffer.position, v, 3*sizeof(float));
    VertexBuffer.nvertex++;
  }
  else
    glVertex3d (x, y, z);
}

static void vertex_buffer_glVertex2d (double x, double y)
{
  if (VertexBuffer.position) {
    if (VertexBuffer.dim < 2)
      VertexBuffer.dim = 2;
    float v[4] = {x, y, 0, 1.};
    vector_multiply (v, VertexBuffer.modelview);
    array_append (VertexBuffer.position, v, 3*sizeof(float));
    VertexBuffer.nvertex++;
  }
  else
    glVertex3d (x, y, 0.);
}
# 415 "/home/spencer/basilisk/src/view.h" 2






# 1 "./draw.h" 1
# 1 "/home/spencer/basilisk/src/draw.h"




# 1 "./fractions.h" 1
# 1 "/home/spencer/basilisk/src/fractions.h"
# 12 "/home/spencer/basilisk/src/fractions.h"
# 1 "./geometry.h" 1
# 1 "/home/spencer/basilisk/src/geometry.h"
# 35 "/home/spencer/basilisk/src/geometry.h"
double line_alpha (double c, coord n)
{
  double alpha, n1, n2;

  n1 = fabs (n.x); n2 = fabs (n.y);
  if (n1 > n2)
    do { double __tmp = n1; n1 = n2; n2 = __tmp; } while(0);

  c = clamp (c, 0., 1.);
  double v1 = n1/2.;
  if (c <= v1/n2)
    alpha = sqrt (2.*c*n1*n2);
  else if (c <= 1. - v1/n2)
    alpha = c*n2 + v1;
  else
    alpha = n1 + n2 - sqrt (2.*n1*n2*(1. - c));

  if (n.x < 0.)
    alpha += n.x;
  if (n.y < 0.)
    alpha += n.y;

  return alpha - (n.x + n.y)/2.;
}
# 163 "/home/spencer/basilisk/src/geometry.h"
double line_area (double nx, double ny, double alpha)
{
  double a, v, area;

  alpha += (nx + ny)/2.;
  if (nx < 0.) {
    alpha -= nx;
    nx = - nx;
  }
  if (ny < 0.) {
    alpha -= ny;
    ny = - ny;
  }

  if (alpha <= 0.)
    return 0.;

  if (alpha >= nx + ny)
    return 1.;

  if (nx < 1e-10)
    area = alpha/ny;
  else if (ny < 1e-10)
    area = alpha/nx;
  else {
    v = sq(alpha);

    a = alpha - nx;
    if (a > 0.)
      v -= a*a;

    a = alpha - ny;
    if (a > 0.)
      v -= a*a;

    area = v/(2.*nx*ny);
  }

  return clamp (area, 0., 1.);
}
# 267 "/home/spencer/basilisk/src/geometry.h"
double rectangle_fraction (coord n, double alpha, coord a, coord b)
{
  coord n1;
   {
    alpha -= n.x*(b.x + a.x)/2.;
    n1.x = n.x*(b.x - a.x);
  } 
#line 270
{
    alpha -= n.y*(b.y + a.y)/2.;
    n1.y = n.y*(b.y - a.y);
  }
  return line_area(n1.x, n1.y, alpha);
}
# 292 "/home/spencer/basilisk/src/geometry.h"
int facets (coord n, double alpha, coord p[2])
{
  int i = 0;
  for (double s = -0.5; s <= 0.5; s += 1.)
    {
      if (fabs (n.y) > 1e-4 && i < 2) {
 double a = (alpha - s*n.x)/n.y;
 if (a >= -0.5 && a <= 0.5) {
   p[i].x = s;
   p[i++].y = a;
 }
      }
      
#line 297
if (fabs (n.x) > 1e-4 && i < 2) {
 double a = (alpha - s*n.y)/n.x;
 if (a >= -0.5 && a <= 0.5) {
   p[i].y = s;
   p[i++].x = a;
 }
      }}
  return i;
}
# 382 "/home/spencer/basilisk/src/geometry.h"
double line_length_center (coord m, double alpha, coord * p)
{
  alpha += (m.x + m.y)/2.;

  coord n = m;
  
    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }
    
#line 388
if (n.y < 0.) {
      alpha -= n.y;
      n.y = - n.y;
    }

  p->x = p->y = p->z = 0.;

  if (alpha <= 0. || alpha >= n.x + n.y)
    return 0.;

  
    if (n.x < 1e-4) {
      p->x = 0.;
      p->y = (m.y < 0. ? 1. - alpha : alpha) - 0.5;
      return 1.;
    }
    
#line 399
if (n.y < 1e-4) {
      p->y = 0.;
      p->x = (m.x < 0. ? 1. - alpha : alpha) - 0.5;
      return 1.;
    }

  if (alpha >= n.x) {
    p->x += 1.;
    p->y += (alpha - n.x)/n.y;
  }
  else
    p->x += alpha/n.x;

  double ax = p->x, ay = p->y;
  if (alpha >= n.y) {
    p->y += 1.;
    ay -= 1.;
    p->x += (alpha - n.y)/n.x;
    ax -= (alpha - n.y)/n.x;
  }
  else {
    p->y += alpha/n.y;
    ay -= alpha/n.y;
  }

   {
    p->x /= 2.;
    p->x = clamp (p->x, 0., 1.);
    if (m.x < 0.)
      p->x = 1. - p->x;
    p->x -= 0.5;
  } 
#line 424
{
    p->y /= 2.;
    p->y = clamp (p->y, 0., 1.);
    if (m.y < 0.)
      p->y = 1. - p->y;
    p->y -= 0.5;
  }

  return sqrt (ax*ax + ay*ay);
}
# 512 "/home/spencer/basilisk/src/geometry.h"
void line_center (coord m, double alpha, double a, coord * p)
{
  alpha += (m.x + m.y)/2.;

  coord n = m;
  
    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }
    
#line 518
if (n.y < 0.) {
      alpha -= n.y;
      n.y = - n.y;
    }

  p->z = 0.;
  if (alpha <= 0.) {
    p->x = p->y = -0.5;
    return;
  }

  if (alpha >= n.x + n.y) {
    p->x = p->y = 0.;
    return;
  }

  
    if (n.x < 1e-4) {
      p->x = 0.;
      p->y = sign(m.y)*(a/2. - 0.5);
      return;
    }
    
#line 535
if (n.y < 1e-4) {
      p->y = 0.;
      p->x = sign(m.x)*(a/2. - 0.5);
      return;
    }

  p->x = p->y = cube(alpha);

   {
    double b = alpha - n.x;
    if (b > 0.) {
      p->x -= sq(b)*(alpha + 2.*n.x);
      p->y -= cube(b);
    }
  } 
#line 543
{
    double b = alpha - n.y;
    if (b > 0.) {
      p->y -= sq(b)*(alpha + 2.*n.y);
      p->x -= cube(b);
    }
  }

   {
    p->x /= 6.*sq(n.x)*n.y*a;
    p->x = sign(m.x)*(p->x - 0.5);
  } 
#line 551
{
    p->y /= 6.*sq(n.y)*n.x*a;
    p->y = sign(m.y)*(p->y - 0.5);
  }
}
# 13 "/home/spencer/basilisk/src/fractions.h" 2





# 1 "./myc2d.h" 1
# 1 "/home/spencer/basilisk/src/myc2d.h"





coord mycs (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  int ix;
  double c_t,c_b,c_r,c_l;
  double mx0,my0,mx1,my1,mm1,mm2;


  c_t = val(c,-1,1,0) + val(c,0,1,0) + val(c,1,1,0);
  c_b = val(c,-1,-1,0) + val(c,0,-1,0) + val(c,1,-1,0);
  c_r = val(c,1,-1,0) + val(c,1,0,0) + val(c,1,1,0);
  c_l = val(c,-1,-1,0) + val(c,-1,0,0) + val(c,-1,1,0);



  mx0 = 0.5*(c_l - c_r);
  my0 = 0.5*(c_b - c_t);


  if (fabs(mx0) <= fabs(my0)) {
    my0 = my0 > 0. ? 1. : -1.;
    ix = 1;
  }
  else {
    mx0 = mx0 > 0. ? 1. : -1.;
    ix = 0;
  }


  mm1 = val(c,-1,-1,0) + 2.0*val(c,-1,0,0) + val(c,-1,1,0);
  mm2 = val(c,1,-1,0) + 2.0*val(c,1,0,0) + val(c,1,1,0);
  mx1 = mm1 - mm2 + 1e-30;
  mm1 = val(c,-1,-1,0) + 2.0*val(c,0,-1,0) + val(c,1,-1,0);
  mm2 = val(c,-1,1,0) + 2.0*val(c,0,1,0) + val(c,1,1,0);
  my1 = mm1 - mm2 + 1e-30;


  if (ix) {
    mm1 = fabs(my1);
    mm1 = fabs(mx1)/mm1;
    if (mm1 > fabs(mx0)) {
      mx0 = mx1;
      my0 = my1;
    }
  }
  else {
    mm1 = fabs(mx1);
    mm1 = fabs(my1)/mm1;
    if (mm1 > fabs(my0)) {
      mx0 = mx1;
      my0 = my1;
    }
  }



  mm1 = fabs(mx0) + fabs(my0);
  coord n = {mx0/mm1, my0/mm1};

  return n;
}
# 13 "/home/spencer/basilisk/src/fractions.h" 2





# 1 "./myc2d.h" 1
# 1 "/home/spencer/basilisk/src/myc2d.h"





static void _stencil_mycs (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;   
  
  
   


_stencil_val(c,-1,1,0); _stencil_val(c,0,1,0); _stencil_val(c,1,1,0); 


     
#line 14
_stencil_val(c,-1,-1,0); _stencil_val(c,0,-1,0); _stencil_val(c,1,-1,0); 
     _stencil_val(c,1,-1,0); _stencil_val(c,1,0,0); _stencil_val(c,1,1,0); 
     _stencil_val(c,-1,-1,0); _stencil_val(c,-1,0,0); _stencil_val(c,-1,1,0);
            
      
   
            
     
       
    



     
      





_stencil_val(c,-1,-1,0);_stencil_val(c,-1,0,0); _stencil_val(c,-1,1,0); 


     
  


      
#line 35
_stencil_val(c,1,-1,0);_stencil_val(c,1,0,0); _stencil_val(c,1,1,0);  
     
        _stencil_val(c,-1,-1,0);_stencil_val(c,0,-1,0); _stencil_val(c,1,-1,0);
       _stencil_val(c,-1,1,0);_stencil_val(c,0,1,0); _stencil_val(c,1,1,0);    
        
        
    
      
       
       
   
    
      
       


  
        
        
    
      
       
       
   



      
  

  
#line 64
return ;
}
# 19 "/home/spencer/basilisk/src/fractions.h" 2
# 40 "/home/spencer/basilisk/src/fractions.h"
void fraction_refine (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;





  double cc = val(c,0,0,0);
  if (cc <= 0. || cc >= 1.)
    {foreach_child()
      val(c,0,0,0) = cc;end_foreach_child()}
  else {




    coord n = mycs (point, c);
    double alpha = line_alpha (cc, n);






    {foreach_child() {
      static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
      coord nc;
      
 nc.x = child.x*n.x;
 
#line 68
nc.y = child.y*n.y;
      val(c,0,0,0) = rectangle_fraction (nc, alpha, a, b);
    }end_foreach_child()}
  }
}











static void alpha_refine (Point point, scalar alpha)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  vector n = _attribute[alpha.i].n;
  double alphac = 2.*val(alpha,0,0,0);
  coord m;
  
    m.x = val(n.x,0,0,0);
    
#line 90
m.y = val(n.y,0,0,0);
  {foreach_child() {
    val(alpha,0,0,0) = alphac;
    
      val(alpha,0,0,0) -= child.x*m.x/2.;
      
#line 94
val(alpha,0,0,0) -= child.y*m.y/2.;
  }end_foreach_child()}
}
# 120 "/home/spencer/basilisk/src/fractions.h"
     
void fractions (scalar Phi, scalar c,
  vector s, double val)
{tracing("fractions","/home/spencer/basilisk/src/fractions.h",0);

  vector   as=(s).x.i>0?(s):new_face_vector("as");
# 136 "/home/spencer/basilisk/src/fractions.h"
  vector p;
  p.x = as.y; p.y = as.x;
# 146 "/home/spencer/basilisk/src/fractions.h"
  foreach_face_stencil(1,{(NonLocal[]){{"p","vector",(void *)&p,NULL,0},{"Phi","scalar",(void *)&Phi,NULL,0},{"val","double",(void *)&val,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 146 \"/home/spencer/basilisk/src/fractions.h\"\n{is_face_y(){ {\n\n\n\n\n\n    if ((val(Phi,0,0,0) - val)*(val(Phi,1,0,0) - val) < 0.) {\n\n\n\n\n\n\n      val_out_(p.x,0,0,0) = (val(Phi,0,0,0) - val)/(val(Phi,0,0,0) - val(Phi,1,0,0));\n      if (val(Phi,0,0,0) < val)\n val_out_(p.x,0,0,0) = 1. - val(p.x,0,0,0);\n    }\n// # 171 \"/home/spencer/basilisk/src/fractions.h\"\n    else\n      val_out_(p.x,0,0,0) = (val(Phi,0,0,0) > val || val(Phi,1,0,0) > val);\n  }}end_is_face_y()\n// #line 146\nis_face_x(){ {\n\n\n\n\n\n    if ((val(Phi,0,0,0) - val)*(val(Phi,0,1,0) - val) < 0.) {\n\n\n\n\n\n\n      val_out_(p.y,0,0,0) = (val(Phi,0,0,0) - val)/(val(Phi,0,0,0) - val(Phi,0,1,0));\n      if (val(Phi,0,0,0) < val)\n val_out_(p.y,0,0,0) = 1. - val(p.y,0,0,0);\n    }\n// # 171 \"/home/spencer/basilisk/src/fractions.h\"\n    else\n      val_out_(p.y,0,0,0) = (val(Phi,0,0,0) > val || val(Phi,0,1,0) > val);\n  }}end_is_face_x()}")}){_stencil_is_face_y(){ {





_stencil_val(Phi,0,0,0);_stencil_val(Phi,1,0,0);{ {






_stencil_val(Phi,0,0,0);_stencil_val(Phi,0,0,0); _stencil_val(Phi,1,0,0);






      
#line 159
_stencil_val_a(p.x,0,0,0);
_stencil_val(Phi,0,0,0);
 { _stencil_val(p.x,0,0,0);_stencil_val_a(p.x,0,0,0);   }     
         
    
#line 162
}
      








{_stencil_val(Phi,0,0,0); _stencil_val(Phi,1,0,0);_stencil_val_a(p.x,0,0,0);       }}





           
# 171 "/home/spencer/basilisk/src/fractions.h"
    
  
}}end__stencil_is_face_y()
#line 146
_stencil_is_face_x(){ {





_stencil_val(Phi,0,0,0);_stencil_val(Phi,0,1,0);{ {






_stencil_val(Phi,0,0,0);_stencil_val(Phi,0,0,0); _stencil_val(Phi,0,1,0);






      
#line 159
_stencil_val_a(p.y,0,0,0);
_stencil_val(Phi,0,0,0);
 { _stencil_val(p.y,0,0,0);_stencil_val_a(p.y,0,0,0);   }     
         
    
#line 162
}
      








{_stencil_val(Phi,0,0,0); _stencil_val(Phi,0,1,0);_stencil_val_a(p.y,0,0,0);       }}





           
# 171 "/home/spencer/basilisk/src/fractions.h"
    
  
}}end__stencil_is_face_x()}end_foreach_face_stencil();
# 146 "/home/spencer/basilisk/src/fractions.h"
  {foreach_face_generic(){is_face_y(){ {





    if ((val(Phi,0,0,0) - val)*(val(Phi,1,0,0) - val) < 0.) {






      val(p.x,0,0,0) = (val(Phi,0,0,0) - val)/(val(Phi,0,0,0) - val(Phi,1,0,0));
      if (val(Phi,0,0,0) < val)
 val(p.x,0,0,0) = 1. - val(p.x,0,0,0);
    }
# 171 "/home/spencer/basilisk/src/fractions.h"
    else
      val(p.x,0,0,0) = (val(Phi,0,0,0) > val || val(Phi,1,0,0) > val);
  }}end_is_face_y()
#line 146
is_face_x(){ {





    if ((val(Phi,0,0,0) - val)*(val(Phi,0,1,0) - val) < 0.) {






      val(p.y,0,0,0) = (val(Phi,0,0,0) - val)/(val(Phi,0,0,0) - val(Phi,0,1,0));
      if (val(Phi,0,0,0) < val)
 val(p.y,0,0,0) = 1. - val(p.y,0,0,0);
    }
# 171 "/home/spencer/basilisk/src/fractions.h"
    else
      val(p.y,0,0,0) = (val(Phi,0,0,0) > val || val(Phi,0,1,0) > val);
  }}end_is_face_x()}end_foreach_face_generic();}
# 196 "/home/spencer/basilisk/src/fractions.h"
  scalar s_z = c;
  foreach_stencil(1,{(NonLocal[]){{"Phi","scalar",(void *)&Phi,NULL,0},{"s_z","scalar",(void *)&s_z,NULL,0},{"p","vector",(void *)&p,NULL,0},{"val","double",(void *)&val,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_("\n// # 163 \"/home/spencer/basilisk/src/geometry.h\"\nreal line_area (real nx, real ny, real alpha)\n{\n  real a, v, area;\n\n  alpha += (nx + ny)/2.;\n  if (nx < 0.) {\n    alpha -= nx;\n    nx = - nx;\n  }\n  if (ny < 0.) {\n    alpha -= ny;\n    ny = - ny;\n  }\n\n  if (alpha <= 0.)\n    return 0.;\n\n  if (alpha >= nx + ny)\n    return 1.;\n\n  if (nx < 1e-10)\n    area = alpha/ny;\n  else if (ny < 1e-10)\n    area = alpha/nx;\n  else {\n    v = sq(alpha);\n\n    a = alpha - nx;\n    if (a > 0.)\n      v -= a*a;\n\n    a = alpha - ny;\n    if (a > 0.)\n      v -= a*a;\n\n    area = v/(2.*nx*ny);\n  }\n\n  return clamp (area, 0., 1.);\n}"),_("\n\n  \n// #line 199 \"/home/spencer/basilisk/src/fractions.h\"\n{\n// # 231 \"/home/spencer/basilisk/src/fractions.h\"\n    coord n;\n    real nn = 0.;\n     {\n      n.x = val(p.y,0,0,0) - val(p.y,1,0,0);\n      nn += fabs(n.x);\n    } \n// #line 233\n{\n      n.y = val(p.x,0,0,0) - val(p.x,0,1,0);\n      nn += fabs(n.y);\n    }\n\n\n\n\n\n    if (nn == 0.)\n      val_out_(s_z,0,0,0) = val(p.x,0,0,0);\n    else {\n\n\n\n\n\n      \n n.x /= nn;\n \n// #line 251\nn.y /= nn;\n\n\n\n\n\n\n      real alpha = 0., ni = 0.;\n      for (int i = 0; i <= 1; i++)\n {\n   if (val(p.x,0,i,0) > 0. && val(p.x,0,i,0) < 1.) {\n     real a = sign(val(Phi,0,i,0) - val)*(val(p.x,0,i,0) - 0.5);\n     alpha += n.x*a + n.y*(i - 0.5);\n     ni++;\n   }\n   \n// #line 261\nif (val(p.y,i,0,0) > 0. && val(p.y,i,0,0) < 1.) {\n     real a = sign(val(Phi,i,0,0) - val)*(val(p.y,i,0,0) - 0.5);\n     alpha += n.y*a + n.x*(i - 0.5);\n     ni++;\n   }}\n// # 274 \"/home/spencer/basilisk/src/fractions.h\"\n      if (ni == 0)\n val_out_(s_z,0,0,0) = max (val(p.x,0,0,0), val(p.y,0,0,0));\n      else if (ni != 4)\n val_out_(s_z,0,0,0) = line_area (n.x, n.y, alpha/ni);\n      else {\n\n\n\n val_out_(s_z,0,0,0) = 0.;\n\n      }\n    }\n  }")})

  {    
# 231 "/home/spencer/basilisk/src/fractions.h"
    
    
     { 
_stencil_val(p.y,0,0,0); _stencil_val(p.y,1,0,0);  
       
       
    
#line 236
} 
#line 233
{ 
_stencil_val(p.x,0,0,0); _stencil_val(p.x,0,1,0);  
       
       
    
#line 236
}





{
      { _stencil_val(p.x,0,0,0);_stencil_val_a(s_z,0,0,0); } 
{      





      
   






      
      for (int i = 0; i <= 1; i++)
 {
   {_stencil_val(p.x,0,i,0); _stencil_val(p.x,0,i,0); {       
     _stencil_val(p.x,0,i,0);_stencil_val(Phi,0,i,0); 
          
     
   }      }
   
#line 261
{_stencil_val(p.y,i,0,0); _stencil_val(p.y,i,0,0); {       
     _stencil_val(p.y,i,0,0);_stencil_val(Phi,i,0,0); 
          
     
   }      }}








{
 {_stencil_val(p.x,0,0,0); _stencil_val(p.y,0,0,0);_stencil_val_a(s_z,0,0,0);   }
{
 {_stencil_val_a(s_z,0,0,0);     } 
{



 _stencil_val_a(s_z,0,0,0);  

      }}}
# 274 "/home/spencer/basilisk/src/fractions.h"
         
          
      
    







}}





       
    
  
#line 286
}end_foreach_stencil();
  {
#line 197
foreach()

  {
# 231 "/home/spencer/basilisk/src/fractions.h"
    coord n;
    double nn = 0.;
     {
      n.x = val(p.y,0,0,0) - val(p.y,1,0,0);
      nn += fabs(n.x);
    } 
#line 233
{
      n.y = val(p.x,0,0,0) - val(p.x,0,1,0);
      nn += fabs(n.y);
    }





    if (nn == 0.)
      val(s_z,0,0,0) = val(p.x,0,0,0);
    else {





      
 n.x /= nn;
 
#line 251
n.y /= nn;






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 {
   if (val(p.x,0,i,0) > 0. && val(p.x,0,i,0) < 1.) {
     double a = sign(val(Phi,0,i,0) - val)*(val(p.x,0,i,0) - 0.5);
     alpha += n.x*a + n.y*(i - 0.5);
     ni++;
   }
   
#line 261
if (val(p.y,i,0,0) > 0. && val(p.y,i,0,0) < 1.) {
     double a = sign(val(Phi,i,0,0) - val)*(val(p.y,i,0,0) - 0.5);
     alpha += n.y*a + n.x*(i - 0.5);
     ni++;
   }}
# 274 "/home/spencer/basilisk/src/fractions.h"
      if (ni == 0)
 val(s_z,0,0,0) = max (val(p.x,0,0,0), val(p.y,0,0,0));
      else if (ni != 4)
 val(s_z,0,0,0) = line_area (n.x, n.y, alpha/ni);
      else {



 val(s_z,0,0,0) = 0.;

      }
    }
  }end_foreach();}if((s).x.i<=0)delete((scalar*)((vector[]){as,{{-1},{-1}}}));
# 351 "/home/spencer/basilisk/src/fractions.h"
end_tracing("fractions","/home/spencer/basilisk/src/fractions.h",0);}
# 395 "/home/spencer/basilisk/src/fractions.h"
coord youngs_normal (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  coord n;
  double nn = 0.;
  if (!(2 == 2)) qassert ("/home/spencer/basilisk/src/fractions.h", 0, "dimension == 2");
   {
    n.x = (val(c,-1,1,0) + 2.*val(c,-1,0,0) + val(c,-1,-1,0) -
    val(c,+1,1,0) - 2.*val(c,+1,0,0) - val(c,+1,-1,0));
    nn += fabs(n.x);
  } 
#line 400
{
    n.y = (val(c,1,-1,0) + 2.*val(c,0,-1,0) + val(c,-1,-1,0) -
    val(c,1,+1,0) - 2.*val(c,0,+1,0) - val(c,-1,+1,0));
    nn += fabs(n.y);
  }

  if (nn > 0.)
    {
      n.x /= nn;
      
#line 408
n.y /= nn;}
  else
    n.x = 1.;
  return n;
}





coord facet_normal (Point point, scalar c, vector s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (s.x.i >= 0) {
    coord n;
    double nn = 0.;
     {
      n.x = val(s.x,0,0,0) - val(s.x,1,0,0);
      nn += fabs(n.x);
    } 
#line 423
{
      n.y = val(s.y,0,0,0) - val(s.y,0,1,0);
      nn += fabs(n.y);
    }
    if (nn > 0.)
      {
 n.x /= nn;
 
#line 429
n.y /= nn;}
    else
      {
 n.x = 1./2;
 
#line 432
n.y = 1./2;}
    return n;
  }
  return mycs (point, c);
}






#line 418
static void _stencil_facet_normal (Point point, scalar c, vector s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (s.x.i >= 0) {    
    
    
     { 
_stencil_val(s.x,0,0,0); _stencil_val(s.x,1,0,0);  
       
       
    
#line 426
} 
#line 423
{ 
_stencil_val(s.y,0,0,0); _stencil_val(s.y,0,1,0);  
       
       
    
#line 426
}
      
   
       
    
       
  
    return ;
  } 
_stencil_mycs (point, c);
  
#line 435
return;
}
# 445 "/home/spencer/basilisk/src/fractions.h"
     
void reconstruction (const scalar c, vector n, scalar alpha)
{tracing("reconstruction","/home/spencer/basilisk/src/fractions.h",0);
  foreach_stencil(1,{(NonLocal[]){{"n","vector",(void *)&n,NULL,0},{"alpha","scalar",(void *)&alpha,NULL,0},{"c","scalar",(void *)&c,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_("\n// # 415 \"/home/spencer/basilisk/src/view.h\" 2\n\n\n\n\n\n\n// # 1 \"./draw.h\" 1\n// # 1 \"/home/spencer/basilisk/src/draw.h\"\n\n\n\n\n// # 1 \"./fractions.h\" 1\n// # 1 \"/home/spencer/basilisk/src/fractions.h\"\n// # 12 \"/home/spencer/basilisk/src/fractions.h\"\n// # 1 \"./geometry.h\" 1\n// # 1 \"/home/spencer/basilisk/src/geometry.h\"\n// # 35 \"/home/spencer/basilisk/src/geometry.h\"\nreal line_alpha (real c, coord n)\n{\n  real alpha, n1, n2;\n\n  n1 = fabs (n.x); n2 = fabs (n.y);\n  if (n1 > n2)\n    do { real __tmp = n1; n1 = n2; n2 = __tmp; } while(0);\n\n  c = clamp (c, 0., 1.);\n  real v1 = n1/2.;\n  if (c <= v1/n2)\n    alpha = sqrt (2.*c*n1*n2);\n  else if (c <= 1. - v1/n2)\n    alpha = c*n2 + v1;\n  else\n    alpha = n1 + n2 - sqrt (2.*n1*n2*(1. - c));\n\n  if (n.x < 0.)\n    alpha += n.x;\n  if (n.y < 0.)\n    alpha += n.y;\n\n  return alpha - (n.x + n.y)/2.;\n}\n// # 13 \"/home/spencer/basilisk/src/fractions.h\" 2\n\n\n\n\n\n// # 1 \"./myc2d.h\" 1\n// # 1 \"/home/spencer/basilisk/src/myc2d.h\"\n\n\n\n\n\ncoord mycs (Point point, scalar c)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int ix;\n  real c_t,c_b,c_r,c_l;\n  real mx0,my0,mx1,my1,mm1,mm2;\n\n\n  c_t = val(c,-1,1,0) + val(c,0,1,0) + val(c,1,1,0);\n  c_b = val(c,-1,-1,0) + val(c,0,-1,0) + val(c,1,-1,0);\n  c_r = val(c,1,-1,0) + val(c,1,0,0) + val(c,1,1,0);\n  c_l = val(c,-1,-1,0) + val(c,-1,0,0) + val(c,-1,1,0);\n\n\n\n  mx0 = 0.5*(c_l - c_r);\n  my0 = 0.5*(c_b - c_t);\n\n\n  if (fabs(mx0) <= fabs(my0)) {\n    my0 = my0 > 0. ? 1. : -1.;\n    ix = 1;\n  }\n  else {\n    mx0 = mx0 > 0. ? 1. : -1.;\n    ix = 0;\n  }\n\n\n  mm1 = val(c,-1,-1,0) + 2.0*val(c,-1,0,0) + val(c,-1,1,0);\n  mm2 = val(c,1,-1,0) + 2.0*val(c,1,0,0) + val(c,1,1,0);\n  mx1 = mm1 - mm2 + 1e-30;\n  mm1 = val(c,-1,-1,0) + 2.0*val(c,0,-1,0) + val(c,1,-1,0);\n  mm2 = val(c,-1,1,0) + 2.0*val(c,0,1,0) + val(c,1,1,0);\n  my1 = mm1 - mm2 + 1e-30;\n\n\n  if (ix) {\n    mm1 = fabs(my1);\n    mm1 = fabs(mx1)/mm1;\n    if (mm1 > fabs(mx0)) {\n      mx0 = mx1;\n      my0 = my1;\n    }\n  }\n  else {\n    mm1 = fabs(mx1);\n    mm1 = fabs(my1)/mm1;\n    if (mm1 > fabs(my0)) {\n      mx0 = mx1;\n      my0 = my1;\n    }\n  }\n\n\n\n  mm1 = fabs(mx0) + fabs(my0);\n  coord n = {mx0/mm1, my0/mm1};\n\n  return n;\n}"),_(" \n// #line 448 \"/home/spencer/basilisk/src/fractions.h\"\n{\n\n\n\n\n\n    if (val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) {\n      val_out_(alpha,0,0,0) = 0.;\n      \n val_out_(n.x,0,0,0) = 0.;\n \n// #line 457\nval_out_(n.y,0,0,0) = 0.;\n    }\n    else {\n\n\n\n\n\n\n      coord m = mycs (point, c);\n      \n val_out_(n.x,0,0,0) = m.x;\n \n// #line 468\nval_out_(n.y,0,0,0) = m.y;\n      val_out_(alpha,0,0,0) = line_alpha (val(c,0,0,0), m);\n    }\n  }")}) {





_stencil_val(c,0,0,0); _stencil_val(c,0,0,0);{ {
      _stencil_val_a(alpha,0,0,0);  
      
 {_stencil_val_a(n.x,0,0,0);  }
 
#line 457
{_stencil_val_a(n.y,0,0,0);  }
    } 
{  






       _stencil_mycs (point, c);
      
 {_stencil_val_a(n.x,0,0,0);  }
 
#line 468
{_stencil_val_a(n.y,0,0,0);  }
_stencil_val(c,0,0,0);
      
#line 469
_stencil_val_a(alpha,0,0,0);    
    }}





          
    
  
#line 471
}end_foreach_stencil();
  {
#line 448
foreach() {





    if (val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) {
      val(alpha,0,0,0) = 0.;
      
 val(n.x,0,0,0) = 0.;
 
#line 457
val(n.y,0,0,0) = 0.;
    }
    else {






      coord m = mycs (point, c);
      
 val(n.x,0,0,0) = m.x;
 
#line 468
val(n.y,0,0,0) = m.y;
      val(alpha,0,0,0) = line_alpha (val(c,0,0,0), m);
    }
  }end_foreach();}
# 480 "/home/spencer/basilisk/src/fractions.h"
  
    _attribute[n.x.i].refine = _attribute[n.x.i].prolongation = refine_injection;
    
#line 481
_attribute[n.y.i].refine = _attribute[n.y.i].prolongation = refine_injection;




  _attribute[alpha.i].n = n;
  _attribute[alpha.i].refine = _attribute[alpha.i].prolongation = alpha_refine;

end_tracing("reconstruction","/home/spencer/basilisk/src/fractions.h",0);}
# 509 "/home/spencer/basilisk/src/fractions.h"
     
void output_facets (scalar c, FILE * fp, vector s)
{tracing("output_facets","/home/spencer/basilisk/src/fractions.h",0);
  foreach_stencil(1,{(NonLocal[]){{"fp","not implemented yet",(void *)fp,NULL,1},{"s","vector",(void *)&s,NULL,0},{"c","scalar",(void *)&c,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_("\n// # 292 \"/home/spencer/basilisk/src/geometry.h\"\nint facets (coord n, real alpha, coord p[2])\n{\n  int i = 0;\n  for (real s = -0.5; s <= 0.5; s += 1.)\n    {\n      if (fabs (n.y) > 1e-4 && i < 2) {\n real a = (alpha - s*n.x)/n.y;\n if (a >= -0.5 && a <= 0.5) {\n   p[i].x = s;\n   p[i++].y = a;\n }\n      }\n      \n// #line 297\nif (fabs (n.x) > 1e-4 && i < 2) {\n real a = (alpha - s*n.y)/n.x;\n if (a >= -0.5 && a <= 0.5) {\n   p[i].y = s;\n   p[i++].x = a;\n }\n      }}\n  return i;\n}\n// # 415 \"/home/spencer/basilisk/src/view.h\" 2\n\n\n\n\n\n\n// # 1 \"./draw.h\" 1\n// # 1 \"/home/spencer/basilisk/src/draw.h\"\n\n\n\n\n// # 1 \"./fractions.h\" 1\n// # 1 \"/home/spencer/basilisk/src/fractions.h\"\n// # 12 \"/home/spencer/basilisk/src/fractions.h\"\n// # 1 \"./geometry.h\" 1\n// # 1 \"/home/spencer/basilisk/src/geometry.h\"\n// # 35 \"/home/spencer/basilisk/src/geometry.h\"\nreal line_alpha (real c, coord n)\n{\n  real alpha, n1, n2;\n\n  n1 = fabs (n.x); n2 = fabs (n.y);\n  if (n1 > n2)\n    do { real __tmp = n1; n1 = n2; n2 = __tmp; } while(0);\n\n  c = clamp (c, 0., 1.);\n  real v1 = n1/2.;\n  if (c <= v1/n2)\n    alpha = sqrt (2.*c*n1*n2);\n  else if (c <= 1. - v1/n2)\n    alpha = c*n2 + v1;\n  else\n    alpha = n1 + n2 - sqrt (2.*n1*n2*(1. - c));\n\n  if (n.x < 0.)\n    alpha += n.x;\n  if (n.y < 0.)\n    alpha += n.y;\n\n  return alpha - (n.x + n.y)/2.;\n}\n// # 13 \"/home/spencer/basilisk/src/fractions.h\" 2\n\n\n\n\n\n// # 1 \"./myc2d.h\" 1\n// # 1 \"/home/spencer/basilisk/src/myc2d.h\"\n\n\n\n\n\ncoord mycs (Point point, scalar c)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int ix;\n  real c_t,c_b,c_r,c_l;\n  real mx0,my0,mx1,my1,mm1,mm2;\n\n\n  c_t = val(c,-1,1,0) + val(c,0,1,0) + val(c,1,1,0);\n  c_b = val(c,-1,-1,0) + val(c,0,-1,0) + val(c,1,-1,0);\n  c_r = val(c,1,-1,0) + val(c,1,0,0) + val(c,1,1,0);\n  c_l = val(c,-1,-1,0) + val(c,-1,0,0) + val(c,-1,1,0);\n\n\n\n  mx0 = 0.5*(c_l - c_r);\n  my0 = 0.5*(c_b - c_t);\n\n\n  if (fabs(mx0) <= fabs(my0)) {\n    my0 = my0 > 0. ? 1. : -1.;\n    ix = 1;\n  }\n  else {\n    mx0 = mx0 > 0. ? 1. : -1.;\n    ix = 0;\n  }\n\n\n  mm1 = val(c,-1,-1,0) + 2.0*val(c,-1,0,0) + val(c,-1,1,0);\n  mm2 = val(c,1,-1,0) + 2.0*val(c,1,0,0) + val(c,1,1,0);\n  mx1 = mm1 - mm2 + 1e-30;\n  mm1 = val(c,-1,-1,0) + 2.0*val(c,0,-1,0) + val(c,1,-1,0);\n  mm2 = val(c,-1,1,0) + 2.0*val(c,0,1,0) + val(c,1,1,0);\n  my1 = mm1 - mm2 + 1e-30;\n\n\n  if (ix) {\n    mm1 = fabs(my1);\n    mm1 = fabs(mx1)/mm1;\n    if (mm1 > fabs(mx0)) {\n      mx0 = mx1;\n      my0 = my1;\n    }\n  }\n  else {\n    mm1 = fabs(mx1);\n    mm1 = fabs(my1)/mm1;\n    if (mm1 > fabs(my0)) {\n      mx0 = mx1;\n      my0 = my1;\n    }\n  }\n\n\n\n  mm1 = fabs(mx0) + fabs(my0);\n  coord n = {mx0/mm1, my0/mm1};\n\n  return n;\n}\n\n\n\n\n\n\n// #line 418 \"/home/spencer/basilisk/src/fractions.h\"\ncoord facet_normal (Point point, scalar c, vector s)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  if (s.x.i >= 0) {\n    coord n;\n    real nn = 0.;\n     {\n      n.x = val(s.x,0,0,0) - val(s.x,1,0,0);\n      nn += fabs(n.x);\n    } \n// #line 423\n{\n      n.y = val(s.y,0,0,0) - val(s.y,0,1,0);\n      nn += fabs(n.y);\n    }\n    if (nn > 0.)\n      {\n n.x /= nn;\n \n// #line 429\nn.y /= nn;}\n    else\n      {\n n.x = 1./2;\n \n// #line 432\nn.y = 1./2;}\n    return n;\n  }\n  return mycs (point, c);\n}"),_("\n    \n// #line 513 \"/home/spencer/basilisk/src/fractions.h\"\nif (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {\n      coord n = facet_normal (point, c, s);\n      real alpha = line_alpha (val(c,0,0,0), n);\n\n\n\n      coord segment[2];\n      if (facets (n, alpha, segment) == 2)\n fprintf (fp, \"%g %g\\n%g %g\\n\\n\",\n   x + segment[0].x*Delta, y + segment[0].y*Delta,\n   x + segment[1].x*Delta, y + segment[1].y*Delta);\n// # 533 \"/home/spencer/basilisk/src/fractions.h\"\n    }")})
    {_stencil_val(c,0,0,0); _stencil_val(c,0,0,0); {  
       _stencil_facet_normal (point, c, s);     
      _stencil_val(c,0,0,0); 



            
           
        
     
 
# 533 "/home/spencer/basilisk/src/fractions.h"
    }        }end_foreach_stencil();
  {
#line 512
foreach()
    if (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {
      coord n = facet_normal (point, c, s);
      double alpha = line_alpha (val(c,0,0,0), n);



      coord segment[2];
      if (facets (n, alpha, segment) == 2)
 fprintf (fp, "%g %g\n%g %g\n\n",
   x + segment[0].x*Delta, y + segment[0].y*Delta,
   x + segment[1].x*Delta, y + segment[1].y*Delta);
# 533 "/home/spencer/basilisk/src/fractions.h"
    }end_foreach();}

  fflush (fp);
end_tracing("output_facets","/home/spencer/basilisk/src/fractions.h",0);}







     
double interface_area (scalar c)
{tracing("interface_area","/home/spencer/basilisk/src/fractions.h",0);
  double area = 0.;
  foreach_stencil (1,{(NonLocal[]){{"area","double",(void *)&area,NULL,0,'+'},{"c","scalar",(void *)&c,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_("\n// # 382 \"/home/spencer/basilisk/src/geometry.h\"\nreal line_length_center (coord m, real alpha, coord * p)\n{\n  alpha += (m.x + m.y)/2.;\n\n  coord n = m;\n  \n    if (n.x < 0.) {\n      alpha -= n.x;\n      n.x = - n.x;\n    }\n    \n// #line 388\nif (n.y < 0.) {\n      alpha -= n.y;\n      n.y = - n.y;\n    }\n\n  p->x = p->y = p->z = 0.;\n\n  if (alpha <= 0. || alpha >= n.x + n.y)\n    return 0.;\n\n  \n    if (n.x < 1e-4) {\n      p->x = 0.;\n      p->y = (m.y < 0. ? 1. - alpha : alpha) - 0.5;\n      return 1.;\n    }\n    \n// #line 399\nif (n.y < 1e-4) {\n      p->y = 0.;\n      p->x = (m.x < 0. ? 1. - alpha : alpha) - 0.5;\n      return 1.;\n    }\n\n  if (alpha >= n.x) {\n    p->x += 1.;\n    p->y += (alpha - n.x)/n.y;\n  }\n  else\n    p->x += alpha/n.x;\n\n  real ax = p->x, ay = p->y;\n  if (alpha >= n.y) {\n    p->y += 1.;\n    ay -= 1.;\n    p->x += (alpha - n.y)/n.x;\n    ax -= (alpha - n.y)/n.x;\n  }\n  else {\n    p->y += alpha/n.y;\n    ay -= alpha/n.y;\n  }\n\n   {\n    p->x /= 2.;\n    p->x = clamp (p->x, 0., 1.);\n    if (m.x < 0.)\n      p->x = 1. - p->x;\n    p->x -= 0.5;\n  } \n// #line 424\n{\n    p->y /= 2.;\n    p->y = clamp (p->y, 0., 1.);\n    if (m.y < 0.)\n      p->y = 1. - p->y;\n    p->y -= 0.5;\n  }\n\n  return sqrt (ax*ax + ay*ay);\n}\n// # 415 \"/home/spencer/basilisk/src/view.h\" 2\n\n\n\n\n\n\n// # 1 \"./draw.h\" 1\n// # 1 \"/home/spencer/basilisk/src/draw.h\"\n\n\n\n\n// # 1 \"./fractions.h\" 1\n// # 1 \"/home/spencer/basilisk/src/fractions.h\"\n// # 12 \"/home/spencer/basilisk/src/fractions.h\"\n// # 1 \"./geometry.h\" 1\n// # 1 \"/home/spencer/basilisk/src/geometry.h\"\n// # 35 \"/home/spencer/basilisk/src/geometry.h\"\nreal line_alpha (real c, coord n)\n{\n  real alpha, n1, n2;\n\n  n1 = fabs (n.x); n2 = fabs (n.y);\n  if (n1 > n2)\n    do { real __tmp = n1; n1 = n2; n2 = __tmp; } while(0);\n\n  c = clamp (c, 0., 1.);\n  real v1 = n1/2.;\n  if (c <= v1/n2)\n    alpha = sqrt (2.*c*n1*n2);\n  else if (c <= 1. - v1/n2)\n    alpha = c*n2 + v1;\n  else\n    alpha = n1 + n2 - sqrt (2.*n1*n2*(1. - c));\n\n  if (n.x < 0.)\n    alpha += n.x;\n  if (n.y < 0.)\n    alpha += n.y;\n\n  return alpha - (n.x + n.y)/2.;\n}\n// # 13 \"/home/spencer/basilisk/src/fractions.h\" 2\n\n\n\n\n\n// # 1 \"./myc2d.h\" 1\n// # 1 \"/home/spencer/basilisk/src/myc2d.h\"\n\n\n\n\n\ncoord mycs (Point point, scalar c)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int ix;\n  real c_t,c_b,c_r,c_l;\n  real mx0,my0,mx1,my1,mm1,mm2;\n\n\n  c_t = val(c,-1,1,0) + val(c,0,1,0) + val(c,1,1,0);\n  c_b = val(c,-1,-1,0) + val(c,0,-1,0) + val(c,1,-1,0);\n  c_r = val(c,1,-1,0) + val(c,1,0,0) + val(c,1,1,0);\n  c_l = val(c,-1,-1,0) + val(c,-1,0,0) + val(c,-1,1,0);\n\n\n\n  mx0 = 0.5*(c_l - c_r);\n  my0 = 0.5*(c_b - c_t);\n\n\n  if (fabs(mx0) <= fabs(my0)) {\n    my0 = my0 > 0. ? 1. : -1.;\n    ix = 1;\n  }\n  else {\n    mx0 = mx0 > 0. ? 1. : -1.;\n    ix = 0;\n  }\n\n\n  mm1 = val(c,-1,-1,0) + 2.0*val(c,-1,0,0) + val(c,-1,1,0);\n  mm2 = val(c,1,-1,0) + 2.0*val(c,1,0,0) + val(c,1,1,0);\n  mx1 = mm1 - mm2 + 1e-30;\n  mm1 = val(c,-1,-1,0) + 2.0*val(c,0,-1,0) + val(c,1,-1,0);\n  mm2 = val(c,-1,1,0) + 2.0*val(c,0,1,0) + val(c,1,1,0);\n  my1 = mm1 - mm2 + 1e-30;\n\n\n  if (ix) {\n    mm1 = fabs(my1);\n    mm1 = fabs(mx1)/mm1;\n    if (mm1 > fabs(mx0)) {\n      mx0 = mx1;\n      my0 = my1;\n    }\n  }\n  else {\n    mm1 = fabs(mx1);\n    mm1 = fabs(my1)/mm1;\n    if (mm1 > fabs(my0)) {\n      mx0 = mx1;\n      my0 = my1;\n    }\n  }\n\n\n\n  mm1 = fabs(mx0) + fabs(my0);\n  coord n = {mx0/mm1, my0/mm1};\n\n  return n;\n}"),_("\n    \n// #line 549 \"/home/spencer/basilisk/src/fractions.h\"\nif (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {\n      coord n = mycs (point, c), p;\n      real alpha = line_alpha (val(c,0,0,0), n);\n      area += pow(Delta, 2 - 1)*line_length_center(n,alpha,&p);\n    }")})
    {_stencil_val(c,0,0,0); _stencil_val(c,0,0,0); {   
       _stencil_mycs (point, c);     
      _stencil_val(c,0,0,0); 
          
    }        }end_foreach_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(+:area)){
#line 548
foreach ()
    if (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {
      coord n = mycs (point, c), p;
      double alpha = line_alpha (val(c,0,0,0), n);
      area += pow(Delta, 2 - 1)*line_length_center(n,alpha,&p);
    }end_foreach();mpi_all_reduce_array(&area,double,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
  
#line 554
{end_tracing("interface_area","/home/spencer/basilisk/src/fractions.h",0);return area;}
end_tracing("interface_area","/home/spencer/basilisk/src/fractions.h",0);}
# 6 "/home/spencer/basilisk/src/draw.h" 2
# 1 "./gl/font.h" 1
# 1 "/home/spencer/basilisk/src/gl/font.h"
# 27 "/home/spencer/basilisk/src/gl/font.h"
# 1 "/home/spencer/basilisk/src/ast/std/stdio.h" 1
@include <stdio.h>
# 28 "/home/spencer/basilisk/src/gl/font.h" 2
# 1 "./gl/og_font.h" 1
# 1 "/home/spencer/basilisk/src/gl/og_font.h"




# 1 "./gl/tinygl.h" 1
# 6 "/home/spencer/basilisk/src/gl/og_font.h" 2

typedef struct tagSOG_StrokeVertex SOG_StrokeVertex;
struct tagSOG_StrokeVertex
{
    GLfloat X, Y;
};

typedef struct tagSOG_StrokeStrip SOG_StrokeStrip;
struct tagSOG_StrokeStrip
{
    int Number;
    const SOG_StrokeVertex *Vertices;
};

typedef struct tagSOG_StrokeChar SOG_StrokeChar;
struct tagSOG_StrokeChar
{
    GLfloat Right;
    int Number;
    const SOG_StrokeStrip* Strips;
};

typedef struct tagSOG_StrokeFont SOG_StrokeFont;
struct tagSOG_StrokeFont
{
    char *Name;
    int Quantity;
    GLfloat Height;
    const SOG_StrokeChar **Characters;
};
# 29 "/home/spencer/basilisk/src/gl/font.h" 2
# 39 "/home/spencer/basilisk/src/gl/font.h"
extern SOG_StrokeFont ogStrokeMonoRoman;
# 48 "/home/spencer/basilisk/src/gl/font.h"
static SOG_StrokeFont *oghStrokeByID( void *font )
{


    if( font == ((void *)0x0001) )
        return &ogStrokeMonoRoman;

    fprintf (ferr, "stroke font %p not found", font );
    return 0;
}
# 83 "/home/spencer/basilisk/src/gl/font.h"
void gl_StrokeCharacter( int character )
{
    void *fontID = ((void *)0x0001);
    const SOG_StrokeChar *schar;
    const SOG_StrokeStrip *strip;
    int i, j;
    SOG_StrokeFont *font = oghStrokeByID( fontID );

    if( !font ||
        ( 1 > character ) ||
        ( font->Quantity < character ) )
        return;

    schar = font->Characters[ character ];
    if( schar )
    {
        strip = schar->Strips;

        for( i = 0; i < schar->Number; i++, strip++ )
        {
            vertex_buffer_glBegin( 0x0003 );
            for( j = 0; j < strip->Number; j++ )
                vertex_buffer_glVertex2d( strip->Vertices[ j ].X, strip->Vertices[ j ].Y );
            vertex_buffer_glEnd( );
        }
        glTranslatef( schar->Right, 0.0, 0.0 );
    }
}
# 147 "/home/spencer/basilisk/src/gl/font.h"
void gl_StrokeString( const char *string )
{
    void *fontID = ((void *)0x0001);
    int i, j;
    float length = 0.0;
    SOG_StrokeFont *font = oghStrokeByID( fontID );
    unsigned char c;

    if( font && string )





        while(( c = *string++ ))
       if( c < font->Quantity ) {
                if( c == '\n' )
                {
                    glTranslatef ( -length, -( float )( font->Height ), 0.0 );
                    length = 0.0;
                }
                else
                {
                    const SOG_StrokeChar *schar =
                        font->Characters[ c ];
                    if( schar )
                    {
                        const SOG_StrokeStrip *strip = schar->Strips;

                        for( i = 0; i < schar->Number; i++, strip++ )
                        {
                            vertex_buffer_glBegin( 0x0003 );

                            for( j = 0; j < strip->Number; j++ )
                                vertex_buffer_glVertex2d( strip->Vertices[ j ].X,
                                            strip->Vertices[ j ].Y);

                            vertex_buffer_glEnd( );
                        }

                        length += schar->Right;
                        glTranslatef( schar->Right, 0.0, 0.0 );
                    }
                }
     }
}
# 226 "/home/spencer/basilisk/src/gl/font.h"
float gl_StrokeWidth( int character )
{
    void *fontID = ((void *)0x0001);
    float ret = 0;
    SOG_StrokeFont *font = oghStrokeByID( fontID );

    if( font &&
        ( 0 < character ) &&
        ( font->Quantity > character ) )
    {
        const SOG_StrokeChar *schar = font->Characters[ character ];
        if( schar )
            ret = schar->Right;
    }

    return ret;
}
# 269 "/home/spencer/basilisk/src/gl/font.h"
float gl_StrokeLength( const char *string )
{
    void *fontID = ((void *)0x0001);
    unsigned char c;
    float length = 0.0;
    float this_line_length = 0.0;
    SOG_StrokeFont *font = oghStrokeByID( fontID );

    if( font && string )
        while(( c = *string++ ))
            if( c < font->Quantity )
            {
                if( c == '\n' )
                {
                    if( length < this_line_length )
                        length = this_line_length;
                    this_line_length = 0.0;
                }
                else
                {
                    const SOG_StrokeChar *schar =
                        font->Characters[ c ];
                    if( schar )
                        this_line_length += schar->Right;
                }
            }

    if( length < this_line_length )
        length = this_line_length;
    return length;
}
# 321 "/home/spencer/basilisk/src/gl/font.h"
GLfloat gl_StrokeHeight()
{
    void *fontID = ((void *)0x0001);
    GLfloat ret = 0;
    SOG_StrokeFont *font = oghStrokeByID( fontID );

    if( font )
        ret = font->Height;

    return ret;
}
# 7 "/home/spencer/basilisk/src/draw.h" 2




void clear()
{
  bview * view = get_view();
  if (view->active)
    view->active = false;
  draw();
}
# 49 "/home/spencer/basilisk/src/draw.h"
void view (float tx, float ty,
    float fov,
    float quat[4],
    float sx, float sy, float sz,
    unsigned width, unsigned height, unsigned samples,
    float bg[3],
    float theta, float phi, float psi,
    bool relative,
    float tz, float near, float far,
    float res,
    char * camera,
    MapFunc map,
    int cache,
    float p1x, float p1y, float p2x, float p2y,
    bview * view1)
{
  bview * v = view1 ? view1 : get_view();
  if (fov) {
    if (relative)
      v->fov += (0.1 + 3.*v->fov)*fov;
    else
      v->fov = fov;
    v->fov = clamp(v->fov,0.01,100.);
  }
  for (int i = 0; i < 4; i++)
    if (quat[i]) {
      for (int j = 0; j < 4; j++)
 v->quat[j] = quat[j];
      break;
    }
  v->tx = relative ? v->tx + tx*0.02*(0.01 + 3.*v->fov) : tx;
  v->ty = relative ? v->ty + ty*0.02*(0.01 + 3.*v->fov) : ty;
  v->sx = sx;
  v->sy = sy;
  v->sz = sz;
  if (bg[0] || bg[1] || bg[2])
    for (int i = 0; i < 3; i++)
      v->bg[i] = bg[i];

  if (camera) {
    v->gfsview = false;
    if (strlen(camera) >= 4 &&
 !strcmp (&camera[strlen(camera) - 4], ".gfv")) {
      FILE * fp = fopen (camera, "r");
      if (!fp) {
 perror (camera);
 exit (1);
      }
      char s[81];
      float q[4], fov;
      int nq = 0, nf = 0;
      while (fgets (s, 81, fp) && (!nq || !nf)) {
 if (!nq)
   nq = sscanf (s, "  q0 = %f q1 = %f q2 = %f q3 = %f",
         &q[0], &q[1], &q[2], &q[3]);
 if (!nf)
   nf = sscanf (s, "  fov = %f", &fov);
      }
      if (nq != 4 || nf != 1) {
 fprintf (ferr, "%s: not a valid gfv file\n", camera);
 exit (1);
      }
      for (int j = 0; j < 4; j++)
 v->quat[j] = q[j];
      v->fov = fov;
      v->gfsview = true;
    }
    else if (!strcmp (camera, "left"))
      gl_axis_to_quat ((float[]){0,1,0}, - 3.14159265358979/2., v->quat);
    else if (!strcmp (camera, "right"))
      gl_axis_to_quat ((float[]){0,1,0}, 3.14159265358979/2., v->quat);
    else if (!strcmp (camera, "top"))
      gl_axis_to_quat ((float[]){1,0,0}, - 3.14159265358979/2., v->quat);
    else if (!strcmp (camera, "bottom"))
      gl_axis_to_quat ((float[]){1,0,0}, 3.14159265358979/2., v->quat);
    else if (!strcmp (camera, "front"))
      gl_axis_to_quat ((float[]){0,0,1}, 0., v->quat);
    else if (!strcmp (camera, "back"))
      gl_axis_to_quat ((float[]){0,1,0}, 3.14159265358979, v->quat);
    else if (!strcmp (camera, "iso")) {
      gl_axis_to_quat ((float[]){0,1,0}, 3.14159265358979/4., v->quat);
      float q[4];
      gl_axis_to_quat ((float[]){1,0,0}, - 3.14159265358979/4., q);
      gl_add_quats(q, v->quat, v->quat);
    }
    else {
      fprintf (ferr, "view(): unknown camera '%s'\n", camera);
      exit (1);
    }
  }
  else if (theta || phi || psi) {
    v->gfsview = false;
    float q[4];
    gl_axis_to_quat ((float[]){1,0,0}, - phi, q);
    if (relative) {
      float q1[4];
      gl_axis_to_quat ((float[]){0,1,0}, theta, q1);
      gl_add_quats(q, q1, q1);
      float q2[4];
      gl_axis_to_quat ((float[]){0,0,1}, psi, q2);
      gl_add_quats(q1, q2, q2);
      gl_add_quats(q2, v->quat, v->quat);
    }
    else {
      gl_axis_to_quat ((float[]){0,1,0}, theta, v->quat);
      gl_add_quats(q, v->quat, v->quat);
      gl_axis_to_quat ((float[]){0,0,1}, psi, q);
      gl_add_quats(q, v->quat, v->quat);
    }
  }

  if (map)
    v->map = map;

  if (p1x || p1y || p2x || p2y) {
    float q[4];
    gl_trackball(q, p1x, p1y, p2x, p2y);
    gl_add_quats (q, v->quat, v->quat);
  }

  if (far > near) {
    v->tz = tz;
    v->far = far;
    v->near = near;
  }

  if (res)
    v->res = res;

  if ((width && width != v->width) ||
      (height && height != v->height) ||
      (samples && samples != v->samples)) {
    v->width = v->width/v->samples;
    v->height = v->height/v->samples;
    if (width) v->width = width;
    if (height) v->height = height;
    if (samples) v->samples = samples;
    v->width *= v->samples;
    v->height *= v->samples;
    framebuffer_destroy (v->fb);


    disable_fpe (FE_DIVBYZERO|FE_INVALID);

    v->fb = framebuffer_new (v->width, v->height);
    init_gl();

    enable_fpe (FE_DIVBYZERO|FE_INVALID);
  }

  if (cache > 0) {
    v->cache = pcalloc (1, sizeof (cexpr),__func__,__FILE__,0);
    v->maxlen = cache;
  }

  clear();
}







void begin_translate (float x, float y, float z)
{
  bview * view = draw();
  glMatrixMode (0x1700);
  glPushMatrix();
  glTranslatef (x, y, z);
  gl_get_frustum (&view->frustum);
}

void end_translate()
{
  bview * view = draw();
  glMatrixMode (0x1700);
  glPopMatrix();
  gl_get_frustum (&view->frustum);
}
# 238 "/home/spencer/basilisk/src/draw.h"
void begin_mirror (coord n, double alpha)
{
  bview * view = draw();
  glMatrixMode (0x1700);
  glPushMatrix();
  normalize (&n);
  GLfloat s[16], t[16];
  s[0] = 1. - 2.*n.x*n.x;
  s[1] = - 2.*n.x*n.y; s[2] = - 2.*n.x*n.z;
  s[3] = 0.;
  s[4] = s[1];
  s[5] = 1. - 2.*n.y*n.y; s[6] = - 2.*n.y*n.z;
  s[7] = 0.;
  s[8] = s[2]; s[9] = s[6]; s[10] = 1. - 2.*n.z*n.z;
  s[11] = 0.;
  s[12] = 0.; s[13] = 0.; s[14] = 0.;
  s[15] = 1.;

  t[0] = 1.; t[1] = 0.; t[2] = 0.; t[3] = 0.;
  t[4] = 0.; t[5] = 1.; t[6] = 0.; t[7] = 0.;
  t[8] = 0.; t[9] = 0.; t[10] = 1.; t[11] = 0.;
  t[12] = - 2.*n.x*alpha;
  t[13] = - 2.*n.y*alpha;
  t[14] = - 2.*n.z*alpha;
  t[15] = 1.;
  matrix_multiply (s, t);
  glMultMatrixf (s);
  gl_get_frustum (&view->frustum);
  view->reversed = !view->reversed;
}

void end_mirror() {
  end_translate();
  bview * view = draw();
  view->reversed = !view->reversed;
}







static void mapped_position (bview * view, coord * p, double * r)
{
  double x = p->x, y = p->y, z = p->z, rm = 0.;
  view->map (p);
  for (int i = -1; i <= 1; i += 2)
    for (int j = -1; j <= 1; j += 2)
      for (int k = -1; k <= 1; k += 2) {
 coord q = {x + i**r, y + j**r, z + k**r};
 view->map (&q);
 double pq = sq(p->x - q.x) + sq(p->y - q.y) + sq(p->z - q.z);
 if (pq > rm)
   rm = pq;
      }
  *r = sqrt (rm);
}

@def foreach_visible(view)
foreach_cell() {

  double _r = Delta*0.71;



  coord _p = {x, y, z};
  if ((view)->map)
    mapped_position (view, &_p, &_r);
  if (VertexBuffer.visible &&
      !sphere_in_frustum (_p.x, _p.y, _p.z, _r, &(view)->frustum))
    continue;
  if (is_leaf(cell) ||
      (VertexBuffer.visible &&
       sphere_diameter (_p.x, _p.y, _p.z, _r/L0, &(view)->frustum)
       < (view)->res)) {
    if (is_active(cell) && is_local(cell)) {
@
@def end_foreach_visible()
    }
    continue;
  }
}
end_foreach_cell();
@
# 373 "/home/spencer/basilisk/src/draw.h"
static bool _reversed = false;

static void begin_draw_lines (bview * view, float color[3], float lw)
{
  glMatrixMode (0x1701);
  glPushMatrix();
  glTranslatef (0., 0., view->lc*view->fov/24.);
  vertex_buffer_glColor3f (color[0], color[1], color[2]);
  glLineWidth (view->samples*(lw > 0. ? lw : 1.));
  _reversed = view->reversed;
  view->reversed = false;
}

static void end_draw_lines()
{
  glMatrixMode (0x1701);
  glPopMatrix();
  bview * view = draw();
  view->reversed = _reversed;
}

static inline double interp (Point point, coord p, scalar col) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  return interpolate_linear (point, col
,
        
#line 396
x + p.x*Delta, y + p.y*Delta, z + p.z*Delta);
}

static double evaluate_expression (Point point, Node * n)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (!(n)) qassert ("/home/spencer/basilisk/src/draw.h", 0, "n");
  switch (n->type) {
  case '1': return n->d.value;
  case '+': return (evaluate_expression (point, n->e[0]) +
      evaluate_expression(point, n->e[1]));
  case '-': return (evaluate_expression (point, n->e[0]) -
      evaluate_expression(point, n->e[1]));
  case '*': return (evaluate_expression (point, n->e[0]) *
      evaluate_expression(point, n->e[1]));
  case '/': return (evaluate_expression (point, n->e[0]) /
      evaluate_expression(point, n->e[1]));
  case '^': return pow (evaluate_expression (point, n->e[0]),
   evaluate_expression(point, n->e[1]));
  case '>': return (evaluate_expression (point, n->e[0]) >
      evaluate_expression(point, n->e[1]));
  case '<': return (evaluate_expression (point, n->e[0]) <
      evaluate_expression(point, n->e[1]));
  case 'L': return (evaluate_expression (point, n->e[0]) <=
      evaluate_expression(point, n->e[1]));
  case 'G': return (evaluate_expression (point, n->e[0]) >=
      evaluate_expression(point, n->e[1]));
  case '=': return (evaluate_expression (point, n->e[0]) ==
      evaluate_expression(point, n->e[1]));
  case 'i': return (evaluate_expression (point, n->e[0]) !=
      evaluate_expression(point, n->e[1]));
  case 'O': return (evaluate_expression (point, n->e[0]) ||
      evaluate_expression(point, n->e[1]));
  case 'A': return (evaluate_expression (point, n->e[0]) &&
      evaluate_expression(point, n->e[1]));
  case '?': return (evaluate_expression (point, n->e[0]) ?
      evaluate_expression(point, n->e[1]) :
      evaluate_expression(point, n->e[2]));
  case 'm': return - evaluate_expression (point, n->e[0]);
  case 'f': return n->d.func (evaluate_expression (point, n->e[0]));
  case 'v': {
    scalar s = {n->s};
    int k[3] = {0,0,0};
    for (int i = 0; i < 3; i++)
      if (n->e[i])
 k[i] = evaluate_expression (point, n->e[i]);
    return val(s,k[0],k[1],k[2]);
  }
  case 'D': return Delta;
  case 'x': return x;
  case 'y': return y;
  case 'z': return z;
  default:
    fprintf (ferr, "unknown operation type '%c'\n", n->type);
    if (!(false)) qassert ("/home/spencer/basilisk/src/draw.h", 0, "false");
  }
  return undefined;
}


#line 399
static void _stencil_evaluate_expression (Point point, Node * n)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES; 
      
  switch (n->type) {
  case '1': return ;
  case '+': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 404
return  
;}
  case '-': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 406
return  
;}
  case '*': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 408
return  
;}
  case '/': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 410
return  
;}
  case '^': {_stencil_evaluate_expression (point, n->e[0]);
   _stencil_evaluate_expression(point, n->e[1]);
#line 412
return  
;}
  case '>': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 414
return  
;}
  case '<': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 416
return  
;}
  case 'L': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 418
return  
;}
  case 'G': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 420
return  
;}
  case '=': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 422
return  
;}
  case 'i': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 424
return  
;}
  case 'O': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 426
return  
;}
  case 'A': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 428
return  
;}
  case '?': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
      _stencil_evaluate_expression(point, n->e[2]);
#line 430
return   

;}
  case 'm': { _stencil_evaluate_expression (point, n->e[0]);return ;}
  case 'f': {_stencil_evaluate_expression (point, n->e[0]);return  ;}
  case 'v': {
    scalar s = {n->s};   
    
    for (int i = 0; i < 3; i++)
      if (n->e[i])
 { _stencil_evaluate_expression (point, n->e[i]); } 
_stencil_val(s,o_stencil,o_stencil,o_stencil);
    
#line 441
return;
  }
  case 'D': return ;
  case 'x': return ;
  case 'y': return ;
  case 'z': return ; 
     
    
        
  }
  return ;
}

static bool assemble_node (Node * n)
{
  if (n->type == 'v') {
    char * id = n->d.id;
    scalar s = lookup_field (id);
    if (s.i >= 0)
      n->s = s.i;
    else {
      n->s = -1;
      if (!strcmp (id, "Delta"))
 reset_node_type (n, 'D');
      else if (!strcmp (id, "x"))
 reset_node_type (n, 'x');
      else if (!strcmp (id, "y"))
 reset_node_type (n, 'y');
      else if (!strcmp (id, "z"))
 reset_node_type (n, 'z');
      else {
 typedef struct { char * name; double val; } Constant;
 static Constant constants[] = {
   {"pi", 3.14159265358979 },
   {"nodata", 1e30 },
   {"HUGE", 1e30 },
   { NULL },
 };
 Constant * p = constants;
 while (p->name) {
   if (!strcmp (p->name, id)) {
     reset_node_type (n, '1');
     n->d.value = p->val;
     break;
   }
   p++;
 }
 if (n->type == 'v') {
   fprintf (ferr, "unknown identifier '%s'\n", id);
   return false;
 }
      }
    }
  }
  for (int i = 0; i < 3; i++)
    if (n->e[i] && !assemble_node (n->e[i]))
      return false;
  return true;
}

static scalar compile_expression (char * expr, bool * isexpr)
{
  *isexpr = false;
  if (!expr)
    return (scalar){-1};

  bview * view = get_view();
  scalar s;
  if (view->cache && (s = get_cexpr (view->cache, expr)).i >= 0)
    return s;

  Node * node = parse_node (expr);
  if (node == NULL) {
    fprintf (ferr, "'%s': syntax error\n", expr);
    return (scalar){-1};
  }
  if (!assemble_node (node)) {
    free_node (node);
    return (scalar){-1};
  }
  if (node->type == 'v' && node->e[0] == NULL) {
    scalar s = {node->s};
    if (_attribute[s.i].block > 0) {
      free_node (node);
      return s;
    }
  }
  s = new_scalar("s");
  pfree (_attribute[s.i].name,__func__,__FILE__,0);
  _attribute[s.i].name = pstrdup (expr,__func__,__FILE__,0);
  foreach_stencil(1,{(NonLocal[]){{"node","not implemented yet",(void *)node,NULL,1},{"undefined","double",(void *)&undefined,NULL,0},{"ferr","not implemented yet",(void *)ferr,NULL,1},{"s","scalar",(void *)&s,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_("\n@ def not_mpi_compatible()\ndo {\n  if (npe() > 1) {\n    fprintf (ferr, \"%s() is not compatible with MPI (yet)\\n\", __func__);\n    exit (1);\n  }\n} while(0)\n@\n@ define system(command) (pid() == 0 ? system(command) : 0)\n@else\n@ define qstderr() stderr\n@ define qstdout() stdout\n@ define ferr stderr\n@ define fout stdout\n@ define not_mpi_compatible()\n@endif\n\n\n\n\n// #line 93 \"/home/spencer/basilisk/src/common.h\"\nstatic inline void qassert (const char * file, int line, const char * cond) {\n  fprintf (ferr, \"%s:%d: Assertion `%s' failed.\\n\", file, line, cond);\n  abort();\n}\n\n\n// #line 399 \"/home/spencer/basilisk/src/draw.h\"\nstatic real evaluate_expression (Point point, Node * n)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  if (!(n)) qassert (\"/home/spencer/basilisk/src/draw.h\", 0, \"n\");\n  switch (n->type) {\n  case '1': return n->d.value;\n  case '+': return (evaluate_expression (point, n->e[0]) +\n      evaluate_expression(point, n->e[1]));\n  case '-': return (evaluate_expression (point, n->e[0]) -\n      evaluate_expression(point, n->e[1]));\n  case '*': return (evaluate_expression (point, n->e[0]) *\n      evaluate_expression(point, n->e[1]));\n  case '/': return (evaluate_expression (point, n->e[0]) /\n      evaluate_expression(point, n->e[1]));\n  case '^': return pow (evaluate_expression (point, n->e[0]),\n   evaluate_expression(point, n->e[1]));\n  case '>': return (evaluate_expression (point, n->e[0]) >\n      evaluate_expression(point, n->e[1]));\n  case '<': return (evaluate_expression (point, n->e[0]) <\n      evaluate_expression(point, n->e[1]));\n  case 'L': return (evaluate_expression (point, n->e[0]) <=\n      evaluate_expression(point, n->e[1]));\n  case 'G': return (evaluate_expression (point, n->e[0]) >=\n      evaluate_expression(point, n->e[1]));\n  case '=': return (evaluate_expression (point, n->e[0]) ==\n      evaluate_expression(point, n->e[1]));\n  case 'i': return (evaluate_expression (point, n->e[0]) !=\n      evaluate_expression(point, n->e[1]));\n  case 'O': return (evaluate_expression (point, n->e[0]) ||\n      evaluate_expression(point, n->e[1]));\n  case 'A': return (evaluate_expression (point, n->e[0]) &&\n      evaluate_expression(point, n->e[1]));\n  case '?': return (evaluate_expression (point, n->e[0]) ?\n      evaluate_expression(point, n->e[1]) :\n      evaluate_expression(point, n->e[2]));\n  case 'm': return - evaluate_expression (point, n->e[0]);\n  case 'f': return n->d.func (evaluate_expression (point, n->e[0]));\n  case 'v': {\n    scalar s = {n->s};\n    int k[3] = {0,0,0};\n    for (int i = 0; i < 3; i++)\n      if (n->e[i])\n k[i] = evaluate_expression (point, n->e[i]);\n    return val(s,k[0],k[1],k[2]);\n  }\n  case 'D': return Delta;\n  case 'x': return x;\n  case 'y': return y;\n  case 'z': return z;\n  default:\n    fprintf (ferr, \"unknown operation type '%c'\\n\", n->type);\n    if (!(false)) qassert (\"/home/spencer/basilisk/src/draw.h\", 0, \"false\");\n  }\n  return undefined;\n}"),_("\n    \n// #line 532 \"/home/spencer/basilisk/src/draw.h\"\nval_out_(s,0,0,0) = evaluate_expression (point, node);")})
    { _stencil_evaluate_expression (point, node);_stencil_val_a(s,0,0,0); }end_foreach_stencil();
  {
#line 531
foreach()
    val(s,0,0,0) = evaluate_expression (point, node);end_foreach();}
  restriction (((scalar[]){s,{-1}}));
  free_node (node);

  if (view->cache)
    view->cache = add_cexpr (view->cache, view->maxlen, expr, s);
  else
    *isexpr = true;
  return s;
}
# 604 "/home/spencer/basilisk/src/draw.h"
static void begin_colorized (float fc[3], bool constant_color,
        double cmap[127][3], bool use_texture)
{

  if (use_texture) {
    GLfloat texture[3*256];
    for (int i = 0; i < 256; i++) {
      Color j = colormap_color (cmap, i/255., 0, 1);
      texture[3*i] = j.r/255.;
      texture[3*i + 1] = j.g/255.;
      texture[3*i + 2] = j.b/255.;
    }
    glTexImage1D (0x0DE0, 0, 0x1907, 256,0, 0x1907, 0x1406, texture);
    glTexParameteri (0x0DE0, 0x2801, 0x2601);
    glTexParameteri (0x0DE0, 0x2800, 0x2601);
    glTexParameteri (0x0DE0, 0x2802, 0x812F);
    glTexParameteri (0x0DE0, 0x2803, 0x812F);
    glEnable (0x0DE0);
  }
  if (constant_color)
    vertex_buffer_glColor3f (fc[0], fc[1], fc[2]);
}

static void end_colorized() {
  glDisable (0x0DE0);
}
# 656 "/home/spencer/basilisk/src/draw.h"
     
bool colorbar (Colormap map, float size, float pos[2],
        char * label, double lscale, double min,
        double max, bool horizontal, bool border,
        bool mid, float lc[3], float lw, float fsize,
        char * format, int levels)
{tracing("colorbar","/home/spencer/basilisk/src/draw.h",0);
  bview * view = draw();
  glDisable (0x0B50);
  glMatrixMode (0x1701);
  glPushMatrix();
  glLoadIdentity();
  glMatrixMode (0x1700);
  glPushMatrix();
  glLoadIdentity();

  float fheight = gl_StrokeHeight();
  if (!size)
    size = 15;
  float width = 2./size;
  if (levels < 1) levels = 1;
  float h = 0, height = 4*width, dh = height/levels;
  glTranslatef (pos[0], pos[1], 0);


  double cmap [127][3];
  (* map) (cmap);
  vertex_buffer_glBegin(0x0007);
  for (int i = 0; i < levels; i++) {
    Color c = colormap_color (cmap, (float)i/(levels - 1), 0, 1);
    vertex_buffer_glColor3f (c.r/255., c.g/255., c.b/255.);
    if (horizontal) {
      vertex_buffer_glVertex2d (h + dh, 0);
      vertex_buffer_glVertex2d (h + dh, width);
      vertex_buffer_glVertex2d (h, width);
      vertex_buffer_glVertex2d (h, 0);
    } else {
      vertex_buffer_glVertex2d (0, h + dh);
      vertex_buffer_glVertex2d (width, h + dh);
      vertex_buffer_glVertex2d (width, h);
      vertex_buffer_glVertex2d (0, h);
    }
    h += dh;
    view->ni++;
  }
  vertex_buffer_glEnd();
  glLineWidth (view->samples*(lw > 0. ? lw : 1.));
  vertex_buffer_glColor3f (lc[0], lc[1], lc[2]);


  if (border) {
    vertex_buffer_glBegin (0x0002);
    vertex_buffer_glVertex2d (0, 0);
    if (horizontal) {
      vertex_buffer_glVertex2d (0, width);
      vertex_buffer_glVertex2d (height, width);
      vertex_buffer_glVertex2d (height, 0);
    } else {
      vertex_buffer_glVertex2d (width, 0);
      vertex_buffer_glVertex2d (width, height);
      vertex_buffer_glVertex2d (0, height);
    }
    vertex_buffer_glEnd();
  }


  float fwidth = gl_StrokeWidth ('1');
  if (!fsize)
    fsize = 20;
  float hscale = 2./(fsize*fwidth), vscale = hscale*view->width/view->height;
  char str[99];
  vertex_buffer_glColor3f (lc[0], lc[1], lc[2]);
  if (horizontal)
    glTranslatef (0, -(fheight/(view->height)), 0);
  else
    glTranslatef (width, -(fheight/(3*view->height)), 0);
  glScalef (hscale, vscale, 1.);
  sprintf (str, format, min);
  if (min > -1e30) {
    glPushMatrix();
    if (horizontal)
      glTranslatef (-fwidth*(strlen(str) - 1)/2, 0, 0);
    glScalef (lscale, lscale, 1.);
    gl_StrokeString (str);
    glPopMatrix();
  }
  if (horizontal)
    glTranslatef (height/hscale,0, 0);
  else
    glTranslatef (0, height/vscale, 0);
  sprintf (str, format, max);
  if (max < 1e30) {
    glPushMatrix();
    if (horizontal)
      glTranslatef (-fwidth*(strlen(str) - 1)/2, 0, 0);
    glScalef (lscale, lscale, 1.);
    gl_StrokeString (str);
    glPopMatrix();
  }

  if (mid) {
    sprintf (str, format, (min + max)/2);
    glPushMatrix();
    if (horizontal)
      glTranslatef (-height/(2*hscale) - fwidth*(strlen(str) - 1)/2,0, 0);
    else
      glTranslatef (0, -height/(2*vscale), 0);
    glScalef (lscale, lscale, 1.);
    gl_StrokeString (str);
    glPopMatrix();
  }

  if (horizontal)
    glTranslatef (-height/(2*hscale) - lscale*fwidth*(strlen(label) - 1)/2, width/vscale, 0);
  else
    glTranslatef (-width/hscale, 0, 0);

  glScalef (lscale, lscale, 1.);
  glTranslatef (0, fheight, 0);
  gl_StrokeString (label);

  glMatrixMode (0x1700);
  glPopMatrix();
  glMatrixMode (0x1701);
  glPopMatrix();
  {end_tracing("colorbar","/home/spencer/basilisk/src/draw.h",0);return true;}
end_tracing("colorbar","/home/spencer/basilisk/src/draw.h",0);}
# 801 "/home/spencer/basilisk/src/draw.h"
static bool cfilter (Point point, scalar c, double cmin)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  double cmin1 = 4.*cmin;
  if (val(c,0,0,0) <= cmin) {
    
      if (val(c,1,0,0) >= 1. - cmin1 || val(c,-1,0,0) >= 1. - cmin1)
 return true;
      
#line 806
if (val(c,0,1,0) >= 1. - cmin1 || val(c,0,-1,0) >= 1. - cmin1)
 return true;
    return false;
  }
  if (val(c,0,0,0) >= 1. - cmin) {
    
      if (val(c,1,0,0) <= cmin1 || val(c,-1,0,0) <= cmin1)
 return true;
      
#line 812
if (val(c,0,1,0) <= cmin1 || val(c,0,-1,0) <= cmin1)
 return true;
    return false;
  }
  int n = 0;
  double min = 1e30, max = - 1e30;
  {foreach_neighbor(1) {
    if (val(c,0,0,0) > cmin && val(c,0,0,0) < 1. - cmin && ++n >= (1 << 2))
      return true;
    if (val(c,0,0,0) > max) max = val(c,0,0,0);
    if (val(c,0,0,0) < min) min = val(c,0,0,0);
  }end_foreach_neighbor()}
  return max - min > 0.5;
}
# 801 "/home/spencer/basilisk/src/draw.h"
static void _stencil_cfilter (Point point, scalar c,_stencil_undefined * cmin)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;   
  
_stencil_val(c,0,0,0); {
    
      {_stencil_val(c,1,0,0); _stencil_val(c,-1,0,0); 
           }
      
#line 806
{_stencil_val(c,0,1,0); _stencil_val(c,0,-1,0); 
           } 
    
  }
_stencil_val(c,0,0,0); {
    
      {_stencil_val(c,1,0,0); _stencil_val(c,-1,0,0); 
       }
      
#line 812
{_stencil_val(c,0,1,0); _stencil_val(c,0,-1,0); 
       } 
    
  }          
     
       
  
  
  {
#line 818
foreach_neighbor(1) {
_stencil_val(c,0,0,0); _stencil_val(c,0,0,0);

_stencil_val(c,0,0,0); { _stencil_val(c,0,0,0); }
_stencil_val(c,0,0,0); { _stencil_val(c,0,0,0); } 
      
                  
       
       
  
#line 823
}end_foreach_neighbor()}
  return     ;
}

static void glvertex3d (bview * view, double x, double y, double z) {
  if (view->map) {
    coord p = {x, y, z};
    view->map (&p);
    vertex_buffer_glVertex3d (p.x, p.y, p.z);
  }
  else
    vertex_buffer_glVertex3d (x, y, z);
}


static void glvertex2d (bview * view, double x, double y) {
  if (view->map) {
    coord p = {x, y, 0.};
    view->map (&p);
    vertex_buffer_glVertex2d (p.x, p.y);
  }
  else
    vertex_buffer_glVertex2d (x, y);
}

static void glvertex_normal3d (bview * view, Point point, vector n,
          double xp, double yp, double zp)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  coord v = {(xp - x)/Delta, (yp - y)/Delta}, np;
  
    np.x = - interp (point, v, n.x);
    
#line 853
np.y = - interp (point, v, n.y);
  vertex_buffer_glNormal3d (np.x, np.y, 1.);
  glvertex3d (view, xp, yp, zp);
}
# 885 "/home/spencer/basilisk/src/draw.h"
     
bool draw_vof (char * c, char * s, bool edges,
        double larger, int filled,
        char * color,
        double min, double max, double spread,
        bool linear,
        Colormap map,
        float fc[3], float lc[3], float lw,
        bool expr,
        bool cbar, float size, float pos[2], char * label, double lscale, bool horizontal, bool border, bool mid, float fsize, char * format, int levels)
{tracing("draw_vof","/home/spencer/basilisk/src/draw.h",0);
  scalar d = lookup_field (c);
  if (d.i < 0) {
    fprintf (ferr, "draw_vof(): no field named '%s'\n", c);
    {end_tracing("draw_vof","/home/spencer/basilisk/src/draw.h",0);return false;}
  }
  vector fs = lookup_vector (s);

  scalar col = {-1}; if (color && strcmp (color, "level")) { col = compile_expression (color, &expr); if (col.i < 0) {end_tracing("draw_vof","/home/spencer/basilisk/src/draw.h",0);return false;} } double cmap[127][3]; if (color) { if (min == 0 && max == 0) { if (col.i < 0) min = 0, max = depth(); else { stats s = statsf (col); double avg = s.sum/s.volume; if (spread < 0.) min = s.min, max = s.max; else { if (!spread) spread = 5.; min = avg - spread*s.stddev; max = avg + spread*s.stddev; } } } if (!map) map = jet; (* map) (cmap); } if ((2 > 2 || linear) && !fc[0] && !fc[1] && !fc[2]) fc[0] = fc[1] = fc[2] = 1.; if (cbar) colorbar (map, size, pos, label, lscale, min, max, horizontal, border, mid, lc, lw, fsize, format, levels);;

  double cmin = 1e-3;



  void (* prolongation) (Point, scalar) = _attribute[d.i].prolongation;
  if (prolongation != fraction_refine) {
    _attribute[d.i].prolongation = fraction_refine;
    _attribute[d.i].dirty = true;
  }


  bview * view = draw();

  if (filled) {
    vertex_buffer_glColor3f (fc[0], fc[1], fc[2]);
    vertex_buffer_glNormal3d (0, 0, view->reversed ? -1 : 1);
    foreach_visible_stencil (1,{(NonLocal[]){{"fs","vector",(void *)&fs,NULL,0},{"baseblock","scalar",(void *)baseblock,NULL,1},{"all","scalar",(void *)all,NULL,1},{"max","double",(void *)&max,NULL,0},{"free_solver_funcs","Array",(void *)free_solver_funcs,NULL,1},{"pmtrace","not implemented yet",(void *)&pmtrace,NULL,0},{"pmfuncn","int",(void *)&pmfuncn,NULL,0},{"pmfuncs","pmfunc",(void *)pmfuncs,NULL,1},{"ferr","not implemented yet",(void *)ferr,NULL,1},{"_view","not implemented yet",(void *)_view,NULL,1},{"VertexBuffer","not implemented yet",(void *)&VertexBuffer,NULL,0},{"d","scalar",(void *)&d,NULL,0},{"filled","int",(void *)&filled,NULL,0},{"view","not implemented yet",(void *)view,NULL,1},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_("\n// # 415 \"/home/spencer/basilisk/src/view.h\" 2\n\n\n\n\n\n\n// # 1 \"./draw.h\" 1\n// # 1 \"/home/spencer/basilisk/src/draw.h\"\n\n\n\n\n// # 1 \"./fractions.h\" 1\n// # 1 \"/home/spencer/basilisk/src/fractions.h\"\n// # 12 \"/home/spencer/basilisk/src/fractions.h\"\n// # 1 \"./geometry.h\" 1\n// # 1 \"/home/spencer/basilisk/src/geometry.h\"\n// # 35 \"/home/spencer/basilisk/src/geometry.h\"\nreal line_alpha (real c, coord n)\n{\n  real alpha, n1, n2;\n\n  n1 = fabs (n.x); n2 = fabs (n.y);\n  if (n1 > n2)\n    do { real __tmp = n1; n1 = n2; n2 = __tmp; } while(0);\n\n  c = clamp (c, 0., 1.);\n  real v1 = n1/2.;\n  if (c <= v1/n2)\n    alpha = sqrt (2.*c*n1*n2);\n  else if (c <= 1. - v1/n2)\n    alpha = c*n2 + v1;\n  else\n    alpha = n1 + n2 - sqrt (2.*n1*n2*(1. - c));\n\n  if (n.x < 0.)\n    alpha += n.x;\n  if (n.y < 0.)\n    alpha += n.y;\n\n  return alpha - (n.x + n.y)/2.;\n}\n// # 13 \"/home/spencer/basilisk/src/fractions.h\" 2\n\n\n\n\n\n// # 1 \"./myc2d.h\" 1\n// # 1 \"/home/spencer/basilisk/src/myc2d.h\"\n\n\n\n\n\ncoord mycs (Point point, scalar c)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int ix;\n  real c_t,c_b,c_r,c_l;\n  real mx0,my0,mx1,my1,mm1,mm2;\n\n\n  c_t = val(c,-1,1,0) + val(c,0,1,0) + val(c,1,1,0);\n  c_b = val(c,-1,-1,0) + val(c,0,-1,0) + val(c,1,-1,0);\n  c_r = val(c,1,-1,0) + val(c,1,0,0) + val(c,1,1,0);\n  c_l = val(c,-1,-1,0) + val(c,-1,0,0) + val(c,-1,1,0);\n\n\n\n  mx0 = 0.5*(c_l - c_r);\n  my0 = 0.5*(c_b - c_t);\n\n\n  if (fabs(mx0) <= fabs(my0)) {\n    my0 = my0 > 0. ? 1. : -1.;\n    ix = 1;\n  }\n  else {\n    mx0 = mx0 > 0. ? 1. : -1.;\n    ix = 0;\n  }\n\n\n  mm1 = val(c,-1,-1,0) + 2.0*val(c,-1,0,0) + val(c,-1,1,0);\n  mm2 = val(c,1,-1,0) + 2.0*val(c,1,0,0) + val(c,1,1,0);\n  mx1 = mm1 - mm2 + 1e-30;\n  mm1 = val(c,-1,-1,0) + 2.0*val(c,0,-1,0) + val(c,1,-1,0);\n  mm2 = val(c,-1,1,0) + 2.0*val(c,0,1,0) + val(c,1,1,0);\n  my1 = mm1 - mm2 + 1e-30;\n\n\n  if (ix) {\n    mm1 = fabs(my1);\n    mm1 = fabs(mx1)/mm1;\n    if (mm1 > fabs(mx0)) {\n      mx0 = mx1;\n      my0 = my1;\n    }\n  }\n  else {\n    mm1 = fabs(mx1);\n    mm1 = fabs(my1)/mm1;\n    if (mm1 > fabs(my0)) {\n      mx0 = mx1;\n      my0 = my1;\n    }\n  }\n\n\n\n  mm1 = fabs(mx0) + fabs(my0);\n  coord n = {mx0/mm1, my0/mm1};\n\n  return n;\n}\n\n\n\n\n\n\n// #line 418 \"/home/spencer/basilisk/src/fractions.h\"\ncoord facet_normal (Point point, scalar c, vector s)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  if (s.x.i >= 0) {\n    coord n;\n    real nn = 0.;\n     {\n      n.x = val(s.x,0,0,0) - val(s.x,1,0,0);\n      nn += fabs(n.x);\n    } \n// #line 423\n{\n      n.y = val(s.y,0,0,0) - val(s.y,0,1,0);\n      nn += fabs(n.y);\n    }\n    if (nn > 0.)\n      {\n n.x /= nn;\n \n// #line 429\nn.y /= nn;}\n    else\n      {\n n.x = 1./2;\n \n// #line 432\nn.y = 1./2;}\n    return n;\n  }\n  return mycs (point, c);\n}\n\n\n// #line 32 \"/home/spencer/basilisk/src/vertexbuffer.h\"\nstatic void vertex_buffer_push_index (unsigned int i)\n{\n  i -= VertexBuffer.vertex;\n  array_append (VertexBuffer.index, &i, sizeof(unsigned int));\n}\n\n\n// #line 112 \"/home/spencer/basilisk/src/vertexbuffer.h\"\nstatic void vertex_buffer_glEnd()\n{\n  if (VertexBuffer.index) {\n    int type = -1;\n    switch (VertexBuffer.state) {\n\n    case 0x0002:\n      for (int i = VertexBuffer.line_loop; i < VertexBuffer.nvertex - 1; i++) {\n vertex_buffer_push_index (i);\n vertex_buffer_push_index (i + 1);\n      }\n      vertex_buffer_push_index (VertexBuffer.nvertex - 1);\n      vertex_buffer_push_index (VertexBuffer.line_loop);\n      type = 0;\n      break;\n\n    case 0x0001:\n      for (int i = VertexBuffer.lines; i < VertexBuffer.nvertex; i += 2) {\n vertex_buffer_push_index (i);\n vertex_buffer_push_index (i + 1);\n      }\n      type = 0;\n      break;\n\n    case 0x0003:\n      for (int i = VertexBuffer.line_strip; i < VertexBuffer.nvertex - 1; i++) {\n vertex_buffer_push_index (i);\n vertex_buffer_push_index (i + 1);\n      }\n      type = 0;\n      break;\n\n    case 0x0007:\n      for (int i = VertexBuffer.quads; i < VertexBuffer.nvertex; i += 4)\n for (int j = 1; j <= 2; j++) {\n   vertex_buffer_push_index (i);\n   vertex_buffer_push_index (i + j);\n   vertex_buffer_push_index (i + j + 1);\n }\n      type = 1;\n      break;\n\n    case 0x0009:\n      for (int j = 1; j <= VertexBuffer.nvertex - VertexBuffer.polygon - 2;\n    j++) {\n vertex_buffer_push_index (VertexBuffer.polygon);\n vertex_buffer_push_index (VertexBuffer.polygon + j);\n vertex_buffer_push_index (VertexBuffer.polygon + j + 1);\n      }\n      type = 1;\n      break;\n\n    case 0x0006:\n      for (int i = VertexBuffer.fan + 1; i < VertexBuffer.nvertex - 1; i++) {\n vertex_buffer_push_index (VertexBuffer.fan);\n vertex_buffer_push_index (i);\n vertex_buffer_push_index (i + 1);\n      }\n      type = 1;\n      break;\n\n    default:\n      break;\n    }\n    VertexBuffer.state = 0;\n    if (VertexBuffer.type >= 0 && type >= 0) {\n\n      if (!(VertexBuffer.type == type)) qassert (\"/home/spencer/basilisk/src/vertexbuffer.h\", 0, \"VertexBuffer.type == type\");\n    }\n    else\n      VertexBuffer.type = type;\n  }\n  else\n    glEnd();\n}\n\n\n// #line 222 \"/home/spencer/basilisk/src/vertexbuffer.h\"\nstatic void vertex_buffer_glVertex2d (real x, real y)\n{\n  if (VertexBuffer.position) {\n    if (VertexBuffer.dim < 2)\n      VertexBuffer.dim = 2;\n    real v[4] = {x, y, 0, 1.};\n    vector_multiply (v, VertexBuffer.modelview);\n    array_append (VertexBuffer.position, v, 3*sizeof(real));\n    VertexBuffer.nvertex++;\n  }\n  else\n    glVertex3d (x, y, 0.);\n}\n\n\n\n// #line 838 \"/home/spencer/basilisk/src/draw.h\"\nstatic void glvertex2d (bview * view, real x, real y) {\n  if (view->map) {\n    coord p = {x, y, 0.};\n    view->map (&p);\n    vertex_buffer_glVertex2d (p.x, p.y);\n  }\n  else\n    vertex_buffer_glVertex2d (x, y);\n}\n\n\n// #line 200 \"/home/spencer/basilisk/src/common.h\"\nstatic void * pmfunc_free (void * ptr, char c)\n{\n  if (!ptr)\n    return ptr;\n  pmdata * d = (pmdata *) (((char *)ptr) - sizeof(pmdata));\n  if (d->id < 1 || d->id > pmfuncn) {\n    fputs (\"*** MTRACE: ERROR!: corrupted free()\", ferr);\n    if (d->size == 0)\n      fputs (\", possible double free()\", ferr);\n    else\n      fputs (\", not traced?\", ferr);\n    fputs (\", aborting...\\n\", ferr);\n    abort();\n    return ptr;\n  }\n  else\n  OMP (omp critical)\n  {\n    pmfunc * f = &pmfuncs[d->id - 1];\n    if (f->total < d->size) {\n      fprintf (ferr, \"*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\\n\",\n        f->total, d->size);\n      abort();\n    }\n    else\n      f->total -= d->size;\n    if (pmtrace.total < d->size) {\n      fprintf (ferr, \"*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\\n\",\n        pmtrace.total, d->size);\n      abort();\n    }\n    else {\n      pmtrace.total -= d->size;\n      pmtrace.overhead -= sizeof(pmdata);\n    }\n    d->id = 0;\n    d->size = 0;\n    pmfunc_trace (f, c);\n  }\n  return d;\n}\n\n\n// #line 256 \"/home/spencer/basilisk/src/common.h\"\nstatic void * prealloc (void * ptr, size_t size,\n   const char * func, const char * file, int line)\n{\n  return pmfunc_alloc ((pmdata *) sysrealloc (pmfunc_free(ptr, '<'),\n           sizeof(pmdata) + size),\n         size, func, file, line, '>');\n}\n\n\n// #line 394 \"/home/spencer/basilisk/src/common.h\"\nvoid array_append (Array * a, void * elem, size_t size)\n{\n  if (a->len + size >= a->max) {\n    a->max += max (size, 4096);\n    a->p = prealloc (a->p, a->max,__func__,__FILE__,0);\n  }\n  memcpy (((char *)a->p) + a->len, elem, size);\n  a->len += size;\n}\n\n\n// #line 380 \"/home/spencer/basilisk/src/common.h\"\nArray * array_new()\n{\n  Array * a = ((Array *) pmalloc ((1)*sizeof(Array),__func__,__FILE__,0));\n  a->p = NULL;\n  a->max = a->len = 0;\n  return a;\n}\n\n\n// #line 1410 \"/home/spencer/basilisk/src/common.h\"\nvoid free_solver_func_add (free_solver_func func)\n{\n  if (!free_solver_funcs)\n    free_solver_funcs = array_new();\n  array_append (free_solver_funcs, &func, sizeof(free_solver_func));\n}\n\n\n// #line 154 \"/home/spencer/basilisk/src/common.h\"\nstatic void pmfunc_trace (pmfunc * f, char c)\n{\n  if (pmtrace.fp)\n    fprintf (pmtrace.fp, \"%c %d %ld %ld %ld\",\n      c, f->id, pmtrace.nr, pmtrace.total, f->total);\n@if 1\n  if (pmtrace.nr % 1 == 0) {\n    struct rusage usage;\n    getrusage (RUSAGE_SELF, &usage);\n    if (pmtrace.fp)\n      fprintf (pmtrace.fp, \" %ld\", usage.ru_maxrss*1024);\n    if (!pmtrace.nr)\n      pmtrace.startrss = usage.ru_maxrss;\n    if (usage.ru_maxrss > pmtrace.maxrss)\n      pmtrace.maxrss = usage.ru_maxrss;\n  }\n@endif\n  if (pmtrace.fp)\n    fputc ('\\n', pmtrace.fp);\n  pmtrace.nr++;\n}\n\n\n// #line 135 \"/home/spencer/basilisk/src/common.h\"\nstatic int pmfunc_index (const char * func, const char * file, int line)\n{\n  pmfunc * p = pmfuncs;\n  for (int i = 0; i < pmfuncn; i++, p++)\n    if (p->line == line && !strcmp(func, p->func) && !strcmp(file, p->file))\n      return p->id;\n  pmfuncn++;\n  pmfuncs = (pmfunc *) sysrealloc (pmfuncs, pmfuncn*sizeof(pmfunc));\n  p = &pmfuncs[pmfuncn - 1];\n  memset (p, 0, sizeof(pmfunc));\n  p->func = systrdup(func);\n  p->file = systrdup(file);\n  p->line = line;\n  p->id = pmfuncn;\n  if (pmtrace.fp)\n    fprintf (pmtrace.fp, \"@ %d %s %s %d\\n\", pmfuncn, func, file, line);\n  return pmfuncn;\n}\n@ def not_mpi_compatible()\ndo {\n  if (npe() > 1) {\n    fprintf (ferr, \"%s() is not compatible with MPI (yet)\\n\", __func__);\n    exit (1);\n  }\n} while(0)\n@\n@ define system(command) (pid() == 0 ? system(command) : 0)\n@else\n@ define qstderr() stderr\n@ define qstdout() stdout\n@ define ferr stderr\n@ define fout stdout\n@ define not_mpi_compatible()\n@endif\n\n\n\n\n// #line 93 \"/home/spencer/basilisk/src/common.h\"\nstatic inline void qassert (const char * file, int line, const char * cond) {\n  fprintf (ferr, \"%s:%d: Assertion `%s' failed.\\n\", file, line, cond);\n  abort();\n}\n\n\n// #line 176 \"/home/spencer/basilisk/src/common.h\"\nstatic void * pmfunc_alloc (pmdata * d, size_t size,\n       const char * func, const char * file, int line,\n       char c)\n{\n  if (!(d != NULL)) qassert (\"/home/spencer/basilisk/src/common.h\", 0, \"d != NULL\");\n  OMP (omp critical)\n  {\n    d->id = pmfunc_index(func, file, line);\n    d->size = size;\n    pmfunc * f = &pmfuncs[d->id - 1];\n    f->total += size;\n    if (f->total > f->max)\n      f->max = f->total;\n    pmtrace.total += size;\n    pmtrace.overhead += sizeof(pmdata);\n    if (pmtrace.total > pmtrace.max) {\n      pmtrace.max = pmtrace.total;\n      pmtrace.maxoverhead = pmtrace.overhead;\n    }\n    pmfunc_trace (f, c);\n  }\n  return ((char *)d) + sizeof(pmdata);\n}\n\n\n// #line 242 \"/home/spencer/basilisk/src/common.h\"\nstatic void * pmalloc (size_t size,\n         const char * func, const char * file, int line)\n{\n  return pmfunc_alloc ((pmdata *) sysmalloc (sizeof(pmdata) + size),\n         size, func, file, line, '+');\n}\n\n\n// #line 249 \"/home/spencer/basilisk/src/common.h\"\nstatic void * pcalloc (size_t nmemb, size_t size,\n         const char * func, const char * file, int line)\n{\n  void * p = pmalloc (nmemb*size, func, file, line);\n  return memset (p, 0, nmemb*size);\n}\n\n\n\n\n\n// #line 171 \"/home/spencer/basilisk/src/view.h\"\nbview * bview_new()\n{\n  bview * p = ((bview *) pcalloc (1, sizeof(bview),__func__,__FILE__,0));\n\n  p->tx = p->ty = 0;\n  p->sx = p->sy = p->sz = 1.;\n  p->quat[0] = p->quat[1] = p->quat[2] = 0; p->quat[3] = 1;\n  p->fov = 24.;\n  gl_trackball (p->quat, 0.0, 0.0, 0.0, 0.0);\n\n\n  p->bg[0] = 1; p->bg[1] = 1; p->bg[2] = 1;\n\n\n\n  p->res = 1.;\n  p->lc = 0.004;\n\n  p->samples = 4;\n  p->width = 600*p->samples, p->height = 600*p->samples;\n\n\n  disable_fpe (FE_DIVBYZERO|FE_INVALID);\n\n  p->fb = framebuffer_new (p->width, p->height);\n\n  init_gl();\n  p->active = false;\n\n  enable_fpe (FE_DIVBYZERO|FE_INVALID);\n\n  return p;\n}\n\n\n// #line 232 \"/home/spencer/basilisk/src/view.h\"\nbview * get_view() {\n  if (!_view) {\n    _view = bview_new();\n    free_solver_func_add (destroy_view);\n  }\n  return _view;\n}\n\n\n// #line 61 \"/home/spencer/basilisk/src/vertexbuffer.h\"\nstatic void vertex_buffer_glBegin (unsigned int state)\n{\n  if (VertexBuffer.index) {\n\n    glGetFloatv (0x0BA6, VertexBuffer.modelview);\n\n    bview * view = get_view();\n\n    real q[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,\n      - view->tx, - view->ty, 3, 1 };\n    matrix_multiply (q, VertexBuffer.modelview);\n    for (int i = 0; i < 16; i++)\n      VertexBuffer.modelview[i] = q[i];\n\n    gl_build_rotmatrix ((real (*)[4])q, view->quat);\n    do { real __tmp = q[1]; q[1] = q[4]; q[4] = __tmp; } while(0);\n    do { real __tmp = q[2]; q[2] = q[8]; q[8] = __tmp; } while(0);\n    do { real __tmp = q[6]; q[6] = q[9]; q[9] = __tmp; } while(0);\n    matrix_multiply (q, VertexBuffer.modelview);\n    for (int i = 0; i < 16; i++)\n      VertexBuffer.modelview[i] = q[i];\n\n    VertexBuffer.state = state;\n    switch (state) {\n    case 0x0002:\n      VertexBuffer.line_loop = VertexBuffer.nvertex;\n      break;\n    case 0x0001:\n      VertexBuffer.lines = VertexBuffer.nvertex;\n      break;\n    case 0x0003:\n      VertexBuffer.line_strip = VertexBuffer.nvertex;\n      break;\n    case 0x0007:\n      VertexBuffer.quads = VertexBuffer.nvertex;\n      break;\n    case 0x0009:\n      VertexBuffer.polygon = VertexBuffer.nvertex;\n      break;\n    case 0x0006:\n      VertexBuffer.fan = VertexBuffer.nvertex;\n      break;\n    default:\n      fprintf (ferr, \"glBegin (%d) not implemented yet\\n\", state);\n      break;\n    }\n  }\n  else\n    glBegin (state);\n}"),_(" \n// #line 921 \"/home/spencer/basilisk/src/draw.h\"\n{\n      if ((filled > 0 && val(d,0,0,0) >= 1.) || (filled < 0 && val(d,0,0,0) <= 0.)) {\n vertex_buffer_glBegin (0x0007);\n glvertex2d (view, x - Delta_x/2., y - Delta_y/2.);\n glvertex2d (view, x + Delta_x/2., y - Delta_y/2.);\n glvertex2d (view, x + Delta_x/2., y + Delta_y/2.);\n glvertex2d (view, x - Delta_x/2., y + Delta_y/2.);\n vertex_buffer_glEnd();\n view->ni++;\n      }\n      else if (val(d,0,0,0) > 0. && val(d,0,0,0) < 1.) {\n coord n = facet_normal (point, d, fs), r = {1.,1.};\n if (filled < 0)\n   {\n     n.x = - n.x;\n     \n// #line 935\nn.y = - n.y;}\n real alpha = line_alpha (filled < 0. ? 1. - val(d,0,0,0) : val(d,0,0,0), n);\n alpha += (n.x + n.y)/2.;\n \n   if (n.x < 0.) alpha -= n.x, n.x = - n.x, r.x = - 1.;\n   \n// #line 939\nif (n.y < 0.) alpha -= n.y, n.y = - n.y, r.y = - 1.;\n coord v[5];\n int nv = 0;\n if (alpha >= 0. && alpha <= n.x) {\n   v[nv].x = alpha/n.x, v[nv++].y = 0.;\n   if (alpha <= n.y)\n     v[nv].x = 0., v[nv++].y = alpha/n.y;\n   else if (alpha >= n.y && alpha - n.y <= n.x) {\n     v[nv].x = (alpha - n.y)/n.x, v[nv++].y = 1.;\n     v[nv].x = 0., v[nv++].y = 1.;\n   }\n   v[nv].x = 0., v[nv++].y = 0.;\n }\n else if (alpha >= n.x && alpha - n.x <= n.y) {\n   v[nv].x = 1., v[nv++].y = (alpha - n.x)/n.y;\n   if (alpha >= n.y && alpha - n.y <= n.x) {\n     v[nv].x = (alpha - n.y)/n.x, v[nv++].y = 1.;\n     v[nv].x = 0., v[nv++].y = 1.;\n   }\n   else if (alpha <= n.y)\n     v[nv].x = 0., v[nv++].y = alpha/n.y;\n   v[nv].x = 0., v[nv++].y = 0.;\n   v[nv].x = 1., v[nv++].y = 0.;\n }\n vertex_buffer_glBegin (0x0009);\n if (r.x*r.y < 0.)\n   for (int i = nv - 1; i >= 0; i--)\n     glvertex2d (view, x + r.x*(v[i].x - 0.5)*Delta,\n   y + r.y*(v[i].y - 0.5)*Delta);\n else\n   for (int i = 0; i < nv; i++)\n     glvertex2d (view, x + r.x*(v[i].x - 0.5)*Delta,\n   y + r.y*(v[i].y - 0.5)*Delta);\n vertex_buffer_glEnd ();\n view->ni++;\n      }\n    }")}) { 
_stencil_val(d,0,0,0); _stencil_val(d,0,0,0);{                             
 
 
 
 
 
 
 
        
{_stencil_val(d,0,0,0); _stencil_val(d,0,0,0); {     
  _stencil_facet_normal (point, d, fs);              
 
   
        
  _stencil_val(d,0,0,0); _stencil_val(d,0,0,0);    
     
            
      
 
   
         
       
       
              
            
         
    
       
  
        
              
            
         
    
             
         
        
         
       
       
             
       
         
     
     
 
         
       
         
     
 
 
      }      }}
                   
      
    
#line 975
}end_foreach_visible_stencil();
    {
#line 921
foreach_visible (view) {
      if ((filled > 0 && val(d,0,0,0) >= 1.) || (filled < 0 && val(d,0,0,0) <= 0.)) {
 vertex_buffer_glBegin (0x0007);
 glvertex2d (view, x - Delta_x/2., y - Delta_y/2.);
 glvertex2d (view, x + Delta_x/2., y - Delta_y/2.);
 glvertex2d (view, x + Delta_x/2., y + Delta_y/2.);
 glvertex2d (view, x - Delta_x/2., y + Delta_y/2.);
 vertex_buffer_glEnd();
 view->ni++;
      }
      else if (val(d,0,0,0) > 0. && val(d,0,0,0) < 1.) {
 coord n = facet_normal (point, d, fs), r = {1.,1.};
 if (filled < 0)
   {
     n.x = - n.x;
     
#line 935
n.y = - n.y;}
 double alpha = line_alpha (filled < 0. ? 1. - val(d,0,0,0) : val(d,0,0,0), n);
 alpha += (n.x + n.y)/2.;
 
   if (n.x < 0.) alpha -= n.x, n.x = - n.x, r.x = - 1.;
   
#line 939
if (n.y < 0.) alpha -= n.y, n.y = - n.y, r.y = - 1.;
 coord v[5];
 int nv = 0;
 if (alpha >= 0. && alpha <= n.x) {
   v[nv].x = alpha/n.x, v[nv++].y = 0.;
   if (alpha <= n.y)
     v[nv].x = 0., v[nv++].y = alpha/n.y;
   else if (alpha >= n.y && alpha - n.y <= n.x) {
     v[nv].x = (alpha - n.y)/n.x, v[nv++].y = 1.;
     v[nv].x = 0., v[nv++].y = 1.;
   }
   v[nv].x = 0., v[nv++].y = 0.;
 }
 else if (alpha >= n.x && alpha - n.x <= n.y) {
   v[nv].x = 1., v[nv++].y = (alpha - n.x)/n.y;
   if (alpha >= n.y && alpha - n.y <= n.x) {
     v[nv].x = (alpha - n.y)/n.x, v[nv++].y = 1.;
     v[nv].x = 0., v[nv++].y = 1.;
   }
   else if (alpha <= n.y)
     v[nv].x = 0., v[nv++].y = alpha/n.y;
   v[nv].x = 0., v[nv++].y = 0.;
   v[nv].x = 1., v[nv++].y = 0.;
 }
 vertex_buffer_glBegin (0x0009);
 if (r.x*r.y < 0.)
   for (int i = nv - 1; i >= 0; i--)
     glvertex2d (view, x + r.x*(v[i].x - 0.5)*Delta,
   y + r.y*(v[i].y - 0.5)*Delta);
 else
   for (int i = 0; i < nv; i++)
     glvertex2d (view, x + r.x*(v[i].x - 0.5)*Delta,
   y + r.y*(v[i].y - 0.5)*Delta);
 vertex_buffer_glEnd ();
 view->ni++;
      }
    }end_foreach_visible();}
  }
  else
    {begin_draw_lines (view, lc, lw); {
      vertex_buffer_glBegin (0x0001);
      foreach_visible_stencil (1,{(NonLocal[]){{"pmtrace","not implemented yet",(void *)&pmtrace,NULL,0},{"pmfuncn","int",(void *)&pmfuncn,NULL,0},{"pmfuncs","pmfunc",(void *)pmfuncs,NULL,1},{"ferr","not implemented yet",(void *)ferr,NULL,1},{"max","double",(void *)&max,NULL,0},{"VertexBuffer","not implemented yet",(void *)&VertexBuffer,NULL,0},{"fs","vector",(void *)&fs,NULL,0},{"cmin","double",(void *)&cmin,NULL,0},{"d","scalar",(void *)&d,NULL,0},{"view","not implemented yet",(void *)view,NULL,1},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_("\n\n\n// #line 135 \"/home/spencer/basilisk/src/common.h\"\nstatic int pmfunc_index (const char * func, const char * file, int line)\n{\n  pmfunc * p = pmfuncs;\n  for (int i = 0; i < pmfuncn; i++, p++)\n    if (p->line == line && !strcmp(func, p->func) && !strcmp(file, p->file))\n      return p->id;\n  pmfuncn++;\n  pmfuncs = (pmfunc *) sysrealloc (pmfuncs, pmfuncn*sizeof(pmfunc));\n  p = &pmfuncs[pmfuncn - 1];\n  memset (p, 0, sizeof(pmfunc));\n  p->func = systrdup(func);\n  p->file = systrdup(file);\n  p->line = line;\n  p->id = pmfuncn;\n  if (pmtrace.fp)\n    fprintf (pmtrace.fp, \"@ %d %s %s %d\\n\", pmfuncn, func, file, line);\n  return pmfuncn;\n}\n@ def not_mpi_compatible()\ndo {\n  if (npe() > 1) {\n    fprintf (ferr, \"%s() is not compatible with MPI (yet)\\n\", __func__);\n    exit (1);\n  }\n} while(0)\n@\n@ define system(command) (pid() == 0 ? system(command) : 0)\n@else\n@ define qstderr() stderr\n@ define qstdout() stdout\n@ define ferr stderr\n@ define fout stdout\n@ define not_mpi_compatible()\n@endif\n\n\n\n\n// #line 93 \"/home/spencer/basilisk/src/common.h\"\nstatic inline void qassert (const char * file, int line, const char * cond) {\n  fprintf (ferr, \"%s:%d: Assertion `%s' failed.\\n\", file, line, cond);\n  abort();\n}\n\n\n// #line 176 \"/home/spencer/basilisk/src/common.h\"\nstatic void * pmfunc_alloc (pmdata * d, size_t size,\n       const char * func, const char * file, int line,\n       char c)\n{\n  if (!(d != NULL)) qassert (\"/home/spencer/basilisk/src/common.h\", 0, \"d != NULL\");\n  OMP (omp critical)\n  {\n    d->id = pmfunc_index(func, file, line);\n    d->size = size;\n    pmfunc * f = &pmfuncs[d->id - 1];\n    f->total += size;\n    if (f->total > f->max)\n      f->max = f->total;\n    pmtrace.total += size;\n    pmtrace.overhead += sizeof(pmdata);\n    if (pmtrace.total > pmtrace.max) {\n      pmtrace.max = pmtrace.total;\n      pmtrace.maxoverhead = pmtrace.overhead;\n    }\n    pmfunc_trace (f, c);\n  }\n  return ((char *)d) + sizeof(pmdata);\n}\n\n\n// #line 154 \"/home/spencer/basilisk/src/common.h\"\nstatic void pmfunc_trace (pmfunc * f, char c)\n{\n  if (pmtrace.fp)\n    fprintf (pmtrace.fp, \"%c %d %ld %ld %ld\",\n      c, f->id, pmtrace.nr, pmtrace.total, f->total);\n@if 1\n  if (pmtrace.nr % 1 == 0) {\n    struct rusage usage;\n    getrusage (RUSAGE_SELF, &usage);\n    if (pmtrace.fp)\n      fprintf (pmtrace.fp, \" %ld\", usage.ru_maxrss*1024);\n    if (!pmtrace.nr)\n      pmtrace.startrss = usage.ru_maxrss;\n    if (usage.ru_maxrss > pmtrace.maxrss)\n      pmtrace.maxrss = usage.ru_maxrss;\n  }\n@endif\n  if (pmtrace.fp)\n    fputc ('\\n', pmtrace.fp);\n  pmtrace.nr++;\n}\n\n\n// #line 200 \"/home/spencer/basilisk/src/common.h\"\nstatic void * pmfunc_free (void * ptr, char c)\n{\n  if (!ptr)\n    return ptr;\n  pmdata * d = (pmdata *) (((char *)ptr) - sizeof(pmdata));\n  if (d->id < 1 || d->id > pmfuncn) {\n    fputs (\"*** MTRACE: ERROR!: corrupted free()\", ferr);\n    if (d->size == 0)\n      fputs (\", possible double free()\", ferr);\n    else\n      fputs (\", not traced?\", ferr);\n    fputs (\", aborting...\\n\", ferr);\n    abort();\n    return ptr;\n  }\n  else\n  OMP (omp critical)\n  {\n    pmfunc * f = &pmfuncs[d->id - 1];\n    if (f->total < d->size) {\n      fprintf (ferr, \"*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\\n\",\n        f->total, d->size);\n      abort();\n    }\n    else\n      f->total -= d->size;\n    if (pmtrace.total < d->size) {\n      fprintf (ferr, \"*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\\n\",\n        pmtrace.total, d->size);\n      abort();\n    }\n    else {\n      pmtrace.total -= d->size;\n      pmtrace.overhead -= sizeof(pmdata);\n    }\n    d->id = 0;\n    d->size = 0;\n    pmfunc_trace (f, c);\n  }\n  return d;\n}\n\n\n// #line 256 \"/home/spencer/basilisk/src/common.h\"\nstatic void * prealloc (void * ptr, size_t size,\n   const char * func, const char * file, int line)\n{\n  return pmfunc_alloc ((pmdata *) sysrealloc (pmfunc_free(ptr, '<'),\n           sizeof(pmdata) + size),\n         size, func, file, line, '>');\n}\n\n\n// #line 394 \"/home/spencer/basilisk/src/common.h\"\nvoid array_append (Array * a, void * elem, size_t size)\n{\n  if (a->len + size >= a->max) {\n    a->max += max (size, 4096);\n    a->p = prealloc (a->p, a->max,__func__,__FILE__,0);\n  }\n  memcpy (((char *)a->p) + a->len, elem, size);\n  a->len += size;\n}\n\n\n// #line 222 \"/home/spencer/basilisk/src/vertexbuffer.h\"\nstatic void vertex_buffer_glVertex2d (real x, real y)\n{\n  if (VertexBuffer.position) {\n    if (VertexBuffer.dim < 2)\n      VertexBuffer.dim = 2;\n    real v[4] = {x, y, 0, 1.};\n    vector_multiply (v, VertexBuffer.modelview);\n    array_append (VertexBuffer.position, v, 3*sizeof(real));\n    VertexBuffer.nvertex++;\n  }\n  else\n    glVertex3d (x, y, 0.);\n}\n\n\n\n// #line 838 \"/home/spencer/basilisk/src/draw.h\"\nstatic void glvertex2d (bview * view, real x, real y) {\n  if (view->map) {\n    coord p = {x, y, 0.};\n    view->map (&p);\n    vertex_buffer_glVertex2d (p.x, p.y);\n  }\n  else\n    vertex_buffer_glVertex2d (x, y);\n}\n// # 292 \"/home/spencer/basilisk/src/geometry.h\"\nint facets (coord n, real alpha, coord p[2])\n{\n  int i = 0;\n  for (real s = -0.5; s <= 0.5; s += 1.)\n    {\n      if (fabs (n.y) > 1e-4 && i < 2) {\n real a = (alpha - s*n.x)/n.y;\n if (a >= -0.5 && a <= 0.5) {\n   p[i].x = s;\n   p[i++].y = a;\n }\n      }\n      \n// #line 297\nif (fabs (n.x) > 1e-4 && i < 2) {\n real a = (alpha - s*n.y)/n.x;\n if (a >= -0.5 && a <= 0.5) {\n   p[i].y = s;\n   p[i++].x = a;\n }\n      }}\n  return i;\n}\n// # 415 \"/home/spencer/basilisk/src/view.h\" 2\n\n\n\n\n\n\n// # 1 \"./draw.h\" 1\n// # 1 \"/home/spencer/basilisk/src/draw.h\"\n\n\n\n\n// # 1 \"./fractions.h\" 1\n// # 1 \"/home/spencer/basilisk/src/fractions.h\"\n// # 12 \"/home/spencer/basilisk/src/fractions.h\"\n// # 1 \"./geometry.h\" 1\n// # 1 \"/home/spencer/basilisk/src/geometry.h\"\n// # 35 \"/home/spencer/basilisk/src/geometry.h\"\nreal line_alpha (real c, coord n)\n{\n  real alpha, n1, n2;\n\n  n1 = fabs (n.x); n2 = fabs (n.y);\n  if (n1 > n2)\n    do { real __tmp = n1; n1 = n2; n2 = __tmp; } while(0);\n\n  c = clamp (c, 0., 1.);\n  real v1 = n1/2.;\n  if (c <= v1/n2)\n    alpha = sqrt (2.*c*n1*n2);\n  else if (c <= 1. - v1/n2)\n    alpha = c*n2 + v1;\n  else\n    alpha = n1 + n2 - sqrt (2.*n1*n2*(1. - c));\n\n  if (n.x < 0.)\n    alpha += n.x;\n  if (n.y < 0.)\n    alpha += n.y;\n\n  return alpha - (n.x + n.y)/2.;\n}\n// # 13 \"/home/spencer/basilisk/src/fractions.h\" 2\n\n\n\n\n\n// # 1 \"./myc2d.h\" 1\n// # 1 \"/home/spencer/basilisk/src/myc2d.h\"\n\n\n\n\n\ncoord mycs (Point point, scalar c)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  int ix;\n  real c_t,c_b,c_r,c_l;\n  real mx0,my0,mx1,my1,mm1,mm2;\n\n\n  c_t = val(c,-1,1,0) + val(c,0,1,0) + val(c,1,1,0);\n  c_b = val(c,-1,-1,0) + val(c,0,-1,0) + val(c,1,-1,0);\n  c_r = val(c,1,-1,0) + val(c,1,0,0) + val(c,1,1,0);\n  c_l = val(c,-1,-1,0) + val(c,-1,0,0) + val(c,-1,1,0);\n\n\n\n  mx0 = 0.5*(c_l - c_r);\n  my0 = 0.5*(c_b - c_t);\n\n\n  if (fabs(mx0) <= fabs(my0)) {\n    my0 = my0 > 0. ? 1. : -1.;\n    ix = 1;\n  }\n  else {\n    mx0 = mx0 > 0. ? 1. : -1.;\n    ix = 0;\n  }\n\n\n  mm1 = val(c,-1,-1,0) + 2.0*val(c,-1,0,0) + val(c,-1,1,0);\n  mm2 = val(c,1,-1,0) + 2.0*val(c,1,0,0) + val(c,1,1,0);\n  mx1 = mm1 - mm2 + 1e-30;\n  mm1 = val(c,-1,-1,0) + 2.0*val(c,0,-1,0) + val(c,1,-1,0);\n  mm2 = val(c,-1,1,0) + 2.0*val(c,0,1,0) + val(c,1,1,0);\n  my1 = mm1 - mm2 + 1e-30;\n\n\n  if (ix) {\n    mm1 = fabs(my1);\n    mm1 = fabs(mx1)/mm1;\n    if (mm1 > fabs(mx0)) {\n      mx0 = mx1;\n      my0 = my1;\n    }\n  }\n  else {\n    mm1 = fabs(mx1);\n    mm1 = fabs(my1)/mm1;\n    if (mm1 > fabs(my0)) {\n      mx0 = mx1;\n      my0 = my1;\n    }\n  }\n\n\n\n  mm1 = fabs(mx0) + fabs(my0);\n  coord n = {mx0/mm1, my0/mm1};\n\n  return n;\n}\n\n\n\n\n\n\n// #line 418 \"/home/spencer/basilisk/src/fractions.h\"\ncoord facet_normal (Point point, scalar c, vector s)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  if (s.x.i >= 0) {\n    coord n;\n    real nn = 0.;\n     {\n      n.x = val(s.x,0,0,0) - val(s.x,1,0,0);\n      nn += fabs(n.x);\n    } \n// #line 423\n{\n      n.y = val(s.y,0,0,0) - val(s.y,0,1,0);\n      nn += fabs(n.y);\n    }\n    if (nn > 0.)\n      {\n n.x /= nn;\n \n// #line 429\nn.y /= nn;}\n    else\n      {\n n.x = 1./2;\n \n// #line 432\nn.y = 1./2;}\n    return n;\n  }\n  return mycs (point, c);\n}\n// # 801 \"/home/spencer/basilisk/src/draw.h\"\nstatic bool cfilter (Point point, scalar c, real cmin)\n{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  real cmin1 = 4.*cmin;\n  if (val(c,0,0,0) <= cmin) {\n    \n      if (val(c,1,0,0) >= 1. - cmin1 || val(c,-1,0,0) >= 1. - cmin1)\n return true;\n      \n// #line 806\nif (val(c,0,1,0) >= 1. - cmin1 || val(c,0,-1,0) >= 1. - cmin1)\n return true;\n    return false;\n  }\n  if (val(c,0,0,0) >= 1. - cmin) {\n    \n      if (val(c,1,0,0) <= cmin1 || val(c,-1,0,0) <= cmin1)\n return true;\n      \n// #line 812\nif (val(c,0,1,0) <= cmin1 || val(c,0,-1,0) <= cmin1)\n return true;\n    return false;\n  }\n  int n = 0;\n  real min = 1e30, max = - 1e30;\n  {foreach_neighbor(1) {\n    if (val(c,0,0,0) > cmin && val(c,0,0,0) < 1. - cmin && ++n >= (1 << 2))\n      return true;\n    if (val(c,0,0,0) > max) max = val(c,0,0,0);\n    if (val(c,0,0,0) < min) min = val(c,0,0,0);\n  }end_foreach_neighbor()}\n  return max - min > 0.5;\n}"),_("\n \n// #line 981 \"/home/spencer/basilisk/src/draw.h\"\nif (cfilter (point, d, cmin)) {\n   coord n = facet_normal (point, d, fs);\n   real alpha = line_alpha (val(d,0,0,0), n);\n   coord segment[2];\n   if (facets (n, alpha, segment) == 2) {\n     glvertex2d (view, x + segment[0].x*Delta, y + segment[0].y*Delta);\n     glvertex2d (view, x + segment[1].x*Delta, y + segment[1].y*Delta);\n     view->ni++;\n   }\n }")})
 {_stencil_cfilter (point, d, NULL); {  
    _stencil_facet_normal (point, d, fs);     
   _stencil_val(d,0,0,0); 
                 
     
     
     
    
         
 } }end_foreach_visible_stencil();
      {
#line 980
foreach_visible (view)
 if (cfilter (point, d, cmin)) {
   coord n = facet_normal (point, d, fs);
   double alpha = line_alpha (val(d,0,0,0), n);
   coord segment[2];
   if (facets (n, alpha, segment) == 2) {
     glvertex2d (view, x + segment[0].x*Delta, y + segment[0].y*Delta);
     glvertex2d (view, x + segment[1].x*Delta, y + segment[1].y*Delta);
     view->ni++;
   }
 }end_foreach_visible();}
      vertex_buffer_glEnd ();
    }end_draw_lines();}
# 1044 "/home/spencer/basilisk/src/draw.h"
  if (prolongation != fraction_refine) {
    _attribute[d.i].prolongation = prolongation;
    _attribute[d.i].dirty = true;
  }


  if (expr) delete(((scalar[]){col,{-1}}));
  {end_tracing("draw_vof","/home/spencer/basilisk/src/draw.h",0);return true;}
end_tracing("draw_vof","/home/spencer/basilisk/src/draw.h",0);}
# 1063 "/home/spencer/basilisk/src/draw.h"
     
bool isoline (char * phi,
       double val,
       int n,
       bool edges,
       double larger, int filled,
       char * color,
       double min, double max, double spread,
       bool linear,
       Colormap map,
       float fc[3], float lc[3], float lw,
       bool expr,
       bool cbar, float size, float pos[2], char * label, double lscale, bool horizontal, bool border, bool mid, float fsize, char * format, int levels)
{tracing("isoline","/home/spencer/basilisk/src/draw.h",0);

  if (!color) color = phi;
  scalar col = {-1}; if (color && strcmp (color, "level")) { col = compile_expression (color, &expr); if (col.i < 0) {end_tracing("isoline","/home/spencer/basilisk/src/draw.h",0);return false;} } double cmap[127][3]; if (color) { if (min == 0 && max == 0) { if (col.i < 0) min = 0, max = depth(); else { stats s = statsf (col); double avg = s.sum/s.volume; if (spread < 0.) min = s.min, max = s.max; else { if (!spread) spread = 5.; min = avg - spread*s.stddev; max = avg + spread*s.stddev; } } } if (!map) map = jet; (* map) (cmap); } if ((2 > 2 || linear) && !fc[0] && !fc[1] && !fc[2]) fc[0] = fc[1] = fc[2] = 1.; if (cbar) colorbar (map, size, pos, label, lscale, min, max, horizontal, border, mid, lc, lw, fsize, format, levels);;
  scalar fphi = col,  fiso=new_scalar("fiso");
  if (!is_vertex_scalar (col)) {
    fphi = new_scalar("fphi");
    foreach_vertex_stencil(1,{(NonLocal[]){{"col","scalar",(void *)&col,NULL,0},{"fphi","scalar",(void *)&fphi,NULL,0},{"val","double",(void *)&val,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n      \n// #line 1084 \"/home/spencer/basilisk/src/draw.h\"\nval_out_(fphi,0,0,0) = (val(col,0,0,0) + val(col,-1,0,0) + val(col,0,-1,0) + val(col,-1,-1,0))/4.;")})
      {_stencil_val(col,0,0,0); _stencil_val(col,-1,0,0); _stencil_val(col,0,-1,0); _stencil_val(col,-1,-1,0);_stencil_val_a(fphi,0,0,0);     }end_foreach_vertex_stencil();
    {
#line 1083
foreach_vertex()
      val(fphi,0,0,0) = (val(col,0,0,0) + val(col,-1,0,0) + val(col,0,-1,0) + val(col,-1,-1,0))/4.;end_foreach_vertex();}
  }
  vector  siso=new_face_vector("siso");
  if (n < 2) {
    fractions (fphi, fiso, siso, val);
    draw_vof ("fiso", "siso", edges, larger, filled, color, min, max, spread
,
       
#line 1090
linear, map, fc, lc, lw, expr
#line 893
, 
false, 15,( float[2]) {-.95, -.95}, "", 1, false, false, false, 50, "%g", 50
#line 1090
);
  }
  else if (max > min) {
    double dv = (max - min)/(n - 1);
    for (val = min; val <= max; val += dv) {
      fractions (fphi, fiso, siso, val);
      draw_vof ("fiso", "siso", edges, larger, filled, color, min, max, spread
,
  
#line 1097
linear, map, fc, lc, lw, expr
#line 893
, 
false, 15,( float[2]) {-.95, -.95}, "", 1, false, false, false, 50, "%g", 50
#line 1097
);
    }
  }
  if (!is_vertex_scalar (col))
    delete (((scalar[]){fphi,{-1}}));
  if (expr) delete (((scalar[]){col,{-1}}));



  {delete((scalar*)((scalar[]){siso.x,siso.y,fiso,{-1}}));{end_tracing("isoline","/home/spencer/basilisk/src/draw.h",0);return true;}}delete((scalar*)((scalar[]){siso.x,siso.y,fiso,{-1}}));
end_tracing("isoline","/home/spencer/basilisk/src/draw.h",0);}
# 1120 "/home/spencer/basilisk/src/draw.h"
     
bool cells (coord n, double alpha,
     float lc[3], float lw)
{tracing("cells","/home/spencer/basilisk/src/draw.h",0);
  bview * view = draw();
  {begin_draw_lines (view, lc, lw); {

    {foreach_visible (view) {
      vertex_buffer_glBegin (0x0002);
      glvertex2d (view, x - Delta_x/2., y - Delta_y/2.);
      glvertex2d (view, x + Delta_x/2., y - Delta_y/2.);
      glvertex2d (view, x + Delta_x/2., y + Delta_y/2.);
      glvertex2d (view, x - Delta_x/2., y + Delta_y/2.);
      vertex_buffer_glEnd();
      view->ni++;
    }end_foreach_visible();}
# 1149 "/home/spencer/basilisk/src/draw.h"
  }end_draw_lines();}
  {end_tracing("cells","/home/spencer/basilisk/src/draw.h",0);return true;}
end_tracing("cells","/home/spencer/basilisk/src/draw.h",0);}






     
bool vectors (char * u, double scale, float lc[3], float lw)
{tracing("vectors","/home/spencer/basilisk/src/draw.h",0);

  vector fu;
  struct { char x, y, z; } index = {'x', 'y', 'z'};
   {
    char name[80];
    sprintf (name, "%s.%c", u, index.x);
    fu.x = lookup_field (name);
    if (fu.x.i < 0) {
      fprintf (ferr, "vectors(): no field named '%s'\n", name);
      {end_tracing("vectors","/home/spencer/basilisk/src/draw.h",0);return false;}
    }
  } 
#line 1164
{
    char name[80];
    sprintf (name, "%s.%c", u, index.y);
    fu.y = lookup_field (name);
    if (fu.y.i < 0) {
      fprintf (ferr, "vectors(): no field named '%s'\n", name);
      {end_tracing("vectors","/home/spencer/basilisk/src/draw.h",0);return false;}
    }
  }
  bview * view = draw();
  float res = view->res;
  if (view->res < 15*view->samples)
    view->res = 15*view->samples;
  {begin_draw_lines (view, lc, lw); {
    double fscale = (scale ? scale : 1.)*view->res/view->samples;
    vertex_buffer_glBegin (0x0001);
    foreach_visible_stencil (1,{(NonLocal[]){{"pmtrace","not implemented yet",(void *)&pmtrace,NULL,0},{"pmfuncn","int",(void *)&pmfuncn,NULL,0},{"pmfuncs","pmfunc",(void *)pmfuncs,NULL,1},{"ferr","not implemented yet",(void *)ferr,NULL,1},{"VertexBuffer","not implemented yet",(void *)&VertexBuffer,NULL,0},{"fscale","double",(void *)&fscale,NULL,0},{"fu","vector",(void *)&fu,NULL,0},{"view","not implemented yet",(void *)view,NULL,1},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_("\n\n\n// #line 135 \"/home/spencer/basilisk/src/common.h\"\nstatic int pmfunc_index (const char * func, const char * file, int line)\n{\n  pmfunc * p = pmfuncs;\n  for (int i = 0; i < pmfuncn; i++, p++)\n    if (p->line == line && !strcmp(func, p->func) && !strcmp(file, p->file))\n      return p->id;\n  pmfuncn++;\n  pmfuncs = (pmfunc *) sysrealloc (pmfuncs, pmfuncn*sizeof(pmfunc));\n  p = &pmfuncs[pmfuncn - 1];\n  memset (p, 0, sizeof(pmfunc));\n  p->func = systrdup(func);\n  p->file = systrdup(file);\n  p->line = line;\n  p->id = pmfuncn;\n  if (pmtrace.fp)\n    fprintf (pmtrace.fp, \"@ %d %s %s %d\\n\", pmfuncn, func, file, line);\n  return pmfuncn;\n}\n@ def not_mpi_compatible()\ndo {\n  if (npe() > 1) {\n    fprintf (ferr, \"%s() is not compatible with MPI (yet)\\n\", __func__);\n    exit (1);\n  }\n} while(0)\n@\n@ define system(command) (pid() == 0 ? system(command) : 0)\n@else\n@ define qstderr() stderr\n@ define qstdout() stdout\n@ define ferr stderr\n@ define fout stdout\n@ define not_mpi_compatible()\n@endif\n\n\n\n\n// #line 93 \"/home/spencer/basilisk/src/common.h\"\nstatic inline void qassert (const char * file, int line, const char * cond) {\n  fprintf (ferr, \"%s:%d: Assertion `%s' failed.\\n\", file, line, cond);\n  abort();\n}\n\n\n// #line 176 \"/home/spencer/basilisk/src/common.h\"\nstatic void * pmfunc_alloc (pmdata * d, size_t size,\n       const char * func, const char * file, int line,\n       char c)\n{\n  if (!(d != NULL)) qassert (\"/home/spencer/basilisk/src/common.h\", 0, \"d != NULL\");\n  OMP (omp critical)\n  {\n    d->id = pmfunc_index(func, file, line);\n    d->size = size;\n    pmfunc * f = &pmfuncs[d->id - 1];\n    f->total += size;\n    if (f->total > f->max)\n      f->max = f->total;\n    pmtrace.total += size;\n    pmtrace.overhead += sizeof(pmdata);\n    if (pmtrace.total > pmtrace.max) {\n      pmtrace.max = pmtrace.total;\n      pmtrace.maxoverhead = pmtrace.overhead;\n    }\n    pmfunc_trace (f, c);\n  }\n  return ((char *)d) + sizeof(pmdata);\n}\n\n\n// #line 154 \"/home/spencer/basilisk/src/common.h\"\nstatic void pmfunc_trace (pmfunc * f, char c)\n{\n  if (pmtrace.fp)\n    fprintf (pmtrace.fp, \"%c %d %ld %ld %ld\",\n      c, f->id, pmtrace.nr, pmtrace.total, f->total);\n@if 1\n  if (pmtrace.nr % 1 == 0) {\n    struct rusage usage;\n    getrusage (RUSAGE_SELF, &usage);\n    if (pmtrace.fp)\n      fprintf (pmtrace.fp, \" %ld\", usage.ru_maxrss*1024);\n    if (!pmtrace.nr)\n      pmtrace.startrss = usage.ru_maxrss;\n    if (usage.ru_maxrss > pmtrace.maxrss)\n      pmtrace.maxrss = usage.ru_maxrss;\n  }\n@endif\n  if (pmtrace.fp)\n    fputc ('\\n', pmtrace.fp);\n  pmtrace.nr++;\n}\n\n\n// #line 200 \"/home/spencer/basilisk/src/common.h\"\nstatic void * pmfunc_free (void * ptr, char c)\n{\n  if (!ptr)\n    return ptr;\n  pmdata * d = (pmdata *) (((char *)ptr) - sizeof(pmdata));\n  if (d->id < 1 || d->id > pmfuncn) {\n    fputs (\"*** MTRACE: ERROR!: corrupted free()\", ferr);\n    if (d->size == 0)\n      fputs (\", possible double free()\", ferr);\n    else\n      fputs (\", not traced?\", ferr);\n    fputs (\", aborting...\\n\", ferr);\n    abort();\n    return ptr;\n  }\n  else\n  OMP (omp critical)\n  {\n    pmfunc * f = &pmfuncs[d->id - 1];\n    if (f->total < d->size) {\n      fprintf (ferr, \"*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\\n\",\n        f->total, d->size);\n      abort();\n    }\n    else\n      f->total -= d->size;\n    if (pmtrace.total < d->size) {\n      fprintf (ferr, \"*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\\n\",\n        pmtrace.total, d->size);\n      abort();\n    }\n    else {\n      pmtrace.total -= d->size;\n      pmtrace.overhead -= sizeof(pmdata);\n    }\n    d->id = 0;\n    d->size = 0;\n    pmfunc_trace (f, c);\n  }\n  return d;\n}\n\n\n// #line 256 \"/home/spencer/basilisk/src/common.h\"\nstatic void * prealloc (void * ptr, size_t size,\n   const char * func, const char * file, int line)\n{\n  return pmfunc_alloc ((pmdata *) sysrealloc (pmfunc_free(ptr, '<'),\n           sizeof(pmdata) + size),\n         size, func, file, line, '>');\n}\n\n\n// #line 394 \"/home/spencer/basilisk/src/common.h\"\nvoid array_append (Array * a, void * elem, size_t size)\n{\n  if (a->len + size >= a->max) {\n    a->max += max (size, 4096);\n    a->p = prealloc (a->p, a->max,__func__,__FILE__,0);\n  }\n  memcpy (((char *)a->p) + a->len, elem, size);\n  a->len += size;\n}\n\n\n// #line 222 \"/home/spencer/basilisk/src/vertexbuffer.h\"\nstatic void vertex_buffer_glVertex2d (real x, real y)\n{\n  if (VertexBuffer.position) {\n    if (VertexBuffer.dim < 2)\n      VertexBuffer.dim = 2;\n    real v[4] = {x, y, 0, 1.};\n    vector_multiply (v, VertexBuffer.modelview);\n    array_append (VertexBuffer.position, v, 3*sizeof(real));\n    VertexBuffer.nvertex++;\n  }\n  else\n    glVertex3d (x, y, 0.);\n}\n\n\n\n// #line 838 \"/home/spencer/basilisk/src/draw.h\"\nstatic void glvertex2d (bview * view, real x, real y) {\n  if (view->map) {\n    coord p = {x, y, 0.};\n    view->map (&p);\n    vertex_buffer_glVertex2d (p.x, p.y);\n  }\n  else\n    vertex_buffer_glVertex2d (x, y);\n}"),_("\n      \n// #line 1181 \"/home/spencer/basilisk/src/draw.h\"\nif (val(fu.x,0,0,0) != 1e30) {\n coord f = { fscale*val(fu.x,0,0,0), fscale*val(fu.y,0,0,0) };\n glvertex2d (view, x + f.x - (f.x - f.y/2.)/5.,\n      y + f.y - (f.x/2. + f.y)/5.);\n glvertex2d (view, x + f.x, y + f.y);\n glvertex2d (view, x + f.x, y + f.y);\n glvertex2d (view, x + f.x - (f.x + f.y/2.)/5.,\n      y + f.y + (f.x/2. - f.y)/5.);\n glvertex2d (view, x, y);\n glvertex2d (view, x + f.x, y + f.y);\n view->ni++;\n      }")})
      {_stencil_val(fu.x,0,0,0); {      
 _stencil_val(fu.y,0,0,0);_stencil_val(fu.x,0,0,0);                
                                  
              
 
 
 
 
 
 
 
      }   }end_foreach_visible_stencil();
    {
#line 1180
foreach_visible (view)
      if (val(fu.x,0,0,0) != 1e30) {
 coord f = { fscale*val(fu.x,0,0,0), fscale*val(fu.y,0,0,0) };
 glvertex2d (view, x + f.x - (f.x - f.y/2.)/5.,
      y + f.y - (f.x/2. + f.y)/5.);
 glvertex2d (view, x + f.x, y + f.y);
 glvertex2d (view, x + f.x, y + f.y);
 glvertex2d (view, x + f.x - (f.x + f.y/2.)/5.,
      y + f.y + (f.x/2. - f.y)/5.);
 glvertex2d (view, x, y);
 glvertex2d (view, x + f.x, y + f.y);
 view->ni++;
      }end_foreach_visible();}
    vertex_buffer_glEnd();
  }end_draw_lines();}
  view->res = res;



  {end_tracing("vectors","/home/spencer/basilisk/src/draw.h",0);return true;}
end_tracing("vectors","/home/spencer/basilisk/src/draw.h",0);}
# 1220 "/home/spencer/basilisk/src/draw.h"
     
bool squares (char * color,
       char * z,
       double min, double max, double spread,
       bool linear,
       Colormap map,
       float fc[3], float lc[3],
       bool expr,

       coord n,
       double alpha,
       float lw,
       bool cbar, float size, float pos[2], char * label, double lscale, bool horizontal, bool border, bool mid, float fsize, char * format, int levels)
{tracing("squares","/home/spencer/basilisk/src/draw.h",0);

  scalar Z = {-1};
  vector fn;
  bool zexpr = false;
  if (z) {
    Z = compile_expression (z, &zexpr);
    if (Z.i < 0)
      {end_tracing("squares","/home/spencer/basilisk/src/draw.h",0);return false;}
    fn = new_vector("fn");
    foreach_stencil(1,{(NonLocal[]){{"Z","scalar",(void *)&Z,NULL,0},{"fn","vector",(void *)&fn,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n      \n// #line 1244 \"/home/spencer/basilisk/src/draw.h\"\n{\n        val_out_(fn.x,0,0,0) = (val(Z,1,0,0) - val(Z,-1,0,0))/(2.*Delta_x);\n        \n// #line 1245\nval_out_(fn.y,0,0,0) = (val(Z,0,1,0) - val(Z,0,-1,0))/(2.*Delta_y);}")})
      {
        {_stencil_val(Z,1,0,0); _stencil_val(Z,-1,0,0);_stencil_val_a(fn.x,0,0,0);   }
        
#line 1245
{_stencil_val(Z,0,1,0); _stencil_val(Z,0,-1,0);_stencil_val_a(fn.y,0,0,0);   }}end_foreach_stencil();
    {
#line 1243
foreach()
      {
        val(fn.x,0,0,0) = (val(Z,1,0,0) - val(Z,-1,0,0))/(2.*Delta_x);
        
#line 1245
val(fn.y,0,0,0) = (val(Z,0,1,0) - val(Z,0,-1,0))/(2.*Delta_y);}end_foreach();}
  }

  scalar col = {-1}; if (color && strcmp (color, "level")) { col = compile_expression (color, &expr); if (col.i < 0) {end_tracing("squares","/home/spencer/basilisk/src/draw.h",0);return false;} } double cmap[127][3]; if (color) { if (min == 0 && max == 0) { if (col.i < 0) min = 0, max = depth(); else { stats s = statsf (col); double avg = s.sum/s.volume; if (spread < 0.) min = s.min, max = s.max; else { if (!spread) spread = 5.; min = avg - spread*s.stddev; max = avg + spread*s.stddev; } } } if (!map) map = jet; (* map) (cmap); } if ((2 > 2 || linear) && !fc[0] && !fc[1] && !fc[2]) fc[0] = fc[1] = fc[2] = 1.; if (cbar) colorbar (map, size, pos, label, lscale, min, max, horizontal, border, mid, lc, lw, fsize, format, levels);;
  scalar f = col;

  bview * view = draw();
  glShadeModel (0x1D01);
  if (linear) {
    {begin_colorized (fc, !VertexBuffer.color || !color, cmap, !VertexBuffer.color && color && linear && col.i >= 0); {

      if (Z.i < 0) {
 vertex_buffer_glNormal3d (0, 0, view->reversed ? -1 : 1);
 foreach_visible_stencil (1,{(NonLocal[]){{"min","double",(void *)&min,NULL,0},{"cmap","double",(void *)cmap,(int[]){1273,0},0},{"col","scalar",(void *)&col,NULL,0},{"linear","bool",(void *)&linear,NULL,0},{"color","not implemented yet",(void *)color,NULL,1},{"baseblock","scalar",(void *)baseblock,NULL,1},{"all","scalar",(void *)all,NULL,1},{"max","double",(void *)&max,NULL,0},{"free_solver_funcs","Array",(void *)free_solver_funcs,NULL,1},{"pmtrace","not implemented yet",(void *)&pmtrace,NULL,0},{"pmfuncn","int",(void *)&pmfuncn,NULL,0},{"pmfuncs","pmfunc",(void *)pmfuncs,NULL,1},{"ferr","not implemented yet",(void *)ferr,NULL,1},{"_view","not implemented yet",(void *)_view,NULL,1},{"VertexBuffer","not implemented yet",(void *)&VertexBuffer,NULL,0},{"f","scalar",(void *)&f,NULL,0},{"view","not implemented yet",(void *)view,NULL,1},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_("\n\n\n// #line 32 \"/home/spencer/basilisk/src/vertexbuffer.h\"\nstatic void vertex_buffer_push_index (unsigned int i)\n{\n  i -= VertexBuffer.vertex;\n  array_append (VertexBuffer.index, &i, sizeof(unsigned int));\n}\n\n\n// #line 112 \"/home/spencer/basilisk/src/vertexbuffer.h\"\nstatic void vertex_buffer_glEnd()\n{\n  if (VertexBuffer.index) {\n    int type = -1;\n    switch (VertexBuffer.state) {\n\n    case 0x0002:\n      for (int i = VertexBuffer.line_loop; i < VertexBuffer.nvertex - 1; i++) {\n vertex_buffer_push_index (i);\n vertex_buffer_push_index (i + 1);\n      }\n      vertex_buffer_push_index (VertexBuffer.nvertex - 1);\n      vertex_buffer_push_index (VertexBuffer.line_loop);\n      type = 0;\n      break;\n\n    case 0x0001:\n      for (int i = VertexBuffer.lines; i < VertexBuffer.nvertex; i += 2) {\n vertex_buffer_push_index (i);\n vertex_buffer_push_index (i + 1);\n      }\n      type = 0;\n      break;\n\n    case 0x0003:\n      for (int i = VertexBuffer.line_strip; i < VertexBuffer.nvertex - 1; i++) {\n vertex_buffer_push_index (i);\n vertex_buffer_push_index (i + 1);\n      }\n      type = 0;\n      break;\n\n    case 0x0007:\n      for (int i = VertexBuffer.quads; i < VertexBuffer.nvertex; i += 4)\n for (int j = 1; j <= 2; j++) {\n   vertex_buffer_push_index (i);\n   vertex_buffer_push_index (i + j);\n   vertex_buffer_push_index (i + j + 1);\n }\n      type = 1;\n      break;\n\n    case 0x0009:\n      for (int j = 1; j <= VertexBuffer.nvertex - VertexBuffer.polygon - 2;\n    j++) {\n vertex_buffer_push_index (VertexBuffer.polygon);\n vertex_buffer_push_index (VertexBuffer.polygon + j);\n vertex_buffer_push_index (VertexBuffer.polygon + j + 1);\n      }\n      type = 1;\n      break;\n\n    case 0x0006:\n      for (int i = VertexBuffer.fan + 1; i < VertexBuffer.nvertex - 1; i++) {\n vertex_buffer_push_index (VertexBuffer.fan);\n vertex_buffer_push_index (i);\n vertex_buffer_push_index (i + 1);\n      }\n      type = 1;\n      break;\n\n    default:\n      break;\n    }\n    VertexBuffer.state = 0;\n    if (VertexBuffer.type >= 0 && type >= 0) {\n\n      if (!(VertexBuffer.type == type)) qassert (\"/home/spencer/basilisk/src/vertexbuffer.h\", 0, \"VertexBuffer.type == type\");\n    }\n    else\n      VertexBuffer.type = type;\n  }\n  else\n    glEnd();\n}\n\n\n// #line 222 \"/home/spencer/basilisk/src/vertexbuffer.h\"\nstatic void vertex_buffer_glVertex2d (real x, real y)\n{\n  if (VertexBuffer.position) {\n    if (VertexBuffer.dim < 2)\n      VertexBuffer.dim = 2;\n    real v[4] = {x, y, 0, 1.};\n    vector_multiply (v, VertexBuffer.modelview);\n    array_append (VertexBuffer.position, v, 3*sizeof(real));\n    VertexBuffer.nvertex++;\n  }\n  else\n    glVertex3d (x, y, 0.);\n}\n\n\n\n// #line 838 \"/home/spencer/basilisk/src/draw.h\"\nstatic void glvertex2d (bview * view, real x, real y) {\n  if (view->map) {\n    coord p = {x, y, 0.};\n    view->map (&p);\n    vertex_buffer_glVertex2d (p.x, p.y);\n  }\n  else\n    vertex_buffer_glVertex2d (x, y);\n}\n\n\n// #line 188 \"/home/spencer/basilisk/src/vertexbuffer.h\"\nstatic void vertex_buffer_glColor3f (real r, real g, real b)\n{\n  if (VertexBuffer.color) {\n    struct { real x, y, z; } color = {r, g, b};\n    array_append (VertexBuffer.color, &color, 3*sizeof(real));\n  }\n  else\n    glColor3f (r, g, b);\n}\n\n\n// #line 261 \"/home/spencer/basilisk/src/output.h\"\nColor colormap_color (real cmap[127][3],\n        real val, real min, real max)\n{\n  Color c;\n  if (val == 1e30) {\n    c.r = c.g = c.b = 0;\n    return c;\n  }\n  int i;\n  real coef;\n  if (max != min)\n    val = (val - min)/(max - min);\n  else\n    val = 0.;\n  if (val <= 0.) i = 0, coef = 0.;\n  else if (val >= 1.) i = 127 - 2, coef = 1.;\n  else {\n    i = val*(127 - 1);\n    coef = val*(127 - 1) - i;\n  }\n  if (!(i >= 0 && i < 127 - 1)) qassert (\"/home/spencer/basilisk/src/output.h\", 0, \"i >= 0 && i < NCMAP - 1\");\n  unsigned char * c1 = (unsigned char *) &c;\n  for (int j = 0; j < 3; j++)\n    c1[j] = 255*(cmap[i][j]*(1. - coef) + cmap[i + 1][j]*coef);\n  return c;\n}\n\n\n// #line 200 \"/home/spencer/basilisk/src/common.h\"\nstatic void * pmfunc_free (void * ptr, char c)\n{\n  if (!ptr)\n    return ptr;\n  pmdata * d = (pmdata *) (((char *)ptr) - sizeof(pmdata));\n  if (d->id < 1 || d->id > pmfuncn) {\n    fputs (\"*** MTRACE: ERROR!: corrupted free()\", ferr);\n    if (d->size == 0)\n      fputs (\", possible double free()\", ferr);\n    else\n      fputs (\", not traced?\", ferr);\n    fputs (\", aborting...\\n\", ferr);\n    abort();\n    return ptr;\n  }\n  else\n  OMP (omp critical)\n  {\n    pmfunc * f = &pmfuncs[d->id - 1];\n    if (f->total < d->size) {\n      fprintf (ferr, \"*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\\n\",\n        f->total, d->size);\n      abort();\n    }\n    else\n      f->total -= d->size;\n    if (pmtrace.total < d->size) {\n      fprintf (ferr, \"*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\\n\",\n        pmtrace.total, d->size);\n      abort();\n    }\n    else {\n      pmtrace.total -= d->size;\n      pmtrace.overhead -= sizeof(pmdata);\n    }\n    d->id = 0;\n    d->size = 0;\n    pmfunc_trace (f, c);\n  }\n  return d;\n}\n\n\n// #line 256 \"/home/spencer/basilisk/src/common.h\"\nstatic void * prealloc (void * ptr, size_t size,\n   const char * func, const char * file, int line)\n{\n  return pmfunc_alloc ((pmdata *) sysrealloc (pmfunc_free(ptr, '<'),\n           sizeof(pmdata) + size),\n         size, func, file, line, '>');\n}\n\n\n// #line 394 \"/home/spencer/basilisk/src/common.h\"\nvoid array_append (Array * a, void * elem, size_t size)\n{\n  if (a->len + size >= a->max) {\n    a->max += max (size, 4096);\n    a->p = prealloc (a->p, a->max,__func__,__FILE__,0);\n  }\n  memcpy (((char *)a->p) + a->len, elem, size);\n  a->len += size;\n}\n\n\n// #line 380 \"/home/spencer/basilisk/src/common.h\"\nArray * array_new()\n{\n  Array * a = ((Array *) pmalloc ((1)*sizeof(Array),__func__,__FILE__,0));\n  a->p = NULL;\n  a->max = a->len = 0;\n  return a;\n}\n\n\n// #line 1410 \"/home/spencer/basilisk/src/common.h\"\nvoid free_solver_func_add (free_solver_func func)\n{\n  if (!free_solver_funcs)\n    free_solver_funcs = array_new();\n  array_append (free_solver_funcs, &func, sizeof(free_solver_func));\n}\n\n\n// #line 154 \"/home/spencer/basilisk/src/common.h\"\nstatic void pmfunc_trace (pmfunc * f, char c)\n{\n  if (pmtrace.fp)\n    fprintf (pmtrace.fp, \"%c %d %ld %ld %ld\",\n      c, f->id, pmtrace.nr, pmtrace.total, f->total);\n@if 1\n  if (pmtrace.nr % 1 == 0) {\n    struct rusage usage;\n    getrusage (RUSAGE_SELF, &usage);\n    if (pmtrace.fp)\n      fprintf (pmtrace.fp, \" %ld\", usage.ru_maxrss*1024);\n    if (!pmtrace.nr)\n      pmtrace.startrss = usage.ru_maxrss;\n    if (usage.ru_maxrss > pmtrace.maxrss)\n      pmtrace.maxrss = usage.ru_maxrss;\n  }\n@endif\n  if (pmtrace.fp)\n    fputc ('\\n', pmtrace.fp);\n  pmtrace.nr++;\n}\n\n\n// #line 135 \"/home/spencer/basilisk/src/common.h\"\nstatic int pmfunc_index (const char * func, const char * file, int line)\n{\n  pmfunc * p = pmfuncs;\n  for (int i = 0; i < pmfuncn; i++, p++)\n    if (p->line == line && !strcmp(func, p->func) && !strcmp(file, p->file))\n      return p->id;\n  pmfuncn++;\n  pmfuncs = (pmfunc *) sysrealloc (pmfuncs, pmfuncn*sizeof(pmfunc));\n  p = &pmfuncs[pmfuncn - 1];\n  memset (p, 0, sizeof(pmfunc));\n  p->func = systrdup(func);\n  p->file = systrdup(file);\n  p->line = line;\n  p->id = pmfuncn;\n  if (pmtrace.fp)\n    fprintf (pmtrace.fp, \"@ %d %s %s %d\\n\", pmfuncn, func, file, line);\n  return pmfuncn;\n}\n@ def not_mpi_compatible()\ndo {\n  if (npe() > 1) {\n    fprintf (ferr, \"%s() is not compatible with MPI (yet)\\n\", __func__);\n    exit (1);\n  }\n} while(0)\n@\n@ define system(command) (pid() == 0 ? system(command) : 0)\n@else\n@ define qstderr() stderr\n@ define qstdout() stdout\n@ define ferr stderr\n@ define fout stdout\n@ define not_mpi_compatible()\n@endif\n\n\n\n\n// #line 93 \"/home/spencer/basilisk/src/common.h\"\nstatic inline void qassert (const char * file, int line, const char * cond) {\n  fprintf (ferr, \"%s:%d: Assertion `%s' failed.\\n\", file, line, cond);\n  abort();\n}\n\n\n// #line 176 \"/home/spencer/basilisk/src/common.h\"\nstatic void * pmfunc_alloc (pmdata * d, size_t size,\n       const char * func, const char * file, int line,\n       char c)\n{\n  if (!(d != NULL)) qassert (\"/home/spencer/basilisk/src/common.h\", 0, \"d != NULL\");\n  OMP (omp critical)\n  {\n    d->id = pmfunc_index(func, file, line);\n    d->size = size;\n    pmfunc * f = &pmfuncs[d->id - 1];\n    f->total += size;\n    if (f->total > f->max)\n      f->max = f->total;\n    pmtrace.total += size;\n    pmtrace.overhead += sizeof(pmdata);\n    if (pmtrace.total > pmtrace.max) {\n      pmtrace.max = pmtrace.total;\n      pmtrace.maxoverhead = pmtrace.overhead;\n    }\n    pmfunc_trace (f, c);\n  }\n  return ((char *)d) + sizeof(pmdata);\n}\n\n\n// #line 242 \"/home/spencer/basilisk/src/common.h\"\nstatic void * pmalloc (size_t size,\n         const char * func, const char * file, int line)\n{\n  return pmfunc_alloc ((pmdata *) sysmalloc (sizeof(pmdata) + size),\n         size, func, file, line, '+');\n}\n\n\n// #line 249 \"/home/spencer/basilisk/src/common.h\"\nstatic void * pcalloc (size_t nmemb, size_t size,\n         const char * func, const char * file, int line)\n{\n  void * p = pmalloc (nmemb*size, func, file, line);\n  return memset (p, 0, nmemb*size);\n}\n\n\n\n\n\n// #line 171 \"/home/spencer/basilisk/src/view.h\"\nbview * bview_new()\n{\n  bview * p = ((bview *) pcalloc (1, sizeof(bview),__func__,__FILE__,0));\n\n  p->tx = p->ty = 0;\n  p->sx = p->sy = p->sz = 1.;\n  p->quat[0] = p->quat[1] = p->quat[2] = 0; p->quat[3] = 1;\n  p->fov = 24.;\n  gl_trackball (p->quat, 0.0, 0.0, 0.0, 0.0);\n\n\n  p->bg[0] = 1; p->bg[1] = 1; p->bg[2] = 1;\n\n\n\n  p->res = 1.;\n  p->lc = 0.004;\n\n  p->samples = 4;\n  p->width = 600*p->samples, p->height = 600*p->samples;\n\n\n  disable_fpe (FE_DIVBYZERO|FE_INVALID);\n\n  p->fb = framebuffer_new (p->width, p->height);\n\n  init_gl();\n  p->active = false;\n\n  enable_fpe (FE_DIVBYZERO|FE_INVALID);\n\n  return p;\n}\n\n\n// #line 232 \"/home/spencer/basilisk/src/view.h\"\nbview * get_view() {\n  if (!_view) {\n    _view = bview_new();\n    free_solver_func_add (destroy_view);\n  }\n  return _view;\n}\n\n\n// #line 61 \"/home/spencer/basilisk/src/vertexbuffer.h\"\nstatic void vertex_buffer_glBegin (unsigned int state)\n{\n  if (VertexBuffer.index) {\n\n    glGetFloatv (0x0BA6, VertexBuffer.modelview);\n\n    bview * view = get_view();\n\n    real q[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,\n      - view->tx, - view->ty, 3, 1 };\n    matrix_multiply (q, VertexBuffer.modelview);\n    for (int i = 0; i < 16; i++)\n      VertexBuffer.modelview[i] = q[i];\n\n    gl_build_rotmatrix ((real (*)[4])q, view->quat);\n    do { real __tmp = q[1]; q[1] = q[4]; q[4] = __tmp; } while(0);\n    do { real __tmp = q[2]; q[2] = q[8]; q[8] = __tmp; } while(0);\n    do { real __tmp = q[6]; q[6] = q[9]; q[9] = __tmp; } while(0);\n    matrix_multiply (q, VertexBuffer.modelview);\n    for (int i = 0; i < 16; i++)\n      VertexBuffer.modelview[i] = q[i];\n\n    VertexBuffer.state = state;\n    switch (state) {\n    case 0x0002:\n      VertexBuffer.line_loop = VertexBuffer.nvertex;\n      break;\n    case 0x0001:\n      VertexBuffer.lines = VertexBuffer.nvertex;\n      break;\n    case 0x0003:\n      VertexBuffer.line_strip = VertexBuffer.nvertex;\n      break;\n    case 0x0007:\n      VertexBuffer.quads = VertexBuffer.nvertex;\n      break;\n    case 0x0009:\n      VertexBuffer.polygon = VertexBuffer.nvertex;\n      break;\n    case 0x0006:\n      VertexBuffer.fan = VertexBuffer.nvertex;\n      break;\n    default:\n      fprintf (ferr, \"glBegin (%d) not implemented yet\\n\", state);\n      break;\n    }\n  }\n  else\n    glBegin (state);\n}"),_("\n   \n// #line 1259 \"/home/spencer/basilisk/src/draw.h\"\nif (val(f,0,0,0) != 1e30) {\n     vertex_buffer_glBegin (0x0006);\n     if (color && linear && col.i >= 0) { if (VertexBuffer.color) { Color b = colormap_color (cmap, (4.*val(f,0,0,0) + 2.*(val(f,1,0,0) + val(f,-1,0,0) + val(f,0,1,0) + val(f,0,-1,0)) + val(f,-1,-1,0) + val(f,1,1,0) + val(f,-1,1,0) + val(f,1,-1,0))/16., min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { real _v = (4.*val(f,0,0,0) + 2.*(val(f,1,0,0) + val(f,-1,0,0) + val(f,0,1,0) + val(f,0,-1,0)) + val(f,-1,-1,0) + val(f,1,1,0) + val(f,-1,1,0) + val(f,1,-1,0))/16.; if (max > min) glTexCoord1d (clamp(((_v) - min)/(max - min), 0., 1.)); else glTexCoord1d (0.); } }\n\n                                                 ;\n     glvertex2d (view, x, y);\n     if (color && linear && col.i >= 0) { if (VertexBuffer.color) { Color b = colormap_color (cmap, (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,-1,0) + val(f,0,-1,0))/4., min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { real _v = (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,-1,0) + val(f,0,-1,0))/4.; if (max > min) glTexCoord1d (clamp(((_v) - min)/(max - min), 0., 1.)); else glTexCoord1d (0.); } };\n     glvertex2d (view, x - Delta_x/2., y - Delta_y/2.);\n     if (color && linear && col.i >= 0) { if (VertexBuffer.color) { Color b = colormap_color (cmap, (val(f,0,0,0) + val(f,1,0,0) + val(f,1,-1,0) + val(f,0,-1,0))/4., min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { real _v = (val(f,0,0,0) + val(f,1,0,0) + val(f,1,-1,0) + val(f,0,-1,0))/4.; if (max > min) glTexCoord1d (clamp(((_v) - min)/(max - min), 0., 1.)); else glTexCoord1d (0.); } };\n     glvertex2d (view, x + Delta_x/2., y - Delta_y/2.);\n     if (color && linear && col.i >= 0) { if (VertexBuffer.color) { Color b = colormap_color (cmap, (val(f,0,0,0) + val(f,1,0,0) + val(f,1,1,0) + val(f,0,1,0))/4., min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { real _v = (val(f,0,0,0) + val(f,1,0,0) + val(f,1,1,0) + val(f,0,1,0))/4.; if (max > min) glTexCoord1d (clamp(((_v) - min)/(max - min), 0., 1.)); else glTexCoord1d (0.); } };\n     glvertex2d (view, x + Delta_x/2., y + Delta_y/2.);\n     if (color && linear && col.i >= 0) { if (VertexBuffer.color) { Color b = colormap_color (cmap, (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,1,0) + val(f,0,1,0))/4., min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { real _v = (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,1,0) + val(f,0,1,0))/4.; if (max > min) glTexCoord1d (clamp(((_v) - min)/(max - min), 0., 1.)); else glTexCoord1d (0.); } };\n     glvertex2d (view, x - Delta_x/2., y + Delta_y/2.);\n     if (color && linear && col.i >= 0) { if (VertexBuffer.color) { Color b = colormap_color (cmap, (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,-1,0) + val(f,0,-1,0))/4., min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { real _v = (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,-1,0) + val(f,0,-1,0))/4.; if (max > min) glTexCoord1d (clamp(((_v) - min)/(max - min), 0., 1.)); else glTexCoord1d (0.); } };\n     glvertex2d (view, x - Delta_x/2., y - Delta_y/2.);\n     vertex_buffer_glEnd();\n     view->ni++;\n   }")})
   {_stencil_val(f,0,0,0); { 
     
     if (color && linear && col.i >= 0) { if (VertexBuffer.color) {                  _stencil_val(f,1,-1,0); _stencil_val(f,-1,1,0); _stencil_val(f,1,1,0); _stencil_val(f,-1,-1,0); _stencil_val(f,0,-1,0); _stencil_val(f,0,1,0); _stencil_val(f,-1,0,0);_stencil_val(f,1,0,0);_stencil_val(f,0,0,0);     } else {              _stencil_val(f,1,-1,0); _stencil_val(f,-1,1,0); _stencil_val(f,1,1,0); _stencil_val(f,-1,-1,0); _stencil_val(f,0,-1,0); _stencil_val(f,0,1,0); _stencil_val(f,-1,0,0);_stencil_val(f,1,0,0);_stencil_val(f,0,0,0);                } }   

                                                 
     
     if (color && linear && col.i >= 0) { if (VertexBuffer.color) {            _stencil_val(f,0,-1,0); _stencil_val(f,-1,-1,0); _stencil_val(f,-1,0,0);_stencil_val(f,0,0,0);     } else {        _stencil_val(f,0,-1,0); _stencil_val(f,-1,-1,0); _stencil_val(f,-1,0,0);_stencil_val(f,0,0,0);                } }       
     
     if (color && linear && col.i >= 0) { if (VertexBuffer.color) {            _stencil_val(f,0,-1,0); _stencil_val(f,1,-1,0); _stencil_val(f,1,0,0);_stencil_val(f,0,0,0);     } else {        _stencil_val(f,0,-1,0); _stencil_val(f,1,-1,0); _stencil_val(f,1,0,0);_stencil_val(f,0,0,0);                } }       
     
     if (color && linear && col.i >= 0) { if (VertexBuffer.color) {            _stencil_val(f,0,1,0); _stencil_val(f,1,1,0); _stencil_val(f,1,0,0);_stencil_val(f,0,0,0);     } else {        _stencil_val(f,0,1,0); _stencil_val(f,1,1,0); _stencil_val(f,1,0,0);_stencil_val(f,0,0,0);                } }       
     
     if (color && linear && col.i >= 0) { if (VertexBuffer.color) {            _stencil_val(f,0,1,0); _stencil_val(f,-1,1,0); _stencil_val(f,-1,0,0);_stencil_val(f,0,0,0);     } else {        _stencil_val(f,0,1,0); _stencil_val(f,-1,1,0); _stencil_val(f,-1,0,0);_stencil_val(f,0,0,0);                } }       
     
     if (color && linear && col.i >= 0) { if (VertexBuffer.color) {            _stencil_val(f,0,-1,0); _stencil_val(f,-1,-1,0); _stencil_val(f,-1,0,0);_stencil_val(f,0,0,0);     } else {        _stencil_val(f,0,-1,0); _stencil_val(f,-1,-1,0); _stencil_val(f,-1,0,0);_stencil_val(f,0,0,0);                } }       
     
     
     
   }   }end_foreach_visible_stencil();
 {
#line 1258
foreach_visible (view)
   if (val(f,0,0,0) != 1e30) {
     vertex_buffer_glBegin (0x0006);
     if (color && linear && col.i >= 0) { if (VertexBuffer.color) { Color b = colormap_color (cmap, (4.*val(f,0,0,0) + 2.*(val(f,1,0,0) + val(f,-1,0,0) + val(f,0,1,0) + val(f,0,-1,0)) + val(f,-1,-1,0) + val(f,1,1,0) + val(f,-1,1,0) + val(f,1,-1,0))/16., min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (4.*val(f,0,0,0) + 2.*(val(f,1,0,0) + val(f,-1,0,0) + val(f,0,1,0) + val(f,0,-1,0)) + val(f,-1,-1,0) + val(f,1,1,0) + val(f,-1,1,0) + val(f,1,-1,0))/16.; if (max > min) glTexCoord1d (clamp(((_v) - min)/(max - min), 0., 1.)); else glTexCoord1d (0.); } }

                                                 ;
     glvertex2d (view, x, y);
     if (color && linear && col.i >= 0) { if (VertexBuffer.color) { Color b = colormap_color (cmap, (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,-1,0) + val(f,0,-1,0))/4., min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,-1,0) + val(f,0,-1,0))/4.; if (max > min) glTexCoord1d (clamp(((_v) - min)/(max - min), 0., 1.)); else glTexCoord1d (0.); } };
     glvertex2d (view, x - Delta_x/2., y - Delta_y/2.);
     if (color && linear && col.i >= 0) { if (VertexBuffer.color) { Color b = colormap_color (cmap, (val(f,0,0,0) + val(f,1,0,0) + val(f,1,-1,0) + val(f,0,-1,0))/4., min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,1,0,0) + val(f,1,-1,0) + val(f,0,-1,0))/4.; if (max > min) glTexCoord1d (clamp(((_v) - min)/(max - min), 0., 1.)); else glTexCoord1d (0.); } };
     glvertex2d (view, x + Delta_x/2., y - Delta_y/2.);
     if (color && linear && col.i >= 0) { if (VertexBuffer.color) { Color b = colormap_color (cmap, (val(f,0,0,0) + val(f,1,0,0) + val(f,1,1,0) + val(f,0,1,0))/4., min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,1,0,0) + val(f,1,1,0) + val(f,0,1,0))/4.; if (max > min) glTexCoord1d (clamp(((_v) - min)/(max - min), 0., 1.)); else glTexCoord1d (0.); } };
     glvertex2d (view, x + Delta_x/2., y + Delta_y/2.);
     if (color && linear && col.i >= 0) { if (VertexBuffer.color) { Color b = colormap_color (cmap, (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,1,0) + val(f,0,1,0))/4., min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,1,0) + val(f,0,1,0))/4.; if (max > min) glTexCoord1d (clamp(((_v) - min)/(max - min), 0., 1.)); else glTexCoord1d (0.); } };
     glvertex2d (view, x - Delta_x/2., y + Delta_y/2.);
     if (color && linear && col.i >= 0) { if (VertexBuffer.color) { Color b = colormap_color (cmap, (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,-1,0) + val(f,0,-1,0))/4., min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,-1,0) + val(f,0,-1,0))/4.; if (max > min) glTexCoord1d (clamp(((_v) - min)/(max - min), 0., 1.)); else glTexCoord1d (0.); } };
     glvertex2d (view, x - Delta_x/2., y - Delta_y/2.);
     vertex_buffer_glEnd();
     view->ni++;
   }end_foreach_visible();}
      }
      else
 {foreach_leaf()
   if (val(f,0,0,0) != 1e30) {
     vertex_buffer_glBegin (0x0006);
     if (color && linear && col.i >= 0) { if (VertexBuffer.color) { Color b = colormap_color (cmap, (4.*val(f,0,0,0) + 2.*(val(f,1,0,0) + val(f,-1,0,0) + val(f,0,1,0) + val(f,0,-1,0)) + val(f,-1,-1,0) + val(f,1,1,0) + val(f,-1,1,0) + val(f,1,-1,0))/16., min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (4.*val(f,0,0,0) + 2.*(val(f,1,0,0) + val(f,-1,0,0) + val(f,0,1,0) + val(f,0,-1,0)) + val(f,-1,-1,0) + val(f,1,1,0) + val(f,-1,1,0) + val(f,1,-1,0))/16.; if (max > min) glTexCoord1d (clamp(((_v) - min)/(max - min), 0., 1.)); else glTexCoord1d (0.); } }

                                                 ;
     glvertex_normal3d (view, point, fn, x, y, val(Z,0,0,0));
     if (color && linear && col.i >= 0) { if (VertexBuffer.color) { Color b = colormap_color (cmap, (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,-1,0) + val(f,0,-1,0))/4., min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,-1,0) + val(f,0,-1,0))/4.; if (max > min) glTexCoord1d (clamp(((_v) - min)/(max - min), 0., 1.)); else glTexCoord1d (0.); } };
     glvertex_normal3d (view, point, fn, x - Delta_x/2., y - Delta_y/2.,
          (val(Z,0,0,0) + val(Z,-1,0,0) + val(Z,-1,-1,0) + val(Z,0,-1,0))/4.);
     if (color && linear && col.i >= 0) { if (VertexBuffer.color) { Color b = colormap_color (cmap, (val(f,0,0,0) + val(f,1,0,0) + val(f,1,-1,0) + val(f,0,-1,0))/4., min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,1,0,0) + val(f,1,-1,0) + val(f,0,-1,0))/4.; if (max > min) glTexCoord1d (clamp(((_v) - min)/(max - min), 0., 1.)); else glTexCoord1d (0.); } };
     glvertex_normal3d (view, point, fn, x + Delta_x/2., y - Delta_y/2.,
          (val(Z,0,0,0) + val(Z,1,0,0) + val(Z,1,-1,0) + val(Z,0,-1,0))/4.);
     if (color && linear && col.i >= 0) { if (VertexBuffer.color) { Color b = colormap_color (cmap, (val(f,0,0,0) + val(f,1,0,0) + val(f,1,1,0) + val(f,0,1,0))/4., min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,1,0,0) + val(f,1,1,0) + val(f,0,1,0))/4.; if (max > min) glTexCoord1d (clamp(((_v) - min)/(max - min), 0., 1.)); else glTexCoord1d (0.); } };
     glvertex_normal3d (view, point, fn, x + Delta_x/2., y + Delta_y/2.,
          (val(Z,0,0,0) + val(Z,1,0,0) + val(Z,1,1,0) + val(Z,0,1,0))/4.);
     if (color && linear && col.i >= 0) { if (VertexBuffer.color) { Color b = colormap_color (cmap, (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,1,0) + val(f,0,1,0))/4., min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,1,0) + val(f,0,1,0))/4.; if (max > min) glTexCoord1d (clamp(((_v) - min)/(max - min), 0., 1.)); else glTexCoord1d (0.); } };
     glvertex_normal3d (view, point, fn, x - Delta_x/2., y + Delta_y/2.,
          (val(Z,0,0,0) + val(Z,-1,0,0) + val(Z,-1,1,0) + val(Z,0,1,0))/4.);
     if (color && linear && col.i >= 0) { if (VertexBuffer.color) { Color b = colormap_color (cmap, (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,-1,0) + val(f,0,-1,0))/4., min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,-1,0) + val(f,0,-1,0))/4.; if (max > min) glTexCoord1d (clamp(((_v) - min)/(max - min), 0., 1.)); else glTexCoord1d (0.); } };
     glvertex_normal3d (view, point, fn, x - Delta_x/2., y - Delta_y/2.,
          (val(Z,0,0,0) + val(Z,-1,0,0) + val(Z,-1,-1,0) + val(Z,0,-1,0))/4.);
     vertex_buffer_glEnd();
     view->ni++;
   }end_foreach_leaf();}
# 1331 "/home/spencer/basilisk/src/draw.h"
    }end_colorized();}
  }
  else {

    vertex_buffer_glNormal3d (0, 0, view->reversed ? -1 : 1);
    vertex_buffer_glBegin (0x0007);
    foreach_visible_stencil (1,{(NonLocal[]){{"pmtrace","not implemented yet",(void *)&pmtrace,NULL,0},{"pmfuncn","int",(void *)&pmfuncn,NULL,0},{"pmfuncs","pmfunc",(void *)pmfuncs,NULL,1},{"VertexBuffer","not implemented yet",(void *)&VertexBuffer,NULL,0},{"max","double",(void *)&max,NULL,0},{"min","double",(void *)&min,NULL,0},{"cmap","double",(void *)cmap,(int[]){1273,0},0},{"ferr","not implemented yet",(void *)ferr,NULL,1},{"col","scalar",(void *)&col,NULL,0},{"linear","bool",(void *)&linear,NULL,0},{"color","not implemented yet",(void *)color,NULL,1},{"f","scalar",(void *)&f,NULL,0},{"view","not implemented yet",(void *)view,NULL,1},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_("\n\n\n// #line 222 \"/home/spencer/basilisk/src/vertexbuffer.h\"\nstatic void vertex_buffer_glVertex2d (real x, real y)\n{\n  if (VertexBuffer.position) {\n    if (VertexBuffer.dim < 2)\n      VertexBuffer.dim = 2;\n    real v[4] = {x, y, 0, 1.};\n    vector_multiply (v, VertexBuffer.modelview);\n    array_append (VertexBuffer.position, v, 3*sizeof(real));\n    VertexBuffer.nvertex++;\n  }\n  else\n    glVertex3d (x, y, 0.);\n}\n\n\n\n// #line 838 \"/home/spencer/basilisk/src/draw.h\"\nstatic void glvertex2d (bview * view, real x, real y) {\n  if (view->map) {\n    coord p = {x, y, 0.};\n    view->map (&p);\n    vertex_buffer_glVertex2d (p.x, p.y);\n  }\n  else\n    vertex_buffer_glVertex2d (x, y);\n}\n\n\n// #line 135 \"/home/spencer/basilisk/src/common.h\"\nstatic int pmfunc_index (const char * func, const char * file, int line)\n{\n  pmfunc * p = pmfuncs;\n  for (int i = 0; i < pmfuncn; i++, p++)\n    if (p->line == line && !strcmp(func, p->func) && !strcmp(file, p->file))\n      return p->id;\n  pmfuncn++;\n  pmfuncs = (pmfunc *) sysrealloc (pmfuncs, pmfuncn*sizeof(pmfunc));\n  p = &pmfuncs[pmfuncn - 1];\n  memset (p, 0, sizeof(pmfunc));\n  p->func = systrdup(func);\n  p->file = systrdup(file);\n  p->line = line;\n  p->id = pmfuncn;\n  if (pmtrace.fp)\n    fprintf (pmtrace.fp, \"@ %d %s %s %d\\n\", pmfuncn, func, file, line);\n  return pmfuncn;\n}\n\n\n// #line 176 \"/home/spencer/basilisk/src/common.h\"\nstatic void * pmfunc_alloc (pmdata * d, size_t size,\n       const char * func, const char * file, int line,\n       char c)\n{\n  if (!(d != NULL)) qassert (\"/home/spencer/basilisk/src/common.h\", 0, \"d != NULL\");\n  OMP (omp critical)\n  {\n    d->id = pmfunc_index(func, file, line);\n    d->size = size;\n    pmfunc * f = &pmfuncs[d->id - 1];\n    f->total += size;\n    if (f->total > f->max)\n      f->max = f->total;\n    pmtrace.total += size;\n    pmtrace.overhead += sizeof(pmdata);\n    if (pmtrace.total > pmtrace.max) {\n      pmtrace.max = pmtrace.total;\n      pmtrace.maxoverhead = pmtrace.overhead;\n    }\n    pmfunc_trace (f, c);\n  }\n  return ((char *)d) + sizeof(pmdata);\n}\n\n\n// #line 154 \"/home/spencer/basilisk/src/common.h\"\nstatic void pmfunc_trace (pmfunc * f, char c)\n{\n  if (pmtrace.fp)\n    fprintf (pmtrace.fp, \"%c %d %ld %ld %ld\",\n      c, f->id, pmtrace.nr, pmtrace.total, f->total);\n@if 1\n  if (pmtrace.nr % 1 == 0) {\n    struct rusage usage;\n    getrusage (RUSAGE_SELF, &usage);\n    if (pmtrace.fp)\n      fprintf (pmtrace.fp, \" %ld\", usage.ru_maxrss*1024);\n    if (!pmtrace.nr)\n      pmtrace.startrss = usage.ru_maxrss;\n    if (usage.ru_maxrss > pmtrace.maxrss)\n      pmtrace.maxrss = usage.ru_maxrss;\n  }\n@endif\n  if (pmtrace.fp)\n    fputc ('\\n', pmtrace.fp);\n  pmtrace.nr++;\n}\n\n\n// #line 200 \"/home/spencer/basilisk/src/common.h\"\nstatic void * pmfunc_free (void * ptr, char c)\n{\n  if (!ptr)\n    return ptr;\n  pmdata * d = (pmdata *) (((char *)ptr) - sizeof(pmdata));\n  if (d->id < 1 || d->id > pmfuncn) {\n    fputs (\"*** MTRACE: ERROR!: corrupted free()\", ferr);\n    if (d->size == 0)\n      fputs (\", possible double free()\", ferr);\n    else\n      fputs (\", not traced?\", ferr);\n    fputs (\", aborting...\\n\", ferr);\n    abort();\n    return ptr;\n  }\n  else\n  OMP (omp critical)\n  {\n    pmfunc * f = &pmfuncs[d->id - 1];\n    if (f->total < d->size) {\n      fprintf (ferr, \"*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\\n\",\n        f->total, d->size);\n      abort();\n    }\n    else\n      f->total -= d->size;\n    if (pmtrace.total < d->size) {\n      fprintf (ferr, \"*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\\n\",\n        pmtrace.total, d->size);\n      abort();\n    }\n    else {\n      pmtrace.total -= d->size;\n      pmtrace.overhead -= sizeof(pmdata);\n    }\n    d->id = 0;\n    d->size = 0;\n    pmfunc_trace (f, c);\n  }\n  return d;\n}\n\n\n// #line 256 \"/home/spencer/basilisk/src/common.h\"\nstatic void * prealloc (void * ptr, size_t size,\n   const char * func, const char * file, int line)\n{\n  return pmfunc_alloc ((pmdata *) sysrealloc (pmfunc_free(ptr, '<'),\n           sizeof(pmdata) + size),\n         size, func, file, line, '>');\n}\n\n\n// #line 394 \"/home/spencer/basilisk/src/common.h\"\nvoid array_append (Array * a, void * elem, size_t size)\n{\n  if (a->len + size >= a->max) {\n    a->max += max (size, 4096);\n    a->p = prealloc (a->p, a->max,__func__,__FILE__,0);\n  }\n  memcpy (((char *)a->p) + a->len, elem, size);\n  a->len += size;\n}\n\n\n// #line 188 \"/home/spencer/basilisk/src/vertexbuffer.h\"\nstatic void vertex_buffer_glColor3f (real r, real g, real b)\n{\n  if (VertexBuffer.color) {\n    struct { real x, y, z; } color = {r, g, b};\n    array_append (VertexBuffer.color, &color, 3*sizeof(real));\n  }\n  else\n    glColor3f (r, g, b);\n}\n@ def not_mpi_compatible()\ndo {\n  if (npe() > 1) {\n    fprintf (ferr, \"%s() is not compatible with MPI (yet)\\n\", __func__);\n    exit (1);\n  }\n} while(0)\n@\n@ define system(command) (pid() == 0 ? system(command) : 0)\n@else\n@ define qstderr() stderr\n@ define qstdout() stdout\n@ define ferr stderr\n@ define fout stdout\n@ define not_mpi_compatible()\n@endif\n\n\n\n\n// #line 93 \"/home/spencer/basilisk/src/common.h\"\nstatic inline void qassert (const char * file, int line, const char * cond) {\n  fprintf (ferr, \"%s:%d: Assertion `%s' failed.\\n\", file, line, cond);\n  abort();\n}\n\n\n// #line 261 \"/home/spencer/basilisk/src/output.h\"\nColor colormap_color (real cmap[127][3],\n        real val, real min, real max)\n{\n  Color c;\n  if (val == 1e30) {\n    c.r = c.g = c.b = 0;\n    return c;\n  }\n  int i;\n  real coef;\n  if (max != min)\n    val = (val - min)/(max - min);\n  else\n    val = 0.;\n  if (val <= 0.) i = 0, coef = 0.;\n  else if (val >= 1.) i = 127 - 2, coef = 1.;\n  else {\n    i = val*(127 - 1);\n    coef = val*(127 - 1) - i;\n  }\n  if (!(i >= 0 && i < 127 - 1)) qassert (\"/home/spencer/basilisk/src/output.h\", 0, \"i >= 0 && i < NCMAP - 1\");\n  unsigned char * c1 = (unsigned char *) &c;\n  for (int j = 0; j < 3; j++)\n    c1[j] = 255*(cmap[i][j]*(1. - coef) + cmap[i + 1][j]*coef);\n  return c;\n}"),_("\n      \n// #line 1338 \"/home/spencer/basilisk/src/draw.h\"\nif (val(f,0,0,0) != 1e30) {\n if (color && (!linear || col.i < 0)) { Color b = colormap_color (cmap, col.i < 0 ? (real) level : val(col,0,0,0), min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); };\n glvertex2d (view, x - Delta_x/2., y - Delta_y/2.);\n if (color && (!linear || col.i < 0)) { Color b = colormap_color (cmap, col.i < 0 ? (real) level : val(col,0,0,0), min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); };\n glvertex2d (view, x + Delta_x/2., y - Delta_y/2.);\n if (color && (!linear || col.i < 0)) { Color b = colormap_color (cmap, col.i < 0 ? (real) level : val(col,0,0,0), min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); };\n glvertex2d (view, x + Delta_x/2., y + Delta_y/2.);\n if (color && (!linear || col.i < 0)) { Color b = colormap_color (cmap, col.i < 0 ? (real) level : val(col,0,0,0), min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); };\n glvertex2d (view, x - Delta_x/2., y + Delta_y/2.);\n view->ni++;\n      }")})
      {_stencil_val(f,0,0,0); {
 if (color && (!linear || col.i < 0)) {               _stencil_val(col,0,0,0);     }       
 
 if (color && (!linear || col.i < 0)) {               _stencil_val(col,0,0,0);     }       
 
 if (color && (!linear || col.i < 0)) {               _stencil_val(col,0,0,0);     }       
 
 if (color && (!linear || col.i < 0)) {               _stencil_val(col,0,0,0);     }       
 
 
      }   }end_foreach_visible_stencil();
    {
#line 1337
foreach_visible (view)
      if (val(f,0,0,0) != 1e30) {
 if (color && (!linear || col.i < 0)) { Color b = colormap_color (cmap, col.i < 0 ? (double) level : val(col,0,0,0), min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); };
 glvertex2d (view, x - Delta_x/2., y - Delta_y/2.);
 if (color && (!linear || col.i < 0)) { Color b = colormap_color (cmap, col.i < 0 ? (double) level : val(col,0,0,0), min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); };
 glvertex2d (view, x + Delta_x/2., y - Delta_y/2.);
 if (color && (!linear || col.i < 0)) { Color b = colormap_color (cmap, col.i < 0 ? (double) level : val(col,0,0,0), min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); };
 glvertex2d (view, x + Delta_x/2., y + Delta_y/2.);
 if (color && (!linear || col.i < 0)) { Color b = colormap_color (cmap, col.i < 0 ? (double) level : val(col,0,0,0), min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); };
 glvertex2d (view, x - Delta_x/2., y + Delta_y/2.);
 view->ni++;
      }end_foreach_visible();}
    vertex_buffer_glEnd();
# 1367 "/home/spencer/basilisk/src/draw.h"
  }
  if (expr) delete (((scalar[]){col,{-1}}));

  if (zexpr) delete (((scalar[]){Z,{-1}}));
  if (z) delete ((scalar *)((vector[]){fn,{{-1},{-1}}}));

  {end_tracing("squares","/home/spencer/basilisk/src/draw.h",0);return true;}
end_tracing("squares","/home/spencer/basilisk/src/draw.h",0);}
# 1385 "/home/spencer/basilisk/src/draw.h"
     
bool box (bool notics, float lc[3], float lw)
{tracing("box","/home/spencer/basilisk/src/draw.h",0);
  bview * view = draw();
  {begin_draw_lines (view, lc, lw); {

    float height = 0.5*gl_StrokeHeight();
    float width = gl_StrokeWidth ('1'), scale = L0/(60.*width), length;
    float Z1 = 2 == 2 ? 0. : Z0;
    char label[80];

    glMatrixMode (0x1700);

    if (!notics) {
      int nt = 8;
      for (int i = 0; i <= nt; i++) {
 glPushMatrix();
 glTranslatef (X0 + i*L0/nt - height/2.*scale, Y0 - width/3.*scale, Z1);
 glRotatef (-90, 0, 0, 1);
 glScalef (scale, scale, 1.);
 sprintf (label, "%g", X0 + i*L0/nt);
 gl_StrokeString (label);
 glPopMatrix();

 glPushMatrix();
 sprintf (label, "%g", Y0 + i*L0/nt);
 length = gl_StrokeLength (label);
 glTranslatef (X0 - (length + width/3.)*scale,
        Y0 + i*L0/nt - height/2.*scale, Z1);
 glScalef (scale, scale, 1.);
 gl_StrokeString (label);
 glPopMatrix();
# 1429 "/home/spencer/basilisk/src/draw.h"
      }

      glPushMatrix();
      sprintf (label, "%g", X0 + L0/2.);
      length = gl_StrokeLength (label);
      glTranslatef (X0 + L0/2 - height*scale, Y0 - (length + 4.*width)*scale, Z1);
      glScalef (2.*scale, 2.*scale, 1.);
      gl_StrokeString ("X");
      glPopMatrix();


      glPushMatrix();
      sprintf (label, "%g", Y0 + L0/2.);
      length = gl_StrokeLength (label);
      glTranslatef (X0 - (length + 4.*width)*scale,
      Y0 + L0/2. - height*scale, Z1);
      glScalef (2.*scale, 2.*scale, 1.);
      gl_StrokeString ("Y");
      glPopMatrix();
# 1460 "/home/spencer/basilisk/src/draw.h"
    }


    {foreach_level (0) {
      vertex_buffer_glBegin (0x0002);
      glvertex2d (view, x - Delta_x/2., y - Delta_y/2.);
      glvertex2d (view, x + Delta_x/2., y - Delta_y/2.);
      glvertex2d (view, x + Delta_x/2., y + Delta_y/2.);
      glvertex2d (view, x - Delta_x/2., y + Delta_y/2.);
      vertex_buffer_glEnd ();
      view->ni++;
    }end_foreach_level();}
# 1492 "/home/spencer/basilisk/src/draw.h"
  }end_draw_lines();}
  {end_tracing("box","/home/spencer/basilisk/src/draw.h",0);return true;}
end_tracing("box","/home/spencer/basilisk/src/draw.h",0);}
# 1507 "/home/spencer/basilisk/src/draw.h"
     
bool isosurface (char * f,
   double v,

   char * color,
   double min, double max, double spread,
   bool linear,
   Colormap map,
   float fc[3], float lc[3], float lw,
   bool expr,
   bool cbar, float size, float pos[2], char * label, double lscale, bool horizontal, bool border, bool mid, float fsize, char * format, int levels)
{tracing("isosurface","/home/spencer/basilisk/src/draw.h",0);
# 1577 "/home/spencer/basilisk/src/draw.h"
  {end_tracing("isosurface","/home/spencer/basilisk/src/draw.h",0);return true;}
end_tracing("isosurface","/home/spencer/basilisk/src/draw.h",0);}
# 1591 "/home/spencer/basilisk/src/draw.h"
void travelling (double start, double end,
   float tx, float ty, float quat[4], float fov)
{
  static float stx, sty, squat[4], sfov;
  static double told = -1.;
  if (told < start && t >= start) {
    bview * view = get_view();
    stx = view->tx, sty = view->ty, sfov = view->fov;
    for (int i = 0; i < 4; i++)
      squat[i] = view->quat[i];
  }
  if (t >= start && t <= end)
    view ( (!tx ? stx : ((t - start)*(tx) + (end - t)*(stx))/(end - start)), (!ty ? sty : ((t - start)*(ty) + (end - t)*(sty))/(end - start))
, (!fov ? sfov : ((t - start)*(fov) + (end - t)*(sfov))/(end - start))
,
#line 51
(
    
#line 51
float[4]) 
#line 1605
{(!quat[0] ? squat[0] : ((t - start)*(quat[0]) + (end - t)*(squat[0]))/(end - start)), (!quat[1] ? squat[1] : ((t - start)*(quat[1]) + (end - t)*(squat[1]))/(end - start)),
           (!quat[2] ? squat[2] : ((t - start)*(quat[2]) + (end - t)*(squat[2]))/(end - start)), (!quat[3] ? squat[3] : ((t - start)*(quat[3]) + (end - t)*(squat[3]))/(end - start))}
#line 51
, 
1., 1., 1., 
800, 800, 4,
(
    
#line 54
float[3]) {0}, 
0., 0., 0., 
false, 
0., 0., 0., 
0., 
NULL, 
NULL, 
0, 
0., 0., 0., 0., 
NULL
#line 1606
);
  if (told < end && t >= end) {
    bview * view = get_view();
    stx = view->tx, sty = view->ty, sfov = view->fov;
    for (int i = 0; i < 4; i++)
      squat[i] = view->quat[i];
  }
  told = t;
}
# 1631 "/home/spencer/basilisk/src/draw.h"
     
bool draw_string (char * str,
    int pos,
    float size,
    float lc[3], float lw)

{tracing("draw_string","/home/spencer/basilisk/src/draw.h",0);
  bview * view = draw();

  glMatrixMode (0x1701);
  glPushMatrix();
  glLoadIdentity();

  glMatrixMode (0x1700);
  glPushMatrix();
  glLoadIdentity();

  vertex_buffer_glColor3f (lc[0], lc[1], lc[2]);
  glLineWidth (view->samples*(lw > 0. ? lw : 1.));

  float width = gl_StrokeWidth ('1'), height = gl_StrokeHeight();
  if (!size)
    size = 40;
  float hscale = 2./(size*width), vscale = hscale*view->width/view->height;
  float vmargin = width/2.*vscale;
  if (pos == 0)
    glTranslatef (-1., -1. + vmargin, 0.);
  else if (pos == 1)
    glTranslatef (-1., 1. - height*vscale, 0.);
  else if (pos == 2)
    glTranslatef (1. - strlen(str)*width*hscale, 1. - height*vscale, 0.);
  else
    glTranslatef (1. - strlen(str)*width*hscale, -1. + vmargin, 0.);
  glScalef (hscale, vscale, 1.);
  gl_StrokeString (str);

  glMatrixMode (0x1700);
  glPopMatrix();
  glMatrixMode (0x1701);
  glPopMatrix();

  {end_tracing("draw_string","/home/spencer/basilisk/src/draw.h",0);return true;}
end_tracing("draw_string","/home/spencer/basilisk/src/draw.h",0);}




     
bool labels (char * f,
      float lc[3], float lw)
{tracing("labels","/home/spencer/basilisk/src/draw.h",0);

  bool expr = false;
  scalar ff = compile_expression (f, &expr);
  if (ff.i < 0)
    {end_tracing("labels","/home/spencer/basilisk/src/draw.h",0);return false;}
  bview * view = draw();
  float width = gl_StrokeWidth ('1'), height = gl_StrokeHeight();
  float res = view->res;
  if (view->res < 150*view->samples)
    view->res = 150*view->samples;
  {begin_draw_lines (view, lc, lw); {
    glMatrixMode (0x1700);
    foreach_visible_stencil (1,{(NonLocal[]){{"baseblock","scalar",(void *)baseblock,NULL,1},{"all","scalar",(void *)all,NULL,1},{"free_solver_funcs","Array",(void *)free_solver_funcs,NULL,1},{"pmtrace","not implemented yet",(void *)&pmtrace,NULL,0},{"pmfuncn","int",(void *)&pmfuncn,NULL,0},{"pmfuncs","pmfunc",(void *)pmfuncs,NULL,1},{"_view","not implemented yet",(void *)_view,NULL,1},{"VertexBuffer","not implemented yet",(void *)&VertexBuffer,NULL,0},{"ferr","not implemented yet",(void *)ferr,NULL,1},{"ogStrokeMonoRoman","not implemented yet",(void *)&ogStrokeMonoRoman,NULL,0},{"height","float",(void *)&height,NULL,0},{"width","float",(void *)&width,NULL,0},{"ff","scalar",(void *)&ff,NULL,0},{"view","not implemented yet",(void *)view,NULL,1},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_("\n\n\n// #line 32 \"/home/spencer/basilisk/src/vertexbuffer.h\"\nstatic void vertex_buffer_push_index (unsigned int i)\n{\n  i -= VertexBuffer.vertex;\n  array_append (VertexBuffer.index, &i, sizeof(unsigned int));\n}\n\n\n// #line 112 \"/home/spencer/basilisk/src/vertexbuffer.h\"\nstatic void vertex_buffer_glEnd()\n{\n  if (VertexBuffer.index) {\n    int type = -1;\n    switch (VertexBuffer.state) {\n\n    case 0x0002:\n      for (int i = VertexBuffer.line_loop; i < VertexBuffer.nvertex - 1; i++) {\n vertex_buffer_push_index (i);\n vertex_buffer_push_index (i + 1);\n      }\n      vertex_buffer_push_index (VertexBuffer.nvertex - 1);\n      vertex_buffer_push_index (VertexBuffer.line_loop);\n      type = 0;\n      break;\n\n    case 0x0001:\n      for (int i = VertexBuffer.lines; i < VertexBuffer.nvertex; i += 2) {\n vertex_buffer_push_index (i);\n vertex_buffer_push_index (i + 1);\n      }\n      type = 0;\n      break;\n\n    case 0x0003:\n      for (int i = VertexBuffer.line_strip; i < VertexBuffer.nvertex - 1; i++) {\n vertex_buffer_push_index (i);\n vertex_buffer_push_index (i + 1);\n      }\n      type = 0;\n      break;\n\n    case 0x0007:\n      for (int i = VertexBuffer.quads; i < VertexBuffer.nvertex; i += 4)\n for (int j = 1; j <= 2; j++) {\n   vertex_buffer_push_index (i);\n   vertex_buffer_push_index (i + j);\n   vertex_buffer_push_index (i + j + 1);\n }\n      type = 1;\n      break;\n\n    case 0x0009:\n      for (int j = 1; j <= VertexBuffer.nvertex - VertexBuffer.polygon - 2;\n    j++) {\n vertex_buffer_push_index (VertexBuffer.polygon);\n vertex_buffer_push_index (VertexBuffer.polygon + j);\n vertex_buffer_push_index (VertexBuffer.polygon + j + 1);\n      }\n      type = 1;\n      break;\n\n    case 0x0006:\n      for (int i = VertexBuffer.fan + 1; i < VertexBuffer.nvertex - 1; i++) {\n vertex_buffer_push_index (VertexBuffer.fan);\n vertex_buffer_push_index (i);\n vertex_buffer_push_index (i + 1);\n      }\n      type = 1;\n      break;\n\n    default:\n      break;\n    }\n    VertexBuffer.state = 0;\n    if (VertexBuffer.type >= 0 && type >= 0) {\n\n      if (!(VertexBuffer.type == type)) qassert (\"/home/spencer/basilisk/src/vertexbuffer.h\", 0, \"VertexBuffer.type == type\");\n    }\n    else\n      VertexBuffer.type = type;\n  }\n  else\n    glEnd();\n}\n\n\n// #line 222 \"/home/spencer/basilisk/src/vertexbuffer.h\"\nstatic void vertex_buffer_glVertex2d (real x, real y)\n{\n  if (VertexBuffer.position) {\n    if (VertexBuffer.dim < 2)\n      VertexBuffer.dim = 2;\n    real v[4] = {x, y, 0, 1.};\n    vector_multiply (v, VertexBuffer.modelview);\n    array_append (VertexBuffer.position, v, 3*sizeof(real));\n    VertexBuffer.nvertex++;\n  }\n  else\n    glVertex3d (x, y, 0.);\n}\n\n\n// #line 200 \"/home/spencer/basilisk/src/common.h\"\nstatic void * pmfunc_free (void * ptr, char c)\n{\n  if (!ptr)\n    return ptr;\n  pmdata * d = (pmdata *) (((char *)ptr) - sizeof(pmdata));\n  if (d->id < 1 || d->id > pmfuncn) {\n    fputs (\"*** MTRACE: ERROR!: corrupted free()\", ferr);\n    if (d->size == 0)\n      fputs (\", possible double free()\", ferr);\n    else\n      fputs (\", not traced?\", ferr);\n    fputs (\", aborting...\\n\", ferr);\n    abort();\n    return ptr;\n  }\n  else\n  OMP (omp critical)\n  {\n    pmfunc * f = &pmfuncs[d->id - 1];\n    if (f->total < d->size) {\n      fprintf (ferr, \"*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\\n\",\n        f->total, d->size);\n      abort();\n    }\n    else\n      f->total -= d->size;\n    if (pmtrace.total < d->size) {\n      fprintf (ferr, \"*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\\n\",\n        pmtrace.total, d->size);\n      abort();\n    }\n    else {\n      pmtrace.total -= d->size;\n      pmtrace.overhead -= sizeof(pmdata);\n    }\n    d->id = 0;\n    d->size = 0;\n    pmfunc_trace (f, c);\n  }\n  return d;\n}\n\n\n// #line 256 \"/home/spencer/basilisk/src/common.h\"\nstatic void * prealloc (void * ptr, size_t size,\n   const char * func, const char * file, int line)\n{\n  return pmfunc_alloc ((pmdata *) sysrealloc (pmfunc_free(ptr, '<'),\n           sizeof(pmdata) + size),\n         size, func, file, line, '>');\n}\n\n\n// #line 394 \"/home/spencer/basilisk/src/common.h\"\nvoid array_append (Array * a, void * elem, size_t size)\n{\n  if (a->len + size >= a->max) {\n    a->max += max (size, 4096);\n    a->p = prealloc (a->p, a->max,__func__,__FILE__,0);\n  }\n  memcpy (((char *)a->p) + a->len, elem, size);\n  a->len += size;\n}\n\n\n// #line 380 \"/home/spencer/basilisk/src/common.h\"\nArray * array_new()\n{\n  Array * a = ((Array *) pmalloc ((1)*sizeof(Array),__func__,__FILE__,0));\n  a->p = NULL;\n  a->max = a->len = 0;\n  return a;\n}\n\n\n// #line 1410 \"/home/spencer/basilisk/src/common.h\"\nvoid free_solver_func_add (free_solver_func func)\n{\n  if (!free_solver_funcs)\n    free_solver_funcs = array_new();\n  array_append (free_solver_funcs, &func, sizeof(free_solver_func));\n}\n\n\n// #line 154 \"/home/spencer/basilisk/src/common.h\"\nstatic void pmfunc_trace (pmfunc * f, char c)\n{\n  if (pmtrace.fp)\n    fprintf (pmtrace.fp, \"%c %d %ld %ld %ld\",\n      c, f->id, pmtrace.nr, pmtrace.total, f->total);\n@if 1\n  if (pmtrace.nr % 1 == 0) {\n    struct rusage usage;\n    getrusage (RUSAGE_SELF, &usage);\n    if (pmtrace.fp)\n      fprintf (pmtrace.fp, \" %ld\", usage.ru_maxrss*1024);\n    if (!pmtrace.nr)\n      pmtrace.startrss = usage.ru_maxrss;\n    if (usage.ru_maxrss > pmtrace.maxrss)\n      pmtrace.maxrss = usage.ru_maxrss;\n  }\n@endif\n  if (pmtrace.fp)\n    fputc ('\\n', pmtrace.fp);\n  pmtrace.nr++;\n}\n\n\n// #line 135 \"/home/spencer/basilisk/src/common.h\"\nstatic int pmfunc_index (const char * func, const char * file, int line)\n{\n  pmfunc * p = pmfuncs;\n  for (int i = 0; i < pmfuncn; i++, p++)\n    if (p->line == line && !strcmp(func, p->func) && !strcmp(file, p->file))\n      return p->id;\n  pmfuncn++;\n  pmfuncs = (pmfunc *) sysrealloc (pmfuncs, pmfuncn*sizeof(pmfunc));\n  p = &pmfuncs[pmfuncn - 1];\n  memset (p, 0, sizeof(pmfunc));\n  p->func = systrdup(func);\n  p->file = systrdup(file);\n  p->line = line;\n  p->id = pmfuncn;\n  if (pmtrace.fp)\n    fprintf (pmtrace.fp, \"@ %d %s %s %d\\n\", pmfuncn, func, file, line);\n  return pmfuncn;\n}\n@ def not_mpi_compatible()\ndo {\n  if (npe() > 1) {\n    fprintf (ferr, \"%s() is not compatible with MPI (yet)\\n\", __func__);\n    exit (1);\n  }\n} while(0)\n@\n@ define system(command) (pid() == 0 ? system(command) : 0)\n@else\n@ define qstderr() stderr\n@ define qstdout() stdout\n@ define ferr stderr\n@ define fout stdout\n@ define not_mpi_compatible()\n@endif\n\n\n\n\n// #line 93 \"/home/spencer/basilisk/src/common.h\"\nstatic inline void qassert (const char * file, int line, const char * cond) {\n  fprintf (ferr, \"%s:%d: Assertion `%s' failed.\\n\", file, line, cond);\n  abort();\n}\n\n\n// #line 176 \"/home/spencer/basilisk/src/common.h\"\nstatic void * pmfunc_alloc (pmdata * d, size_t size,\n       const char * func, const char * file, int line,\n       char c)\n{\n  if (!(d != NULL)) qassert (\"/home/spencer/basilisk/src/common.h\", 0, \"d != NULL\");\n  OMP (omp critical)\n  {\n    d->id = pmfunc_index(func, file, line);\n    d->size = size;\n    pmfunc * f = &pmfuncs[d->id - 1];\n    f->total += size;\n    if (f->total > f->max)\n      f->max = f->total;\n    pmtrace.total += size;\n    pmtrace.overhead += sizeof(pmdata);\n    if (pmtrace.total > pmtrace.max) {\n      pmtrace.max = pmtrace.total;\n      pmtrace.maxoverhead = pmtrace.overhead;\n    }\n    pmfunc_trace (f, c);\n  }\n  return ((char *)d) + sizeof(pmdata);\n}\n\n\n// #line 242 \"/home/spencer/basilisk/src/common.h\"\nstatic void * pmalloc (size_t size,\n         const char * func, const char * file, int line)\n{\n  return pmfunc_alloc ((pmdata *) sysmalloc (sizeof(pmdata) + size),\n         size, func, file, line, '+');\n}\n\n\n// #line 249 \"/home/spencer/basilisk/src/common.h\"\nstatic void * pcalloc (size_t nmemb, size_t size,\n         const char * func, const char * file, int line)\n{\n  void * p = pmalloc (nmemb*size, func, file, line);\n  return memset (p, 0, nmemb*size);\n}\n\n\n\n\n\n// #line 171 \"/home/spencer/basilisk/src/view.h\"\nbview * bview_new()\n{\n  bview * p = ((bview *) pcalloc (1, sizeof(bview),__func__,__FILE__,0));\n\n  p->tx = p->ty = 0;\n  p->sx = p->sy = p->sz = 1.;\n  p->quat[0] = p->quat[1] = p->quat[2] = 0; p->quat[3] = 1;\n  p->fov = 24.;\n  gl_trackball (p->quat, 0.0, 0.0, 0.0, 0.0);\n\n\n  p->bg[0] = 1; p->bg[1] = 1; p->bg[2] = 1;\n\n\n\n  p->res = 1.;\n  p->lc = 0.004;\n\n  p->samples = 4;\n  p->width = 600*p->samples, p->height = 600*p->samples;\n\n\n  disable_fpe (FE_DIVBYZERO|FE_INVALID);\n\n  p->fb = framebuffer_new (p->width, p->height);\n\n  init_gl();\n  p->active = false;\n\n  enable_fpe (FE_DIVBYZERO|FE_INVALID);\n\n  return p;\n}\n\n\n// #line 232 \"/home/spencer/basilisk/src/view.h\"\nbview * get_view() {\n  if (!_view) {\n    _view = bview_new();\n    free_solver_func_add (destroy_view);\n  }\n  return _view;\n}\n\n\n// #line 61 \"/home/spencer/basilisk/src/vertexbuffer.h\"\nstatic void vertex_buffer_glBegin (unsigned int state)\n{\n  if (VertexBuffer.index) {\n\n    glGetFloatv (0x0BA6, VertexBuffer.modelview);\n\n    bview * view = get_view();\n\n    real q[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,\n      - view->tx, - view->ty, 3, 1 };\n    matrix_multiply (q, VertexBuffer.modelview);\n    for (int i = 0; i < 16; i++)\n      VertexBuffer.modelview[i] = q[i];\n\n    gl_build_rotmatrix ((real (*)[4])q, view->quat);\n    do { real __tmp = q[1]; q[1] = q[4]; q[4] = __tmp; } while(0);\n    do { real __tmp = q[2]; q[2] = q[8]; q[8] = __tmp; } while(0);\n    do { real __tmp = q[6]; q[6] = q[9]; q[9] = __tmp; } while(0);\n    matrix_multiply (q, VertexBuffer.modelview);\n    for (int i = 0; i < 16; i++)\n      VertexBuffer.modelview[i] = q[i];\n\n    VertexBuffer.state = state;\n    switch (state) {\n    case 0x0002:\n      VertexBuffer.line_loop = VertexBuffer.nvertex;\n      break;\n    case 0x0001:\n      VertexBuffer.lines = VertexBuffer.nvertex;\n      break;\n    case 0x0003:\n      VertexBuffer.line_strip = VertexBuffer.nvertex;\n      break;\n    case 0x0007:\n      VertexBuffer.quads = VertexBuffer.nvertex;\n      break;\n    case 0x0009:\n      VertexBuffer.polygon = VertexBuffer.nvertex;\n      break;\n    case 0x0006:\n      VertexBuffer.fan = VertexBuffer.nvertex;\n      break;\n    default:\n      fprintf (ferr, \"glBegin (%d) not implemented yet\\n\", state);\n      break;\n    }\n  }\n  else\n    glBegin (state);\n}\n// # 48 \"/home/spencer/basilisk/src/gl/font.h\"\nstatic SOG_StrokeFont *oghStrokeByID( void *font )\n{\n\n\n    if( font == ((void *)0x0001) )\n        return &ogStrokeMonoRoman;\n\n    fprintf (ferr, \"stroke font %p not found\", font );\n    return 0;\n}\n// # 147 \"/home/spencer/basilisk/src/gl/font.h\"\nvoid gl_StrokeString( const char *string )\n{\n    void *fontID = ((void *)0x0001);\n    int i, j;\n    real length = 0.0;\n    SOG_StrokeFont *font = oghStrokeByID( fontID );\n    unsigned char c;\n\n    if( font && string )\n\n\n\n\n\n        while(( c = *string++ ))\n       if( c < font->Quantity ) {\n                if( c == '\\n' )\n                {\n                    glTranslatef ( -length, -( real )( font->Height ), 0.0 );\n                    length = 0.0;\n                }\n                else\n                {\n                    const SOG_StrokeChar *schar =\n                        font->Characters[ c ];\n                    if( schar )\n                    {\n                        const SOG_StrokeStrip *strip = schar->Strips;\n\n                        for( i = 0; i < schar->Number; i++, strip++ )\n                        {\n                            vertex_buffer_glBegin( 0x0003 );\n\n                            for( j = 0; j < strip->Number; j++ )\n                                vertex_buffer_glVertex2d( strip->Vertices[ j ].X,\n                                            strip->Vertices[ j ].Y);\n\n                            vertex_buffer_glEnd( );\n                        }\n\n                        length += schar->Right;\n                        glTranslatef( schar->Right, 0.0, 0.0 );\n                    }\n                }\n     }\n}"),_("\n      \n// #line 1695 \"/home/spencer/basilisk/src/draw.h\"\nif (val(ff,0,0,0) != 1e30) {\n glPushMatrix();\n char s[80];\n sprintf (s, \"%g\", val(ff,0,0,0));\n real scale = 0.8*Delta_x/(strlen(s)*width);\n glTranslatef (x - 0.4*Delta_x, y - scale*height/3., 0.);\n glScalef (scale, scale, 1.);\n gl_StrokeString (s);\n glPopMatrix();\n      }")})
      {_stencil_val(ff,0,0,0); { 
 
  
_stencil_val(ff,0,0,0);     
 
            
 
 
 
 
      
#line 1704
}   }end_foreach_visible_stencil();
    {
#line 1694
foreach_visible (view)
      if (val(ff,0,0,0) != 1e30) {
 glPushMatrix();
 char s[80];
 sprintf (s, "%g", val(ff,0,0,0));
 float scale = 0.8*Delta_x/(strlen(s)*width);
 glTranslatef (x - 0.4*Delta_x, y - scale*height/3., 0.);
 glScalef (scale, scale, 1.);
 gl_StrokeString (s);
 glPopMatrix();
      }end_foreach_visible();}
  }end_draw_lines();}
  view->res = res;
  if (expr) delete (((scalar[]){ff,{-1}}));
  {end_tracing("labels","/home/spencer/basilisk/src/draw.h",0);return true;}




end_tracing("labels","/home/spencer/basilisk/src/draw.h",0);}







# 1 "./draw_json.h" 1
# 1 "/home/spencer/basilisk/src/draw_json.h"

int _view_json (char * s, int len) {
  int i, len1 = 0;
  i = snprintf (s, len, "  \"view\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"tx\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"ty\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"fov\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"quat\": { \"type\": \"pfloat\", \"cardinality\": 4, \"value\": [%f,%f,%f,%f] }", 0., 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"sx\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"sy\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"sz\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"width\": { \"type\": \"punsigned\", \"cardinality\": 1, \"value\": \"%d\" }", 0);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"height\": { \"type\": \"punsigned\", \"cardinality\": 1, \"value\": \"%d\" }", 0);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"samples\": { \"type\": \"punsigned\", \"cardinality\": 1, \"value\": \"%d\" }", 0);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"bg\": { \"type\": \"pfloat\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"theta\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"phi\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"psi\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"relative\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", 0);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"tz\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"near\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"far\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"res\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"camera\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", "");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"cache\": { \"type\": \"pint\", \"cardinality\": 1, \"value\": \"%d\" }", 0);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"p1x\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"p1y\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"p2x\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"p2y\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}
int _begin_translate_json (char * s, int len) {
  int i, len1 = 0;
  i = snprintf (s, len, "  \"begin_translate\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"x\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"y\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"z\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}
int _begin_mirror_json (char * s, int len) {
  int i, len1 = 0;
  i = snprintf (s, len, "  \"begin_mirror\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"n\": { \"type\": \"pdouble\", \"cardinality\": 3, \"value\": [%lf,%lf,%lf] }", 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"alpha\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}
int _draw_vof_json (char * s, int len) {
  int i, len1 = 0;
  i = snprintf (s, len, "  \"draw_vof\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"c\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", "");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"s\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", "");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"edges\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", 0);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"larger\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"filled\": { \"type\": \"pint\", \"cardinality\": 1, \"value\": \"%d\" }", 0);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"color\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", "");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"min\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"max\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"spread\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"linear\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", 0);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"fc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}
int _isoline_json (char * s, int len) {
  int i, len1 = 0;
  i = snprintf (s, len, "  \"isoline\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"phi\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", "");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"val\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"n\": { \"type\": \"pint\", \"cardinality\": 1, \"value\": \"%d\" }", 0);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"edges\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", 0);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"larger\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"filled\": { \"type\": \"pint\", \"cardinality\": 1, \"value\": \"%d\" }", 0);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"color\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", "");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"min\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"max\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"spread\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"linear\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", 0);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"fc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}
int _cells_json (char * s, int len) {
  int i, len1 = 0;
  i = snprintf (s, len, "  \"cells\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"n\": { \"type\": \"pdouble\", \"cardinality\": 3, \"value\": [%lf,%lf,%lf] }", 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"alpha\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}
int _vectors_json (char * s, int len) {
  int i, len1 = 0;
  i = snprintf (s, len, "  \"vectors\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"u\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", "");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"scale\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}
int _squares_json (char * s, int len) {
  int i, len1 = 0;
  i = snprintf (s, len, "  \"squares\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"color\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", "");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"z\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", "");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"min\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"max\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"spread\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"linear\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", 0);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"fc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"n\": { \"type\": \"pdouble\", \"cardinality\": 3, \"value\": [%lf,%lf,%lf] }", 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"alpha\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}
int _box_json (char * s, int len) {
  int i, len1 = 0;
  i = snprintf (s, len, "  \"box\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"notics\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", 0);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}
int _isosurface_json (char * s, int len) {
  int i, len1 = 0;
  i = snprintf (s, len, "  \"isosurface\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"f\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", "");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"v\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"color\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", "");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"min\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"max\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"spread\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"linear\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", 0);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"fc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}
int _travelling_json (char * s, int len) {
  int i, len1 = 0;
  i = snprintf (s, len, "  \"travelling\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"start\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"end\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"tx\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"ty\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"quat\": { \"type\": \"pfloat\", \"cardinality\": 4, \"value\": [%f,%f,%f,%f] }", 0., 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"fov\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}
int _draw_string_json (char * s, int len) {
  int i, len1 = 0;
  i = snprintf (s, len, "  \"draw_string\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"str\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", "");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"pos\": { \"type\": \"pint\", \"cardinality\": 1, \"value\": \"%d\" }", 0);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"size\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}
int _labels_json (char * s, int len) {
  int i, len1 = 0;
  i = snprintf (s, len, "  \"labels\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"f\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", "");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}
# 1722 "/home/spencer/basilisk/src/draw.h" 2

struct {
  int (* json) (char * s, int len);
} bview_interface[] = {
  { _draw_vof_json },
  { _squares_json },
  { _cells_json },
  { _box_json },

  { _isoline_json },
  { _labels_json },
  { _vectors_json },



  { NULL }
};
# 422 "/home/spencer/basilisk/src/view.h" 2
# 437 "/home/spencer/basilisk/src/view.h"
bool load (FILE * fp, char * file, Array * buf);

static void bview_draw (bview * view)
{
  if (!view->active)
    return;
  view->active = false;
  glFinish ();
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}
# 494 "/home/spencer/basilisk/src/view.h"
     
bool save (char * file, char * format, char * opt,
    FILE * fp,
    float lw,
    int sort, int options,

    bview * view)
{tracing("save","/home/spencer/basilisk/src/view.h",0);
  if (file) {
    char * s = strchr (file, '.'), * dot = s;
    while (s) {
      dot = s;
      s = strchr (s + 1, '.');
    }
    if (dot)
      format = dot + 1;
  }

  if (!view)
    view = get_view();

  if ((!strcmp (format, "png") && which ("convert")) ||
      !strcmp (format, "jpg") ||
      (file && is_animation (file))) {
    bview_draw (view);
    unsigned char * image = (unsigned char *) compose_image (view);
    if (pid() == 0) {
      FILE * fp = open_image (file, opt);
      if (!fp) {
 perror (file);
 {end_tracing("save","/home/spencer/basilisk/src/view.h",0);return false;}
      }
      gl_write_image (fp, image, view->width, view->height, view->samples);
      close_image (file, fp);
    }
    {end_tracing("save","/home/spencer/basilisk/src/view.h",0);return true;}
  }

  if (file && (fp = fopen (file, "w")) == NULL) {
    perror (file);
    {end_tracing("save","/home/spencer/basilisk/src/view.h",0);return false;}
  }
  if (!fp)
    fp = fout;

  if (!strcmp (format, "ppm")) {
    bview_draw (view);
    unsigned char * image = (unsigned char *) compose_image (view);
    if (pid() == 0)
      gl_write_image (fp, image, view->width, view->height, view->samples);
  }

  else if (!strcmp (format, "png")) {
    bview_draw (view);
    unsigned char * image = (unsigned char *) compose_image (view);
    if (pid() == 0)
      gl_write_image_png (fp, image, view->width, view->height, view->samples);
  }

  else if (!strcmp (format, "bv")) {

    fprintf (ferr, "save(): error: the '%s' format is no longer supported\n",
      format);
# 573 "/home/spencer/basilisk/src/view.h"
  }

  else if (!strcmp (format, "gnu") ||
    !strcmp (format, "obj") ||
    !strcmp (format, "kml") ||
    !strcmp (format, "ps") ||
    !strcmp (format, "eps") ||
    !strcmp (format, "tex") ||
    !strcmp (format, "pdf") ||
    !strcmp (format, "svg") ||
    !strcmp (format, "pgf"))
    fprintf (ferr, "save(): error: the '%s' format is no longer supported\n",
      format);

  else {
    fprintf (ferr, "save(): unknown format '%s'\n", format);
    if (file) {
      fclose (fp);
      remove (file);
    }
    {end_tracing("save","/home/spencer/basilisk/src/view.h",0);return false;}
  }

  fflush (fp);
  if (file)
    fclose (fp);

  {end_tracing("save","/home/spencer/basilisk/src/view.h",0);return true;}
end_tracing("save","/home/spencer/basilisk/src/view.h",0);}







static char * remove_blanks (char * line)
{
  while (strchr (" \t", *line)) line++;
  char * s = line, * cur = line;
  bool instring = false;
  while (*s != '\0' && *s != '#') {
    if (*s == '"')
      instring = !instring;
    if (instring || !strchr (" \t", *s))
      *cur++ = *s;
    s++;
  }
  *cur = '\0';
  return line;
}

# 1 "./parse.h" 1
# 1 "/home/spencer/basilisk/src/parse.h"




enum ParamsType { pstring, pint, punsigned, pbool, pfloat, pdouble, pcolormap };

typedef struct {
  char * key;
  enum ParamsType type;
  void * val;
  int n;
} Params;

static bool atobool (char * s)
{
  if (!strcmp (s, "true"))
    return true;
  if (!strcmp (s, "false"))
    return false;
  return atoi (s) != 0;
}

static bool args (Params * p, char * val)
{
  static char * name[] = { "string", "int", "unsigned",
      "bool", "float", "double", "colormap" };
  switch (p->type) {

  case pstring:
    if (val[0] != '"') {
      fprintf (ferr, "expecting a string for '%s' got '%s'\n", p->key, val);
      return false;
    }
    if (val[strlen(val) - 1] != '"') {
      fprintf (ferr, "unterminated quoted string '%s'\n", val);
      return false;
    }
    val[strlen(val) - 1] = '\0';
    char * s = &val[1];
    int nc = 0;
    while (*s != '\0') {
      if (!strchr (" \t\n\r", *s))
 nc++;
      s++;
    }
    *((char **)p->val) = nc > 0 ? &val[1] : NULL;
    break;

  case pcolormap:
    if (!strcmp (val, "jet"))
      *((Colormap *)p->val) = jet;
    else if (!strcmp (val, "cool_warm"))
      *((Colormap *)p->val) = cool_warm;
    else if (!strcmp (val, "gray"))
      *((Colormap *)p->val) = gray;
    else if (!strcmp (val, "randomap"))
      *((Colormap *)p->val) = randomap;
    else {
      fprintf (ferr, "unknown colormap '%s'\n", val);
      return false;
    }
    break;

  case pint: case punsigned: case pbool: case pdouble: case pfloat:
    if (val[0] == '"') {
      fprintf (ferr, "expecting a %s for '%s' got %s\n",
        name[p->type], p->key, val);
      return false;
    }
    if (!p->n) {
      switch (p->type) {
      case pint: *((int *)p->val) = atoi(val); break;
      case punsigned: *((unsigned *)p->val) = atoi(val); break;
      case pbool: *((bool *)p->val) = atobool(val); break;
      case pfloat: *((float *)p->val) = atof(val); break;
      case pdouble: *((double *)p->val) = atof(val); break;
      default: if (!(false)) qassert ("/home/spencer/basilisk/src/parse.h", 0, "false");
      }
    }
    else {
      if (val[0] != '{') {
 fprintf (ferr, "expecting an array for '%s' got %s\n", p->key, val);
 return false;
      }
      val++;
      int i = 0;
      char c = ',';
      while (i < p->n && c != '}') {
 char * s = strchr (val, ',');
 if (!s)
   s = strchr (val, '}');
 if (!s) {
   fprintf (ferr, "expecting an array for '%s' got %s\n", p->key, val);
   return false;
 }
 c = *s;
 *s++ = '\0';
 switch (p->type) {
 case pint: ((int *)p->val)[i++] = atoi (val); break;
 case punsigned: ((unsigned *)p->val)[i++] = atoi (val); break;
 case pbool: ((bool *)p->val)[i++] = atobool (val); break;
 case pfloat: ((float *)p->val)[i++] = atof (val); break;
 case pdouble: ((double *)p->val)[i++] = atof (val); break;
 default: if (!(false)) qassert ("/home/spencer/basilisk/src/parse.h", 0, "false");
 }
 val = s;
      }
      if (c != '}') {
 fprintf (ferr, "expecting '}' for '%s' got %s\n", p->key, val);
 return false;
      }
    }
    break;

  default:
    if (!(false)) qassert ("/home/spencer/basilisk/src/parse.h", 0, "false");
  }
  return true;
}

static char * find_comma (char * s)
{
  int par = 0;
  while (*s != '\0') {
    if (*s == ',' && par == 0) {
      *s = '\0';
      return s + 1;
    }
    if (*s == '{')
      par++;
    else if (*s == '}')
      par--;
    s++;
  }
  return NULL;
}

static char * mystrtok (char * str, const char * delim)
{
  static char * s = NULL;
  char * start = str ? str : s;
  bool string = false;
  s = start;
  while (*s != '\0') {
    if (*s == '"')
      string = !string;
    if (!string && strchr(delim, *s))
      break;
    s++;
  }
  if (*s != '\0')
    *s++ = '\0';
  return start;
}

int parse_params (Params * params)
{
  char * s;
  int i = 0, n = 0;
  Params * p = params;
  while (p->key) p++, n++;
  if (!(s = mystrtok (NULL, ");")) || s[0] == '\n')
    return false;
  while (s && *s != '\0') {
    char * next = find_comma (s), * key = s;
    if ((s = strchr (key, '='))) {
      s[0] = '\0', s++;
      i = -1;
      Params * p = params;
      while (p->key && strcmp(p->key, key)) p++;
      if (!p->key) {
 fprintf (ferr, "unknown key '%s'\n", key);
 return false;
      }
      if (!args (p, s))
 return false;
    }
    else {
      if (i < 0) {
 fprintf (ferr, "anonymous value '%s' after keys\n", key);
 return false;
      }
      if (i >= n) {
 fprintf (ferr, "too many parameters: '%s' %d %d\n", key, i, n);
 return false;
      }
      if (!args (&params[i], key))
 return false;
      i++;
    }
    s = next;
  }
  return true;
}
# 626 "/home/spencer/basilisk/src/view.h" 2

bool process_line (char * line)
{
  if (line[0] == '\0')
    return true;
  char * s = mystrtok (remove_blanks (line), "(");
  if (!s)
    return true;

  if (!strcmp (s, "restore")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file) {
      bview * view = get_view();
      if (view->cache) {
 free_cexpr (view->cache);
 view->cache = pcalloc (1, sizeof (cexpr),__func__,__FILE__,0);
      }
      if (!restore ( file, all
#line 1117 "/home/spencer/basilisk/src/output.h"
, 
NULL
#line 644 "/home/spencer/basilisk/src/view.h"
))
 fprintf (ferr, "could not restore from '%s'\n", file);
      else {
 restriction (all);
 fields_stats();
 clear();
      }
    }
  }

  else if (!strcmp (s, "dump")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    dump ( file
#line 1041 "/home/spencer/basilisk/src/output.h"
, 
all, 
NULL, 
false
#line 657 "/home/spencer/basilisk/src/view.h"
);
  }

  else if (!strcmp (s, "input_gfs")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file) {
      input_gfs ( 
#line 157 "/home/spencer/basilisk/src/input.h"
stdin
#line 664 "/home/spencer/basilisk/src/view.h"
, all, file);
      restriction (all);
      fields_stats();
      clear();
    }
  }

  else if (!strcmp (s, "save")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file)
      save ( file
#line 495
, "ppm", NULL, 
NULL, 
0, 
0, 0, 

NULL
#line 675
);
  }

  else if (!strcmp (s, "load")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file)
      load ( 
#line 437
NULL
#line 682
, file
#line 437
, NULL
#line 682
);
  }






# 1 "./draw_get.h" 1
# 1 "/home/spencer/basilisk/src/draw_get.h"

else if (!strcmp (s, "view")) {
  float tx = 0.;
  float ty = 0.;
  float fov = 0.;
  float quat[4] = {0};
  float sx = 1.;
  float sy = 1.;
  float sz = 1.;
  unsigned width = 800;
  unsigned height = 800;
  unsigned samples = 4;
  float bg[3] = {0};
  float theta = 0.;
  float phi = 0.;
  float psi = 0.;
  bool relative = false;
  float tz = 0.;
  float near = 0.;
  float far = 0.;
  float res = 0.;
  char * camera = NULL;
  MapFunc map = NULL;
  int cache = 0;
  float p1x = 0.;
  float p1y = 0.;
  float p2x = 0.;
  float p2y = 0.;
  bview * view1 = NULL;
  Params params[] = {
    {"tx", pfloat, &tx},
    {"ty", pfloat, &ty},
    {"fov", pfloat, &fov},
    {"quat", pfloat, quat, 4},
    {"sx", pfloat, &sx},
    {"sy", pfloat, &sy},
    {"sz", pfloat, &sz},
    {"width", punsigned, &width},
    {"height", punsigned, &height},
    {"samples", punsigned, &samples},
    {"bg", pfloat, bg, 3},
    {"theta", pfloat, &theta},
    {"phi", pfloat, &phi},
    {"psi", pfloat, &psi},
    {"relative", pbool, &relative},
    {"tz", pfloat, &tz},
    {"near", pfloat, &near},
    {"far", pfloat, &far},
    {"res", pfloat, &res},
    {"camera", pstring, &camera},
    {"cache", pint, &cache},
    {"p1x", pfloat, &p1x},
    {"p1y", pfloat, &p1y},
    {"p2x", pfloat, &p2x},
    {"p2y", pfloat, &p2y},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  view (tx,ty,fov,quat,sx,sy,sz,width,height,samples,bg,theta,phi,psi,relative,tz,near,far,res,camera,map,cache,p1x,p1y,p2x,p2y,view1);
}
else if (!strcmp (s, "begin_translate")) {
  float x = 0;
  float y = 0.;
  float z = 0.;
  Params params[] = {
    {"x", pfloat, &x},
    {"y", pfloat, &y},
    {"z", pfloat, &z},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  begin_translate (x,y,z);
}
else if (!strcmp (s, "begin_mirror")) {
  coord n = {0};
  double alpha = 0.;
  Params params[] = {
    {"n", pdouble, &n, 3},
    {"alpha", pdouble, &alpha},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  begin_mirror (n,alpha);
}
else if (!strcmp (s, "draw_vof")) {
  char * c = NULL;
  char * s = NULL;
  bool edges = false;
  double larger = 0.;
  int filled = 0;
  char * color = NULL;
  double min = 0;
  double max = 0;
  double spread = 0;
  bool linear = false;
  Colormap map = jet;
  float fc[3] = {0};
  float lc[3] = {0};
  float lw = 1.;
  bool expr = false;
  Params params[] = {
    {"c", pstring, &c},
    {"s", pstring, &s},
    {"edges", pbool, &edges},
    {"larger", pdouble, &larger},
    {"filled", pint, &filled},
    {"color", pstring, &color},
    {"min", pdouble, &min},
    {"max", pdouble, &max},
    {"spread", pdouble, &spread},
    {"linear", pbool, &linear},
    {"fc", pfloat, fc, 3},
    {"lc", pfloat, lc, 3},
    {"lw", pfloat, &lw},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  if (!draw_vof (c,s,edges,larger,filled,color,min,max,spread,linear,map,fc,lc,lw,expr
#line 893 "/home/spencer/basilisk/src/draw.h"
, 
false, 15,( float[2]) {-.95, -.95}, "", 1, false, false, false, 50, "%g", 50
#line 122 "/home/spencer/basilisk/src/draw_get.h"
))
    return false;
}
else if (!strcmp (s, "isoline")) {
  char * phi = NULL;
  double val = 0.;
  int n = 1;
  bool edges = false;
  double larger = 0.;
  int filled = 0;
  char * color = NULL;
  double min = 0;
  double max = 0;
  double spread = 0;
  bool linear = false;
  Colormap map = jet;
  float fc[3] = {0};
  float lc[3] = {0};
  float lw = 1.;
  bool expr = false;
  Params params[] = {
    {"phi", pstring, &phi},
    {"val", pdouble, &val},
    {"n", pint, &n},
    {"edges", pbool, &edges},
    {"larger", pdouble, &larger},
    {"filled", pint, &filled},
    {"color", pstring, &color},
    {"min", pdouble, &min},
    {"max", pdouble, &max},
    {"spread", pdouble, &spread},
    {"linear", pbool, &linear},
    {"fc", pfloat, fc, 3},
    {"lc", pfloat, lc, 3},
    {"lw", pfloat, &lw},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  if (!isoline (phi,val,n,edges,larger,filled,color,min,max,spread,linear,map,fc,lc,lw,expr
#line 1074 "/home/spencer/basilisk/src/draw.h"
, 
false, 15,( float[2]) {-.95, -.95}, "", 1, false, false, false, 50, "%g", 50
#line 161 "/home/spencer/basilisk/src/draw_get.h"
))
    return false;
}
else if (!strcmp (s, "cells")) {
  coord n = {0,0,1};
  double alpha = 0.;
  float lc[3] = {0};
  float lw = 1.;
  Params params[] = {
    {"n", pdouble, &n, 3},
    {"alpha", pdouble, &alpha},
    {"lc", pfloat, lc, 3},
    {"lw", pfloat, &lw},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  if (!cells (n,alpha,lc,lw))
    return false;
}
else if (!strcmp (s, "vectors")) {
  char * u = NULL;
  double scale = 1;
  float lc[3] = {0};
  float lw = 1.;
  Params params[] = {
    {"u", pstring, &u},
    {"scale", pdouble, &scale},
    {"lc", pfloat, lc, 3},
    {"lw", pfloat, &lw},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  if (!vectors (u,scale,lc,lw))
    return false;
}
else if (!strcmp (s, "squares")) {
  char * color = NULL;
  char * z = NULL;
  double min = 0;
  double max = 0;
  double spread = 0;
  bool linear = false;
  Colormap map = jet;
  float fc[3] = {0};
  float lc[3] = {0};
  bool expr = false;
  coord n = {0,0,1};
  double alpha = 0;
  Params params[] = {
    {"color", pstring, &color},
    {"z", pstring, &z},
    {"min", pdouble, &min},
    {"max", pdouble, &max},
    {"spread", pdouble, &spread},
    {"linear", pbool, &linear},
    {"fc", pfloat, fc, 3},
    {"lc", pfloat, lc, 3},
    {"n", pdouble, &n, 3},
    {"alpha", pdouble, &alpha},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  if (!squares (color,z,min,max,spread,linear,map,fc,lc,expr,n,alpha
#line 1230 "/home/spencer/basilisk/src/draw.h"
, 
1, 
false, 15,( float[2]) {-.95, -.95}, "", 1, false, false, false, 50, "%g", 50
#line 226 "/home/spencer/basilisk/src/draw_get.h"
))
    return false;
}
else if (!strcmp (s, "box")) {
  bool notics = false;
  float lc[3] = {0};
  float lw = 1.;
  Params params[] = {
    {"notics", pbool, &notics},
    {"lc", pfloat, lc, 3},
    {"lw", pfloat, &lw},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  if (!box (notics,lc,lw))
    return false;
}
else if (!strcmp (s, "isosurface")) {
  char * f = NULL;
  double v = 0.;
  char * color = NULL;
  double min = 0;
  double max = 0;
  double spread = 0;
  bool linear = false;
  Colormap map = jet;
  float fc[3] = {0};
  float lc[3] = {0};
  float lw = 1;
  bool expr = false;
  Params params[] = {
    {"f", pstring, &f},
    {"v", pdouble, &v},
    {"color", pstring, &color},
    {"min", pdouble, &min},
    {"max", pdouble, &max},
    {"spread", pdouble, &spread},
    {"linear", pbool, &linear},
    {"fc", pfloat, fc, 3},
    {"lc", pfloat, lc, 3},
    {"lw", pfloat, &lw},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  if (!isosurface (f,v,color,min,max,spread,linear,map,fc,lc,lw,expr
#line 1516 "/home/spencer/basilisk/src/draw.h"
, 
false, 15,( float[2]) {-.95, -.95}, "", 1, false, false, false, 50, "%g", 50
#line 272 "/home/spencer/basilisk/src/draw_get.h"
))
    return false;
}
else if (!strcmp (s, "travelling")) {
  double start = 0;
  double end = 0;
  float tx = 0;
  float ty = 0;
  float quat[4] = {0};
  float fov = 0;
  Params params[] = {
    {"start", pdouble, &start},
    {"end", pdouble, &end},
    {"tx", pfloat, &tx},
    {"ty", pfloat, &ty},
    {"quat", pfloat, quat, 4},
    {"fov", pfloat, &fov},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  travelling (start,end,tx,ty,quat,fov);
}
else if (!strcmp (s, "draw_string")) {
  char * str = NULL;
  int pos = 0;
  float size = 40;
  float lc[3] = {0};
  float lw = 1;
  Params params[] = {
    {"str", pstring, &str},
    {"pos", pint, &pos},
    {"size", pfloat, &size},
    {"lc", pfloat, lc, 3},
    {"lw", pfloat, &lw},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  if (!draw_string (str,pos,size,lc,lw))
    return false;
}
else if (!strcmp (s, "labels")) {
  char * f = NULL;
  float lc[3] = {0};
  float lw = 1;
  Params params[] = {
    {"f", pstring, &f},
    {"lc", pfloat, lc, 3},
    {"lw", pfloat, &lw},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  if (!labels (f,lc,lw))
    return false;
}
# 691 "/home/spencer/basilisk/src/view.h" 2

  else if (!strcmp (s, "end_mirror"))
    end_mirror();

  else if (!strcmp (s, "end_translate"))
    end_translate();

  else if (!strcmp (s, "clear"))
    clear();

  else if (s[0] != '\n' && s[0] != '\0')
    fprintf (ferr, "load(): syntax error: '%s'\n", s);

  return true;
}

bool load (FILE * fp, char * file, Array * buf)
{
  if (file) {
    fp = fopen (file, "r");
    if (!fp) {
      perror (file);
      return false;
    }
  }

  if (fp) {
    char line[256];
    while (fgets (line, 256, fp) && process_line (line));
  }
  else if (buf) {
    int i = 0;
    char * s = (char *) buf->p;
    while (i < buf->len) {
      char * start = s;
      while (i < buf->len && *s != '\n')
 s++, i++;
      if (*s == '\n' && ++s > start) {
 char line[s - start + 1];
 strncpy (line, start, s - start);
 line[s - start] = '\0';
 process_line (line);
      }
    }
  }
  return true;
}
# 91 "/home/spencer/basilisk/src/display.h" 2
# 1 "./khash.h" 1
# 1 "/home/spencer/basilisk/src/khash.h"
# 128 "/home/spencer/basilisk/src/khash.h"
# 1 "/home/spencer/basilisk/src/ast/std/stdlib.h" 1
@include <stdlib.h>
# 129 "/home/spencer/basilisk/src/khash.h" 2
# 1 "/home/spencer/basilisk/src/ast/std/string.h" 1
@include <string.h>
# 130 "/home/spencer/basilisk/src/khash.h" 2
# 1 "/home/spencer/basilisk/src/ast/std/stdint.h" 1
@include <stdint.h>
# 131 "/home/spencer/basilisk/src/khash.h" 2

typedef uint32_t khint32_t;
typedef uint64_t khint64_t;
# 151 "/home/spencer/basilisk/src/khash.h"
typedef khint32_t khint_t;
typedef khint_t khiter_t;
# 181 "/home/spencer/basilisk/src/khash.h"
static const double __ac_HASH_UPPER = 0.77;
# 384 "/home/spencer/basilisk/src/khash.h"
static inline khint_t __ac_X31_hash_string(const char *s)
{
 khint_t h = (khint_t)*s;
 if (h) for (++s ; *s; ++s) h = (h << 5) - h + (khint_t)*s;
 return h;
}
# 401 "/home/spencer/basilisk/src/khash.h"
static inline khint_t __ac_Wang_hash(khint_t key)
{
    key += ~(key << 15);
    key ^= (key >> 10);
    key += (key << 3);
    key ^= (key >> 6);
    key += ~(key << 11);
    key ^= (key >> 16);
    return key;
}
# 600 "/home/spencer/basilisk/src/khash.h"
typedef const char *kh_cstr_t;
# 92 "/home/spencer/basilisk/src/display.h" 2

typedef struct {
  int fd;
  int iter;
} DisplayClient;

typedef struct kh_strhash_s { khint_t n_buckets, size, n_occupied, upper_bound; khint32_t *flags; kh_cstr_t *keys; DisplayClient * *vals; } kh_strhash_t; static inline __attribute__ ((__unused__)) kh_strhash_t *kh_init_strhash(void) { return (kh_strhash_t*)pcalloc(1,sizeof(kh_strhash_t),__func__,__FILE__,0); } static inline __attribute__ ((__unused__)) void kh_destroy_strhash(kh_strhash_t *h) { if (h) { pfree((void *)h->keys,__func__,__FILE__,0); pfree(h->flags,__func__,__FILE__,0); pfree((void *)h->vals,__func__,__FILE__,0); pfree(h,__func__,__FILE__,0); } } static inline __attribute__ ((__unused__)) void kh_clear_strhash(kh_strhash_t *h) { if (h && h->flags) { memset(h->flags, 0xaa, ((h->n_buckets) < 16? 1 : (h->n_buckets)>>4) * sizeof(khint32_t)); h->size = h->n_occupied = 0; } } static inline __attribute__ ((__unused__)) khint_t kh_get_strhash(const kh_strhash_t *h, kh_cstr_t key) { if (h->n_buckets) { khint_t k, i, last, mask, step = 0; mask = h->n_buckets - 1; k = __ac_X31_hash_string(key); i = k & mask; last = i; while (!((h->flags[i>>4]>>((i&0xfU)<<1))&2) && (((h->flags[i>>4]>>((i&0xfU)<<1))&1) || !(strcmp(h->keys[i], key) == 0))) { i = (i + (++step)) & mask; if (i == last) return h->n_buckets; } return ((h->flags[i>>4]>>((i&0xfU)<<1))&3)? h->n_buckets : i; } else return 0; } static inline __attribute__ ((__unused__)) int kh_resize_strhash(kh_strhash_t *h, khint_t new_n_buckets) { khint32_t *new_flags = 0; khint_t j = 1; { (--(new_n_buckets), (new_n_buckets)|=(new_n_buckets)>>1, (new_n_buckets)|=(new_n_buckets)>>2, (new_n_buckets)|=(new_n_buckets)>>4, (new_n_buckets)|=(new_n_buckets)>>8, (new_n_buckets)|=(new_n_buckets)>>16, ++(new_n_buckets)); if (new_n_buckets < 4) new_n_buckets = 4; if (h->size >= (khint_t)(new_n_buckets * __ac_HASH_UPPER + 0.5)) j = 0; else { new_flags = (khint32_t*)pmalloc(((new_n_buckets) < 16? 1 : (new_n_buckets)>>4) * sizeof(khint32_t),__func__,__FILE__,0); if (!new_flags) return -1; memset(new_flags, 0xaa, ((new_n_buckets) < 16? 1 : (new_n_buckets)>>4) * sizeof(khint32_t)); if (h->n_buckets < new_n_buckets) { kh_cstr_t *new_keys = (kh_cstr_t*)prealloc((void *)h->keys,new_n_buckets * sizeof(kh_cstr_t),__func__,__FILE__,0); if (!new_keys) { pfree(new_flags,__func__,__FILE__,0); return -1; } h->keys = new_keys; if (1) { DisplayClient * *new_vals = (DisplayClient **)prealloc((void *)h->vals,new_n_buckets * sizeof(DisplayClient *),__func__,__FILE__,0); if (!new_vals) { pfree(new_flags,__func__,__FILE__,0); return -1; } h->vals = new_vals; } } } } if (j) { for (j = 0; j != h->n_buckets; ++j) { if (((h->flags[j>>4]>>((j&0xfU)<<1))&3) == 0) { kh_cstr_t key = h->keys[j]; DisplayClient * val; khint_t new_mask; new_mask = new_n_buckets - 1; if (1) val = h->vals[j]; (h->flags[j>>4]|=1ul<<((j&0xfU)<<1)); while (1) { khint_t k, i, step = 0; k = __ac_X31_hash_string(key); i = k & new_mask; while (!((new_flags[i>>4]>>((i&0xfU)<<1))&2)) i = (i + (++step)) & new_mask; (new_flags[i>>4]&=~(2ul<<((i&0xfU)<<1))); if (i < h->n_buckets && ((h->flags[i>>4]>>((i&0xfU)<<1))&3) == 0) { { kh_cstr_t tmp = h->keys[i]; h->keys[i] = key; key = tmp; } if (1) { DisplayClient * tmp = h->vals[i]; h->vals[i] = val; val = tmp; } (h->flags[i>>4]|=1ul<<((i&0xfU)<<1)); } else { h->keys[i] = key; if (1) h->vals[i] = val; break; } } } } if (h->n_buckets > new_n_buckets) { h->keys = (kh_cstr_t*)prealloc((void *)h->keys,new_n_buckets * sizeof(kh_cstr_t),__func__,__FILE__,0); if (1) h->vals = (DisplayClient **)prealloc((void *)h->vals,new_n_buckets * sizeof(DisplayClient *),__func__,__FILE__,0); } pfree(h->flags,__func__,__FILE__,0); h->flags = new_flags; h->n_buckets = new_n_buckets; h->n_occupied = h->size; h->upper_bound = (khint_t)(h->n_buckets * __ac_HASH_UPPER + 0.5); } return 0; } static inline __attribute__ ((__unused__)) khint_t kh_put_strhash(kh_strhash_t *h, kh_cstr_t key, int *ret) { khint_t x; if (h->n_occupied >= h->upper_bound) { if (h->n_buckets > (h->size<<1)) { if (kh_resize_strhash(h, h->n_buckets - 1) < 0) { *ret = -1; return h->n_buckets; } } else if (kh_resize_strhash(h, h->n_buckets + 1) < 0) { *ret = -1; return h->n_buckets; } } { khint_t k, i, site, last, mask = h->n_buckets - 1, step = 0; x = site = h->n_buckets; k = __ac_X31_hash_string(key); i = k & mask; if (((h->flags[i>>4]>>((i&0xfU)<<1))&2)) x = i; else { last = i; while (!((h->flags[i>>4]>>((i&0xfU)<<1))&2) && (((h->flags[i>>4]>>((i&0xfU)<<1))&1) || !(strcmp(h->keys[i], key) == 0))) { if (((h->flags[i>>4]>>((i&0xfU)<<1))&1)) site = i; i = (i + (++step)) & mask; if (i == last) { x = site; break; } } if (x == h->n_buckets) { if (((h->flags[i>>4]>>((i&0xfU)<<1))&2) && site != h->n_buckets) x = site; else x = i; } } } if (((h->flags[x>>4]>>((x&0xfU)<<1))&2)) { h->keys[x] = key; (h->flags[x>>4]&=~(3ul<<((x&0xfU)<<1))); ++h->size; ++h->n_occupied; *ret = 1; } else if (((h->flags[x>>4]>>((x&0xfU)<<1))&1)) { h->keys[x] = key; (h->flags[x>>4]&=~(3ul<<((x&0xfU)<<1))); ++h->size; *ret = 2; } else *ret = 0; return x; } static inline __attribute__ ((__unused__)) void kh_del_strhash(kh_strhash_t *h, khint_t x) { if (x != h->n_buckets && !((h->flags[x>>4]>>((x&0xfU)<<1))&3)) { (h->flags[x>>4]|=1ul<<((x&0xfU)<<1)); --h->size; } }

static struct {
  kh_strhash_t * objects;
  int sock, port;
  char * error;
  Array * controls;
} Display = { .sock = -1 };

static void display_display()
{
  ;
  for (khiter_t k = (khint_t)(0); k != ((Display.objects)->n_buckets);
       ++k)
    if ((!(((Display.objects)->flags[(k)>>4]>>(((k)&0xfU)<<1))&3))) {
      ;
      DisplayClient * client = ((Display.objects)->vals[k]);
      while (client->fd >= 0) {
 ;
 client++;
      }
      ;
    }
  ;
}

static char * read_file_into_buffer (FILE * fp)
{
  if (fseek (fp, 0L, SEEK_END) < 0)
    return NULL;
  long bufsize = ftell (fp);
  if (bufsize <= 0)
    return NULL;
  char * buf = pmalloc (sizeof(char)*(bufsize + 1),__func__,__FILE__,0);
  if (fseek (fp, 0L, SEEK_SET) < 0) {
    pfree (buf,__func__,__FILE__,0);
    return NULL;
  }
  size_t newLen = fread (buf, sizeof(char), bufsize, fp);
  buf[newLen] = '\0';
  return buf;
}

static void display_command (const char * command)
{
  ;

  vertex_buffer_setup();
  VertexBuffer.vertex = 0;


  int bak = -1;
  if (pid() == 0) {
    fflush (ferr);
    bak = dup (2);
    FILE * fp = tmpfile();
    dup2 (fileno (fp), 2);
    fclose (fp);
  }

  char * line = pstrdup (command,__func__,__FILE__,0);
  bool status = process_line (line);
  pfree (line,__func__,__FILE__,0);

  pfree (Display.error,__func__,__FILE__,0);
  Display.error = NULL;
  if (status) {
    if (VertexBuffer.type < 0)
      VertexBuffer.type = 0;
  }
  else {
    if (pid() == 0) {
      fflush (ferr);
      FILE * fp = fdopen (2, "r");
      Display.error = read_file_into_buffer (fp);
      int len = Display.error ? strlen (Display.error) : 0;
      if (len > 0 && Display.error[len - 1] == '\n')
 Display.error[len - 1] = '\0';
      fclose(fp);
    }
    else
      Display.error = pstrdup ("error (on slave process)",__func__,__FILE__,0);
    VertexBuffer.type = - 1;
  }

  if (pid() == 0) {

    dup2 (bak, 2);
    close (bak);
  }

  if (VertexBuffer.type < 0)
    ;
  else {
    if (pid() > 0) {
      if (VertexBuffer.normal->len < VertexBuffer.position->len)
 VertexBuffer.normal->len = 0;
      if (VertexBuffer.color->len < VertexBuffer.position->len)
 VertexBuffer.color->len = 0;
    }
   



                                               ;
  }
}

static int ws_send_array (int fd, Array * a, int status, int type,
     unsigned int * shift)
{
@if _MPI
  if (pid() == 0) {
    void * p = NULL;
    long len;
    for (int pe = 0; pe < npe(); pe++) {
      if (pe == 0)
 p = a->p, len = a->len;
      else {
 MPI_Status status;
 MPI_Recv (&len, 1, MPI_LONG, pe, 22, MPI_COMM_WORLD, &status);
 if (len > 0) {
   p = pmalloc (len,__func__,__FILE__,0);
   MPI_Recv (p, len, MPI_BYTE, pe, 23, MPI_COMM_WORLD, &status);
 }
 else
   p = NULL;
      }
      if (type == 0)
 shift[pe] = (pe > 0 ? shift[pe - 1] : 0) + len/(3*sizeof(float));
      else if (type == 1 && pe > 0)
 for (unsigned int i = 0; i < len/sizeof(unsigned int); i++)
   ((unsigned int *) p)[i] += shift[pe - 1];
      if (status >= 0 && len > 0 && ws_send (fd, p, len) < len)
 status = -1;
      if (pe > 0)
 pfree (p,__func__,__FILE__,0);
    }
  }
  else {
    MPI_Send (&a->len, 1, MPI_LONG, 0, 22, MPI_COMM_WORLD);
    if (a->len > 0)
      MPI_Send (a->p, a->len, MPI_BYTE, 0, 23, MPI_COMM_WORLD);
  }
@else
  if (status >= 0 && a->len > 0 && ws_send (fd, a->p, a->len) < a->len)
    status = -1;
@endif
  return status;
}

static int display_send (const char * command, int fd)
{
  int status = 0;

  ;

  unsigned int commandlen = strlen (command);
  unsigned int errorlen = Display.error ? strlen (Display.error) : 0;

  int paddedlen = 4*ceil(commandlen/4.);
  size_t len = 2*sizeof(unsigned int) + paddedlen;

  unsigned int lens[] = {VertexBuffer.position->len,
    VertexBuffer.normal->len,
    VertexBuffer.color->len,
    VertexBuffer.index->len}, glens[4];
  int type = VertexBuffer.type, gtype;
@if _MPI
  MPI_Reduce (lens, glens, 4, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Allreduce (&type, &gtype, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
@else
  for (int i = 0; i < 4; i++) glens[i] = lens[i];
  gtype = type;
@endif

  if (gtype < 0)
    len += errorlen;
  else
    len += 2*sizeof(int) +
      4*sizeof (unsigned int) + glens[0] + glens[1] + glens[2] + glens[3];

  if (pid() == 0) {
    if (ws_sendframe_init (fd, len, false, 2) < 0 ||
 ws_send (fd, &commandlen, sizeof(unsigned int)) < sizeof(unsigned int) ||
 ws_send (fd, &errorlen, sizeof(unsigned int)) < sizeof(unsigned int) ||
 ws_send (fd, command, commandlen) < commandlen)
      status = -1;


    for (int i = 0; i < paddedlen - commandlen && status >= 0; i++) {
      char c = '\0';
      if (ws_send (fd, &c, 1) < 1)
 status = -1;
    }
  }

  if (gtype < 0) {
    if (pid() == 0 && status >= 0 &&
 ws_send (fd, Display.error, errorlen) < errorlen)
      status = -1;
  }
  else {
    if (pid() == 0 && status >= 0 &&
 (ws_send (fd, &VertexBuffer.dim, sizeof(int)) < sizeof(int) ||
  ws_send (fd, &gtype, sizeof(int)) < sizeof(int) ||
  ws_send (fd, glens, 4*sizeof (unsigned int)) < 4*sizeof (unsigned int)))
      status = -1;
    unsigned int * shift = pmalloc (sizeof(unsigned int)*npe(),__func__,__FILE__,0);
    status = ws_send_array (fd, VertexBuffer.position, status, 0, shift);
    status = ws_send_array (fd, VertexBuffer.normal, status, -1, shift);
    status = ws_send_array (fd, VertexBuffer.color, status, -1, shift);
    status = ws_send_array (fd, VertexBuffer.index, status, 1, shift);
    pfree (shift,__func__,__FILE__,0);
  }

@if _MPI
  MPI_Bcast (&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
@endif

  return status;
}

static void display_add (const char * command, int fd)
{
  ;
  khiter_t k = kh_get_strhash(Display.objects, command);
  if (k == ((Display.objects)->n_buckets)) {
    int ret;
    k = kh_put_strhash(Display.objects, pstrdup (command,__func__,__FILE__,0), &ret);
    DisplayClient * client = pmalloc (sizeof(DisplayClient),__func__,__FILE__,0);
    client->fd = -1;
    ((Display.objects)->vals[k]) = client;
  }
  DisplayClient * clients = ((Display.objects)->vals[k]), * client = clients;
  int len = 0;
  while (client->fd >= 0) {
    if (client->fd == fd)
      fd = -1;
    client++, len++;
  }
  if (fd >= 0) {
    ((Display.objects)->vals[k]) = clients =
      prealloc (clients, (len + 2)*sizeof (DisplayClient),__func__,__FILE__,0);
    clients[len].fd = fd;
    clients[len].iter = -1;
    clients[len + 1].fd = -1;
  }
  display_display();
}




static void array_remove (khiter_t k, int fd)
{
  DisplayClient * clients = ((Display.objects)->vals[k]), * client = clients;
  int i = -1, len = 0;
  while (client->fd >= 0) {
    if (client->fd == fd) {
      if (i != -1)

                                        ;
      i = len;
    }
    client++, len++;
  }
  if (i < 0)
   
                                    ;
  else if (len == 1) {
    pfree ((void *) ((Display.objects)->keys[k]),__func__,__FILE__,0);
    pfree ((void *) ((Display.objects)->vals[k]),__func__,__FILE__,0);
    kh_del_strhash(Display.objects, k);
  }
  else
    for (int j = i; j < len; j++)
      clients[j] = clients[j + 1];
}

static void display_remove (const char * command, int fd)
{
  ;
  khiter_t k = kh_get_strhash(Display.objects, command);
  if (k == ((Display.objects)->n_buckets))
   
                ;
  else
    array_remove (k, fd);
  display_display();
}

typedef struct {
  char * name, * tooltip;
  void * ptr;
  double min, max;
  int size;
} DisplayControl;

static char * display_control_json()
{
  char * build = pmalloc (4096,__func__,__FILE__,0);
  int len = 0;
  len += snprintf (build + len, 4096, "#{"), build = prealloc (build, len + 4096,__func__,__FILE__,0);
  char sep = ' ';
  DisplayControl * d = Display.controls->p;
  for (int i = 0; i < Display.controls->len/sizeof(DisplayControl); i++, d++) {
    len += snprintf (build + len, 4096, "%c\n  \"%s\": { ", sep, d->name), build = prealloc (build, len + 4096,__func__,__FILE__,0); sep = ',';
    if (d->tooltip)
      len += snprintf (build + len, 4096, "\"tooltip\": \"%s\", ", d->tooltip), build = prealloc (build, len + 4096,__func__,__FILE__,0);
    switch (d->size) {
    case 4:
      len += snprintf (build + len, 4096, "\"type\": \"int\", \"value\": %d, \"min\": %g, \"max\": %g", *((int *)d->ptr), d->min, d->max), build = prealloc (build, len + 4096,__func__,__FILE__,0)
                                     ;
      break;
    case 8:
      len += snprintf (build + len, 4096, "\"type\": \"double\", \"value\": %g, " "\"min\": %g, \"max\": %g", *((double *)d->ptr), d->min, d->max), build = prealloc (build, len + 4096,__func__,__FILE__,0)

                                        ;
      break;
    default:
      if (!(false)) qassert ("/home/spencer/basilisk/src/display.h", 0, "false");
    }
    len += snprintf (build + len, 4096, " }"), build = prealloc (build, len + 4096,__func__,__FILE__,0);
  }
  len += snprintf (build + len, 4096, "}"), build = prealloc (build, len + 4096,__func__,__FILE__,0);
  return build;
}

static DisplayControl * display_control_lookup (const char * name)
{
  DisplayControl * d = Display.controls->p;
  for (int i = 0; i < Display.controls->len/sizeof(DisplayControl); i++, d++)
      if (!strcmp (d->name, name))
 return d;
  return NULL;
}

void display_control_internal (void * ptr,
          double min, double max,
          char * name, char * tooltip,
          char * ptr_name, int size)
{
  DisplayControl d;
  if (!name)
    name = ptr_name;

  if (display_control_lookup (name))
    return;

  d.name = pstrdup (name,__func__,__FILE__,0);
  d.tooltip = tooltip ? pstrdup (tooltip,__func__,__FILE__,0) : NULL;
  d.ptr = ptr;
  d.size = size;
  d.min = min, d.max = max;
  array_append (Display.controls, &d, sizeof (DisplayControl));

  if (pid() == 0) {
    char * controls = display_control_json();
    ws_sendframe_txt (0, controls, true);
    pfree (controls,__func__,__FILE__,0);
  }
}





static void display_control_update (const char * command, int fd)
{
  char * s = pstrdup (command,__func__,__FILE__,0), * s1 = strchr (command, ':');
  *s1++ = '\0';
  DisplayControl * d = display_control_lookup (command);
  if (d == NULL)
   
                ;
  else {
    ;
    double val = atof(s1);
    if (d->max > d->min)
      val = clamp (val, d->min, d->max);
    switch (d->size) {
    case 4: *((int *)d->ptr) = val; break;
    case 8: *((double *)d->ptr) = val; break;
    default: if (!(false)) qassert ("/home/spencer/basilisk/src/display.h", 0, "false");
    }

    if (pid() == 0) {
      char * controls = display_control_json();
      ws_sendframe_txt (- fd, controls, true);
      pfree (controls,__func__,__FILE__,0);
    }
  }
  pfree (s,__func__,__FILE__,0);
}

static char * bview_interface_json()
{
  char * build = pmalloc (4096,__func__,__FILE__,0);
  int len = 0;

  len += snprintf (build + len, 4096, "{\n"), build = prealloc (build, len + 4096,__func__,__FILE__,0);

  int i = 0;
  while (bview_interface[i].json) {
    len += snprintf (build + len, 4096, "%s", i ? ",\n" : ""), build = prealloc (build, len + 4096,__func__,__FILE__,0);
    len += bview_interface[i].json (build + len, 4096);
    build = prealloc (build, len + 4096,__func__,__FILE__,0);
    len += snprintf (build + len, 4096, "\n"), build = prealloc (build, len + 4096,__func__,__FILE__,0);
    i++;
  }
  len += snprintf (build + len, 4096, "}"), build = prealloc (build, len + 4096,__func__,__FILE__,0);
  return build;
}

void display_onclose (int fd)
{
  ;
  for (khiter_t k = (khint_t)(0); k != ((Display.objects)->n_buckets);
       ++k)
    if ((!(((Display.objects)->flags[(k)>>4]>>(((k)&0xfU)<<1))&3)))
      array_remove (k, fd);
  display_display();
}

void display_onmessage (int fd, const char * msg, size_t size, int type)
{
  if (type == 1) {
    if (!msg)
      fprintf (ferr, "error receiving data on websocket\n");
    else switch (msg[0]) {
      case '+': display_add (msg + 1, fd); break;
      case '-': display_remove (msg + 1, fd); break;
      case '#': display_control_update (msg + 1, fd); break;
      default: fprintf (ferr,
   "display_onmessage: error: unknown message type '%s'\n",
   msg);
 break;
      }
  }
  else
    fprintf (ferr, "display_onmessage: error: unexpected message type '%d'\n",
      type);
}

void display_onopen (int fd)
{
  char * interface = bview_interface_json();
  char * controls = display_control_json();
  int status = 0;

  if (pid() == 0)
    if (ws_sendframe_txt (fd, interface, false) < 0 ||
 ws_sendframe_txt (fd, controls, false) < 0 ||
 (display_defaults && ws_sendframe_txt (fd, display_defaults, false) < 0))
      status = -1;
  pfree (interface,__func__,__FILE__,0);
  pfree (controls,__func__,__FILE__,0);

@if _MPI
  MPI_Bcast (&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
@endif

  ;

  if (status < 0) {
    display_onclose (fd);
    if (pid() == 0)
      close (fd);
  }
}

static void display_update (int i)
{
  for (khiter_t k = (khint_t)(0); k != ((Display.objects)->n_buckets);
       ++k)
    if ((!(((Display.objects)->flags[(k)>>4]>>(((k)&0xfU)<<1))&3))) {
      DisplayClient * client = ((Display.objects)->vals[k]);
      while (client->fd >= 0) {
 if (client->iter < i)
   break;
 client++;
      }
      if (client->fd >= 0) {
 const char * command = ((Display.objects)->keys[k]);
 display_command (command);
 client = ((Display.objects)->vals[k]);
 while (client->fd >= 0) {
   if (client->iter < i) {
       client->iter = i;
       if (display_send (command, client->fd) < 0) {
  ;
  if (pid() == 0)
    close (client->fd);
  display_onclose (client->fd);
  if (!(!(((Display.objects)->flags[(k)>>4]>>(((k)&0xfU)<<1))&3)))
    break;
       }
       else
  client++;
   }
   else
     client++;
 }
 vertex_buffer_free();
 if (Display.error && (!(((Display.objects)->flags[(k)>>4]>>(((k)&0xfU)<<1))&3))) {
   pfree ((void *) ((Display.objects)->keys[k]),__func__,__FILE__,0);
   pfree ((void *) ((Display.objects)->vals[k]),__func__,__FILE__,0);
   kh_del_strhash(Display.objects, k);
 }
      }
    }
}

@if _MPI
static Array * pack_messages (struct ws_message * messages)
{
  struct ws_message * msg = messages;
  Array * packed = array_new();
  while (msg && msg->fd >= 0) {
    array_append (packed, msg, sizeof(struct ws_message));
    array_append (packed, msg->msg, msg->size);
    msg++;
  }
  return packed;
}

static struct ws_message * unpack_messages (Array * packed)
{
  Array * array = array_new();
  char * p = packed->p;
  while (p - (char *) packed->p < packed->len) {
    struct ws_message * msg = (struct ws_message *) p;
    msg->msg = sysmalloc (msg->size + 1);
    p += sizeof(struct ws_message);
    memcpy (msg->msg, p, msg->size);
    msg->msg[msg->size] = '\0';
    array_append (array, msg, sizeof(struct ws_message));
    p += msg->size;
  }
  struct ws_message msg = {-1};
  array_append (array, &msg, sizeof(struct ws_message));
  struct ws_message * messages = array->p;
  pfree (array,__func__,__FILE__,0);
  return messages;
}
@endif

int display_poll (int timeout)
{
  struct ws_message * messages = NULL, * msg;
  int nmsg = 0;
  if (pid() == 0) {
    msg = messages = ws_socket_poll (Display.sock, timeout);
    while (msg && msg->fd >= 0) msg++, nmsg++;
  }

@if _MPI
  Array * packed;
  if (pid() == 0)
    packed = pack_messages (messages);
  else
    packed = pcalloc (1, sizeof (Array),__func__,__FILE__,0);
  MPI_Bcast (&packed->len, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  if (packed->len > 0) {
    if (pid() != 0)
      packed->p = pmalloc (packed->len,__func__,__FILE__,0);
    MPI_Bcast (packed->p, packed->len, MPI_BYTE, 0, MPI_COMM_WORLD);
    if (pid() != 0)
      messages = unpack_messages (packed);
  }
  array_free (packed);
@endif

  msg = messages;
  nmsg = 0;
  while (msg && msg->fd >= 0) {
    switch (msg->type) {

    case 0:
      display_onopen (msg->fd); break;

    case 1: case 2:
      display_onmessage (msg->fd, msg->msg, msg->size, msg->type);
      break;

    case 8:
      display_onclose (msg->fd); break;

    default:
      if (!(false)) qassert ("/home/spencer/basilisk/src/display.h", 0, "false");

    }
    sysfree (msg->msg);
    msg++, nmsg++;
  }
  sysfree (messages);
  return nmsg;
}

void display_url (FILE * fp)
{
  char hostname[1024];
  hostname[1023] = '\0';
  gethostname (hostname, 1023);
  struct hostent * h = gethostbyname (hostname);
  if (!h)
    fprintf (ferr,
      "display_url(): warning: gethostbyname(\"%s\") returned NULL\n",
      hostname);
  fprintf (fp, "http://basilisk.fr/three.js/editor/index.html" "?ws://%s:%d", h ? h->h_name : "127.0.0.1",
    Display.port);
}

int display_usage = 20;



int display_play = 1 < 0 ? -1 : 0;




static void display_destroy()
{
  for (khiter_t k = (khint_t)(0); k != ((Display.objects)->n_buckets);
       ++k)
    if ((!(((Display.objects)->flags[(k)>>4]>>(((k)&0xfU)<<1))&3))) {
      pfree ((void *) ((Display.objects)->keys[k]),__func__,__FILE__,0);
      pfree ((void *) ((Display.objects)->vals[k]),__func__,__FILE__,0);
    }
  kh_destroy_strhash(Display.objects);

  DisplayControl * d = Display.controls->p;
  for (int i = 0; i < Display.controls->len/sizeof(DisplayControl); i++, d++) {
    pfree (d->name,__func__,__FILE__,0);
    pfree (d->tooltip,__func__,__FILE__,0);
  }
  array_free (Display.controls);

  if (pid() == 0) {
    remove ("display.html");


    close (Display.sock);
  }
}

            
void display_init()
{
  if (pid() == 0) {
    const char * port = "7100:7200";
    if (!strchr (port, ':'))
      Display.sock = ws_socket_open (atoi (port));
    else {
      char * s = pstrdup (port,__func__,__FILE__,0);
      char * s1 = strchr (s, ':');
      *s1++ = '\0';
      int pmax = atoi(s1);
      Display.port = atoi(s);
      while ((Display.sock = ws_socket_open (Display.port)) < 0 &&
      Display.port < pmax)
 Display.port++;
      pfree (s,__func__,__FILE__,0);
    }
    if (Display.sock < 0) {
      char s[80];
      sprintf (s, "display(): could not open port '%s'", port);
      perror (s);
      exit (1);
    }

    FILE * fp = fopen ("display.html", "w");
    if (!fp)
      perror ("display.html");
    else {
      fputs ("<head><meta http-equiv=\"refresh\" content=\"0;URL=", fp);
      display_url (fp);
      fputs ("\"></head>\n", fp);
      fclose (fp);
    }
  }

  Display.objects = kh_init_strhash();
  Display.controls = array_new();

  free_solver_func_add (display_destroy);




  display_control_internal (&(display_play), -1, 1, "Run/Pause"
#line 438
, NULL
#line 790
, "display_play", sizeof(display_play));


  display_control_internal (&(display_usage), 0, 50, "Display %", "maximum % of runtime used by display", "display_usage", sizeof(display_usage))
                                            ;

}

static int refresh_display_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}


#line 798
      static int refresh_display(const int i,const double t,Event *_ev){tracing("refresh_display","/home/spencer/basilisk/src/display.h",0);
{
  do {
    if (display_play)
      display_update (i);
    if (display_poll (display_play ? - 1 : 0))
      display_update (i);
  } while (display_play < 0);

  static timer global_timer = {0};
  static double poll_elapsed = 0.;
  int refresh = (poll_elapsed <=
   display_usage/100.*timer_elapsed (global_timer));
@if _MPI
  MPI_Bcast (&refresh, 1, MPI_INT, 0, MPI_COMM_WORLD);
@endif
  if (refresh) {
    global_timer = timer_start();
    display_update (i);
    poll_elapsed = timer_elapsed (global_timer);
  }
}{end_tracing("refresh_display","/home/spencer/basilisk/src/display.h",0);return 0;}end_tracing("refresh_display","/home/spencer/basilisk/src/display.h",0);}
# 434 "/home/spencer/basilisk/src/utils.h" 2
# 12 "/home/spencer/basilisk/src/run.h" 2

     
void run (void)
{tracing("run","/home/spencer/basilisk/src/run.h",0);
  iter = 0, t = 0., dt = 1.;
  init_grid (N);

  perf.nc = perf.tnc = 0;
  perf.gt = timer_start();
  while (events (true)) {





    update_perf();
    iter = inext, t = tnext;
  }




  timer_print (perf.gt, iter, perf.tnc);

  free_grid();
end_tracing("run","/home/spencer/basilisk/src/run.h",0);}




static int defaults_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0)!=0;*ip=i;*tp=t;return ret;}





#line 42
      static int defaults(const int i,const double t,Event *_ev){tracing("defaults","/home/spencer/basilisk/src/run.h",0); {
  display ("box();"
#line 1425 "/home/spencer/basilisk/src/common.h"
, false
#line 43 "/home/spencer/basilisk/src/run.h"
);
}{end_tracing("defaults","/home/spencer/basilisk/src/run.h",0);return 0;}end_tracing("defaults","/home/spencer/basilisk/src/run.h",0);}





static int cleanup_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(t = TEND_EVENT)!=0;*ip=i;*tp=t;return ret;}






#line 50
      static int cleanup(const int i,const double t,Event *_ev){tracing("cleanup","/home/spencer/basilisk/src/run.h",0); {
  display ("", true);
}{end_tracing("cleanup","/home/spencer/basilisk/src/run.h",0);return 0;}end_tracing("cleanup","/home/spencer/basilisk/src/run.h",0);}
# 28 "/home/spencer/basilisk/src/navier-stokes/centered.h" 2
# 1 "./timestep.h" 1
# 1 "/home/spencer/basilisk/src/timestep.h"

double timestep (const vector u, double dtmax)
{
  static double previous = 0.;
  dtmax /= CFL;
  if(!is_constant(cm)){
  
#line 6
foreach_face_stencil(1,{(NonLocal[]){{"dtmax","double",(void *)&dtmax,NULL,0,'m'},{"cm","scalar",(void *)&cm,NULL,0},{"u","vector",(void *)&u,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 6 \"/home/spencer/basilisk/src/timestep.h\"\n{is_face_x(){\n    if (val(u.x,0,0,0) != 0.) {\n      real dt = Delta/fabs(val(u.x,0,0,0));\n\n\n\n\n      dt *= val(cm,0,0,0);\n\n      if (dt < dtmax) dtmax = dt;\n    }}end_is_face_x()\n// #line 6\nis_face_y(){\n    if (val(u.y,0,0,0) != 0.) {\n      real dt = Delta/fabs(val(u.y,0,0,0));\n\n\n\n\n      dt *= val(cm,0,0,0);\n\n      if (dt < dtmax) dtmax = dt;\n    }}end_is_face_y()}")}){_stencil_is_face_x(){
    {_stencil_val(u.x,0,0,0); {   
      _stencil_val(u.x,0,0,0); 




_stencil_val(cm,0,0,0);   




       

         
    
#line 16
}   }}end__stencil_is_face_x()
#line 6
_stencil_is_face_y(){
    {_stencil_val(u.y,0,0,0); {   
      _stencil_val(u.y,0,0,0); 




_stencil_val(cm,0,0,0);   




       

         
    
#line 16
}   }}end__stencil_is_face_y()}end_foreach_face_stencil();
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(min:dtmax)){
#line 6
foreach_face_generic(){is_face_x(){
    if (val(u.x,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.x,0,0,0));




      dt *= val(cm,0,0,0);

      if (dt < dtmax) dtmax = dt;
    }}end_is_face_x()
#line 6
is_face_y(){
    if (val(u.y,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.y,0,0,0));




      dt *= val(cm,0,0,0);

      if (dt < dtmax) dtmax = dt;
    }}end_is_face_y()}end_foreach_face_generic();mpi_all_reduce_array(&dtmax,double,MPI_MIN,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 16
}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#line 6
foreach_face_stencil(1,{(NonLocal[]){{"dtmax","double",(void *)&dtmax,NULL,0,'m'},{"_const_cm","double",(void *)&_const_cm,NULL,0},{"u","vector",(void *)&u,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 6 \"/home/spencer/basilisk/src/timestep.h\"\n{is_face_x(){\n    if (val(u.x,0,0,0) != 0.) {\n      real dt = Delta/fabs(val(u.x,0,0,0));\n\n\n\n\n      dt *= _const_cm;\n\n      if (dt < dtmax) dtmax = dt;\n    }}end_is_face_x()\n// #line 6\nis_face_y(){\n    if (val(u.y,0,0,0) != 0.) {\n      real dt = Delta/fabs(val(u.y,0,0,0));\n\n\n\n\n      dt *= _const_cm;\n\n      if (dt < dtmax) dtmax = dt;\n    }}end_is_face_y()}")}){_stencil_is_face_x(){
    {_stencil_val(u.x,0,0,0); {   
      _stencil_val(u.x,0,0,0);




;   




       

         
    
#line 16
}   }}end__stencil_is_face_x()
#line 6
_stencil_is_face_y(){
    {_stencil_val(u.y,0,0,0); {   
      _stencil_val(u.y,0,0,0);




;   




       

         
    
#line 16
}   }}end__stencil_is_face_y()}end_foreach_face_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(min:dtmax)){
#line 6
foreach_face_generic(){is_face_x(){
    if (val(u.x,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.x,0,0,0));




      dt *= _const_cm;

      if (dt < dtmax) dtmax = dt;
    }}end_is_face_x()
#line 6
is_face_y(){
    if (val(u.y,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.y,0,0,0));




      dt *= _const_cm;

      if (dt < dtmax) dtmax = dt;
    }}end_is_face_y()}end_foreach_face_generic();mpi_all_reduce_array(&dtmax,double,MPI_MIN,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 16
}
  dtmax *= CFL;
  if (dtmax > previous)
    dtmax = (previous + 0.1*dtmax)/1.1;
  previous = dtmax;
  return dtmax;
}
# 29 "/home/spencer/basilisk/src/navier-stokes/centered.h" 2
# 1 "./bcg.h" 1
# 1 "/home/spencer/basilisk/src/bcg.h"
# 11 "/home/spencer/basilisk/src/bcg.h"
void tracer_fluxes (scalar f,
      vector uf,
      vector flux,
      double dt,
              scalar src)
{





  vector  g=new_vector("g");
  gradients (((scalar[]){f,{-1}}),((vector[]) {g,{{-1},{-1}}}));




  if(!is_constant(fm.x) && !is_constant(src)){




  
#line 28
foreach_face_stencil(1,{(NonLocal[]){{"flux","vector",(void *)&flux,NULL,0},{"g","vector",(void *)&g,NULL,0},{"src","scalar",(void *)&src,NULL,0},{"f","scalar",(void *)&f,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 28 \"/home/spencer/basilisk/src/bcg.h\"\n{is_face_x(){ {\n\n\n\n\n\n\n\n    real un = dt*val(uf.x,0,0,0)/(val(fm.x,0,0,0)*Delta + 0.), s = sign(un);\n    int i = -(s + 1.)/2.;\n    real f2 = val(f,i,0,0) + (val(src,0,0,0) + val(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;\n\n\n\n\n\n    if (val(fm.y,i,0,0) && val(fm.y,i,1,0)) {\n      real vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(val(fm.y,i,0,0) + val(fm.y,i,1,0));\n      real fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);\n      f2 -= dt*vn*fyy/(2.*Delta);\n    }\n// # 58 \"/home/spencer/basilisk/src/bcg.h\"\n    val_out_(flux.x,0,0,0) = f2*val(uf.x,0,0,0);\n  }}end_is_face_x()\n// #line 28\nis_face_y(){ {\n\n\n\n\n\n\n\n    real un = dt*val(uf.y,0,0,0)/(val(fm.y,0,0,0)*Delta + 0.), s = sign(un);\n    int i = -(s + 1.)/2.;\n    real f2 = val(f,0,i,0) + (val(src,0,0,0) + val(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;\n\n\n\n\n\n    if (val(fm.x,0,i,0) && val(fm.x,1,i,0)) {\n      real vn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(val(fm.x,0,i,0) + val(fm.x,1,i,0));\n      real fyy = vn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);\n      f2 -= dt*vn*fyy/(2.*Delta);\n    }\n// # 58 \"/home/spencer/basilisk/src/bcg.h\"\n    val_out_(flux.y,0,0,0) = f2*val(uf.y,0,0,0);\n  }}end_is_face_y()}")}){_stencil_is_face_x(){ {        







    _stencil_val(fm.x,0,0,0);_stencil_val(uf.x,0,0,0);              
    
    _stencil_val(g.x,o_stencil,0,0); _stencil_val(src,-1,0,0);_stencil_val(src,0,0,0);_stencil_val(f, o_stencil,0,0);





_stencil_val(fm.y,o_stencil,0,0); _stencil_val(fm.y,o_stencil,1,0); {     
       _stencil_val(fm.y,o_stencil,1,0);_stencil_val(fm.y,o_stencil,0,0); _stencil_val(uf.y,o_stencil,1,0);_stencil_val(uf.y,o_stencil,0,0);         
       _stencil_val(f,o_stencil,-1,0);_stencil_val(f, o_stencil,0,0);_stencil_val(f, o_stencil,0,0); _stencil_val(f,o_stencil,1,0);
        
    }









_stencil_val(uf.x,0,0,0);





      
# 58 "/home/spencer/basilisk/src/bcg.h"
    _stencil_val_a(flux.x,0,0,0);  
  }}end__stencil_is_face_x()
#line 28
_stencil_is_face_y(){ {        







    _stencil_val(fm.y,0,0,0);_stencil_val(uf.y,0,0,0);              
    
    _stencil_val(g.y,0,o_stencil,0); _stencil_val(src,0,-1,0);_stencil_val(src,0,0,0);_stencil_val(f,0, o_stencil,0);





_stencil_val(fm.x,0,o_stencil,0); _stencil_val(fm.x,1,o_stencil,0); {     
       _stencil_val(fm.x,1,o_stencil,0);_stencil_val(fm.x,0,o_stencil,0); _stencil_val(uf.x,1,o_stencil,0);_stencil_val(uf.x,0,o_stencil,0);         
       _stencil_val(f,-1,o_stencil,0);_stencil_val(f,0, o_stencil,0);_stencil_val(f,0, o_stencil,0); _stencil_val(f,1,o_stencil,0);
        
    }









_stencil_val(uf.y,0,0,0);





      
# 58 "/home/spencer/basilisk/src/bcg.h"
    _stencil_val_a(flux.y,0,0,0);  
  }}end__stencil_is_face_y()}end_foreach_face_stencil();{
#line 28
foreach_face_generic(){is_face_x(){ {







    double un = dt*val(uf.x,0,0,0)/(val(fm.x,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val(src,0,0,0) + val(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





    if (val(fm.y,i,0,0) && val(fm.y,i,1,0)) {
      double vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(val(fm.y,i,0,0) + val(fm.y,i,1,0));
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
# 58 "/home/spencer/basilisk/src/bcg.h"
    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  }}end_is_face_x()
#line 28
is_face_y(){ {







    double un = dt*val(uf.y,0,0,0)/(val(fm.y,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val(src,0,0,0) + val(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





    if (val(fm.x,0,i,0) && val(fm.x,1,i,0)) {
      double vn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(val(fm.x,0,i,0) + val(fm.x,1,i,0));
      double fyy = vn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
# 58 "/home/spencer/basilisk/src/bcg.h"
    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  }}end_is_face_y()}end_foreach_face_generic();}}else if(is_constant(fm.x) && !is_constant(src)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);




  
#line 28
foreach_face_stencil(1,{(NonLocal[]){{"flux","vector",(void *)&flux,NULL,0},{"g","vector",(void *)&g,NULL,0},{"src","scalar",(void *)&src,NULL,0},{"f","scalar",(void *)&f,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 28 \"/home/spencer/basilisk/src/bcg.h\"\n{is_face_x(){ {\n\n\n\n\n\n\n\n    real un = dt*val(uf.x,0,0,0)/(_const_fm.x*Delta + 0.), s = sign(un);\n    int i = -(s + 1.)/2.;\n    real f2 = val(f,i,0,0) + (val(src,0,0,0) + val(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;\n\n\n\n\n\n    if (_const_fm.y && _const_fm.y) {\n      real vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(_const_fm.y + _const_fm.y);\n      real fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);\n      f2 -= dt*vn*fyy/(2.*Delta);\n    }\n// # 58 \"/home/spencer/basilisk/src/bcg.h\"\n    val_out_(flux.x,0,0,0) = f2*val(uf.x,0,0,0);\n  }}end_is_face_x()\n// #line 28\nis_face_y(){ {\n\n\n\n\n\n\n\n    real un = dt*val(uf.y,0,0,0)/(_const_fm.y*Delta + 0.), s = sign(un);\n    int i = -(s + 1.)/2.;\n    real f2 = val(f,0,i,0) + (val(src,0,0,0) + val(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;\n\n\n\n\n\n    if (_const_fm.x && _const_fm.x) {\n      real vn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(_const_fm.x + _const_fm.x);\n      real fyy = vn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);\n      f2 -= dt*vn*fyy/(2.*Delta);\n    }\n// # 58 \"/home/spencer/basilisk/src/bcg.h\"\n    val_out_(flux.y,0,0,0) = f2*val(uf.y,0,0,0);\n  }}end_is_face_y()}")}){_stencil_is_face_x(){ {







;_stencil_val(uf.x,0,0,0);              
    
    _stencil_val(g.x,o_stencil,0,0); _stencil_val(src,-1,0,0);_stencil_val(src,0,0,0);_stencil_val(f, o_stencil,0,0);





;; {
;; _stencil_val(uf.y,o_stencil,1,0);_stencil_val(uf.y,o_stencil,0,0);         
       _stencil_val(f,o_stencil,-1,0);_stencil_val(f, o_stencil,0,0);_stencil_val(f, o_stencil,0,0); _stencil_val(f,o_stencil,1,0);
        
    }









_stencil_val(uf.x,0,0,0);





      
# 58 "/home/spencer/basilisk/src/bcg.h"
    _stencil_val_a(flux.x,0,0,0);  
  }}end__stencil_is_face_x()
#line 28
_stencil_is_face_y(){ {







;_stencil_val(uf.y,0,0,0);              
    
    _stencil_val(g.y,0,o_stencil,0); _stencil_val(src,0,-1,0);_stencil_val(src,0,0,0);_stencil_val(f,0, o_stencil,0);





;; {
;; _stencil_val(uf.x,1,o_stencil,0);_stencil_val(uf.x,0,o_stencil,0);         
       _stencil_val(f,-1,o_stencil,0);_stencil_val(f,0, o_stencil,0);_stencil_val(f,0, o_stencil,0); _stencil_val(f,1,o_stencil,0);
        
    }









_stencil_val(uf.y,0,0,0);





      
# 58 "/home/spencer/basilisk/src/bcg.h"
    _stencil_val_a(flux.y,0,0,0);  
  }}end__stencil_is_face_y()}end_foreach_face_stencil();




  {
#line 28
foreach_face_generic(){is_face_x(){ {







    double un = dt*val(uf.x,0,0,0)/(_const_fm.x*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val(src,0,0,0) + val(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





    if (_const_fm.y && _const_fm.y) {
      double vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(_const_fm.y + _const_fm.y);
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
# 58 "/home/spencer/basilisk/src/bcg.h"
    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  }}end_is_face_x()
#line 28
is_face_y(){ {







    double un = dt*val(uf.y,0,0,0)/(_const_fm.y*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val(src,0,0,0) + val(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





    if (_const_fm.x && _const_fm.x) {
      double vn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(_const_fm.x + _const_fm.x);
      double fyy = vn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
# 58 "/home/spencer/basilisk/src/bcg.h"
    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  }}end_is_face_y()}end_foreach_face_generic();}}else if(!is_constant(fm.x) && is_constant(src)){double _const_src=_constant[src.i-_NVARMAX];NOT_UNUSED(_const_src);




  
#line 28
foreach_face_stencil(1,{(NonLocal[]){{"flux","vector",(void *)&flux,NULL,0},{"g","vector",(void *)&g,NULL,0},{"_const_src","double",(void *)&_const_src,NULL,0},{"f","scalar",(void *)&f,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 28 \"/home/spencer/basilisk/src/bcg.h\"\n{is_face_x(){ {\n\n\n\n\n\n\n\n    real un = dt*val(uf.x,0,0,0)/(val(fm.x,0,0,0)*Delta + 0.), s = sign(un);\n    int i = -(s + 1.)/2.;\n    real f2 = val(f,i,0,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;\n\n\n\n\n\n    if (val(fm.y,i,0,0) && val(fm.y,i,1,0)) {\n      real vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(val(fm.y,i,0,0) + val(fm.y,i,1,0));\n      real fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);\n      f2 -= dt*vn*fyy/(2.*Delta);\n    }\n// # 58 \"/home/spencer/basilisk/src/bcg.h\"\n    val_out_(flux.x,0,0,0) = f2*val(uf.x,0,0,0);\n  }}end_is_face_x()\n// #line 28\nis_face_y(){ {\n\n\n\n\n\n\n\n    real un = dt*val(uf.y,0,0,0)/(val(fm.y,0,0,0)*Delta + 0.), s = sign(un);\n    int i = -(s + 1.)/2.;\n    real f2 = val(f,0,i,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;\n\n\n\n\n\n    if (val(fm.x,0,i,0) && val(fm.x,1,i,0)) {\n      real vn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(val(fm.x,0,i,0) + val(fm.x,1,i,0));\n      real fyy = vn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);\n      f2 -= dt*vn*fyy/(2.*Delta);\n    }\n// # 58 \"/home/spencer/basilisk/src/bcg.h\"\n    val_out_(flux.y,0,0,0) = f2*val(uf.y,0,0,0);\n  }}end_is_face_y()}")}){_stencil_is_face_x(){ {        







    _stencil_val(fm.x,0,0,0);_stencil_val(uf.x,0,0,0);              
    
    _stencil_val(g.x,o_stencil,0,0);;;_stencil_val(f, o_stencil,0,0);





_stencil_val(fm.y,o_stencil,0,0); _stencil_val(fm.y,o_stencil,1,0); {     
       _stencil_val(fm.y,o_stencil,1,0);_stencil_val(fm.y,o_stencil,0,0); _stencil_val(uf.y,o_stencil,1,0);_stencil_val(uf.y,o_stencil,0,0);         
       _stencil_val(f,o_stencil,-1,0);_stencil_val(f, o_stencil,0,0);_stencil_val(f, o_stencil,0,0); _stencil_val(f,o_stencil,1,0);
        
    }









_stencil_val(uf.x,0,0,0);





      
# 58 "/home/spencer/basilisk/src/bcg.h"
    _stencil_val_a(flux.x,0,0,0);  
  }}end__stencil_is_face_x()
#line 28
_stencil_is_face_y(){ {        







    _stencil_val(fm.y,0,0,0);_stencil_val(uf.y,0,0,0);              
    
    _stencil_val(g.y,0,o_stencil,0);;;_stencil_val(f,0, o_stencil,0);





_stencil_val(fm.x,0,o_stencil,0); _stencil_val(fm.x,1,o_stencil,0); {     
       _stencil_val(fm.x,1,o_stencil,0);_stencil_val(fm.x,0,o_stencil,0); _stencil_val(uf.x,1,o_stencil,0);_stencil_val(uf.x,0,o_stencil,0);         
       _stencil_val(f,-1,o_stencil,0);_stencil_val(f,0, o_stencil,0);_stencil_val(f,0, o_stencil,0); _stencil_val(f,1,o_stencil,0);
        
    }









_stencil_val(uf.y,0,0,0);





      
# 58 "/home/spencer/basilisk/src/bcg.h"
    _stencil_val_a(flux.y,0,0,0);  
  }}end__stencil_is_face_y()}end_foreach_face_stencil();




  {
#line 28
foreach_face_generic(){is_face_x(){ {







    double un = dt*val(uf.x,0,0,0)/(val(fm.x,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





    if (val(fm.y,i,0,0) && val(fm.y,i,1,0)) {
      double vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(val(fm.y,i,0,0) + val(fm.y,i,1,0));
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
# 58 "/home/spencer/basilisk/src/bcg.h"
    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  }}end_is_face_x()
#line 28
is_face_y(){ {







    double un = dt*val(uf.y,0,0,0)/(val(fm.y,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





    if (val(fm.x,0,i,0) && val(fm.x,1,i,0)) {
      double vn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(val(fm.x,0,i,0) + val(fm.x,1,i,0));
      double fyy = vn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
# 58 "/home/spencer/basilisk/src/bcg.h"
    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  }}end_is_face_y()}end_foreach_face_generic();}}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_src=_constant[src.i-_NVARMAX];NOT_UNUSED(_const_src);




  
#line 28
foreach_face_stencil(1,{(NonLocal[]){{"flux","vector",(void *)&flux,NULL,0},{"g","vector",(void *)&g,NULL,0},{"_const_src","double",(void *)&_const_src,NULL,0},{"f","scalar",(void *)&f,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 28 \"/home/spencer/basilisk/src/bcg.h\"\n{is_face_x(){ {\n\n\n\n\n\n\n\n    real un = dt*val(uf.x,0,0,0)/(_const_fm.x*Delta + 0.), s = sign(un);\n    int i = -(s + 1.)/2.;\n    real f2 = val(f,i,0,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;\n\n\n\n\n\n    if (_const_fm.y && _const_fm.y) {\n      real vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(_const_fm.y + _const_fm.y);\n      real fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);\n      f2 -= dt*vn*fyy/(2.*Delta);\n    }\n// # 58 \"/home/spencer/basilisk/src/bcg.h\"\n    val_out_(flux.x,0,0,0) = f2*val(uf.x,0,0,0);\n  }}end_is_face_x()\n// #line 28\nis_face_y(){ {\n\n\n\n\n\n\n\n    real un = dt*val(uf.y,0,0,0)/(_const_fm.y*Delta + 0.), s = sign(un);\n    int i = -(s + 1.)/2.;\n    real f2 = val(f,0,i,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;\n\n\n\n\n\n    if (_const_fm.x && _const_fm.x) {\n      real vn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(_const_fm.x + _const_fm.x);\n      real fyy = vn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);\n      f2 -= dt*vn*fyy/(2.*Delta);\n    }\n// # 58 \"/home/spencer/basilisk/src/bcg.h\"\n    val_out_(flux.y,0,0,0) = f2*val(uf.y,0,0,0);\n  }}end_is_face_y()}")}){_stencil_is_face_x(){ {







;_stencil_val(uf.x,0,0,0);              
    
    _stencil_val(g.x,o_stencil,0,0);;;_stencil_val(f, o_stencil,0,0);





;; {
;; _stencil_val(uf.y,o_stencil,1,0);_stencil_val(uf.y,o_stencil,0,0);         
       _stencil_val(f,o_stencil,-1,0);_stencil_val(f, o_stencil,0,0);_stencil_val(f, o_stencil,0,0); _stencil_val(f,o_stencil,1,0);
        
    }









_stencil_val(uf.x,0,0,0);





      
# 58 "/home/spencer/basilisk/src/bcg.h"
    _stencil_val_a(flux.x,0,0,0);  
  }}end__stencil_is_face_x()
#line 28
_stencil_is_face_y(){ {







;_stencil_val(uf.y,0,0,0);              
    
    _stencil_val(g.y,0,o_stencil,0);;;_stencil_val(f,0, o_stencil,0);





;; {
;; _stencil_val(uf.x,1,o_stencil,0);_stencil_val(uf.x,0,o_stencil,0);         
       _stencil_val(f,-1,o_stencil,0);_stencil_val(f,0, o_stencil,0);_stencil_val(f,0, o_stencil,0); _stencil_val(f,1,o_stencil,0);
        
    }









_stencil_val(uf.y,0,0,0);





      
# 58 "/home/spencer/basilisk/src/bcg.h"
    _stencil_val_a(flux.y,0,0,0);  
  }}end__stencil_is_face_y()}end_foreach_face_stencil();




  {
#line 28
foreach_face_generic(){is_face_x(){ {







    double un = dt*val(uf.x,0,0,0)/(_const_fm.x*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





    if (_const_fm.y && _const_fm.y) {
      double vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(_const_fm.y + _const_fm.y);
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
# 58 "/home/spencer/basilisk/src/bcg.h"
    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  }}end_is_face_x()
#line 28
is_face_y(){ {







    double un = dt*val(uf.y,0,0,0)/(_const_fm.y*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





    if (_const_fm.x && _const_fm.x) {
      double vn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(_const_fm.x + _const_fm.x);
      double fyy = vn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
# 58 "/home/spencer/basilisk/src/bcg.h"
    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  }}end_is_face_y()}end_foreach_face_generic();}}delete((scalar*)((vector[]){g,{{-1},{-1}}}));
}






void advection (scalar * tracers, vector u, double dt,
  scalar * src)
{




  scalar * psrc = src;
  if (!src)
    {scalar*_i=(scalar*)( tracers);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
      const scalar zero = new_const_scalar("zero",9, 0.);
      src = list_append (src, zero);
    }}}
  if (!(list_len (tracers) == list_len (src))) qassert ("/home/spencer/basilisk/src/bcg.h", 0, "list_len (tracers) == list_len (src)");

  scalar f, source;
  {scalar*_i0=src;scalar*_i1= tracers;if(_i0)for(source=*_i0,f=*_i1;_i0->i>= 0;source=*++_i0,f=*++_i1){ {
    vector  flux=new_face_vector("flux");
    tracer_fluxes (f, u, flux, dt, source);

    if(!is_constant(cm)){

    
#line 87
foreach_stencil(1,{(NonLocal[]){{"cm","scalar",(void *)&cm,NULL,0},{"flux","vector",(void *)&flux,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"f","scalar",(void *)&f,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n      \n// #line 88 \"/home/spencer/basilisk/src/bcg.h\"\n{\n        val_out_(f,0,0,0) += dt*(val(flux.x,0,0,0) - val(flux.x,1,0,0))/(Delta*val(cm,0,0,0));\n        \n// #line 89\nval_out_(f,0,0,0) += dt*(val(flux.y,0,0,0) - val(flux.y,0,1,0))/(Delta*val(cm,0,0,0));}")})
      {
        {_stencil_val(flux.x,0,0,0); _stencil_val(flux.x,1,0,0);_stencil_val(cm,0,0,0);_stencil_val_r(f,0,0,0);   }
        
#line 89
{_stencil_val(flux.y,0,0,0); _stencil_val(flux.y,0,1,0);_stencil_val(cm,0,0,0);_stencil_val_r(f,0,0,0);   }}end_foreach_stencil();{
#line 87
foreach()
      {
        val(f,0,0,0) += dt*(val(flux.x,0,0,0) - val(flux.x,1,0,0))/(Delta*val(cm,0,0,0));
        
#line 89
val(f,0,0,0) += dt*(val(flux.y,0,0,0) - val(flux.y,0,1,0))/(Delta*val(cm,0,0,0));}end_foreach();}}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);

    
#line 87
foreach_stencil(1,{(NonLocal[]){{"_const_cm","double",(void *)&_const_cm,NULL,0},{"flux","vector",(void *)&flux,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"f","scalar",(void *)&f,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n      \n// #line 88 \"/home/spencer/basilisk/src/bcg.h\"\n{\n        val_out_(f,0,0,0) += dt*(val(flux.x,0,0,0) - val(flux.x,1,0,0))/(Delta*_const_cm);\n        \n// #line 89\nval_out_(f,0,0,0) += dt*(val(flux.y,0,0,0) - val(flux.y,0,1,0))/(Delta*_const_cm);}")})
      {
        {_stencil_val(flux.x,0,0,0); _stencil_val(flux.x,1,0,0);;_stencil_val_r(f,0,0,0);   }
        
#line 89
{_stencil_val(flux.y,0,0,0); _stencil_val(flux.y,0,1,0);;_stencil_val_r(f,0,0,0);   }}end_foreach_stencil();

    {
#line 87
foreach()
      {
        val(f,0,0,0) += dt*(val(flux.x,0,0,0) - val(flux.x,1,0,0))/(Delta*_const_cm);
        
#line 89
val(f,0,0,0) += dt*(val(flux.y,0,0,0) - val(flux.y,0,1,0))/(Delta*_const_cm);}end_foreach();}}delete((scalar*)((vector[]){flux,{{-1},{-1}}}));



  }}}

  if (!psrc)
    pfree (src,__func__,__FILE__,0);
}
# 30 "/home/spencer/basilisk/src/navier-stokes/centered.h" 2



# 1 "./viscosity.h" 1
# 1 "/home/spencer/basilisk/src/viscosity.h"
# 50 "/home/spencer/basilisk/src/viscosity.h"
# 1 "./poisson.h" 1
# 1 "/home/spencer/basilisk/src/poisson.h"
# 32 "/home/spencer/basilisk/src/poisson.h"
void mg_cycle (scalar * a, scalar * res, scalar * da,
        void (* relax) (scalar * da, scalar * res,
          int depth, void * data),
        void * data,
        int nrelax, int minlevel, int maxlevel)
{




  restriction (res);





  minlevel = min (minlevel, maxlevel);
  for (int l = minlevel; l <= maxlevel; l++) {




    if (l == minlevel)
      {foreach_level_or_leaf (l)
 {scalar*_i=(scalar*)( da);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
  
     val(s,0,0,0) = 0.;}}end_foreach_level_or_leaf();}





    else
      {foreach_level (l)
 {scalar*_i=(scalar*)( da);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
  
     val(s,0,0,0) = bilinear (point, s);}}end_foreach_level();}





    boundary_level (da, l);
    for (int i = 0; i < nrelax; i++) {
      relax (da, res, l, data);
      boundary_level (da, l);
    }
  }




  foreach_stencil(1,{(NonLocal[]){{"da","scalar",(void *)da,NULL,1},{"a","scalar",(void *)a,NULL,1},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_(" \n// #line 84 \"/home/spencer/basilisk/src/poisson.h\"\n{\n    scalar s, ds;\n    {forin2 (s, ds , a, da)\n     \n val_out_(s,0,0,0) += val(ds,0,0,0); endforin2()}\n  }")}) {
    scalar s, ds;
    {scalar*_i0= da;scalar*_i1= a;if(_i0)for(ds=*_i0,s=*_i1;_i0->i>= 0;ds=*++_i0,s=*++_i1){
     
 { _stencil_val(ds,0,0,0);_stencil_val_r(s,0,0,0); }}}
  }end_foreach_stencil();




  {
#line 84
foreach() {
    scalar s, ds;
    {scalar*_i0= da;scalar*_i1= a;if(_i0)for(ds=*_i0,s=*_i1;_i0->i>= 0;ds=*++_i0,s=*++_i1){
     
 val(s,0,0,0) += val(ds,0,0,0);}}
  }end_foreach();}
}
# 102 "/home/spencer/basilisk/src/poisson.h"
int NITERMAX = 100, NITERMIN = 1;
double TOLERANCE = 1e-3;




typedef struct {
  int i;
  double resb, resa;
  double sum;
  int nrelax;
  int minlevel;
} mgstats;
# 125 "/home/spencer/basilisk/src/poisson.h"
mgstats mg_solve (scalar * a, scalar * b,
    double (* residual) (scalar * a, scalar * b, scalar * res,
           void * data),
    void (* relax) (scalar * da, scalar * res, int depth,
      void * data),
    void * data,
    int nrelax,
    scalar * res,
    int minlevel,
    double tolerance)
{





  scalar * da = list_clone (a), * pres = res;
  if (!res)
    res = list_clone (b);






  for (int b = 0; b < nboundary; b++)
    {scalar*_i=(scalar*)( da);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b];}}




  mgstats s = {0};
  double sum = 0.;
  scalar rhs = b[0];
  foreach_stencil (1,{(NonLocal[]){{"rhs","scalar",(void *)&rhs,NULL,0},{"sum","double",(void *)&sum,NULL,0,'+'},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n    \n// #line 161 \"/home/spencer/basilisk/src/poisson.h\"\nsum += val(rhs,0,0,0);")})
    { _stencil_val(rhs,0,0,0); }end_foreach_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(+:sum)){
#line 160
foreach ()
    sum += val(rhs,0,0,0);end_foreach();mpi_all_reduce_array(&sum,double,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
  
#line 162
s.sum = sum;
  s.nrelax = nrelax > 0 ? nrelax : 4;




  double resb;
  resb = s.resb = s.resa = (* residual) (a, b, res, data);






  for (s.i = 0;
       s.i < NITERMAX && (s.i < NITERMIN || s.resa > tolerance);
       s.i++) {
    mg_cycle (a, res, da, relax, data,
       s.nrelax,
       minlevel,
       grid->maxdepth);
    s.resa = (* residual) (a, b, res, data);
# 192 "/home/spencer/basilisk/src/poisson.h"
    if (s.resa > tolerance) {
      if (resb/s.resa < 1.2 && s.nrelax < 100)
 s.nrelax++;
      else if (resb/s.resa > 10 && s.nrelax > 2)
 s.nrelax--;
    }







    resb = s.resa;
  }
  s.minlevel = minlevel;




  if (s.resa > tolerance) {
    scalar v = a[0];
    fprintf (ferr,
      "WARNING: convergence for %s not reached after %d iterations\n"
      "  res: %g sum: %g nrelax: %d tolerance: %g\n", _attribute[v.i].name,
      s.i, s.resa, s.sum, s.nrelax, tolerance), fflush (ferr);
  }




  if (!pres)
    delete (res), pfree (res,__func__,__FILE__,0);
  delete (da), pfree (da,__func__,__FILE__,0);

  return s;
}
# 251 "/home/spencer/basilisk/src/poisson.h"
struct Poisson {
  scalar a, b;
          vector alpha;
          scalar lambda;
  double tolerance;
  int nrelax, minlevel;
  scalar * res;



};





static void relax (scalar * al, scalar * bl, int l, void * data)
{
  scalar a = al[0], b = bl[0];
  struct Poisson * p = (struct Poisson *) data;
          vector alpha = p->alpha;
          scalar lambda = p->lambda;
# 289 "/home/spencer/basilisk/src/poisson.h"
  scalar c = a;






  if(!is_constant(lambda) && !is_constant(alpha.x)){{foreach_level_or_leaf (l) {
    double n = - sq(Delta)*val(b,0,0,0), d = - val(lambda,0,0,0)*sq(Delta);
     {
      n += val(alpha.x,1,0,0)*val(a,1,0,0) + val(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val(alpha.x,1,0,0) + val(alpha.x,0,0,0);
    } 
#line 298
{
      n += val(alpha.y,0,1,0)*val(a,0,1,0) + val(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val(alpha.y,0,1,0) + val(alpha.y,0,0,0);
    }
# 312 "/home/spencer/basilisk/src/poisson.h"
      val(c,0,0,0) = n/d;
  }end_foreach_level_or_leaf();}}else if(is_constant(lambda) && !is_constant(alpha.x)){double _const_lambda=_constant[lambda.i-_NVARMAX];NOT_UNUSED(_const_lambda);






  {
#line 296
foreach_level_or_leaf (l) {
    double n = - sq(Delta)*val(b,0,0,0), d = - _const_lambda*sq(Delta);
     {
      n += val(alpha.x,1,0,0)*val(a,1,0,0) + val(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val(alpha.x,1,0,0) + val(alpha.x,0,0,0);
    } 
#line 298
{
      n += val(alpha.y,0,1,0)*val(a,0,1,0) + val(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val(alpha.y,0,1,0) + val(alpha.y,0,0,0);
    }
# 312 "/home/spencer/basilisk/src/poisson.h"
      val(c,0,0,0) = n/d;
  }end_foreach_level_or_leaf();}}else if(!is_constant(lambda) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);






  {
#line 296
foreach_level_or_leaf (l) {
    double n = - sq(Delta)*val(b,0,0,0), d = - val(lambda,0,0,0)*sq(Delta);
     {
      n += _const_alpha.x*val(a,1,0,0) + _const_alpha.x*val(a,-1,0,0);
      d += _const_alpha.x + _const_alpha.x;
    } 
#line 298
{
      n += _const_alpha.y*val(a,0,1,0) + _const_alpha.y*val(a,0,-1,0);
      d += _const_alpha.y + _const_alpha.y;
    }
# 312 "/home/spencer/basilisk/src/poisson.h"
      val(c,0,0,0) = n/d;
  }end_foreach_level_or_leaf();}}else {double _const_lambda=_constant[lambda.i-_NVARMAX];NOT_UNUSED(_const_lambda);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);






  {
#line 296
foreach_level_or_leaf (l) {
    double n = - sq(Delta)*val(b,0,0,0), d = - _const_lambda*sq(Delta);
     {
      n += _const_alpha.x*val(a,1,0,0) + _const_alpha.x*val(a,-1,0,0);
      d += _const_alpha.x + _const_alpha.x;
    } 
#line 298
{
      n += _const_alpha.y*val(a,0,1,0) + _const_alpha.y*val(a,0,-1,0);
      d += _const_alpha.y + _const_alpha.y;
    }
# 312 "/home/spencer/basilisk/src/poisson.h"
      val(c,0,0,0) = n/d;
  }end_foreach_level_or_leaf();}}
# 331 "/home/spencer/basilisk/src/poisson.h"
}






static double residual (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a = al[0], b = bl[0], res = resl[0];
  struct Poisson * p = (struct Poisson *) data;
          vector alpha = p->alpha;
          scalar lambda = p->lambda;
  double maxres = 0.;


  vector  g=new_face_vector("g");
  if(!is_constant(alpha.x)){
  
#line 348
foreach_face_stencil(1,{(NonLocal[]){{"a","scalar",(void *)&a,NULL,0},{"alpha","vector",(void *)&alpha,NULL,0},{"g","vector",(void *)&g,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 348 \"/home/spencer/basilisk/src/poisson.h\"\n{is_face_x(){\n    val_out_(g.x,0,0,0) = val(alpha.x,0,0,0)*((val(a,0,0,0) - val(a,0 -1,0,0))/Delta);}end_is_face_x()\n// #line 348\nis_face_y(){\n    val_out_(g.y,0,0,0) = val(alpha.y,0,0,0)*((val(a,0,0,0) - val(a,0,0 -1,0))/Delta);}end_is_face_y()}")}){_stencil_is_face_x(){
    { _stencil_val(alpha.x,0,0,0);_stencil_val(a,0,0,0); _stencil_val(a,0 -1,0,0);_stencil_val_a(g.x,0,0,0);  }}end__stencil_is_face_x()
#line 348
_stencil_is_face_y(){
    { _stencil_val(alpha.y,0,0,0);_stencil_val(a,0,0,0); _stencil_val(a,0,0 -1,0);_stencil_val_a(g.y,0,0,0);  }}end__stencil_is_face_y()}end_foreach_face_stencil();{
#line 348
foreach_face_generic(){is_face_x(){
    val(g.x,0,0,0) = val(alpha.x,0,0,0)*((val(a,0,0,0) - val(a,0 -1,0,0))/Delta);}end_is_face_x()
#line 348
is_face_y(){
    val(g.y,0,0,0) = val(alpha.y,0,0,0)*((val(a,0,0,0) - val(a,0,0 -1,0))/Delta);}end_is_face_y()}end_foreach_face_generic();}}else {_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  
#line 348
foreach_face_stencil(1,{(NonLocal[]){{"a","scalar",(void *)&a,NULL,0},{"_const_alpha","_coord",(void *)&_const_alpha,NULL,0},{"g","vector",(void *)&g,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 348 \"/home/spencer/basilisk/src/poisson.h\"\n{is_face_x(){\n    val_out_(g.x,0,0,0) = _const_alpha.x*((val(a,0,0,0) - val(a,0 -1,0,0))/Delta);}end_is_face_x()\n// #line 348\nis_face_y(){\n    val_out_(g.y,0,0,0) = _const_alpha.y*((val(a,0,0,0) - val(a,0,0 -1,0))/Delta);}end_is_face_y()}")}){_stencil_is_face_x(){
    {;_stencil_val(a,0,0,0); _stencil_val(a,0 -1,0,0);_stencil_val_a(g.x,0,0,0);  }}end__stencil_is_face_x()
#line 348
_stencil_is_face_y(){
    {;_stencil_val(a,0,0,0); _stencil_val(a,0,0 -1,0);_stencil_val_a(g.y,0,0,0);  }}end__stencil_is_face_y()}end_foreach_face_stencil();
  {
#line 348
foreach_face_generic(){is_face_x(){
    val(g.x,0,0,0) = _const_alpha.x*((val(a,0,0,0) - val(a,0 -1,0,0))/Delta);}end_is_face_x()
#line 348
is_face_y(){
    val(g.y,0,0,0) = _const_alpha.y*((val(a,0,0,0) - val(a,0,0 -1,0))/Delta);}end_is_face_y()}end_foreach_face_generic();}}
  if(!is_constant(lambda)){
  
#line 350
foreach_stencil (1,{(NonLocal[]){{"maxres","double",(void *)&maxres,NULL,0,'M'},{"g","vector",(void *)&g,NULL,0},{"a","scalar",(void *)&a,NULL,0},{"lambda","scalar",(void *)&lambda,NULL,0},{"b","scalar",(void *)&b,NULL,0},{"res","scalar",(void *)&res,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_(" \n// #line 350 \"/home/spencer/basilisk/src/poisson.h\"\n{\n    val_out_(res,0,0,0) = val(b,0,0,0) - val(lambda,0,0,0)*val(a,0,0,0);\n    \n      val_out_(res,0,0,0) -= (val(g.x,1,0,0) - val(g.x,0,0,0))/Delta;\n      \n// #line 353\nval_out_(res,0,0,0) -= (val(g.y,0,1,0) - val(g.y,0,0,0))/Delta;\n\n\n\n\n\n\n    if (fabs (val(res,0,0,0)) > maxres)\n      maxres = fabs (val(res,0,0,0));\n  }")}) { 
_stencil_val(b,0,0,0); _stencil_val(lambda,0,0,0);_stencil_val(a,0,0,0);
    
#line 351
_stencil_val_a(res,0,0,0);  
    
      {_stencil_val(g.x,1,0,0); _stencil_val(g.x,0,0,0);_stencil_val_r(res,0,0,0);   }
      
#line 353
{_stencil_val(g.y,0,1,0); _stencil_val(g.y,0,0,0);_stencil_val_r(res,0,0,0);   }






_stencil_val(res,0,0,0);
      {_stencil_val(res,0,0,0);   }






        
  
#line 362
}end_foreach_stencil();
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 350
foreach () {
    val(res,0,0,0) = val(b,0,0,0) - val(lambda,0,0,0)*val(a,0,0,0);
    
      val(res,0,0,0) -= (val(g.x,1,0,0) - val(g.x,0,0,0))/Delta;
      
#line 353
val(res,0,0,0) -= (val(g.y,0,1,0) - val(g.y,0,0,0))/Delta;






    if (fabs (val(res,0,0,0)) > maxres)
      maxres = fabs (val(res,0,0,0));
  }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 362
}else {double _const_lambda=_constant[lambda.i-_NVARMAX];NOT_UNUSED(_const_lambda);
  
#line 350
foreach_stencil (1,{(NonLocal[]){{"maxres","double",(void *)&maxres,NULL,0,'M'},{"g","vector",(void *)&g,NULL,0},{"a","scalar",(void *)&a,NULL,0},{"_const_lambda","double",(void *)&_const_lambda,NULL,0},{"b","scalar",(void *)&b,NULL,0},{"res","scalar",(void *)&res,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_(" \n// #line 350 \"/home/spencer/basilisk/src/poisson.h\"\n{\n    val_out_(res,0,0,0) = val(b,0,0,0) - _const_lambda*val(a,0,0,0);\n    \n      val_out_(res,0,0,0) -= (val(g.x,1,0,0) - val(g.x,0,0,0))/Delta;\n      \n// #line 353\nval_out_(res,0,0,0) -= (val(g.y,0,1,0) - val(g.y,0,0,0))/Delta;\n\n\n\n\n\n\n    if (fabs (val(res,0,0,0)) > maxres)\n      maxres = fabs (val(res,0,0,0));\n  }")}) { 
_stencil_val(b,0,0,0);;_stencil_val(a,0,0,0);
    
#line 351
_stencil_val_a(res,0,0,0);  
    
      {_stencil_val(g.x,1,0,0); _stencil_val(g.x,0,0,0);_stencil_val_r(res,0,0,0);   }
      
#line 353
{_stencil_val(g.y,0,1,0); _stencil_val(g.y,0,0,0);_stencil_val_r(res,0,0,0);   }






_stencil_val(res,0,0,0);
      {_stencil_val(res,0,0,0);   }






        
  
#line 362
}end_foreach_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 350
foreach () {
    val(res,0,0,0) = val(b,0,0,0) - _const_lambda*val(a,0,0,0);
    
      val(res,0,0,0) -= (val(g.x,1,0,0) - val(g.x,0,0,0))/Delta;
      
#line 353
val(res,0,0,0) -= (val(g.y,0,1,0) - val(g.y,0,0,0))/Delta;






    if (fabs (val(res,0,0,0)) > maxres)
      maxres = fabs (val(res,0,0,0));
  }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 362
}
# 380 "/home/spencer/basilisk/src/poisson.h"
  {delete((scalar*)((vector[]){g,{{-1},{-1}}}));return maxres;}delete((scalar*)((vector[]){g,{{-1},{-1}}}));
}
# 392 "/home/spencer/basilisk/src/poisson.h"
mgstats poisson (scalar a, scalar b,
           vector alpha,
           scalar lambda,
   double tolerance,
   int nrelax,
   int minlevel,
   scalar * res,
   double (* flux) (Point, scalar, vector, double *))
{






  if (alpha.x.i < 0)
    alpha = unityf;
  if (lambda.i < 0) {
    const scalar zeroc = new_const_scalar("zeroc",10, 0.);
    lambda = zeroc;
  }




  restriction (((scalar[]){alpha.x,alpha.y,lambda,{-1}}));





  double defaultol = TOLERANCE;
  if (tolerance)
    TOLERANCE = tolerance;

  struct Poisson p = {a, b, alpha, lambda, tolerance, nrelax, minlevel, res };






  mgstats s = mg_solve ((
#line 125
scalar *
#line 434
)((scalar[]){a,{-1}}),( 
#line 125
scalar *
#line 434
)((scalar[]) {b,{-1}}), residual, relax, &p
,
   
#line 435
nrelax, res, max(1, minlevel)
#line 133
, 
TOLERANCE
#line 435
);




  if (tolerance)
    TOLERANCE = defaultol;

  return s;
}
# 463 "/home/spencer/basilisk/src/poisson.h"
     
mgstats project (vector uf, scalar p,
           vector alpha,
   double dt,
   int nrelax)
{tracing("project","/home/spencer/basilisk/src/poisson.h",0);






  scalar  div=new_scalar("div");
  foreach_stencil(1,{(NonLocal[]){{"dt","double",(void *)&dt,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"div","scalar",(void *)&div,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_(" \n// #line 476 \"/home/spencer/basilisk/src/poisson.h\"\n{\n    val_out_(div,0,0,0) = 0.;\n    \n      val_out_(div,0,0,0) += val(uf.x,1,0,0) - val(uf.x,0,0,0);\n      \n// #line 479\nval_out_(div,0,0,0) += val(uf.y,0,1,0) - val(uf.y,0,0,0);\n    val_out_(div,0,0,0) /= dt*Delta;\n  }")}) {
    _stencil_val_a(div,0,0,0);  
    
      { _stencil_val(uf.x,1,0,0); _stencil_val(uf.x,0,0,0);_stencil_val_r(div,0,0,0);  }
      
#line 479
{ _stencil_val(uf.y,0,1,0); _stencil_val(uf.y,0,0,0);_stencil_val_r(div,0,0,0);  }
    _stencil_val_r(div,0,0,0);  
  }end_foreach_stencil();
  {
#line 476
foreach() {
    val(div,0,0,0) = 0.;
    
      val(div,0,0,0) += val(uf.x,1,0,0) - val(uf.x,0,0,0);
      
#line 479
val(div,0,0,0) += val(uf.y,0,1,0) - val(uf.y,0,0,0);
    val(div,0,0,0) /= dt*Delta;
  }end_foreach();}
# 492 "/home/spencer/basilisk/src/poisson.h"
  mgstats mgp = poisson (p, div, alpha
#line 393
,
( scalar) {-1}
#line 493
, TOLERANCE/sq(dt), nrelax
#line 396
, 
0, 
NULL, 
NULL
#line 493
);




  if(!is_constant(alpha.x)){




  
#line 498
foreach_face_stencil(1,{(NonLocal[]){{"p","scalar",(void *)&p,NULL,0},{"alpha","vector",(void *)&alpha,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 498 \"/home/spencer/basilisk/src/poisson.h\"\n{is_face_x(){\n    val_out_(uf.x,0,0,0) -= dt*val(alpha.x,0,0,0)*((val(p,0,0,0) - val(p,0 -1,0,0))/Delta);}end_is_face_x()\n// #line 498\nis_face_y(){\n    val_out_(uf.y,0,0,0) -= dt*val(alpha.y,0,0,0)*((val(p,0,0,0) - val(p,0,0 -1,0))/Delta);}end_is_face_y()}")}){_stencil_is_face_x(){
    {_stencil_val(alpha.x,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0 -1,0,0);_stencil_val_r(uf.x,0,0,0);   }}end__stencil_is_face_x()
#line 498
_stencil_is_face_y(){
    {_stencil_val(alpha.y,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0,0 -1,0);_stencil_val_r(uf.y,0,0,0);   }}end__stencil_is_face_y()}end_foreach_face_stencil();{
#line 498
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) -= dt*val(alpha.x,0,0,0)*((val(p,0,0,0) - val(p,0 -1,0,0))/Delta);}end_is_face_x()
#line 498
is_face_y(){
    val(uf.y,0,0,0) -= dt*val(alpha.y,0,0,0)*((val(p,0,0,0) - val(p,0,0 -1,0))/Delta);}end_is_face_y()}end_foreach_face_generic();}}else {_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);




  
#line 498
foreach_face_stencil(1,{(NonLocal[]){{"p","scalar",(void *)&p,NULL,0},{"_const_alpha","_coord",(void *)&_const_alpha,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 498 \"/home/spencer/basilisk/src/poisson.h\"\n{is_face_x(){\n    val_out_(uf.x,0,0,0) -= dt*_const_alpha.x*((val(p,0,0,0) - val(p,0 -1,0,0))/Delta);}end_is_face_x()\n// #line 498\nis_face_y(){\n    val_out_(uf.y,0,0,0) -= dt*_const_alpha.y*((val(p,0,0,0) - val(p,0,0 -1,0))/Delta);}end_is_face_y()}")}){_stencil_is_face_x(){
    {;_stencil_val(p,0,0,0); _stencil_val(p,0 -1,0,0);_stencil_val_r(uf.x,0,0,0);   }}end__stencil_is_face_x()
#line 498
_stencil_is_face_y(){
    {;_stencil_val(p,0,0,0); _stencil_val(p,0,0 -1,0);_stencil_val_r(uf.y,0,0,0);   }}end__stencil_is_face_y()}end_foreach_face_stencil();




  {
#line 498
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) -= dt*_const_alpha.x*((val(p,0,0,0) - val(p,0 -1,0,0))/Delta);}end_is_face_x()
#line 498
is_face_y(){
    val(uf.y,0,0,0) -= dt*_const_alpha.y*((val(p,0,0,0) - val(p,0,0 -1,0))/Delta);}end_is_face_y()}end_foreach_face_generic();}}

  {delete((scalar*)((scalar[]){div,{-1}}));{end_tracing("project","/home/spencer/basilisk/src/poisson.h",0);return mgp;}}delete((scalar*)((scalar[]){div,{-1}}));
end_tracing("project","/home/spencer/basilisk/src/poisson.h",0);}
# 51 "/home/spencer/basilisk/src/viscosity.h" 2

struct Viscosity {
  vector mu;
  scalar rho;
  double dt;
};
# 135 "/home/spencer/basilisk/src/viscosity.h"
static void relax_viscosity (scalar * a, scalar * b, int l, void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
          vector mu = p->mu;
          scalar rho = p->rho;
  double dt = p->dt;
  vector u = (*((vector *)&(a[0]))), r = (*((vector *)&(b[0])));




  vector w = u;


  if(!is_constant(rho) && !is_constant(mu.x)){{foreach_level_or_leaf (l) {
    
      val(w.x,0,0,0) = (dt/val(rho,0,0,0)*(2.*val(mu.x,1,0,0)*val(u.x,1,0,0) + 2.*val(mu.x,0,0,0)*val(u.x,-1,0,0)

      + val(mu.y,0,1,0)*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - val(mu.y,0,0,0)*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)
# 168 "/home/spencer/basilisk/src/viscosity.h"
      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.}).x + dt/val(rho,0,0,0)*(2.*val(mu.x,1,0,0) + 2.*val(mu.x,0,0,0)

          + val(mu.y,0,1,0) + val(mu.y,0,0,0)




        ));
      
#line 151
val(w.y,0,0,0) = (dt/val(rho,0,0,0)*(2.*val(mu.y,0,1,0)*val(u.y,0,1,0) + 2.*val(mu.y,0,0,0)*val(u.y,0,-1,0)

      + val(mu.x,1,0,0)*(val(u.y,1,0,0) +
     (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
     (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - val(mu.x,0,0,0)*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)
# 168 "/home/spencer/basilisk/src/viscosity.h"
      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.}).y + dt/val(rho,0,0,0)*(2.*val(mu.y,0,1,0) + 2.*val(mu.y,0,0,0)

          + val(mu.x,1,0,0) + val(mu.x,0,0,0)




        ));
  }end_foreach_level_or_leaf();}}else if(is_constant(rho) && !is_constant(mu.x)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);


  {
#line 149
foreach_level_or_leaf (l) {
    
      val(w.x,0,0,0) = (dt/_const_rho*(2.*val(mu.x,1,0,0)*val(u.x,1,0,0) + 2.*val(mu.x,0,0,0)*val(u.x,-1,0,0)

      + val(mu.y,0,1,0)*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - val(mu.y,0,0,0)*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)
# 168 "/home/spencer/basilisk/src/viscosity.h"
      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.}).x + dt/_const_rho*(2.*val(mu.x,1,0,0) + 2.*val(mu.x,0,0,0)

          + val(mu.y,0,1,0) + val(mu.y,0,0,0)




        ));
      
#line 151
val(w.y,0,0,0) = (dt/_const_rho*(2.*val(mu.y,0,1,0)*val(u.y,0,1,0) + 2.*val(mu.y,0,0,0)*val(u.y,0,-1,0)

      + val(mu.x,1,0,0)*(val(u.y,1,0,0) +
     (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
     (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - val(mu.x,0,0,0)*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)
# 168 "/home/spencer/basilisk/src/viscosity.h"
      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.}).y + dt/_const_rho*(2.*val(mu.y,0,1,0) + 2.*val(mu.y,0,0,0)

          + val(mu.x,1,0,0) + val(mu.x,0,0,0)




        ));
  }end_foreach_level_or_leaf();}}else if(!is_constant(rho) && is_constant(mu.x)){_coord _const_mu={_constant[mu.x.i-_NVARMAX],_constant[mu.y.i-_NVARMAX]};NOT_UNUSED(_const_mu);


  {
#line 149
foreach_level_or_leaf (l) {
    
      val(w.x,0,0,0) = (dt/val(rho,0,0,0)*(2.*_const_mu.x*val(u.x,1,0,0) + 2.*_const_mu.x*val(u.x,-1,0,0)

      + _const_mu.y*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - _const_mu.y*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)
# 168 "/home/spencer/basilisk/src/viscosity.h"
      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.}).x + dt/val(rho,0,0,0)*(2.*_const_mu.x + 2.*_const_mu.x

          + _const_mu.y + _const_mu.y




        ));
      
#line 151
val(w.y,0,0,0) = (dt/val(rho,0,0,0)*(2.*_const_mu.y*val(u.y,0,1,0) + 2.*_const_mu.y*val(u.y,0,-1,0)

      + _const_mu.x*(val(u.y,1,0,0) +
     (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
     (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - _const_mu.x*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)
# 168 "/home/spencer/basilisk/src/viscosity.h"
      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.}).y + dt/val(rho,0,0,0)*(2.*_const_mu.y + 2.*_const_mu.y

          + _const_mu.x + _const_mu.x




        ));
  }end_foreach_level_or_leaf();}}else {double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);_coord _const_mu={_constant[mu.x.i-_NVARMAX],_constant[mu.y.i-_NVARMAX]};NOT_UNUSED(_const_mu);


  {
#line 149
foreach_level_or_leaf (l) {
    
      val(w.x,0,0,0) = (dt/_const_rho*(2.*_const_mu.x*val(u.x,1,0,0) + 2.*_const_mu.x*val(u.x,-1,0,0)

      + _const_mu.y*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - _const_mu.y*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)
# 168 "/home/spencer/basilisk/src/viscosity.h"
      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.}).x + dt/_const_rho*(2.*_const_mu.x + 2.*_const_mu.x

          + _const_mu.y + _const_mu.y




        ));
      
#line 151
val(w.y,0,0,0) = (dt/_const_rho*(2.*_const_mu.y*val(u.y,0,1,0) + 2.*_const_mu.y*val(u.y,0,-1,0)

      + _const_mu.x*(val(u.y,1,0,0) +
     (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
     (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - _const_mu.x*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)
# 168 "/home/spencer/basilisk/src/viscosity.h"
      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.}).y + dt/_const_rho*(2.*_const_mu.y + 2.*_const_mu.y

          + _const_mu.x + _const_mu.x




        ));
  }end_foreach_level_or_leaf();}}
# 195 "/home/spencer/basilisk/src/viscosity.h"
}
# 204 "/home/spencer/basilisk/src/viscosity.h"
static double residual_viscosity (scalar * a, scalar * b, scalar * resl,
      void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
          vector mu = p->mu;
          scalar rho = p->rho;
  double dt = p->dt;
  vector u = (*((vector *)&(a[0]))), r = (*((vector *)&(b[0]))), res = (*((vector *)&(resl[0])));
  double maxres = 0.;
# 221 "/home/spencer/basilisk/src/viscosity.h"
  boundary_internal ((scalar *)((vector[]){u,{{-1},{-1}}}), "/home/spencer/basilisk/src/viscosity.h", 0);

   {
    vector  taux=new_face_vector("taux");
    if(!is_constant(mu.x)){
    
#line 225
foreach_face_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"mu","vector",(void *)&mu,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 225 \"/home/spencer/basilisk/src/viscosity.h\"\nis_face_x(){\n      val_out_(taux.x,0,0,0) = 2.*val(mu.x,0,0,0)*(val(u.x,0,0,0) - val(u.x,-1,0,0))/Delta;}end_is_face_x()")})_stencil_is_face_x(){
      {_stencil_val(mu.x,0,0,0);_stencil_val(u.x,0,0,0); _stencil_val(u.x,-1,0,0);_stencil_val_a(taux.x,0,0,0);   }}end__stencil_is_face_x()end_foreach_face_stencil();{
#line 225
foreach_face_generic()is_face_x(){
      val(taux.x,0,0,0) = 2.*val(mu.x,0,0,0)*(val(u.x,0,0,0) - val(u.x,-1,0,0))/Delta;}end_is_face_x()end_foreach_face_generic();}}else {_coord _const_mu={_constant[mu.x.i-_NVARMAX],_constant[mu.y.i-_NVARMAX]};NOT_UNUSED(_const_mu);
    
#line 225
foreach_face_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"_const_mu","_coord",(void *)&_const_mu,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 225 \"/home/spencer/basilisk/src/viscosity.h\"\nis_face_x(){\n      val_out_(taux.x,0,0,0) = 2.*_const_mu.x*(val(u.x,0,0,0) - val(u.x,-1,0,0))/Delta;}end_is_face_x()")})_stencil_is_face_x(){
      {;_stencil_val(u.x,0,0,0); _stencil_val(u.x,-1,0,0);_stencil_val_a(taux.x,0,0,0);   }}end__stencil_is_face_x()end_foreach_face_stencil();
    {
#line 225
foreach_face_generic()is_face_x(){
      val(taux.x,0,0,0) = 2.*_const_mu.x*(val(u.x,0,0,0) - val(u.x,-1,0,0))/Delta;}end_is_face_x()end_foreach_face_generic();}}

      if(!is_constant(mu.x)){

      
#line 228
foreach_face_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"mu","vector",(void *)&mu,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 228 \"/home/spencer/basilisk/src/viscosity.h\"\nis_face_y(){\n val_out_(taux.y,0,0,0) = val(mu.y,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,-1,0) +\n      (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -\n      (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)/Delta;}end_is_face_y()")})_stencil_is_face_y(){
 { _stencil_val(mu.y,0,0,0);_stencil_val(u.x,0,0,0); _stencil_val(u.x,0,-1,0);
_stencil_val(u.y,1,-1,0); _stencil_val(u.y,1,0,0);
_stencil_val(u.y,-1,-1,0); _stencil_val(u.y,-1,0,0);
#line 229
_stencil_val_a(taux.y,0,0,0);    
        
      }}end__stencil_is_face_y()end_foreach_face_stencil();{
#line 228
foreach_face_generic()is_face_y(){
 val(taux.y,0,0,0) = val(mu.y,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,-1,0) +
      (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
      (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)/Delta;}end_is_face_y()end_foreach_face_generic();}}else {_coord _const_mu={_constant[mu.x.i-_NVARMAX],_constant[mu.y.i-_NVARMAX]};NOT_UNUSED(_const_mu);

      
#line 228
foreach_face_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"_const_mu","_coord",(void *)&_const_mu,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 228 \"/home/spencer/basilisk/src/viscosity.h\"\nis_face_y(){\n val_out_(taux.y,0,0,0) = _const_mu.y*(val(u.x,0,0,0) - val(u.x,0,-1,0) +\n      (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -\n      (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)/Delta;}end_is_face_y()")})_stencil_is_face_y(){
 {;_stencil_val(u.x,0,0,0); _stencil_val(u.x,0,-1,0);
_stencil_val(u.y,1,-1,0); _stencil_val(u.y,1,0,0);
_stencil_val(u.y,-1,-1,0); _stencil_val(u.y,-1,0,0);
#line 229
_stencil_val_a(taux.y,0,0,0);    
        
      }}end__stencil_is_face_y()end_foreach_face_stencil();

      {
#line 228
foreach_face_generic()is_face_y(){
 val(taux.y,0,0,0) = _const_mu.y*(val(u.x,0,0,0) - val(u.x,0,-1,0) +
      (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
      (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)/Delta;}end_is_face_y()end_foreach_face_generic();}}







    if(!is_constant(rho)){







    
#line 239
foreach_stencil (1,{(NonLocal[]){{"maxres","double",(void *)&maxres,NULL,0,'M'},{"rho","scalar",(void *)&rho,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"u","vector",(void *)&u,NULL,0},{"r","vector",(void *)&r,NULL,0},{"res","vector",(void *)&res,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_(" \n// #line 239 \"/home/spencer/basilisk/src/viscosity.h\"\n{\n      real d = 0.;\n      \n d += val(taux.x,1,0,0) - val(taux.x,0,0,0);\n \n// #line 242\nd += val(taux.y,0,1,0) - val(taux.y,0,0,0);\n      val_out_(res.x,0,0,0) = val(r.x,0,0,0) - ((coord){1.,1.}).x*val(u.x,0,0,0) + dt/val(rho,0,0,0)*d/Delta;\n      if (fabs (val(res.x,0,0,0)) > maxres)\n maxres = fabs (val(res.x,0,0,0));\n    }")}) {   
      
      
 { _stencil_val(taux.x,1,0,0); _stencil_val(taux.x,0,0,0);  }
 
#line 242
{ _stencil_val(taux.y,0,1,0); _stencil_val(taux.y,0,0,0);  } 
_stencil_val(r.x,0,0,0);_stencil_val(u.x,0,0,0);_stencil_val(rho,0,0,0);
      
#line 243
_stencil_val_a(res.x,0,0,0);
_stencil_val(res.x,0,0,0);
 {_stencil_val(res.x,0,0,0);   }     
          
    
#line 246
}end_foreach_stencil();
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 239
foreach () {
      double d = 0.;
      
 d += val(taux.x,1,0,0) - val(taux.x,0,0,0);
 
#line 242
d += val(taux.y,0,1,0) - val(taux.y,0,0,0);
      val(res.x,0,0,0) = val(r.x,0,0,0) - ((coord){1.,1.}).x*val(u.x,0,0,0) + dt/val(rho,0,0,0)*d/Delta;
      if (fabs (val(res.x,0,0,0)) > maxres)
 maxres = fabs (val(res.x,0,0,0));
    }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 246
}else {double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);







    
#line 239
foreach_stencil (1,{(NonLocal[]){{"maxres","double",(void *)&maxres,NULL,0,'M'},{"_const_rho","double",(void *)&_const_rho,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"u","vector",(void *)&u,NULL,0},{"r","vector",(void *)&r,NULL,0},{"res","vector",(void *)&res,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_(" \n// #line 239 \"/home/spencer/basilisk/src/viscosity.h\"\n{\n      real d = 0.;\n      \n d += val(taux.x,1,0,0) - val(taux.x,0,0,0);\n \n// #line 242\nd += val(taux.y,0,1,0) - val(taux.y,0,0,0);\n      val_out_(res.x,0,0,0) = val(r.x,0,0,0) - ((coord){1.,1.}).x*val(u.x,0,0,0) + dt/_const_rho*d/Delta;\n      if (fabs (val(res.x,0,0,0)) > maxres)\n maxres = fabs (val(res.x,0,0,0));\n    }")}) {   
      
      
 { _stencil_val(taux.x,1,0,0); _stencil_val(taux.x,0,0,0);  }
 
#line 242
{ _stencil_val(taux.y,0,1,0); _stencil_val(taux.y,0,0,0);  } 
_stencil_val(r.x,0,0,0);_stencil_val(u.x,0,0,0);;
      
#line 243
_stencil_val_a(res.x,0,0,0);
_stencil_val(res.x,0,0,0);
 {_stencil_val(res.x,0,0,0);   }     
          
    
#line 246
}end_foreach_stencil();







    
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 239
foreach () {
      double d = 0.;
      
 d += val(taux.x,1,0,0) - val(taux.x,0,0,0);
 
#line 242
d += val(taux.y,0,1,0) - val(taux.y,0,0,0);
      val(res.x,0,0,0) = val(r.x,0,0,0) - ((coord){1.,1.}).x*val(u.x,0,0,0) + dt/_const_rho*d/Delta;
      if (fabs (val(res.x,0,0,0)) > maxres)
 maxres = fabs (val(res.x,0,0,0));
    }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 246
}delete((scalar*)((vector[]){taux,{{-1},{-1}}}));
  } 
#line 223
{
    vector  taux=new_face_vector("taux");
    if(!is_constant(mu.y)){
    
#line 225
foreach_face_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"mu","vector",(void *)&mu,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 225 \"/home/spencer/basilisk/src/viscosity.h\"\nis_face_y(){\n      val_out_(taux.y,0,0,0) = 2.*val(mu.y,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,-1,0))/Delta;}end_is_face_y()")})_stencil_is_face_y(){
      {_stencil_val(mu.y,0,0,0);_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,-1,0);_stencil_val_a(taux.y,0,0,0);   }}end__stencil_is_face_y()end_foreach_face_stencil();{
#line 225
foreach_face_generic()is_face_y(){
      val(taux.y,0,0,0) = 2.*val(mu.y,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,-1,0))/Delta;}end_is_face_y()end_foreach_face_generic();}}else {_coord _const_mu={_constant[mu.y.i-_NVARMAX],_constant[mu.x.i-_NVARMAX]};NOT_UNUSED(_const_mu);
    
#line 225
foreach_face_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"_const_mu","_coord",(void *)&_const_mu,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 225 \"/home/spencer/basilisk/src/viscosity.h\"\nis_face_y(){\n      val_out_(taux.y,0,0,0) = 2.*_const_mu.y*(val(u.y,0,0,0) - val(u.y,0,-1,0))/Delta;}end_is_face_y()")})_stencil_is_face_y(){
      {;_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,-1,0);_stencil_val_a(taux.y,0,0,0);   }}end__stencil_is_face_y()end_foreach_face_stencil();
    {
#line 225
foreach_face_generic()is_face_y(){
      val(taux.y,0,0,0) = 2.*_const_mu.y*(val(u.y,0,0,0) - val(u.y,0,-1,0))/Delta;}end_is_face_y()end_foreach_face_generic();}}

      if(!is_constant(mu.y)){

      
#line 228
foreach_face_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"mu","vector",(void *)&mu,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 228 \"/home/spencer/basilisk/src/viscosity.h\"\nis_face_x(){\n val_out_(taux.x,0,0,0) = val(mu.x,0,0,0)*(val(u.y,0,0,0) - val(u.y,-1,0,0) +\n      (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -\n      (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)/Delta;}end_is_face_x()")})_stencil_is_face_x(){
 { _stencil_val(mu.x,0,0,0);_stencil_val(u.y,0,0,0); _stencil_val(u.y,-1,0,0);
_stencil_val(u.x,-1,1,0); _stencil_val(u.x,0,1,0);
_stencil_val(u.x,-1,-1,0); _stencil_val(u.x,0,-1,0);
#line 229
_stencil_val_a(taux.x,0,0,0);    
        
      }}end__stencil_is_face_x()end_foreach_face_stencil();{
#line 228
foreach_face_generic()is_face_x(){
 val(taux.x,0,0,0) = val(mu.x,0,0,0)*(val(u.y,0,0,0) - val(u.y,-1,0,0) +
      (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
      (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)/Delta;}end_is_face_x()end_foreach_face_generic();}}else {_coord _const_mu={_constant[mu.y.i-_NVARMAX],_constant[mu.x.i-_NVARMAX]};NOT_UNUSED(_const_mu);

      
#line 228
foreach_face_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"_const_mu","_coord",(void *)&_const_mu,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 228 \"/home/spencer/basilisk/src/viscosity.h\"\nis_face_x(){\n val_out_(taux.x,0,0,0) = _const_mu.x*(val(u.y,0,0,0) - val(u.y,-1,0,0) +\n      (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -\n      (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)/Delta;}end_is_face_x()")})_stencil_is_face_x(){
 {;_stencil_val(u.y,0,0,0); _stencil_val(u.y,-1,0,0);
_stencil_val(u.x,-1,1,0); _stencil_val(u.x,0,1,0);
_stencil_val(u.x,-1,-1,0); _stencil_val(u.x,0,-1,0);
#line 229
_stencil_val_a(taux.x,0,0,0);    
        
      }}end__stencil_is_face_x()end_foreach_face_stencil();

      {
#line 228
foreach_face_generic()is_face_x(){
 val(taux.x,0,0,0) = _const_mu.x*(val(u.y,0,0,0) - val(u.y,-1,0,0) +
      (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
      (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)/Delta;}end_is_face_x()end_foreach_face_generic();}}







    if(!is_constant(rho)){







    
#line 239
foreach_stencil (1,{(NonLocal[]){{"maxres","double",(void *)&maxres,NULL,0,'M'},{"rho","scalar",(void *)&rho,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"u","vector",(void *)&u,NULL,0},{"r","vector",(void *)&r,NULL,0},{"res","vector",(void *)&res,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_(" \n// #line 239 \"/home/spencer/basilisk/src/viscosity.h\"\n{\n      real d = 0.;\n      \n d += val(taux.y,0,1,0) - val(taux.y,0,0,0);\n \n// #line 242\nd += val(taux.x,1,0,0) - val(taux.x,0,0,0);\n      val_out_(res.y,0,0,0) = val(r.y,0,0,0) - ((coord){1.,1.}).y*val(u.y,0,0,0) + dt/val(rho,0,0,0)*d/Delta;\n      if (fabs (val(res.y,0,0,0)) > maxres)\n maxres = fabs (val(res.y,0,0,0));\n    }")}) {   
      
      
 { _stencil_val(taux.y,0,1,0); _stencil_val(taux.y,0,0,0);  }
 
#line 242
{ _stencil_val(taux.x,1,0,0); _stencil_val(taux.x,0,0,0);  } 
_stencil_val(r.y,0,0,0);_stencil_val(u.y,0,0,0);_stencil_val(rho,0,0,0);
      
#line 243
_stencil_val_a(res.y,0,0,0);
_stencil_val(res.y,0,0,0);
 {_stencil_val(res.y,0,0,0);   }     
          
    
#line 246
}end_foreach_stencil();
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 239
foreach () {
      double d = 0.;
      
 d += val(taux.y,0,1,0) - val(taux.y,0,0,0);
 
#line 242
d += val(taux.x,1,0,0) - val(taux.x,0,0,0);
      val(res.y,0,0,0) = val(r.y,0,0,0) - ((coord){1.,1.}).y*val(u.y,0,0,0) + dt/val(rho,0,0,0)*d/Delta;
      if (fabs (val(res.y,0,0,0)) > maxres)
 maxres = fabs (val(res.y,0,0,0));
    }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 246
}else {double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);







    
#line 239
foreach_stencil (1,{(NonLocal[]){{"maxres","double",(void *)&maxres,NULL,0,'M'},{"_const_rho","double",(void *)&_const_rho,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"u","vector",(void *)&u,NULL,0},{"r","vector",(void *)&r,NULL,0},{"res","vector",(void *)&res,NULL,0},{"taux","vector",(void *)&taux,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_(" \n// #line 239 \"/home/spencer/basilisk/src/viscosity.h\"\n{\n      real d = 0.;\n      \n d += val(taux.y,0,1,0) - val(taux.y,0,0,0);\n \n// #line 242\nd += val(taux.x,1,0,0) - val(taux.x,0,0,0);\n      val_out_(res.y,0,0,0) = val(r.y,0,0,0) - ((coord){1.,1.}).y*val(u.y,0,0,0) + dt/_const_rho*d/Delta;\n      if (fabs (val(res.y,0,0,0)) > maxres)\n maxres = fabs (val(res.y,0,0,0));\n    }")}) {   
      
      
 { _stencil_val(taux.y,0,1,0); _stencil_val(taux.y,0,0,0);  }
 
#line 242
{ _stencil_val(taux.x,1,0,0); _stencil_val(taux.x,0,0,0);  } 
_stencil_val(r.y,0,0,0);_stencil_val(u.y,0,0,0);;
      
#line 243
_stencil_val_a(res.y,0,0,0);
_stencil_val(res.y,0,0,0);
 {_stencil_val(res.y,0,0,0);   }     
          
    
#line 246
}end_foreach_stencil();







    
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 239
foreach () {
      double d = 0.;
      
 d += val(taux.y,0,1,0) - val(taux.y,0,0,0);
 
#line 242
d += val(taux.x,1,0,0) - val(taux.x,0,0,0);
      val(res.y,0,0,0) = val(r.y,0,0,0) - ((coord){1.,1.}).y*val(u.y,0,0,0) + dt/_const_rho*d/Delta;
      if (fabs (val(res.y,0,0,0)) > maxres)
 maxres = fabs (val(res.y,0,0,0));
    }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 246
}delete((scalar*)((vector[]){taux,{{-1},{-1}}}));
  }
# 276 "/home/spencer/basilisk/src/viscosity.h"
  return maxres;
}
# 288 "/home/spencer/basilisk/src/viscosity.h"
     
mgstats viscosity (vector u, vector mu, scalar rho, double dt,
     int nrelax, scalar * res)
{tracing("viscosity","/home/spencer/basilisk/src/viscosity.h",0);





  vector  r=new_vector("r");
  foreach_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"r","vector",(void *)&r,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n    \n// #line 299 \"/home/spencer/basilisk/src/viscosity.h\"\n{\n      val_out_(r.x,0,0,0) = val(u.x,0,0,0);\n      \n// #line 300\nval_out_(r.y,0,0,0) = val(u.y,0,0,0);}")})
    {
      { _stencil_val(u.x,0,0,0);_stencil_val_a(r.x,0,0,0); }
      
#line 300
{ _stencil_val(u.y,0,0,0);_stencil_val_a(r.y,0,0,0); }}end_foreach_stencil();
  {
#line 298
foreach()
    {
      val(r.x,0,0,0) = val(u.x,0,0,0);
      
#line 300
val(r.y,0,0,0) = val(u.y,0,0,0);}end_foreach();}




  restriction (((scalar[]){mu.x,mu.y,rho,{-1}}));
  struct Viscosity p = { mu, rho, dt };
  { mgstats _ret= mg_solve ((scalar *)((vector[]){u,{{-1},{-1}}}), (scalar *)((vector[]){r,{{-1},{-1}}})
,
     
#line 308
residual_viscosity, relax_viscosity, &p, nrelax, res
#line 132 "/home/spencer/basilisk/src/poisson.h"
, 
0, 
TOLERANCE
#line 308 "/home/spencer/basilisk/src/viscosity.h"
);delete((scalar*)((vector[]){r,{{-1},{-1}}}));{end_tracing("viscosity","/home/spencer/basilisk/src/viscosity.h",0);return _ret;}}delete((scalar*)((vector[]){r,{{-1},{-1}}}));
end_tracing("viscosity","/home/spencer/basilisk/src/viscosity.h",0);}
# 318 "/home/spencer/basilisk/src/viscosity.h"
     
mgstats viscosity_explicit (vector u, vector mu, scalar rho, double dt)
{tracing("viscosity_explicit","/home/spencer/basilisk/src/viscosity.h",0);
  vector  r=new_vector("r");
  mgstats mg = {0};
  struct Viscosity p = { mu, rho, dt };
  mg.resb = residual_viscosity ((scalar *)((vector[]){u,{{-1},{-1}}}), (scalar *)((vector[]){u,{{-1},{-1}}}), (scalar *)((vector[]){r,{{-1},{-1}}}), &p);
  foreach_stencil(1,{(NonLocal[]){{"r","vector",(void *)&r,NULL,0},{"u","vector",(void *)&u,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n    \n// #line 326 \"/home/spencer/basilisk/src/viscosity.h\"\n{\n      val_out_(u.x,0,0,0) += val(r.x,0,0,0);\n      \n// #line 327\nval_out_(u.y,0,0,0) += val(r.y,0,0,0);}")})
    {
      { _stencil_val(r.x,0,0,0);_stencil_val_r(u.x,0,0,0); }
      
#line 327
{ _stencil_val(r.y,0,0,0);_stencil_val_r(u.y,0,0,0); }}end_foreach_stencil();
  {
#line 325
foreach()
    {
      val(u.x,0,0,0) += val(r.x,0,0,0);
      
#line 327
val(u.y,0,0,0) += val(r.y,0,0,0);}end_foreach();}
  {delete((scalar*)((vector[]){r,{{-1},{-1}}}));{end_tracing("viscosity_explicit","/home/spencer/basilisk/src/viscosity.h",0);return mg;}}delete((scalar*)((vector[]){r,{{-1},{-1}}}));
end_tracing("viscosity_explicit","/home/spencer/basilisk/src/viscosity.h",0);}
# 34 "/home/spencer/basilisk/src/navier-stokes/centered.h" 2
# 44 "/home/spencer/basilisk/src/navier-stokes/centered.h"
scalar  p={0};
vector  u={{1},{2}},  g={{3},{4}};
scalar  pf={5};
vector  uf={{6},{7}};
# 70 "/home/spencer/basilisk/src/navier-stokes/centered.h"
        vector mu = {{_NVARMAX+0},{_NVARMAX+1}}, a = {{_NVARMAX+0},{_NVARMAX+1}}, alpha = {{_NVARMAX+2},{_NVARMAX+3}};
        scalar rho = {_NVARMAX+4};
mgstats mgp = {0}, mgpf = {0}, mgu = {0};
bool stokes = false;
#line 91
static double _boundary0(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(a.x,1,0,0)*val(fm.x,1,0,0)/val(alpha.x,1,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_a.x*val(fm.x,1,0,0)/val(alpha.x,1,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(a.x,1,0,0)*_const_fm.x/val(alpha.x,1,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_a.x*_const_fm.x/val(alpha.x,1,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(a.x,1,0,0)*val(fm.x,1,0,0)/_const_alpha.x), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_a.x*val(fm.x,1,0,0)/_const_alpha.x), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(a.x,1,0,0)*_const_fm.x/_const_alpha.x), point, neighbor, _s, data));}}}}else {_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_a.x*_const_fm.x/_const_alpha.x), point, neighbor, _s, data));}}}}}static double _boundary0_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(a.x,1,0,0)*val(fm.x,1,0,0)/val(alpha.x,1,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_a.x*val(fm.x,1,0,0)/val(alpha.x,1,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(a.x,1,0,0)*_const_fm.x/val(alpha.x,1,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_a.x*_const_fm.x/val(alpha.x,1,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(a.x,1,0,0)*val(fm.x,1,0,0)/_const_alpha.x), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_a.x*val(fm.x,1,0,0)/_const_alpha.x), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(a.x,1,0,0)*_const_fm.x/_const_alpha.x), point, neighbor, _s, data));}}}}else {_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_a.x*_const_fm.x/_const_alpha.x), point, neighbor, _s, data));}}}}}
static double _boundary1(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(a.x,0,0,0)*val(fm.x,0,0,0)/val(alpha.x,0,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_a.x*val(fm.x,0,0,0)/val(alpha.x,0,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(a.x,0,0,0)*_const_fm.x/val(alpha.x,0,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_a.x*_const_fm.x/val(alpha.x,0,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(a.x,0,0,0)*val(fm.x,0,0,0)/_const_alpha.x), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_a.x*val(fm.x,0,0,0)/_const_alpha.x), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(a.x,0,0,0)*_const_fm.x/_const_alpha.x), point, neighbor, _s, data));}}}}else {_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_a.x*_const_fm.x/_const_alpha.x), point, neighbor, _s, data));}}}}}static double _boundary1_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(a.x,0,0,0)*val(fm.x,0,0,0)/val(alpha.x,0,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_a.x*val(fm.x,0,0,0)/val(alpha.x,0,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(a.x,0,0,0)*_const_fm.x/val(alpha.x,0,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_a.x*_const_fm.x/val(alpha.x,0,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(a.x,0,0,0)*val(fm.x,0,0,0)/_const_alpha.x), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_a.x*val(fm.x,0,0,0)/_const_alpha.x), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(a.x,0,0,0)*_const_fm.x/_const_alpha.x), point, neighbor, _s, data));}}}}else {_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_a.x*_const_fm.x/_const_alpha.x), point, neighbor, _s, data));}}}}}








static double _boundary2(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(a.y,0,1,0)*val(fm.y,0,1,0)/val(alpha.y,0,1,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_a.y*val(fm.y,0,1,0)/val(alpha.y,0,1,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(a.y,0,1,0)*_const_fm.y/val(alpha.y,0,1,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_a.y*_const_fm.y/val(alpha.y,0,1,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(a.y,0,1,0)*val(fm.y,0,1,0)/_const_alpha.y), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_a.y*val(fm.y,0,1,0)/_const_alpha.y), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((val(a.y,0,1,0)*_const_fm.y/_const_alpha.y), point, neighbor, _s, data));}}}}else {_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann((_const_a.y*_const_fm.y/_const_alpha.y), point, neighbor, _s, data));}}}}}static double _boundary2_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(a.y,0,1,0)*val(fm.y,0,1,0)/val(alpha.y,0,1,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_a.y*val(fm.y,0,1,0)/val(alpha.y,0,1,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(a.y,0,1,0)*_const_fm.y/val(alpha.y,0,1,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_a.y*_const_fm.y/val(alpha.y,0,1,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(a.y,0,1,0)*val(fm.y,0,1,0)/_const_alpha.y), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_a.y*val(fm.y,0,1,0)/_const_alpha.y), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((val(a.y,0,1,0)*_const_fm.y/_const_alpha.y), point, neighbor, _s, data));}}}}else {_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_a.y*_const_fm.y/_const_alpha.y), point, neighbor, _s, data));}}}}}
static double _boundary3(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(a.y,0,0,0)*val(fm.y,0,0,0)/val(alpha.y,0,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_a.y*val(fm.y,0,0,0)/val(alpha.y,0,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(a.y,0,0,0)*_const_fm.y/val(alpha.y,0,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_a.y*_const_fm.y/val(alpha.y,0,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(a.y,0,0,0)*val(fm.y,0,0,0)/_const_alpha.y), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_a.y*val(fm.y,0,0,0)/_const_alpha.y), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (val(a.y,0,0,0)*_const_fm.y/_const_alpha.y), point, neighbor, _s, data));}}}}else {_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(- (_const_a.y*_const_fm.y/_const_alpha.y), point, neighbor, _s, data));}}}}}static double _boundary3_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(a.y,0,0,0)*val(fm.y,0,0,0)/val(alpha.y,0,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_a.y*val(fm.y,0,0,0)/val(alpha.y,0,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(a.y,0,0,0)*_const_fm.y/val(alpha.y,0,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_a.y*_const_fm.y/val(alpha.y,0,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(a.y,0,0,0)*val(fm.y,0,0,0)/_const_alpha.y), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_a.y*val(fm.y,0,0,0)/_const_alpha.y), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(a.y,0,0,0)*_const_fm.y/_const_alpha.y), point, neighbor, _s, data));}}}}else {_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_a.y*_const_fm.y/_const_alpha.y), point, neighbor, _s, data));}}}}}
#line 126
static int defaults_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0)!=0;*ip=i;*tp=t;return ret;}
# 126 "/home/spencer/basilisk/src/navier-stokes/centered.h"
      static int defaults_0(const int i,const double t,Event *_ev){tracing("defaults_0","/home/spencer/basilisk/src/navier-stokes/centered.h",0);
{

  CFL = 0.8;




  _attribute[p.i].nodump = _attribute[pf.i].nodump = true;





  if (alpha.x.i == unityf.x.i) {
    alpha = fm;
    rho = cm;
  }
  else if (!is_constant(alpha.x)) {
    vector alphav = alpha;
    if(!is_constant(fm.x)){
    
#line 146
foreach_face_stencil(1,{(NonLocal[]){{"fm","vector",(void *)&fm,NULL,0},{"alphav","vector",(void *)&alphav,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 146 \"/home/spencer/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n      val_out_(alphav.x,0,0,0) = val(fm.x,0,0,0);}end_is_face_x()\n// #line 146\nis_face_y(){\n      val_out_(alphav.y,0,0,0) = val(fm.y,0,0,0);}end_is_face_y()}")}){_stencil_is_face_x(){
      { _stencil_val(fm.x,0,0,0);_stencil_val_a(alphav.x,0,0,0); }}end__stencil_is_face_x()
#line 146
_stencil_is_face_y(){
      { _stencil_val(fm.y,0,0,0);_stencil_val_a(alphav.y,0,0,0); }}end__stencil_is_face_y()}end_foreach_face_stencil();{
#line 146
foreach_face_generic(){is_face_x(){
      val(alphav.x,0,0,0) = val(fm.x,0,0,0);}end_is_face_x()
#line 146
is_face_y(){
      val(alphav.y,0,0,0) = val(fm.y,0,0,0);}end_is_face_y()}end_foreach_face_generic();}}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
    
#line 146
foreach_face_stencil(1,{(NonLocal[]){{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"alphav","vector",(void *)&alphav,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 146 \"/home/spencer/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n      val_out_(alphav.x,0,0,0) = _const_fm.x;}end_is_face_x()\n// #line 146\nis_face_y(){\n      val_out_(alphav.y,0,0,0) = _const_fm.y;}end_is_face_y()}")}){_stencil_is_face_x(){
      {;_stencil_val_a(alphav.x,0,0,0); }}end__stencil_is_face_x()
#line 146
_stencil_is_face_y(){
      {;_stencil_val_a(alphav.y,0,0,0); }}end__stencil_is_face_y()}end_foreach_face_stencil();
    {
#line 146
foreach_face_generic(){is_face_x(){
      val(alphav.x,0,0,0) = _const_fm.x;}end_is_face_x()
#line 146
is_face_y(){
      val(alphav.y,0,0,0) = _const_fm.y;}end_is_face_y()}end_foreach_face_generic();}}
  }






  _attribute[uf.x.i].refine = refine_face_solenoidal;
# 178 "/home/spencer/basilisk/src/navier-stokes/centered.h"
  foreach_stencil(1,{(NonLocal[]){{"t","double",(void *)&t,NULL,0},{"u","vector",(void *)&u,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n    \n// #line 179 \"/home/spencer/basilisk/src/navier-stokes/centered.h\"\n{\n      dimensional (val(u.x,0,0,0) == Delta/t);\n      \n// #line 180\ndimensional (val(u.y,0,0,0) == Delta/t);}")})
    {
      {_stencil_val(u.x,0,0,0);   }
      
#line 180
{_stencil_val(u.y,0,0,0);   }}end_foreach_stencil();
# 178 "/home/spencer/basilisk/src/navier-stokes/centered.h"
  {foreach()
    {
      dimensional (val(u.x,0,0,0) == Delta/t);
      
#line 180
dimensional (val(u.y,0,0,0) == Delta/t);}end_foreach();}
}{end_tracing("defaults_0","/home/spencer/basilisk/src/navier-stokes/centered.h",0);return 0;}end_tracing("defaults_0","/home/spencer/basilisk/src/navier-stokes/centered.h",0);}





static int default_display_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0)!=0;*ip=i;*tp=t;return ret;}






#line 187
      static int default_display(const int i,const double t,Event *_ev){tracing("default_display","/home/spencer/basilisk/src/navier-stokes/centered.h",0);
  display ("squares (color = 'u.x', spread = -1);"
#line 1425 "/home/spencer/basilisk/src/common.h"
, false
#line 188 "/home/spencer/basilisk/src/navier-stokes/centered.h"
);{end_tracing("default_display","/home/spencer/basilisk/src/navier-stokes/centered.h",0);return 0;}end_tracing("default_display","/home/spencer/basilisk/src/navier-stokes/centered.h",0);}





double dtmax;

static int init_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0)!=0;*ip=i;*tp=t;return ret;}


#line 196
      static int init(const int i,const double t,Event *_ev){tracing("init","/home/spencer/basilisk/src/navier-stokes/centered.h",0);
{
  trash (((vector[]){uf,{{-1},{-1}}}));
  if(!is_constant(fm.x)){
  
#line 199
foreach_face_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 199 \"/home/spencer/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(uf.x,0,0,0) = val(fm.x,0,0,0)*((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.);}end_is_face_x()\n// #line 199\nis_face_y(){\n    val_out_(uf.y,0,0,0) = val(fm.y,0,0,0)*((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.);}end_is_face_y()}")}){_stencil_is_face_x(){
    { _stencil_val(fm.x,0,0,0);_stencil_val(u.x,0,0,0); _stencil_val(u.x,0 -1,0,0);_stencil_val_a(uf.x,0,0,0);  }}end__stencil_is_face_x()
#line 199
_stencil_is_face_y(){
    { _stencil_val(fm.y,0,0,0);_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,0 -1,0);_stencil_val_a(uf.y,0,0,0);  }}end__stencil_is_face_y()}end_foreach_face_stencil();{
#line 199
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = val(fm.x,0,0,0)*((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.);}end_is_face_x()
#line 199
is_face_y(){
    val(uf.y,0,0,0) = val(fm.y,0,0,0)*((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.);}end_is_face_y()}end_foreach_face_generic();}}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 199
foreach_face_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 199 \"/home/spencer/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(uf.x,0,0,0) = _const_fm.x*((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.);}end_is_face_x()\n// #line 199\nis_face_y(){\n    val_out_(uf.y,0,0,0) = _const_fm.y*((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.);}end_is_face_y()}")}){_stencil_is_face_x(){
    {;_stencil_val(u.x,0,0,0); _stencil_val(u.x,0 -1,0,0);_stencil_val_a(uf.x,0,0,0);  }}end__stencil_is_face_x()
#line 199
_stencil_is_face_y(){
    {;_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,0 -1,0);_stencil_val_a(uf.y,0,0,0);  }}end__stencil_is_face_y()}end_foreach_face_stencil();
  {
#line 199
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = _const_fm.x*((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.);}end_is_face_x()
#line 199
is_face_y(){
    val(uf.y,0,0,0) = _const_fm.y*((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.);}end_is_face_y()}end_foreach_face_generic();}}




  event ("properties");





  dtmax = DT;
  event ("stability");
}{end_tracing("init","/home/spencer/basilisk/src/navier-stokes/centered.h",0);return 0;}end_tracing("init","/home/spencer/basilisk/src/navier-stokes/centered.h",0);}








static int set_dtmax_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}
# 222 "/home/spencer/basilisk/src/navier-stokes/centered.h"
      static int set_dtmax(const int i,const double t,Event *_ev){tracing("set_dtmax","/home/spencer/basilisk/src/navier-stokes/centered.h",0); dtmax = DT;{end_tracing("set_dtmax","/home/spencer/basilisk/src/navier-stokes/centered.h",0);return 0;}end_tracing("set_dtmax","/home/spencer/basilisk/src/navier-stokes/centered.h",0);}

static int stability_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}


#line 224
      static int stability(const int i,const double t,Event *_ev){tracing("stability","/home/spencer/basilisk/src/navier-stokes/centered.h",0); {
  dt = dtnext (stokes ? dtmax : timestep (uf, dtmax));
}{end_tracing("stability","/home/spencer/basilisk/src/navier-stokes/centered.h",0);return 0;}end_tracing("stability","/home/spencer/basilisk/src/navier-stokes/centered.h",0);}







static int vof_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}








#line 234
static int vof(const int i,const double t,Event *_ev){;return 0;}
static int tracer_advection_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}

#line 235
static int tracer_advection(const int i,const double t,Event *_ev){;return 0;}
static int tracer_diffusion_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}

#line 236
static int tracer_diffusion(const int i,const double t,Event *_ev){;return 0;}






static int properties_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}







#line 243
static int properties(const int i,const double t,Event *_ev){;return 0;}
# 255 "/home/spencer/basilisk/src/navier-stokes/centered.h"
void prediction()
{
  vector du;
   {
    scalar s = new_scalar("s");
    du.x = s;
  } 
#line 258
{
    scalar s = new_scalar("s");
    du.y = s;
  }

  if (_attribute[u.x.i].gradient)
    {
    
#line 264
foreach_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"du","vector",(void *)&du,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n      \n// #line 265 \"/home/spencer/basilisk/src/navier-stokes/centered.h\"\n{ {\n\n\n\n\n\n   val_out_(du.x,0,0,0) = u.x.gradient (val(u.x,-1,0,0), val(u.x,0,0,0), val(u.x,1,0,0))/Delta;\n      } \n// #line 265\n{\n\n\n\n\n\n   val_out_(du.y,0,0,0) = u.y.gradient (val(u.y,0,-1,0), val(u.y,0,0,0), val(u.y,0,1,0))/Delta;\n      }}")})
      { {





_stencil_val(u.x,-1,0,0); _stencil_val(u.x,0,0,0); _stencil_val(u.x,1,0,0);





   
#line 271
_stencil_val_a(du.x,0,0,0);   
      } 
#line 265
{





_stencil_val(u.y,0,-1,0); _stencil_val(u.y,0,0,0); _stencil_val(u.y,0,1,0);





   
#line 271
_stencil_val_a(du.y,0,0,0);   
      }}end_foreach_stencil();{
#line 264
foreach()
      { {





   val(du.x,0,0,0) = _attribute[u.x.i].gradient (val(u.x,-1,0,0), val(u.x,0,0,0), val(u.x,1,0,0))/Delta;
      } 
#line 265
{





   val(du.y,0,0,0) = _attribute[u.y.i].gradient (val(u.y,0,-1,0), val(u.y,0,0,0), val(u.y,0,1,0))/Delta;
      }}end_foreach();}}
  else
    {
    
#line 274
foreach_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"du","vector",(void *)&du,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n      \n// #line 275 \"/home/spencer/basilisk/src/navier-stokes/centered.h\"\n{ {\n\n\n\n\n\n   val_out_(du.x,0,0,0) = (val(u.x,1,0,0) - val(u.x,-1,0,0))/(2.*Delta);\n    } \n// #line 275\n{\n\n\n\n\n\n   val_out_(du.y,0,0,0) = (val(u.y,0,1,0) - val(u.y,0,-1,0))/(2.*Delta);\n    }}")})
      { {





_stencil_val(u.x,1,0,0); _stencil_val(u.x,-1,0,0);





   
#line 281
_stencil_val_a(du.x,0,0,0);   
    } 
#line 275
{





_stencil_val(u.y,0,1,0); _stencil_val(u.y,0,-1,0);





   
#line 281
_stencil_val_a(du.y,0,0,0);   
    }}end_foreach_stencil();{
#line 274
foreach()
      { {





   val(du.x,0,0,0) = (val(u.x,1,0,0) - val(u.x,-1,0,0))/(2.*Delta);
    } 
#line 275
{





   val(du.y,0,0,0) = (val(u.y,0,1,0) - val(u.y,0,-1,0))/(2.*Delta);
    }}end_foreach();}}

  trash (((vector[]){uf,{{-1},{-1}}}));
  if(!is_constant(fm.x)){
  
#line 285
foreach_face_stencil(1,{(NonLocal[]){{"fm","vector",(void *)&fm,NULL,0},{"du","vector",(void *)&du,NULL,0},{"g","vector",(void *)&g,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"u","vector",(void *)&u,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 285 \"/home/spencer/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){ {\n    real un = dt*(val(u.x,0,0,0) + val(u.x,-1,0,0))/(2.*Delta), s = sign(un);\n    int i = -(s + 1.)/2.;\n    val_out_(uf.x,0,0,0) = val(u.x,i,0,0) + (val(g.x,0,0,0) + val(g.x,-1,0,0))*dt/4. + s*(1. - s*un)*val(du.x,i,0,0)*Delta/2.;\n\n    if (val(fm.y,i,0,0) && val(fm.y,i,1,0)) {\n      real fyy = val(u.y,i,0,0) < 0. ? val(u.x,i,1,0) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,-1,0);\n      val_out_(uf.x,0,0,0) -= dt*val(u.y,i,0,0)*fyy/(2.*Delta);\n    }\n\n\n\n\n\n\n\n    val_out_(uf.x,0,0,0) *= val(fm.x,0,0,0);\n  }}end_is_face_x()\n// #line 285\nis_face_y(){ {\n    real un = dt*(val(u.y,0,0,0) + val(u.y,0,-1,0))/(2.*Delta), s = sign(un);\n    int i = -(s + 1.)/2.;\n    val_out_(uf.y,0,0,0) = val(u.y,0,i,0) + (val(g.y,0,0,0) + val(g.y,0,-1,0))*dt/4. + s*(1. - s*un)*val(du.y,0,i,0)*Delta/2.;\n\n    if (val(fm.x,0,i,0) && val(fm.x,1,i,0)) {\n      real fyy = val(u.x,0,i,0) < 0. ? val(u.y,1,i,0) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,-1,i,0);\n      val_out_(uf.y,0,0,0) -= dt*val(u.x,0,i,0)*fyy/(2.*Delta);\n    }\n\n\n\n\n\n\n\n    val_out_(uf.y,0,0,0) *= val(fm.y,0,0,0);\n  }}end_is_face_y()}")}){_stencil_is_face_x(){ {       
     _stencil_val(u.x,-1,0,0);_stencil_val(u.x,0,0,0);     
    
_stencil_val(u.x, o_stencil,0,0);_stencil_val(g.x,0,0,0); _stencil_val(g.x,-1,0,0);_stencil_val(du.x,o_stencil,0,0);
    
#line 288
_stencil_val_a(uf.x,0,0,0);

_stencil_val(fm.y,o_stencil,0,0); _stencil_val(fm.y,o_stencil,1,0); {        
       _stencil_val(u.x,o_stencil,-1,0);_stencil_val(u.x, o_stencil,0,0);_stencil_val(u.x, o_stencil,0,0); _stencil_val(u.x,o_stencil,1,0);_stencil_val(u.y, o_stencil,0,0);
_stencil_val(u.y,o_stencil,0,0);
      
#line 292
_stencil_val_r(uf.x,0,0,0);  
    } 







_stencil_val(fm.x,0,0,0);        

      







    
#line 301
_stencil_val_r(uf.x,0,0,0); 
  }}end__stencil_is_face_x()
#line 285
_stencil_is_face_y(){ {       
     _stencil_val(u.y,0,-1,0);_stencil_val(u.y,0,0,0);     
    
_stencil_val(u.y,0, o_stencil,0);_stencil_val(g.y,0,0,0); _stencil_val(g.y,0,-1,0);_stencil_val(du.y,0,o_stencil,0);
    
#line 288
_stencil_val_a(uf.y,0,0,0);

_stencil_val(fm.x,0,o_stencil,0); _stencil_val(fm.x,1,o_stencil,0); {        
       _stencil_val(u.y,-1,o_stencil,0);_stencil_val(u.y,0, o_stencil,0);_stencil_val(u.y,0, o_stencil,0); _stencil_val(u.y,1,o_stencil,0);_stencil_val(u.x,0, o_stencil,0);
_stencil_val(u.x,0,o_stencil,0);
      
#line 292
_stencil_val_r(uf.y,0,0,0);  
    } 







_stencil_val(fm.y,0,0,0);        

      







    
#line 301
_stencil_val_r(uf.y,0,0,0); 
  }}end__stencil_is_face_y()}end_foreach_face_stencil();{
#line 285
foreach_face_generic(){is_face_x(){ {
    double un = dt*(val(u.x,0,0,0) + val(u.x,-1,0,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.x,0,0,0) = val(u.x,i,0,0) + (val(g.x,0,0,0) + val(g.x,-1,0,0))*dt/4. + s*(1. - s*un)*val(du.x,i,0,0)*Delta/2.;

    if (val(fm.y,i,0,0) && val(fm.y,i,1,0)) {
      double fyy = val(u.y,i,0,0) < 0. ? val(u.x,i,1,0) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,-1,0);
      val(uf.x,0,0,0) -= dt*val(u.y,i,0,0)*fyy/(2.*Delta);
    }







    val(uf.x,0,0,0) *= val(fm.x,0,0,0);
  }}end_is_face_x()
#line 285
is_face_y(){ {
    double un = dt*(val(u.y,0,0,0) + val(u.y,0,-1,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.y,0,0,0) = val(u.y,0,i,0) + (val(g.y,0,0,0) + val(g.y,0,-1,0))*dt/4. + s*(1. - s*un)*val(du.y,0,i,0)*Delta/2.;

    if (val(fm.x,0,i,0) && val(fm.x,1,i,0)) {
      double fyy = val(u.x,0,i,0) < 0. ? val(u.y,1,i,0) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,-1,i,0);
      val(uf.y,0,0,0) -= dt*val(u.x,0,i,0)*fyy/(2.*Delta);
    }







    val(uf.y,0,0,0) *= val(fm.y,0,0,0);
  }}end_is_face_y()}end_foreach_face_generic();}}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 285
foreach_face_stencil(1,{(NonLocal[]){{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"du","vector",(void *)&du,NULL,0},{"g","vector",(void *)&g,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"u","vector",(void *)&u,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 285 \"/home/spencer/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){ {\n    real un = dt*(val(u.x,0,0,0) + val(u.x,-1,0,0))/(2.*Delta), s = sign(un);\n    int i = -(s + 1.)/2.;\n    val_out_(uf.x,0,0,0) = val(u.x,i,0,0) + (val(g.x,0,0,0) + val(g.x,-1,0,0))*dt/4. + s*(1. - s*un)*val(du.x,i,0,0)*Delta/2.;\n\n    if (_const_fm.y && _const_fm.y) {\n      real fyy = val(u.y,i,0,0) < 0. ? val(u.x,i,1,0) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,-1,0);\n      val_out_(uf.x,0,0,0) -= dt*val(u.y,i,0,0)*fyy/(2.*Delta);\n    }\n\n\n\n\n\n\n\n    val_out_(uf.x,0,0,0) *= _const_fm.x;\n  }}end_is_face_x()\n// #line 285\nis_face_y(){ {\n    real un = dt*(val(u.y,0,0,0) + val(u.y,0,-1,0))/(2.*Delta), s = sign(un);\n    int i = -(s + 1.)/2.;\n    val_out_(uf.y,0,0,0) = val(u.y,0,i,0) + (val(g.y,0,0,0) + val(g.y,0,-1,0))*dt/4. + s*(1. - s*un)*val(du.y,0,i,0)*Delta/2.;\n\n    if (_const_fm.x && _const_fm.x) {\n      real fyy = val(u.x,0,i,0) < 0. ? val(u.y,1,i,0) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,-1,i,0);\n      val_out_(uf.y,0,0,0) -= dt*val(u.x,0,i,0)*fyy/(2.*Delta);\n    }\n\n\n\n\n\n\n\n    val_out_(uf.y,0,0,0) *= _const_fm.y;\n  }}end_is_face_y()}")}){_stencil_is_face_x(){ {       
     _stencil_val(u.x,-1,0,0);_stencil_val(u.x,0,0,0);     
    
_stencil_val(u.x, o_stencil,0,0);_stencil_val(g.x,0,0,0); _stencil_val(g.x,-1,0,0);_stencil_val(du.x,o_stencil,0,0);
    
#line 288
_stencil_val_a(uf.x,0,0,0);

;; {        
       _stencil_val(u.x,o_stencil,-1,0);_stencil_val(u.x, o_stencil,0,0);_stencil_val(u.x, o_stencil,0,0); _stencil_val(u.x,o_stencil,1,0);_stencil_val(u.y, o_stencil,0,0);
_stencil_val(u.y,o_stencil,0,0);
      
#line 292
_stencil_val_r(uf.x,0,0,0);  
    }







;        

      







    
#line 301
_stencil_val_r(uf.x,0,0,0); 
  }}end__stencil_is_face_x()
#line 285
_stencil_is_face_y(){ {       
     _stencil_val(u.y,0,-1,0);_stencil_val(u.y,0,0,0);     
    
_stencil_val(u.y,0, o_stencil,0);_stencil_val(g.y,0,0,0); _stencil_val(g.y,0,-1,0);_stencil_val(du.y,0,o_stencil,0);
    
#line 288
_stencil_val_a(uf.y,0,0,0);

;; {        
       _stencil_val(u.y,-1,o_stencil,0);_stencil_val(u.y,0, o_stencil,0);_stencil_val(u.y,0, o_stencil,0); _stencil_val(u.y,1,o_stencil,0);_stencil_val(u.x,0, o_stencil,0);
_stencil_val(u.x,0,o_stencil,0);
      
#line 292
_stencil_val_r(uf.y,0,0,0);  
    }







;        

      







    
#line 301
_stencil_val_r(uf.y,0,0,0); 
  }}end__stencil_is_face_y()}end_foreach_face_stencil();
  {
#line 285
foreach_face_generic(){is_face_x(){ {
    double un = dt*(val(u.x,0,0,0) + val(u.x,-1,0,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.x,0,0,0) = val(u.x,i,0,0) + (val(g.x,0,0,0) + val(g.x,-1,0,0))*dt/4. + s*(1. - s*un)*val(du.x,i,0,0)*Delta/2.;

    if (_const_fm.y && _const_fm.y) {
      double fyy = val(u.y,i,0,0) < 0. ? val(u.x,i,1,0) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,-1,0);
      val(uf.x,0,0,0) -= dt*val(u.y,i,0,0)*fyy/(2.*Delta);
    }







    val(uf.x,0,0,0) *= _const_fm.x;
  }}end_is_face_x()
#line 285
is_face_y(){ {
    double un = dt*(val(u.y,0,0,0) + val(u.y,0,-1,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.y,0,0,0) = val(u.y,0,i,0) + (val(g.y,0,0,0) + val(g.y,0,-1,0))*dt/4. + s*(1. - s*un)*val(du.y,0,i,0)*Delta/2.;

    if (_const_fm.x && _const_fm.x) {
      double fyy = val(u.x,0,i,0) < 0. ? val(u.y,1,i,0) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,-1,i,0);
      val(uf.y,0,0,0) -= dt*val(u.x,0,i,0)*fyy/(2.*Delta);
    }







    val(uf.y,0,0,0) *= _const_fm.y;
  }}end_is_face_y()}end_foreach_face_generic();}}

  delete ((scalar *)((vector[]){du,{{-1},{-1}}}));
}
#line 316
static int advection_term_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}
# 316 "/home/spencer/basilisk/src/navier-stokes/centered.h"
      static int advection_term(const int i,const double t,Event *_ev){tracing("advection_term","/home/spencer/basilisk/src/navier-stokes/centered.h",0);
{
  if (!stokes) {
    prediction();
    mgpf = project (uf, pf, alpha, dt/2., mgpf.nrelax);
    advection ((scalar *)((vector[]){u,{{-1},{-1}}}), uf, dt, (scalar *)((vector[]){g,{{-1},{-1}}}));
  }
}{end_tracing("advection_term","/home/spencer/basilisk/src/navier-stokes/centered.h",0);return 0;}end_tracing("advection_term","/home/spencer/basilisk/src/navier-stokes/centered.h",0);}







static void correction (double dt)
{
  foreach_stencil(1,{(NonLocal[]){{"g","vector",(void *)&g,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"u","vector",(void *)&u,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n    \n// #line 334 \"/home/spencer/basilisk/src/navier-stokes/centered.h\"\n{\n      val_out_(u.x,0,0,0) += dt*val(g.x,0,0,0);\n      \n// #line 335\nval_out_(u.y,0,0,0) += dt*val(g.y,0,0,0);}")})
    {
      {_stencil_val(g.x,0,0,0);_stencil_val_r(u.x,0,0,0);  }
      
#line 335
{_stencil_val(g.y,0,0,0);_stencil_val_r(u.y,0,0,0);  }}end_foreach_stencil();
  {
#line 333
foreach()
    {
      val(u.x,0,0,0) += dt*val(g.x,0,0,0);
      
#line 335
val(u.y,0,0,0) += dt*val(g.y,0,0,0);}end_foreach();}
}








static int viscous_term_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}
# 345 "/home/spencer/basilisk/src/navier-stokes/centered.h"
      static int viscous_term(const int i,const double t,Event *_ev){tracing("viscous_term","/home/spencer/basilisk/src/navier-stokes/centered.h",0);
{
  if (constant(mu.x) != 0.) {
    correction (dt);
    mgu = viscosity (u, mu, rho, dt, mgu.nrelax
#line 290 "/home/spencer/basilisk/src/viscosity.h"
, NULL
#line 349 "/home/spencer/basilisk/src/navier-stokes/centered.h"
);
    correction (-dt);
  }




  if (!is_constant(a.x)) {
    vector af = a;
    trash (((vector[]){af,{{-1},{-1}}}));
    foreach_face_stencil(1,{(NonLocal[]){{"af","vector",(void *)&af,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 359 \"/home/spencer/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n      val_out_(af.x,0,0,0) = 0.;}end_is_face_x()\n// #line 359\nis_face_y(){\n      val_out_(af.y,0,0,0) = 0.;}end_is_face_y()}")}){_stencil_is_face_x(){
      {_stencil_val_a(af.x,0,0,0);  }}end__stencil_is_face_x()
#line 359
_stencil_is_face_y(){
      {_stencil_val_a(af.y,0,0,0);  }}end__stencil_is_face_y()}end_foreach_face_stencil();
    {
#line 359
foreach_face_generic(){is_face_x(){
      val(af.x,0,0,0) = 0.;}end_is_face_x()
#line 359
is_face_y(){
      val(af.y,0,0,0) = 0.;}end_is_face_y()}end_foreach_face_generic();}
  }
}{end_tracing("viscous_term","/home/spencer/basilisk/src/navier-stokes/centered.h",0);return 0;}end_tracing("viscous_term","/home/spencer/basilisk/src/navier-stokes/centered.h",0);}
#line 381
static int acceleration_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}
# 381 "/home/spencer/basilisk/src/navier-stokes/centered.h"
      static int acceleration(const int i,const double t,Event *_ev){tracing("acceleration","/home/spencer/basilisk/src/navier-stokes/centered.h",0);
{
  trash (((vector[]){uf,{{-1},{-1}}}));
  if(!is_constant(fm.x) && !is_constant(a.x)){
  
#line 384
foreach_face_stencil(1,{(NonLocal[]){{"a","vector",(void *)&a,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"u","vector",(void *)&u,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 384 \"/home/spencer/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(uf.x,0,0,0) = val(fm.x,0,0,0)*(((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.) + dt*val(a.x,0,0,0));}end_is_face_x()\n// #line 384\nis_face_y(){\n    val_out_(uf.y,0,0,0) = val(fm.y,0,0,0)*(((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.) + dt*val(a.y,0,0,0));}end_is_face_y()}")}){_stencil_is_face_x(){
    { _stencil_val(fm.x,0,0,0);_stencil_val(u.x,0,0,0); _stencil_val(u.x,0 -1,0,0);_stencil_val(a.x,0,0,0);_stencil_val_a(uf.x,0,0,0);    }}end__stencil_is_face_x()
#line 384
_stencil_is_face_y(){
    { _stencil_val(fm.y,0,0,0);_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,0 -1,0);_stencil_val(a.y,0,0,0);_stencil_val_a(uf.y,0,0,0);    }}end__stencil_is_face_y()}end_foreach_face_stencil();{
#line 384
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = val(fm.x,0,0,0)*(((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.) + dt*val(a.x,0,0,0));}end_is_face_x()
#line 384
is_face_y(){
    val(uf.y,0,0,0) = val(fm.y,0,0,0)*(((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.) + dt*val(a.y,0,0,0));}end_is_face_y()}end_foreach_face_generic();}}else if(is_constant(fm.x) && !is_constant(a.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 384
foreach_face_stencil(1,{(NonLocal[]){{"a","vector",(void *)&a,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"u","vector",(void *)&u,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 384 \"/home/spencer/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(uf.x,0,0,0) = _const_fm.x*(((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.) + dt*val(a.x,0,0,0));}end_is_face_x()\n// #line 384\nis_face_y(){\n    val_out_(uf.y,0,0,0) = _const_fm.y*(((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.) + dt*val(a.y,0,0,0));}end_is_face_y()}")}){_stencil_is_face_x(){
    {;_stencil_val(u.x,0,0,0); _stencil_val(u.x,0 -1,0,0);_stencil_val(a.x,0,0,0);_stencil_val_a(uf.x,0,0,0);    }}end__stencil_is_face_x()
#line 384
_stencil_is_face_y(){
    {;_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,0 -1,0);_stencil_val(a.y,0,0,0);_stencil_val_a(uf.y,0,0,0);    }}end__stencil_is_face_y()}end_foreach_face_stencil();
  {
#line 384
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = _const_fm.x*(((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.) + dt*val(a.x,0,0,0));}end_is_face_x()
#line 384
is_face_y(){
    val(uf.y,0,0,0) = _const_fm.y*(((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.) + dt*val(a.y,0,0,0));}end_is_face_y()}end_foreach_face_generic();}}else if(!is_constant(fm.x) && is_constant(a.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);
  
#line 384
foreach_face_stencil(1,{(NonLocal[]){{"_const_a","_coord",(void *)&_const_a,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"u","vector",(void *)&u,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 384 \"/home/spencer/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(uf.x,0,0,0) = val(fm.x,0,0,0)*(((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.) + dt*_const_a.x);}end_is_face_x()\n// #line 384\nis_face_y(){\n    val_out_(uf.y,0,0,0) = val(fm.y,0,0,0)*(((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.) + dt*_const_a.y);}end_is_face_y()}")}){_stencil_is_face_x(){
    { _stencil_val(fm.x,0,0,0);_stencil_val(u.x,0,0,0); _stencil_val(u.x,0 -1,0,0);;_stencil_val_a(uf.x,0,0,0);    }}end__stencil_is_face_x()
#line 384
_stencil_is_face_y(){
    { _stencil_val(fm.y,0,0,0);_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,0 -1,0);;_stencil_val_a(uf.y,0,0,0);    }}end__stencil_is_face_y()}end_foreach_face_stencil();
  {
#line 384
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = val(fm.x,0,0,0)*(((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.) + dt*_const_a.x);}end_is_face_x()
#line 384
is_face_y(){
    val(uf.y,0,0,0) = val(fm.y,0,0,0)*(((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.) + dt*_const_a.y);}end_is_face_y()}end_foreach_face_generic();}}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);
  
#line 384
foreach_face_stencil(1,{(NonLocal[]){{"_const_a","_coord",(void *)&_const_a,NULL,0},{"dt","double",(void *)&dt,NULL,0},{"u","vector",(void *)&u,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"uf","vector",(void *)&uf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 384 \"/home/spencer/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(uf.x,0,0,0) = _const_fm.x*(((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.) + dt*_const_a.x);}end_is_face_x()\n// #line 384\nis_face_y(){\n    val_out_(uf.y,0,0,0) = _const_fm.y*(((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.) + dt*_const_a.y);}end_is_face_y()}")}){_stencil_is_face_x(){
    {;_stencil_val(u.x,0,0,0); _stencil_val(u.x,0 -1,0,0);;_stencil_val_a(uf.x,0,0,0);    }}end__stencil_is_face_x()
#line 384
_stencil_is_face_y(){
    {;_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,0 -1,0);;_stencil_val_a(uf.y,0,0,0);    }}end__stencil_is_face_y()}end_foreach_face_stencil();
  {
#line 384
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = _const_fm.x*(((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.) + dt*_const_a.x);}end_is_face_x()
#line 384
is_face_y(){
    val(uf.y,0,0,0) = _const_fm.y*(((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.) + dt*_const_a.y);}end_is_face_y()}end_foreach_face_generic();}}
}{end_tracing("acceleration","/home/spencer/basilisk/src/navier-stokes/centered.h",0);return 0;}end_tracing("acceleration","/home/spencer/basilisk/src/navier-stokes/centered.h",0);}
# 395 "/home/spencer/basilisk/src/navier-stokes/centered.h"
void centered_gradient (scalar p, vector g)
{





  vector  gf=new_face_vector("gf");
  if(!is_constant(fm.x) && !is_constant(a.x) && !is_constant(alpha.x)){
  
#line 403
foreach_face_stencil(1,{(NonLocal[]){{"p","scalar",(void *)&p,NULL,0},{"alpha","vector",(void *)&alpha,NULL,0},{"a","vector",(void *)&a,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"gf","vector",(void *)&gf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 403 \"/home/spencer/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(gf.x,0,0,0) = val(fm.x,0,0,0)*val(a.x,0,0,0) - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()\n// #line 403\nis_face_y(){\n    val_out_(gf.y,0,0,0) = val(fm.y,0,0,0)*val(a.y,0,0,0) - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}")}){_stencil_is_face_x(){
    { _stencil_val(fm.x,0,0,0);_stencil_val(a.x,0,0,0); _stencil_val(alpha.x,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,-1,0,0);_stencil_val_a(gf.x,0,0,0);   }}end__stencil_is_face_x()
#line 403
_stencil_is_face_y(){
    { _stencil_val(fm.y,0,0,0);_stencil_val(a.y,0,0,0); _stencil_val(alpha.y,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0,-1,0);_stencil_val_a(gf.y,0,0,0);   }}end__stencil_is_face_y()}end_foreach_face_stencil();{
#line 403
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = val(fm.x,0,0,0)*val(a.x,0,0,0) - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 403
is_face_y(){
    val(gf.y,0,0,0) = val(fm.y,0,0,0)*val(a.y,0,0,0) - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}end_foreach_face_generic();}}else if(is_constant(fm.x) && !is_constant(a.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 403
foreach_face_stencil(1,{(NonLocal[]){{"p","scalar",(void *)&p,NULL,0},{"alpha","vector",(void *)&alpha,NULL,0},{"a","vector",(void *)&a,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"gf","vector",(void *)&gf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 403 \"/home/spencer/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(gf.x,0,0,0) = _const_fm.x*val(a.x,0,0,0) - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()\n// #line 403\nis_face_y(){\n    val_out_(gf.y,0,0,0) = _const_fm.y*val(a.y,0,0,0) - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}")}){_stencil_is_face_x(){
    {;_stencil_val(a.x,0,0,0); _stencil_val(alpha.x,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,-1,0,0);_stencil_val_a(gf.x,0,0,0);   }}end__stencil_is_face_x()
#line 403
_stencil_is_face_y(){
    {;_stencil_val(a.y,0,0,0); _stencil_val(alpha.y,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0,-1,0);_stencil_val_a(gf.y,0,0,0);   }}end__stencil_is_face_y()}end_foreach_face_stencil();
  {
#line 403
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = _const_fm.x*val(a.x,0,0,0) - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 403
is_face_y(){
    val(gf.y,0,0,0) = _const_fm.y*val(a.y,0,0,0) - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}end_foreach_face_generic();}}else if(!is_constant(fm.x) && is_constant(a.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);
  
#line 403
foreach_face_stencil(1,{(NonLocal[]){{"p","scalar",(void *)&p,NULL,0},{"alpha","vector",(void *)&alpha,NULL,0},{"_const_a","_coord",(void *)&_const_a,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"gf","vector",(void *)&gf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 403 \"/home/spencer/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(gf.x,0,0,0) = val(fm.x,0,0,0)*_const_a.x - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()\n// #line 403\nis_face_y(){\n    val_out_(gf.y,0,0,0) = val(fm.y,0,0,0)*_const_a.y - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}")}){_stencil_is_face_x(){
    { _stencil_val(fm.x,0,0,0);; _stencil_val(alpha.x,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,-1,0,0);_stencil_val_a(gf.x,0,0,0);   }}end__stencil_is_face_x()
#line 403
_stencil_is_face_y(){
    { _stencil_val(fm.y,0,0,0);; _stencil_val(alpha.y,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0,-1,0);_stencil_val_a(gf.y,0,0,0);   }}end__stencil_is_face_y()}end_foreach_face_stencil();
  {
#line 403
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = val(fm.x,0,0,0)*_const_a.x - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 403
is_face_y(){
    val(gf.y,0,0,0) = val(fm.y,0,0,0)*_const_a.y - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}end_foreach_face_generic();}}else if(is_constant(fm.x) && is_constant(a.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);
  
#line 403
foreach_face_stencil(1,{(NonLocal[]){{"p","scalar",(void *)&p,NULL,0},{"alpha","vector",(void *)&alpha,NULL,0},{"_const_a","_coord",(void *)&_const_a,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"gf","vector",(void *)&gf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 403 \"/home/spencer/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(gf.x,0,0,0) = _const_fm.x*_const_a.x - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()\n// #line 403\nis_face_y(){\n    val_out_(gf.y,0,0,0) = _const_fm.y*_const_a.y - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}")}){_stencil_is_face_x(){
    {;; _stencil_val(alpha.x,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,-1,0,0);_stencil_val_a(gf.x,0,0,0);   }}end__stencil_is_face_x()
#line 403
_stencil_is_face_y(){
    {;; _stencil_val(alpha.y,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0,-1,0);_stencil_val_a(gf.y,0,0,0);   }}end__stencil_is_face_y()}end_foreach_face_stencil();
  {
#line 403
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = _const_fm.x*_const_a.x - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 403
is_face_y(){
    val(gf.y,0,0,0) = _const_fm.y*_const_a.y - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}end_foreach_face_generic();}}else if(!is_constant(fm.x) && !is_constant(a.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  
#line 403
foreach_face_stencil(1,{(NonLocal[]){{"p","scalar",(void *)&p,NULL,0},{"_const_alpha","_coord",(void *)&_const_alpha,NULL,0},{"a","vector",(void *)&a,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"gf","vector",(void *)&gf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 403 \"/home/spencer/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(gf.x,0,0,0) = val(fm.x,0,0,0)*val(a.x,0,0,0) - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()\n// #line 403\nis_face_y(){\n    val_out_(gf.y,0,0,0) = val(fm.y,0,0,0)*val(a.y,0,0,0) - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}")}){_stencil_is_face_x(){
    { _stencil_val(fm.x,0,0,0);_stencil_val(a.x,0,0,0);;_stencil_val(p,0,0,0); _stencil_val(p,-1,0,0);_stencil_val_a(gf.x,0,0,0);   }}end__stencil_is_face_x()
#line 403
_stencil_is_face_y(){
    { _stencil_val(fm.y,0,0,0);_stencil_val(a.y,0,0,0);;_stencil_val(p,0,0,0); _stencil_val(p,0,-1,0);_stencil_val_a(gf.y,0,0,0);   }}end__stencil_is_face_y()}end_foreach_face_stencil();
  {
#line 403
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = val(fm.x,0,0,0)*val(a.x,0,0,0) - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 403
is_face_y(){
    val(gf.y,0,0,0) = val(fm.y,0,0,0)*val(a.y,0,0,0) - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}end_foreach_face_generic();}}else if(is_constant(fm.x) && !is_constant(a.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  
#line 403
foreach_face_stencil(1,{(NonLocal[]){{"p","scalar",(void *)&p,NULL,0},{"_const_alpha","_coord",(void *)&_const_alpha,NULL,0},{"a","vector",(void *)&a,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"gf","vector",(void *)&gf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 403 \"/home/spencer/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(gf.x,0,0,0) = _const_fm.x*val(a.x,0,0,0) - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()\n// #line 403\nis_face_y(){\n    val_out_(gf.y,0,0,0) = _const_fm.y*val(a.y,0,0,0) - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}")}){_stencil_is_face_x(){
    {;_stencil_val(a.x,0,0,0);;_stencil_val(p,0,0,0); _stencil_val(p,-1,0,0);_stencil_val_a(gf.x,0,0,0);   }}end__stencil_is_face_x()
#line 403
_stencil_is_face_y(){
    {;_stencil_val(a.y,0,0,0);;_stencil_val(p,0,0,0); _stencil_val(p,0,-1,0);_stencil_val_a(gf.y,0,0,0);   }}end__stencil_is_face_y()}end_foreach_face_stencil();
  {
#line 403
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = _const_fm.x*val(a.x,0,0,0) - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 403
is_face_y(){
    val(gf.y,0,0,0) = _const_fm.y*val(a.y,0,0,0) - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}end_foreach_face_generic();}}else if(!is_constant(fm.x) && is_constant(a.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  
#line 403
foreach_face_stencil(1,{(NonLocal[]){{"p","scalar",(void *)&p,NULL,0},{"_const_alpha","_coord",(void *)&_const_alpha,NULL,0},{"_const_a","_coord",(void *)&_const_a,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"gf","vector",(void *)&gf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 403 \"/home/spencer/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(gf.x,0,0,0) = val(fm.x,0,0,0)*_const_a.x - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()\n// #line 403\nis_face_y(){\n    val_out_(gf.y,0,0,0) = val(fm.y,0,0,0)*_const_a.y - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}")}){_stencil_is_face_x(){
    { _stencil_val(fm.x,0,0,0);;;_stencil_val(p,0,0,0); _stencil_val(p,-1,0,0);_stencil_val_a(gf.x,0,0,0);   }}end__stencil_is_face_x()
#line 403
_stencil_is_face_y(){
    { _stencil_val(fm.y,0,0,0);;;_stencil_val(p,0,0,0); _stencil_val(p,0,-1,0);_stencil_val_a(gf.y,0,0,0);   }}end__stencil_is_face_y()}end_foreach_face_stencil();
  {
#line 403
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = val(fm.x,0,0,0)*_const_a.x - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 403
is_face_y(){
    val(gf.y,0,0,0) = val(fm.y,0,0,0)*_const_a.y - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}end_foreach_face_generic();}}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  
#line 403
foreach_face_stencil(1,{(NonLocal[]){{"p","scalar",(void *)&p,NULL,0},{"_const_alpha","_coord",(void *)&_const_alpha,NULL,0},{"_const_a","_coord",(void *)&_const_a,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"gf","vector",(void *)&gf,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 403 \"/home/spencer/basilisk/src/navier-stokes/centered.h\"\n{is_face_x(){\n    val_out_(gf.x,0,0,0) = _const_fm.x*_const_a.x - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()\n// #line 403\nis_face_y(){\n    val_out_(gf.y,0,0,0) = _const_fm.y*_const_a.y - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}")}){_stencil_is_face_x(){
    {;;;_stencil_val(p,0,0,0); _stencil_val(p,-1,0,0);_stencil_val_a(gf.x,0,0,0);   }}end__stencil_is_face_x()
#line 403
_stencil_is_face_y(){
    {;;;_stencil_val(p,0,0,0); _stencil_val(p,0,-1,0);_stencil_val_a(gf.y,0,0,0);   }}end__stencil_is_face_y()}end_foreach_face_stencil();
  {
#line 403
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = _const_fm.x*_const_a.x - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 403
is_face_y(){
    val(gf.y,0,0,0) = _const_fm.y*_const_a.y - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}end_foreach_face_generic();}}





  trash (((vector[]){g,{{-1},{-1}}}));
  if(!is_constant(fm.x)){
  
#line 411
foreach_stencil(1,{(NonLocal[]){{"fm","vector",(void *)&fm,NULL,0},{"gf","vector",(void *)&gf,NULL,0},{"g","vector",(void *)&g,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n    \n// #line 412 \"/home/spencer/basilisk/src/navier-stokes/centered.h\"\n{\n      val_out_(g.x,0,0,0) = (val(gf.x,0,0,0) + val(gf.x,1,0,0))/(val(fm.x,0,0,0) + val(fm.x,1,0,0) + 0.);\n      \n// #line 413\nval_out_(g.y,0,0,0) = (val(gf.y,0,0,0) + val(gf.y,0,1,0))/(val(fm.y,0,0,0) + val(fm.y,0,1,0) + 0.);}")})
    {
      {_stencil_val(gf.x,0,0,0); _stencil_val(gf.x,1,0,0);_stencil_val(fm.x,0,0,0); _stencil_val(fm.x,1,0,0);_stencil_val_a(g.x,0,0,0);      }
      
#line 413
{_stencil_val(gf.y,0,0,0); _stencil_val(gf.y,0,1,0);_stencil_val(fm.y,0,0,0); _stencil_val(fm.y,0,1,0);_stencil_val_a(g.y,0,0,0);      }}end_foreach_stencil();{
#line 411
foreach()
    {
      val(g.x,0,0,0) = (val(gf.x,0,0,0) + val(gf.x,1,0,0))/(val(fm.x,0,0,0) + val(fm.x,1,0,0) + 0.);
      
#line 413
val(g.y,0,0,0) = (val(gf.y,0,0,0) + val(gf.y,0,1,0))/(val(fm.y,0,0,0) + val(fm.y,0,1,0) + 0.);}end_foreach();}}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 411
foreach_stencil(1,{(NonLocal[]){{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"gf","vector",(void *)&gf,NULL,0},{"g","vector",(void *)&g,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n    \n// #line 412 \"/home/spencer/basilisk/src/navier-stokes/centered.h\"\n{\n      val_out_(g.x,0,0,0) = (val(gf.x,0,0,0) + val(gf.x,1,0,0))/(_const_fm.x + _const_fm.x + 0.);\n      \n// #line 413\nval_out_(g.y,0,0,0) = (val(gf.y,0,0,0) + val(gf.y,0,1,0))/(_const_fm.y + _const_fm.y + 0.);}")})
    {
      {_stencil_val(gf.x,0,0,0); _stencil_val(gf.x,1,0,0);;;_stencil_val_a(g.x,0,0,0);      }
      
#line 413
{_stencil_val(gf.y,0,0,0); _stencil_val(gf.y,0,1,0);;;_stencil_val_a(g.y,0,0,0);      }}end_foreach_stencil();
  {
#line 411
foreach()
    {
      val(g.x,0,0,0) = (val(gf.x,0,0,0) + val(gf.x,1,0,0))/(_const_fm.x + _const_fm.x + 0.);
      
#line 413
val(g.y,0,0,0) = (val(gf.y,0,0,0) + val(gf.y,0,1,0))/(_const_fm.y + _const_fm.y + 0.);}end_foreach();}}delete((scalar*)((vector[]){gf,{{-1},{-1}}}));
}






static int projection_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}







#line 421
      static int projection(const int i,const double t,Event *_ev){tracing("projection","/home/spencer/basilisk/src/navier-stokes/centered.h",0);
{
  mgp = project (uf, p, alpha, dt, mgp.nrelax);
  centered_gradient (p, g);




  correction (dt);
}{end_tracing("projection","/home/spencer/basilisk/src/navier-stokes/centered.h",0);return 0;}end_tracing("projection","/home/spencer/basilisk/src/navier-stokes/centered.h",0);}





static int end_timestep_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}






#line 436
static int end_timestep(const int i,const double t,Event *_ev){;return 0;}









static int adapt_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}
# 446 "/home/spencer/basilisk/src/navier-stokes/centered.h"
      static int adapt(const int i,const double t,Event *_ev){tracing("adapt","/home/spencer/basilisk/src/navier-stokes/centered.h",0); {






  event ("properties");
}{end_tracing("adapt","/home/spencer/basilisk/src/navier-stokes/centered.h",0);return 0;}end_tracing("adapt","/home/spencer/basilisk/src/navier-stokes/centered.h",0);}
# 2 "cylinder-osc.c" 2
# 1 "../immersed.h" 1
# 1 "./../immersed.h"
extern coord vc;
extern scalar airfoil;
extern vector sf;
vector  aF={{8},{9}};







static int end_timestep_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}








#line 12
      static int end_timestep_0(const int i,const double t,Event *_ev){tracing("end_timestep_0","./../immersed.h",0); {
  foreach_stencil(1,{(NonLocal[]){{"dt","double",(void *)&dt,NULL,0},{"u","vector",(void *)&u,NULL,0},{"vc","coord",(void *)&vc,NULL,0},{"airfoil","scalar",(void *)&airfoil,NULL,0},{"aF","vector",(void *)&aF,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n    \n// #line 14 \"./../immersed.h\"\n{ {\n      val_out_(aF.x,0,0,0) = val(airfoil,0,0,0)*((vc.x)-val(u.x,0,0,0))/dt;\n    } \n// #line 14\n{\n      val_out_(aF.y,0,0,0) = val(airfoil,0,0,0)*((vc.y)-val(u.y,0,0,0))/dt;\n    }}")})
    { { 
_stencil_val(airfoil,0,0,0);_stencil_val(u.x,0,0,0);
      
#line 15
_stencil_val_a(aF.x,0,0,0); 
    } 
#line 14
{ 
_stencil_val(airfoil,0,0,0);_stencil_val(u.y,0,0,0);
      
#line 15
_stencil_val_a(aF.y,0,0,0); 
    }}end_foreach_stencil();
  {
#line 13
foreach()
    { {
      val(aF.x,0,0,0) = val(airfoil,0,0,0)*((vc.x)-val(u.x,0,0,0))/dt;
    } 
#line 14
{
      val(aF.y,0,0,0) = val(airfoil,0,0,0)*((vc.y)-val(u.y,0,0,0))/dt;
    }}end_foreach();}
  correction(-dt);
  foreach_stencil(1,{(NonLocal[]){{"aF","vector",(void *)&aF,NULL,0},{"g","vector",(void *)&g,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n    \n// #line 19 \"./../immersed.h\"\n{\n      val_out_(g.x,0,0,0) += val(aF.x,0,0,0);\n      \n// #line 20\nval_out_(g.y,0,0,0) += val(aF.y,0,0,0);}")})
    {
      { _stencil_val(aF.x,0,0,0);_stencil_val_r(g.x,0,0,0); }
      
#line 20
{ _stencil_val(aF.y,0,0,0);_stencil_val_r(g.y,0,0,0); }}end_foreach_stencil();
  {
#line 18
foreach()
    {
      val(g.x,0,0,0) += val(aF.x,0,0,0);
      
#line 20
val(g.y,0,0,0) += val(aF.y,0,0,0);}end_foreach();}
  correction(dt);
}{end_tracing("end_timestep_0","./../immersed.h",0);return 0;}end_tracing("end_timestep_0","./../immersed.h",0);}
# 33 "./../immersed.h"
void normal_vector_face (Point point, vector fc, coord * nv) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  coord n;
   {
    n.x = val(fc.x,0,0,0) - val(fc.x,1,0,0);
  } 
#line 35
{
    n.y = val(fc.y,0,0,0) - val(fc.y,0,1,0);
  }
  double nn = sqrt(sq(n.x) + sq(n.y));
  if (nn > 0)
    {
      n.x /= nn;
      
#line 41
n.y /= nn;}
  else
    {
      n.x /= 1./2;
      
#line 44
n.y /= 1./2;}
  *nv = n;
}

void normal_vector (Point point, scalar s, coord * nv) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  coord n;
  
      n.x = (val(s,1,0,0) - val(s,-1,0,0))/(2.*Delta);
      
#line 51
n.y = (val(s,0,1,0) - val(s,0,-1,0))/(2.*Delta);
  double mag = 0;
  
    mag += sq(n.x);
    
#line 54
mag += sq(n.y);
  mag = sqrt(mag);
  if (mag > 0)
    {
      n.x /= mag;
      
#line 58
n.y /= mag;}
  *nv = n;

}


#line 48
static void _stencil_normal_vector (Point point, scalar s,_stencil_undefined * nv) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES; 
  
  
      {_stencil_val(s,1,0,0); _stencil_val(s,-1,0,0);   }
      
#line 51
{_stencil_val(s,0,1,0); _stencil_val(s,0,-1,0);   }   
     
   
     
    
  
    
        
    

}
# 82 "./../immersed.h"
void immersed_force (scalar c, vector fc, coord * F) {
  coord Fi = {0, 0};
   if(!is_constant(cm) && !is_constant(fm.x)){
   
#line 84
foreach_stencil(1,{(NonLocal[]){{"aF","vector",(void *)&aF,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"Fi","coord",(void *)&Fi,NULL,0},{"cm","scalar",(void *)&cm,NULL,0},{"c","scalar",(void *)&c,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_(" \n// #line 84 \"./../immersed.h\"\n{\n    real area = val(c,0,0,0)*(sq(Delta)*val(cm,0,0,0));\n    \n      Fi.x += - val(fm.x,0,0,0)*(val(aF.x,0,0,0))*area;\n      \n// #line 87\nFi.y += - val(fm.y,0,0,0)*(val(aF.y,0,0,0))*area;\n   }")}) {  
    _stencil_val(cm,0,0,0); _stencil_val(c,0,0,0);
    
      { _stencil_val(fm.x,0,0,0);_stencil_val(aF.x,0,0,0);  }
      
#line 87
{ _stencil_val(fm.y,0,0,0);_stencil_val(aF.y,0,0,0);  }
   }end_foreach_stencil();{
#line 84
foreach() {
    double area = val(c,0,0,0)*(sq(Delta)*val(cm,0,0,0));
    
      Fi.x += - val(fm.x,0,0,0)*(val(aF.x,0,0,0))*area;
      
#line 87
Fi.y += - val(fm.y,0,0,0)*(val(aF.y,0,0,0))*area;
   }end_foreach();}}else if(is_constant(cm) && !is_constant(fm.x)){double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
   
#line 84
foreach_stencil(1,{(NonLocal[]){{"aF","vector",(void *)&aF,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"Fi","coord",(void *)&Fi,NULL,0},{"_const_cm","double",(void *)&_const_cm,NULL,0},{"c","scalar",(void *)&c,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_(" \n// #line 84 \"./../immersed.h\"\n{\n    real area = val(c,0,0,0)*(sq(Delta)*_const_cm);\n    \n      Fi.x += - val(fm.x,0,0,0)*(val(aF.x,0,0,0))*area;\n      \n// #line 87\nFi.y += - val(fm.y,0,0,0)*(val(aF.y,0,0,0))*area;\n   }")}) {
; _stencil_val(c,0,0,0);
    
      { _stencil_val(fm.x,0,0,0);_stencil_val(aF.x,0,0,0);  }
      
#line 87
{ _stencil_val(fm.y,0,0,0);_stencil_val(aF.y,0,0,0);  }
   }end_foreach_stencil();
   {
#line 84
foreach() {
    double area = val(c,0,0,0)*(sq(Delta)*_const_cm);
    
      Fi.x += - val(fm.x,0,0,0)*(val(aF.x,0,0,0))*area;
      
#line 87
Fi.y += - val(fm.y,0,0,0)*(val(aF.y,0,0,0))*area;
   }end_foreach();}}else if(!is_constant(cm) && is_constant(fm.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
   
#line 84
foreach_stencil(1,{(NonLocal[]){{"aF","vector",(void *)&aF,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"Fi","coord",(void *)&Fi,NULL,0},{"cm","scalar",(void *)&cm,NULL,0},{"c","scalar",(void *)&c,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_(" \n// #line 84 \"./../immersed.h\"\n{\n    real area = val(c,0,0,0)*(sq(Delta)*val(cm,0,0,0));\n    \n      Fi.x += - _const_fm.x*(val(aF.x,0,0,0))*area;\n      \n// #line 87\nFi.y += - _const_fm.y*(val(aF.y,0,0,0))*area;\n   }")}) {  
    _stencil_val(cm,0,0,0); _stencil_val(c,0,0,0);
    
      {;_stencil_val(aF.x,0,0,0);  }
      
#line 87
{;_stencil_val(aF.y,0,0,0);  }
   }end_foreach_stencil();
   {
#line 84
foreach() {
    double area = val(c,0,0,0)*(sq(Delta)*val(cm,0,0,0));
    
      Fi.x += - _const_fm.x*(val(aF.x,0,0,0))*area;
      
#line 87
Fi.y += - _const_fm.y*(val(aF.y,0,0,0))*area;
   }end_foreach();}}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
   
#line 84
foreach_stencil(1,{(NonLocal[]){{"aF","vector",(void *)&aF,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"Fi","coord",(void *)&Fi,NULL,0},{"_const_cm","double",(void *)&_const_cm,NULL,0},{"c","scalar",(void *)&c,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_(" \n// #line 84 \"./../immersed.h\"\n{\n    real area = val(c,0,0,0)*(sq(Delta)*_const_cm);\n    \n      Fi.x += - _const_fm.x*(val(aF.x,0,0,0))*area;\n      \n// #line 87\nFi.y += - _const_fm.y*(val(aF.y,0,0,0))*area;\n   }")}) {
; _stencil_val(c,0,0,0);
    
      {;_stencil_val(aF.x,0,0,0);  }
      
#line 87
{;_stencil_val(aF.y,0,0,0);  }
   }end_foreach_stencil();
   {
#line 84
foreach() {
    double area = val(c,0,0,0)*(sq(Delta)*_const_cm);
    
      Fi.x += - _const_fm.x*(val(aF.x,0,0,0))*area;
      
#line 87
Fi.y += - _const_fm.y*(val(aF.y,0,0,0))*area;
   }end_foreach();}}
  *F = Fi;
}
# 3 "cylinder-osc.c" 2
# 1 "view.h" 1
# 4 "cylinder-osc.c" 2

int LEVEL = 9;

double RE, KC, Amp, freq = 1;
const double U0 = 1;
const double D = 0.5;
const double T_max = 30;

coord ci = {0};
coord vc = {0};
coord xc = {0};

scalar  airfoil={10};
vector  sf={{11},{12}};
vector  muv={{13},{14}};

static double _boundary4(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _dirichlet(0, point, neighbor, _s, data));}}}static double _boundary4_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _dirichlet_homogeneous(0, point, neighbor, _s, data));}}}
static double _boundary5(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _dirichlet(0, point, neighbor, _s, data));}}}static double _boundary5_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _dirichlet_homogeneous(0, point, neighbor, _s, data));}}}
static double _boundary6(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(0, point, neighbor, _s, data));}}}static double _boundary6_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(0, point, neighbor, _s, data));}}}
static double _boundary7(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(0, point, neighbor, _s, data));}}}static double _boundary7_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(0, point, neighbor, _s, data));}}}

static double _boundary8(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _dirichlet(0, point, neighbor, _s, data));}}}static double _boundary8_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _dirichlet_homogeneous(0, point, neighbor, _s, data));}}}
static double _boundary9(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _dirichlet(0, point, neighbor, _s, data));}}}static double _boundary9_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _dirichlet_homogeneous(0, point, neighbor, _s, data));}}}
static double _boundary10(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(0, point, neighbor, _s, data));}}}static double _boundary10_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(0, point, neighbor, _s, data));}}}
static double _boundary11(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(0, point, neighbor, _s, data));}}}static double _boundary11_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(0, point, neighbor, _s, data));}}}

static double _boundary12(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(0, point, neighbor, _s, data));}}}static double _boundary12_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(0, point, neighbor, _s, data));}}}
static double _boundary13(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _dirichlet(0, point, neighbor, _s, data));}}}static double _boundary13_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _dirichlet_homogeneous(0, point, neighbor, _s, data));}}}
static double _boundary14(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _dirichlet(0, point, neighbor, _s, data));}}}static double _boundary14_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _dirichlet_homogeneous(0, point, neighbor, _s, data));}}}

static double _boundary15(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann(0, point, neighbor, _s, data));}}}static double _boundary15_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _neumann_homogeneous(0, point, neighbor, _s, data));}}}
static double _boundary16(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _dirichlet(0, point, neighbor, _s, data));}}}static double _boundary16_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _dirichlet_homogeneous(0, point, neighbor, _s, data));}}}
static double _boundary17(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _dirichlet(0, point, neighbor, _s, data));}}}static double _boundary17_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;{return( _dirichlet_homogeneous(0, point, neighbor, _s, data));}}}

double SDF (double x, double y) {

   return sq (x - ci.x) + sq (y - ci.y) - sq(D/2);
}

int j, T;
int main() {
#line 71
_init_solver();
  
#line 45
L0 = 12.8*D;
  ci.x = L0/2;
  ci.y = L0/2;
  size(L0);
  init_grid (2 << (LEVEL-4));
  mu = muv;
  TOLERANCE = 1.e-5;
  DT = 5e-4;

  T = 0;
  j = 0;
  RE = 100;
  KC = 5;
  run();

  T = 0;
  j = 1;
  RE = 200;
  KC = 10;
  L0 = 25.6*D;
  size(L0);
  LEVEL = 10;
  init_grid (2 << (LEVEL-4));
  ci.x = L0/2;
  ci.y = L0/2;
  run();
free_solver();

#line 71
}

static int init_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(t = 0)!=0;*ip=i;*tp=t;return ret;}


#line 73
      static int init_0(const int i,const double t,Event *_ev){tracing("init_0","cylinder-osc.c",0); {
  Amp = KC*D/(2*M_PI);
  freq = U0/(D*KC);
  foreach_stencil(1,{(NonLocal[]){{"u","vector",(void *)&u,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n    \n// #line 77 \"cylinder-osc.c\"\n{\n      val_out_(u.x,0,0,0) = 0.;\n      \n// #line 78\nval_out_(u.y,0,0,0) = 0.;}")})
    {
      {_stencil_val_a(u.x,0,0,0);  }
      
#line 78
{_stencil_val_a(u.y,0,0,0);  }}end_foreach_stencil();
  {
#line 76
foreach()
    {
      val(u.x,0,0,0) = 0.;
      
#line 78
val(u.y,0,0,0) = 0.;}end_foreach();}
}{end_tracing("init_0","cylinder-osc.c",0);return 0;}end_tracing("init_0","cylinder-osc.c",0);}

scalar  ref={15};
vector  fref={{16},{17}};
static int moving_cylinder_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}

#line 83
      static int moving_cylinder(const int i,const double t,Event *_ev){tracing("moving_cylinder","cylinder-osc.c",0); {
  vc.y = -Amp*2*M_PI*freq*cos(2*M_PI*freq*t);
  xc.y = -Amp*sin(2*M_PI*freq*t);
  do { scalar  phi=new_vertex_scalar("phi"); foreach_vertex_stencil(1,{(NonLocal[]){{"D","double",(void *)&D,NULL,0},{"xc","coord",(void *)&xc,NULL,0},{"ci","coord",(void *)&ci,NULL,0},{"phi","scalar",(void *)&phi,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_(" \n// #line 86 \"cylinder-osc.c\"\nval_out_(phi,0,0,0) = - sq(x - ci.x - xc.x) - sq(y - ci.y - xc.y) + sq(D/2);")}) {_stencil_val_a(phi,0,0,0);               }end_foreach_vertex_stencil(); {foreach_vertex() val(phi,0,0,0) = - sq(x - ci.x - xc.x) - sq(y - ci.y - xc.y) + sq(D/2);end_foreach_vertex();} fractions (phi, airfoil, sf
#line 122 "/home/spencer/basilisk/src/fractions.h"
, 0.
#line 86 "cylinder-osc.c"
);delete((scalar*)((scalar[]){phi,{-1}})); } while(0);
  do { scalar  phi=new_vertex_scalar("phi"); foreach_vertex_stencil(1,{(NonLocal[]){{"D","double",(void *)&D,NULL,0},{"xc","coord",(void *)&xc,NULL,0},{"ci","coord",(void *)&ci,NULL,0},{"phi","scalar",(void *)&phi,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_(" \n// #line 87 \"cylinder-osc.c\"\nval_out_(phi,0,0,0) = - sq(x - ci.x - xc.x) - sq(y - ci.y - xc.y) + sq(D/2);")}) {_stencil_val_a(phi,0,0,0);               }end_foreach_vertex_stencil(); {foreach_vertex() val(phi,0,0,0) = - sq(x - ci.x - xc.x) - sq(y - ci.y - xc.y) + sq(D/2);end_foreach_vertex();} fractions (phi, ref, fref
#line 122 "/home/spencer/basilisk/src/fractions.h"
, 0.
#line 87 "cylinder-osc.c"
);delete((scalar*)((scalar[]){phi,{-1}})); } while(0);

  coord nv;
  foreach_stencil(1,{(NonLocal[]){{"D","double",(void *)&D,NULL,0},{"ci","coord",(void *)&ci,NULL,0},{"nv","coord",(void *)&nv,NULL,0},{"airfoil","scalar",(void *)&airfoil,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_("\n\n\n// #line 38 \"cylinder-osc.c\"\nreal SDF (real x, real y) {\n\n   return sq (x - ci.x) + sq (y - ci.y) - sq(D/2);\n}\n\n\n// #line 48 \"./../immersed.h\"\nvoid normal_vector (Point point, scalar s, coord * nv) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;\n  coord n;\n  \n      n.x = (val(s,1,0,0) - val(s,-1,0,0))/(2.*Delta);\n      \n// #line 51\nn.y = (val(s,0,1,0) - val(s,0,-1,0))/(2.*Delta);\n  real mag = 0;\n  \n    mag += sq(n.x);\n    \n// #line 54\nmag += sq(n.y);\n  mag = sqrt(mag);\n  if (mag > 0)\n    {\n      n.x /= mag;\n      \n// #line 58\nn.y /= mag;}\n  *nv = n;\n\n}"),_(" \n// #line 90 \"cylinder-osc.c\"\n{\n    normal_vector (point, airfoil, &nv);\n    real lambda = fabs(nv.x) + fabs(nv.y);\n    if (lambda != 0.) {\n      real eta = 0.065*(1 - sq(lambda)) + 0.39;\n      real num = SDF (x,y);\n      real den = lambda * eta * sqrt(2) * Delta;\n      val_out_(airfoil,0,0,0) = den != 0? 0.5*(1 - tanh (num/den*0.5)): val(airfoil,0,0,0);\n    }\n  }")}) {
    _stencil_normal_vector (point, airfoil,NULL );     
     
{                    
      
      
       
_stencil_val(airfoil,0,0,0);
      
#line 97
_stencil_val_a(airfoil,0,0,0);        
    }
       
  
#line 99
}end_foreach_stencil();
  {
#line 90
foreach() {
    normal_vector (point, airfoil, &nv);
    double lambda = fabs(nv.x) + fabs(nv.y);
    if (lambda != 0.) {
      double eta = 0.065*(1 - sq(lambda)) + 0.39;
      double num = SDF (x,y);
      double den = lambda * eta * sqrt(2) * Delta;
      val(airfoil,0,0,0) = den != 0? 0.5*(1 - tanh (num/den*0.5)): val(airfoil,0,0,0);
    }
  }end_foreach();}
}{end_tracing("moving_cylinder","cylinder-osc.c",0);return 0;}end_tracing("moving_cylinder","cylinder-osc.c",0);}

static int properties_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}


#line 102
      static int properties_0(const int i,const double t,Event *_ev){tracing("properties_0","cylinder-osc.c",0); {
  if(!is_constant(fm.x)){
  
#line 103
foreach_face_stencil(1,{(NonLocal[]){{"RE","double",(void *)&RE,NULL,0},{"D","double",(void *)&D,NULL,0},{"U0","double",(void *)&U0,NULL,0},{"fm","vector",(void *)&fm,NULL,0},{"muv","vector",(void *)&muv,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 103 \"cylinder-osc.c\"\n{is_face_x(){\n    val_out_(muv.x,0,0,0) = val(fm.x,0,0,0)*(U0)*(D)/(RE);}end_is_face_x()\n// #line 103\nis_face_y(){\n    val_out_(muv.y,0,0,0) = val(fm.y,0,0,0)*(U0)*(D)/(RE);}end_is_face_y()}")}){_stencil_is_face_x(){
    { _stencil_val(fm.x,0,0,0);_stencil_val_a(muv.x,0,0,0); }}end__stencil_is_face_x()
#line 103
_stencil_is_face_y(){
    { _stencil_val(fm.y,0,0,0);_stencil_val_a(muv.y,0,0,0); }}end__stencil_is_face_y()}end_foreach_face_stencil();{
#line 103
foreach_face_generic(){is_face_x(){
    val(muv.x,0,0,0) = val(fm.x,0,0,0)*(U0)*(D)/(RE);}end_is_face_x()
#line 103
is_face_y(){
    val(muv.y,0,0,0) = val(fm.y,0,0,0)*(U0)*(D)/(RE);}end_is_face_y()}end_foreach_face_generic();}}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 103
foreach_face_stencil(1,{(NonLocal[]){{"RE","double",(void *)&RE,NULL,0},{"D","double",(void *)&D,NULL,0},{"U0","double",(void *)&U0,NULL,0},{"_const_fm","_coord",(void *)&_const_fm,NULL,0},{"muv","vector",(void *)&muv,NULL,0},{"N","int",(void *)&N,NULL,0},{"L0","double",(void *)&L0,NULL,0},{"Z0","double",(void *)&Z0,NULL,0},{"Y0","double",(void *)&Y0,NULL,0},{"X0","double",(void *)&X0,NULL,0},{0}},_(""),_("\n// #line 103 \"cylinder-osc.c\"\n{is_face_x(){\n    val_out_(muv.x,0,0,0) = _const_fm.x*(U0)*(D)/(RE);}end_is_face_x()\n// #line 103\nis_face_y(){\n    val_out_(muv.y,0,0,0) = _const_fm.y*(U0)*(D)/(RE);}end_is_face_y()}")}){_stencil_is_face_x(){
    {;_stencil_val_a(muv.x,0,0,0); }}end__stencil_is_face_x()
#line 103
_stencil_is_face_y(){
    {;_stencil_val_a(muv.y,0,0,0); }}end__stencil_is_face_y()}end_foreach_face_stencil();
  {
#line 103
foreach_face_generic(){is_face_x(){
    val(muv.x,0,0,0) = _const_fm.x*(U0)*(D)/(RE);}end_is_face_x()
#line 103
is_face_y(){
    val(muv.y,0,0,0) = _const_fm.y*(U0)*(D)/(RE);}end_is_face_y()}end_foreach_face_generic();}}
   boundary_internal ((scalar *)(scalar *)((vector[]) {muv,{{-1},{-1}}}), "cylinder-osc.c", 0);
}{end_tracing("properties_0","cylinder-osc.c",0);return 0;}end_tracing("properties_0","cylinder-osc.c",0);}

static int period_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(t += 1/freq)!=0;*ip=i;*tp=t;return ret;}


#line 108
      static int period(const int i,const double t,Event *_ev){tracing("period","cylinder-osc.c",0); {
  T++;
}{end_tracing("period","cylinder-osc.c",0);return 0;}end_tracing("period","cylinder-osc.c",0);}

static int logfile_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=( T <= T_max)!=0;*ip=i;*tp=t;return ret;}static int logfile_expr1(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}


#line 112
      static int logfile(const int i,const double t,Event *_ev){tracing("logfile","cylinder-osc.c",0);{
  coord Fp = {0};
  immersed_force (airfoil, sf, &Fp);
  double CD = (Fp.x)/(0.5*sq(U0)*D);
  double CL = (Fp.y)/(0.5*sq(U0)*D);

  fprintf (ferr, "%d %g %d %d %d %d %d %g %g %g %g %g %d\n",
    i, t, j, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax, CD, CL, xc.y, vc.y, dt, T);
}{end_tracing("logfile","cylinder-osc.c",0);return 0;}end_tracing("logfile","cylinder-osc.c",0);}

static int movie_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(t += 0.05)!=0;*ip=i;*tp=t;return ret;}


#line 122
      static int movie(const int i,const double t,Event *_ev){tracing("movie","cylinder-osc.c",0); {
  int FOV = j > 0? 5:10;
  scalar  omega=new_scalar("omega");
  vorticity (u, omega);
  char name[80];
  sprintf (name, "%d-snap-%g.png", j, t);
  FILE * fp1 = fopen(name, "w");
  view ( -0.5, -0.5, FOV
#line 50 "/home/spencer/basilisk/src/draw.h"
,
(
    
#line 51
float[4]) {0}, 
1., 1., 1.
#line 129 "cylinder-osc.c"
, 600, 675
#line 53 "/home/spencer/basilisk/src/draw.h"
, 4,
(
    
#line 54
float[3]) {0}, 
0., 0., 0., 
false, 
0., 0., 0., 
0., 
NULL, 
NULL, 
0, 
0., 0., 0., 0., 
NULL
#line 129 "cylinder-osc.c"
);
  draw_vof ("airfoil", "sf"
#line 886 "/home/spencer/basilisk/src/draw.h"
, false, 
0., 0, 
NULL, 
0, 0, 0, 
false, 
jet,
(
        
#line 892
float[3]) {0},( float[3]) {0}
#line 130 "cylinder-osc.c"
, 3
#line 892 "/home/spencer/basilisk/src/draw.h"
, 
false, 
false, 15,( float[2]) {-.95, -.95}, "", 1, false, false, false, 50, "%g", 50
#line 130 "cylinder-osc.c"
);
  squares ("omega"
#line 1221 "/home/spencer/basilisk/src/draw.h"
, 
NULL, 
0, 0, 0, 
false
#line 131 "cylinder-osc.c"
, cool_warm
#line 1225 "/home/spencer/basilisk/src/draw.h"
,
(
       
#line 1226
float[3]) {0},( float[3]) {0}, 
false,

(

       
#line 1229
coord) {0,0,1}, 
0, 
1, 
false, 15,( float[2]) {-.95, -.95}, "", 1, false, false, false, 50, "%g", 50
#line 131 "cylinder-osc.c"
);
  cells(
#line 1121 "/home/spencer/basilisk/src/draw.h"
(coord) {0,0,1}, 0.,
(
     
#line 1122
float[3]) {0}, 1.
#line 132 "cylinder-osc.c"
);
  save ( 
#line 495 "/home/spencer/basilisk/src/view.h"
NULL, "ppm", NULL
#line 133 "cylinder-osc.c"
, fp1
#line 496 "/home/spencer/basilisk/src/view.h"
, 
0, 
0, 0, 

NULL
#line 133 "cylinder-osc.c"
);
  fclose (fp1);delete((scalar*)((scalar[]){omega,{-1}}));
}{end_tracing("movie","cylinder-osc.c",0);return 0;}end_tracing("movie","cylinder-osc.c",0);}

static int adapt_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}


#line 137
      static int adapt_0(const int i,const double t,Event *_ev){tracing("adapt_0","cylinder-osc.c",0); {
  adapt_wavelet (
#line 138 "/home/spencer/basilisk/src/grid/tree-common.h"
(
#line 173
scalar *
#line 138
)
#line 138 "cylinder-osc.c"
((scalar[]){airfoil,u.x,u.y,{-1}}), (double[]){1.e-2,3e-3,3e-3}
, LEVEL, 2
#line 176 "/home/spencer/basilisk/src/grid/tree-common.h"
, 
all
#line 139 "cylinder-osc.c"
);
}{end_tracing("adapt_0","cylinder-osc.c",0);return 0;}end_tracing("adapt_0","cylinder-osc.c",0);}
# 2 "ast/init_solver.h"

static void _init_solver (void)
{
  void init_solver();
datasize=18*sizeof(double);
  
#line 6
init_solver();
  {
#line 24
quadtree_methods();
#line 796 "/home/spencer/basilisk/src/display.h"
display_init();

    /**
    Placeholder for last event registrations. Must be empty. */

    
#line 12 "ast/init_solver.h"
{
      
      /**
      Placeholder for field allocations. Must be empty. */
    
      {  
#line 42 "/home/spencer/basilisk/src/run.h"
event_register((Event){0,1,defaults,{defaults_expr0},((int *)0),((double *)0),"/home/spencer/basilisk/src/run.h",0,"defaults"});  
#line 126 "/home/spencer/basilisk/src/navier-stokes/centered.h"
event_register((Event){0,1,defaults_0,{defaults_0_expr0},((int *)0),((double *)0),"/home/spencer/basilisk/src/navier-stokes/centered.h",0,"defaults"});  
#line 187
event_register((Event){0,1,default_display,{default_display_expr0},((int *)0),((double *)0),"/home/spencer/basilisk/src/navier-stokes/centered.h",0,"default_display"});  








event_register((Event){0,1,init,{init_expr0},((int *)0),((double *)0),"/home/spencer/basilisk/src/navier-stokes/centered.h",0,"init"});  
#line 73 "cylinder-osc.c"
event_register((Event){0,1,init_0,{init_0_expr0},((int *)0),((double *)0),"cylinder-osc.c",0,"init"});  









event_register((Event){0,1,moving_cylinder,{moving_cylinder_expr0},((int *)0),((double *)0),"cylinder-osc.c",0,"moving_cylinder"});  
#line 108
event_register((Event){0,1,period,{period_expr0},((int *)0),((double *)0),"cylinder-osc.c",0,"period"});  



event_register((Event){0,2,logfile,{logfile_expr0,logfile_expr1},((int *)0),((double *)0),"cylinder-osc.c",0,"logfile"});  









event_register((Event){0,1,movie,{movie_expr0},((int *)0),((double *)0),"cylinder-osc.c",0,"movie"});
	
	/**
	Placeholder for event registrations. Must be empty. */
	
      
#line 22 "ast/init_solver.h"
}
#line 1257 "/home/spencer/basilisk/src/common.h"
init_const_vector((vector){{_NVARMAX+0},{_NVARMAX+1}},"zerof",(double[]){0.,0.,0.});
init_const_vector((vector){{_NVARMAX+2},{_NVARMAX+3}},"unityf",(double[]){1.,1.,1.});
init_const_scalar((scalar){_NVARMAX+4},"unity", 1.);
init_const_scalar((scalar){_NVARMAX+5},"zeroc", 0.);



init_const_vector((vector){{_NVARMAX+6},{_NVARMAX+7}},"unityf0",(double[]){1.,1.,1.});
init_const_scalar((scalar){_NVARMAX+8},"unity0", 1.);  init_scalar((scalar){0},"p");  init_vector((vector){{1},{2}},"u");  init_vector((vector){{3},{4}},"g");  init_scalar((scalar){5},"pf");  init_face_vector((vector){{6},{7}},"uf");  init_vector((vector){{8},{9}},"aF");  init_scalar((scalar){10},"airfoil");  init_face_vector((vector){{11},{12}},"sf");  init_face_vector((vector){{13},{14}},"muv");  init_scalar((scalar){15},"ref");  init_face_vector((vector){{16},{17}},"fref");
    
#line 23 "ast/init_solver.h"
}_attribute[p.i].dirty=1,_attribute[p.i].boundary[right]=_boundary0,_attribute[p.i].boundary_homogeneous[right]=_boundary0_homogeneous;_attribute[p.i].dirty=1,_attribute[p.i].boundary[left]=_boundary1,_attribute[p.i].boundary_homogeneous[left]=_boundary1_homogeneous;_attribute[p.i].dirty=1,_attribute[p.i].boundary[top]=_boundary2,_attribute[p.i].boundary_homogeneous[top]=_boundary2_homogeneous;_attribute[p.i].dirty=1,_attribute[p.i].boundary[bottom]=_boundary3,_attribute[p.i].boundary_homogeneous[bottom]=_boundary3_homogeneous;_attribute[u.x.i].dirty=1,_attribute[u.x.i].boundary[left]=_boundary4,_attribute[u.x.i].boundary_homogeneous[left]=_boundary4_homogeneous;_attribute[u.y.i].dirty=1,_attribute[u.y.i].boundary[left]=_boundary5,_attribute[u.y.i].boundary_homogeneous[left]=_boundary5_homogeneous;_attribute[p.i].dirty=1,_attribute[p.i].boundary[left]=_boundary6,_attribute[p.i].boundary_homogeneous[left]=_boundary6_homogeneous;_attribute[pf.i].dirty=1,_attribute[pf.i].boundary[left]=_boundary7,_attribute[pf.i].boundary_homogeneous[left]=_boundary7_homogeneous;_attribute[u.x.i].dirty=1,_attribute[u.x.i].boundary[right]=_boundary8,_attribute[u.x.i].boundary_homogeneous[right]=_boundary8_homogeneous;_attribute[u.y.i].dirty=1,_attribute[u.y.i].boundary[right]=_boundary9,_attribute[u.y.i].boundary_homogeneous[right]=_boundary9_homogeneous;_attribute[p.i].dirty=1,_attribute[p.i].boundary[right]=_boundary10,_attribute[p.i].boundary_homogeneous[right]=_boundary10_homogeneous;_attribute[pf.i].dirty=1,_attribute[pf.i].boundary[right]=_boundary11,_attribute[pf.i].boundary_homogeneous[right]=_boundary11_homogeneous;_attribute[u.x.i].dirty=1,_attribute[u.x.i].boundary[top]=_boundary12,_attribute[u.x.i].boundary_homogeneous[top]=_boundary12_homogeneous;_attribute[p.i].dirty=1,_attribute[p.i].boundary[top]=_boundary13,_attribute[p.i].boundary_homogeneous[top]=_boundary13_homogeneous;_attribute[pf.i].dirty=1,_attribute[pf.i].boundary[top]=_boundary14,_attribute[pf.i].boundary_homogeneous[top]=_boundary14_homogeneous;_attribute[u.x.i].dirty=1,_attribute[u.x.i].boundary[bottom]=_boundary15,_attribute[u.x.i].boundary_homogeneous[bottom]=_boundary15_homogeneous;_attribute[p.i].dirty=1,_attribute[p.i].boundary[bottom]=_boundary16,_attribute[p.i].boundary_homogeneous[bottom]=_boundary16_homogeneous;_attribute[pf.i].dirty=1,_attribute[pf.i].boundary[bottom]=_boundary17,_attribute[pf.i].boundary_homogeneous[bottom]=_boundary17_homogeneous;  
#line 798 "/home/spencer/basilisk/src/display.h"
event_register((Event){0,1,refresh_display,{refresh_display_expr0},((int *)0),((double *)0),"/home/spencer/basilisk/src/display.h",0,"refresh_display"});  
#line 50 "/home/spencer/basilisk/src/run.h"
event_register((Event){0,1,cleanup,{cleanup_expr0},((int *)0),((double *)0),"/home/spencer/basilisk/src/run.h",0,"cleanup"});  
#line 222 "/home/spencer/basilisk/src/navier-stokes/centered.h"
event_register((Event){0,1,set_dtmax,{set_dtmax_expr0},((int *)0),((double *)0),"/home/spencer/basilisk/src/navier-stokes/centered.h",0,"set_dtmax"});  

event_register((Event){0,1,stability,{stability_expr0},((int *)0),((double *)0),"/home/spencer/basilisk/src/navier-stokes/centered.h",0,"stability"});  









event_register((Event){0,1,vof,{vof_expr0},((int *)0),((double *)0),"/home/spencer/basilisk/src/navier-stokes/centered.h",0,"vof"});  
event_register((Event){0,1,tracer_advection,{tracer_advection_expr0},((int *)0),((double *)0),"/home/spencer/basilisk/src/navier-stokes/centered.h",0,"tracer_advection"});  
event_register((Event){0,1,tracer_diffusion,{tracer_diffusion_expr0},((int *)0),((double *)0),"/home/spencer/basilisk/src/navier-stokes/centered.h",0,"tracer_diffusion"});  






event_register((Event){0,1,properties,{properties_expr0},((int *)0),((double *)0),"/home/spencer/basilisk/src/navier-stokes/centered.h",0,"properties"});  
#line 316
event_register((Event){0,1,advection_term,{advection_term_expr0},((int *)0),((double *)0),"/home/spencer/basilisk/src/navier-stokes/centered.h",0,"advection_term"});  
#line 345
event_register((Event){0,1,viscous_term,{viscous_term_expr0},((int *)0),((double *)0),"/home/spencer/basilisk/src/navier-stokes/centered.h",0,"viscous_term"});  
#line 381
event_register((Event){0,1,acceleration,{acceleration_expr0},((int *)0),((double *)0),"/home/spencer/basilisk/src/navier-stokes/centered.h",0,"acceleration"});  
#line 421
event_register((Event){0,1,projection,{projection_expr0},((int *)0),((double *)0),"/home/spencer/basilisk/src/navier-stokes/centered.h",0,"projection"});  
#line 436
event_register((Event){0,1,end_timestep,{end_timestep_expr0},((int *)0),((double *)0),"/home/spencer/basilisk/src/navier-stokes/centered.h",0,"end_timestep"});  









event_register((Event){0,1,adapt,{adapt_expr0},((int *)0),((double *)0),"/home/spencer/basilisk/src/navier-stokes/centered.h",0,"adapt"});  
#line 12 "./../immersed.h"
event_register((Event){0,1,end_timestep_0,{end_timestep_0_expr0},((int *)0),((double *)0),"./../immersed.h",0,"end_timestep"});  
#line 102 "cylinder-osc.c"
event_register((Event){0,1,properties_0,{properties_0_expr0},((int *)0),((double *)0),"cylinder-osc.c",0,"properties"});  
#line 137
event_register((Event){0,1,adapt_0,{adapt_0_expr0},((int *)0),((double *)0),"cylinder-osc.c",0,"adapt"});
  
#line 24 "ast/init_solver.h"
}
  set_fpe();
}
