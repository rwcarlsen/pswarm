#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdio.h>
#include <math.h>
#ifdef linux
#include <strings.h>
#endif
#include <signal.h>
#include <setjmp.h>

#include "pswarm.h"

double r_objfun(int n, double *x, double *lb, double *ub);
double py_outfcn(int n, int s, int iter, int gbest, struct swarm *pop);

static PyObject *py_objf, *py_outf, *Problem, *Options;

extern struct Options opt;
extern struct Stats stats;
extern int PSwarm(int n, void (*objf)(), double *lb, double *ub, int lincons,
                  double *A, double *b, double **sol, double *f, double *x);

static PyObject *pswarm_py(PyObject *self, PyObject *args);

static int saved_options = 0;
static struct Options opt_backup;
