#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <signal.h>
#include <setjmp.h>
#include <math.h>
#include <sys/types.h>
#ifdef linux
#include <sys/wait.h>
#include <unistd.h>
#include <string.h>
#endif
#include <sys/stat.h>
#include <fcntl.h>

#include "pswarm_main.h"
#include "pswarm.h"

typedef struct { char *msg; } Exit_code;
extern struct Stats stats;

static Exit_code exit_codes[] = {
    {"Normal exit"},                    /* 0 */
    {"Abnormal exit"},                  /* 1 */
    {"Failed to allocate memory"},      /* 2 */
    {"Unable to initalize population"}, /* 3 */
};

extern void *pswarm_malloc(size_t size);

void objfun(int, int, double *, double *, double *, double *);

extern int PSwarm(int, void (*)(), double *, double *, int, double *, double *,
                  double **, double *, double *);
extern struct Options opt;

extern void save_cache_file(int n, int m);
extern void load_cache_file(int n, int m);
extern int read_cesam_cache(int n, double *x, double *age, double *teff,
                            double *lum, double *r);
extern void write_cesam_cache(int n, double *x, double *age, double *teff,
                              double *lum, double *r);

extern void set_problem(double *x, double *lb, double *ub, double *A,
                        double *b);
extern void set_problem_dimension(int *n, int *lincons);

extern void user_init();

static jmp_buf Jb;

void catchfpe(int n) {
  printf("\nFloating point error.\n");
  fflush(stdout);
  longjmp(Jb, 1);
}

/************************************************
                    Main
*************************************************/
int main(int argc, char **argv) {
  int exit_code, i;
  double *sol = NULL;
  double *f = NULL;
  double *lb = NULL;
  double *ub = NULL;
  double *A = NULL;
  double *b = NULL; /* linear constraints defined in a Fortran way */
  int lincons = n_cons();   /* number of linear constraints */

  int n = n_dims();
  double *X0;

  if (!setjmp(Jb)) {
    signal(SIGFPE, catchfpe);

    /* lower and upper bounds on variables */
    lb = pswarm_malloc(n * sizeof(double));

    ub = pswarm_malloc(n * sizeof(double));

    if (!lb || !ub) {
      printf("Unable to allocate memory for variable bounds\n");
      exit(1);
    }

    X0 = pswarm_malloc(n * sizeof(double));
    if (!X0) {
      printf("Unable to allocate memory for initial guess\n");
      exit(1);
    }

    if (lincons) {
      A = pswarm_malloc(lincons * n * sizeof(double));
      b = pswarm_malloc(lincons * sizeof(double));

      memset(A, 0, (lincons)*n * sizeof(double));
      memset(b, 0, (lincons) * sizeof(double));
    } else {
      A = b = NULL;
    }

    set_problem(X0, lb, ub, A, b);

    /* check for finite bound on the variables */
    /*    for(i=0;i<n;i++){
          if(lb[i]<=-Inf || ub[i]>=Inf){
        printf("Variable %i has no finite bounds\nPlease provide finite bounds
       for all variables\n",i);
        exit(1);
          }
        }*/

    user_init();
    //      load_cache_file(n, 4);
    exit_code = PSwarm(n, &objfun, lb, ub, lincons, A, b, &sol, f, X0);
    //      save_cache_file(n, 4);
    if (opt.IPrint >= 0) printf("\n%s\n", exit_codes[exit_code].msg);
  }

  if (lb) free(lb);
  if (ub) free(ub);
  if (A) free(A);
  if (b) free(b);

  return exit_code;
}
