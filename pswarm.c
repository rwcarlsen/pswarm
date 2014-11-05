#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <memory.h>

#include "pattern.h"
#include "pswarm.h"

#ifdef LINEAR
extern void dgemv_();
extern void dpotrf_();
extern double dnrm2_();
extern int mve_presolve(int m, int n, double *A, double *b, int maxiter,
                        double tol, double *x);
extern int mve_solver(int m, int n, double *A, double *b, double *x0,
                      double maxiter, double tol, double *x, double *E2);

void check_feasible_pop(int n, int s, int lincons, struct swarm *pop,
                        double *lb, double *ub, double *A, double *b);
#endif

int feasible_p(int n, double *x, int lincons, double *A, double *b, double *lb,
               double *ub);
void print_pop(int, int, int, struct swarm *);
int init_pop(int, int, int, struct swarm *, double *, double *, double *,
             double *, int, double *);
double projection(double, double, double);
void matlab_write_pop(int n, int gbest, int s, struct swarm *pop, int iter);
void print_best(int n, int gbest, int s, struct swarm *pop, int iter);
void print_array(int n, double *x);
double outfcn(int n, int s, int iter, int gbest, struct swarm *pop);
void *pswarm_malloc(size_t size);

extern void pollstep(int n, int lincons, int pi, void (*objf)(), double *lb,
                     double *ub, double *A, double *b,
                     struct poll_vector **last_sucess);
extern void init_pattern(int);
extern void clean_pattern();

#if SYS_RANDOM != 1
static long rand_seed;

#define SHUFFLE 256 /* size of random array */
#define DBL_MIN 2.2250738585072014e-308

static long rand_seed;
static double randflt(long *);
static double resettable_randflt(long *rand_seed, int reset);
#endif

struct swarm pop;
struct Stats stats;

/* default options */
struct Options opt = {
    42,      /* swarm size */
    0.5,     /* cognitial parameter */
    0.5,     /* social parameter */
    0.5,     /* maximum velocity factor */
    2000,    /* maximum of iterations */
    2000,    /* maximum of function evaluations */
    0.9,     /* initial weight */
    0.4,     /* final weight */
    0.5,     /* max norm 2 for gradient */
    10,      /* bound limit */
    1.0e-5,  /* tolerance for stopping criteria */
    Inf,     /* initial delta -- is computed or user provided */
    5.0,     /* factor for initial delta */
    2,       /* increase delta by a factor of */
    0.5,     /* decrease delta by a factor of */
    N2,      /* type of basis on pattern search */
    0.1,     /* active constraints epsilon */
    10,      /* IPrint - print info each IPrint iterations */
    &outfcn, /* outputfcn - function to call for printing */
    1        /* vectorized call to the objective function */
};

/* Pattern Swarm algorithm */
int PSwarm(int n, void (*objf)(), double *lb, double *ub, int lincons,
           double *A, double *b, double **sol, double *f, double *x) {
  int i, j, iter, gbest, success, actives, iterunsuc = 0, process;
  double *maxv, maxnormv, normtmp, weight, mindelta, normtmp2;
  double *vectorx, *vectorfx;
  char *buff;
  time_t tt;
  double *AlphaMax;

#ifdef LINEAR
  double AlphaMaxLinCons;
  double one = 1.0;
  double none = -1.0;
  double zero = 0.0;
  int oneI = 1;
  char Normal = 'N';

  double *velocity;
  double *tmpm1;
  double *tmpm2;
#endif

  static struct poll_vector *last_success = NULL;

  /* Initial time */
  time(&tt);
  buff = ctime(&tt);
  if (opt.IPrint >= 0) printf("\nInitial time: %s\n", buff);

#if SYS_RANDOM == 1
  /* seed random number with time */
  srand((unsigned int)tt);
#else
  rand_seed = (long)abs(tt);
  resettable_randflt(&rand_seed, 1); /* initialize random number generator */
  /* seed random number with time */
  srand((unsigned int)tt);
#endif

  if (n <= 0) {
    printf("Number of variables must be positive\n");
    return EXIT_KO;
  }

  AlphaMax = pswarm_malloc(n * sizeof(double));
#ifdef LINEAR
  velocity = pswarm_malloc(n * sizeof(double));
  tmpm1 = pswarm_malloc(lincons * sizeof(double));
  tmpm2 = pswarm_malloc(lincons * sizeof(double));
#endif

  if (!ub || !lb) {
    printf("Lower and upper bounds must be defined\n");
    return EXIT_KO;
  }

/* checkup on variables bounds */
/* we will be adding ficticious bounds whenever necessary */
/*  for(i=0;i<n;i++){
if((lb[i]<=-Inf && !x) || (ub[i]>=Inf && !x)){
printf("Not all variables have finite bound and no initial guess given\n");
printf("All variables must have finite simple bounds or an initial guess should
be provided\n");

return EXIT_INITIAL;
}
}*/

/* allocate memory for swarm */

#ifndef LINEAR
  if (lincons > 0 && opt.IPrint >= 0) {
    printf(
        "\n**** Linear constraints defined but PSwarm was compiled without "
        "linear constraints support\n");
    printf("**** Ignoring linear constraints\n\n");
    lincons = 0;
  }
#endif

  pop.x = pswarm_malloc(opt.s * n * sizeof(double));
  pop.v = pswarm_malloc(opt.s * n * sizeof(double));
  pop.y = pswarm_malloc(opt.s * n * sizeof(double));
  pop.fx = pswarm_malloc(opt.s * sizeof(double));
  pop.fy = pswarm_malloc(opt.s * sizeof(double));
  pop.active = pswarm_malloc(opt.s * sizeof(int));

  /* allocate memory for maximum velocity allowed */
  maxv = pswarm_malloc(n * sizeof(double));

  /*	if(!pop.x || !pop.v || !pop.y || !pop.fx || !pop.fy || !pop.active || !
     maxv)
                  return EXIT_MEM; */

  pop.scale = 1.0; /* no scale */
  /* initialize maximum velocity allowed and compute delta. Compute scale for
   * active constraints */
  if (opt.delta >= Inf) { /* User didn't provide delta */
    mindelta = Inf;
    for (j = 0; j < n; j++) {
      if (lb[j] > -Inf && ub[j] < Inf) {
        if (pop.scale < ub[j] - lb[j]) pop.scale = ub[j] - lb[j];
        if (mindelta > (ub[j] - lb[j])) mindelta = (ub[j] - lb[j]);
        maxv[j] = (ub[j] - lb[j]) * opt.maxvfactor;
      } else {
        maxv[j] = Inf;
      }
    }
    if (mindelta >= Inf || mindelta < 2 * sqrt(opt.tol))
      opt.delta = 2 * sqrt(sqrt(opt.tol));
    else
      opt.delta = mindelta / opt.fdelta;
  }

  if (opt.IPrint >= 0) printf("Delta for pattern search: %f\n", opt.delta);

#ifdef LINEAR
  /* update scale with constraints */
  for (i = 0; i < lincons; i++)
    if (pop.scale < fabs(b[i])) pop.scale = fabs(b[i]);
#endif

  if (opt.IPrint >= 0) printf("Population scale: %f\n", pop.scale);

  /* initialize population */
  if (x) {
    if (opt.IPrint >= 0)
      printf("Initial guess provided, including in initial population\n\n");
    if (init_pop(n, opt.s, lincons, &pop, lb, ub, A, b, 1, x)) {
      printf("Unable to initialize population\n");

      free(AlphaMax);
#ifdef LINEAR
      free(velocity);
      free(tmpm1);
      free(tmpm2);
#endif

      return EXIT_INITIAL;
    }
  } else {
    if (init_pop(n, opt.s, lincons, &pop, lb, ub, A, b, 0, NULL)) {
      printf("Unable to initialize population\n");

      free(AlphaMax);
#ifdef LINEAR
      free(velocity);
      free(tmpm1);
      free(tmpm2);
#endif

      return EXIT_INITIAL;
    }
  }

#ifdef LINEAR
// printf("Initial Population check\n");
// check_feasible_pop(n, opt.s, lincons, &pop, lb, ub, A, b);
#endif

  actives = opt.s;

  iter = 0;
  stats.pollsteps = 0;
  stats.sucpollsteps = 0;
  stats.objfunctions = 0;
  gbest = 0;      /* global best */
  maxnormv = Inf; /* don't stop in first iteration */

  init_pattern(n); /* initialize pattern search */

  /* Main cycle */
  while (iter < opt.maxiter && stats.objfunctions < opt.maxf) {
    if (opt.IPrint > 0 && (iter == 0 || iter % opt.IPrint == 0))
      if (opt.outfcn(n, opt.s, iter, gbest, &pop) < 0) {
        printf("User requested to stop\n");
        break;
      }

    if (maxnormv < opt.tol && pop.delta < opt.tol) {
      if (opt.IPrint >= 0)
        printf("\n\nStopping due to velocity and tolerance\n\n");
      break;
    }

    if (actives <= 1 && pop.delta < opt.tol) {
      if (opt.IPrint >= 0)
        printf("\n\nStopping due to single particle and tolerance\n\n");
      break;
    }

    iter++;

    //		check_feasible_pop(n, opt.s, lincons, &pop, lb, ub, A, b);

    success = 0; /* controls if gbest was updated with success */

    if (opt.vectorized) { /* call objf once with all the points */

      vectorx = (double *)pswarm_malloc(
          opt.s * n *
          sizeof(
              double)); /* we could avoid to wast memory on points not active */
      vectorfx = (double *)pswarm_malloc(opt.s * sizeof(double));

      for (j = 0, i = 0; i < opt.s; i++)
        if (pop.active[i] &&
            feasible_p(n, &pop.x[i * n], lincons, A, b, lb, ub)) {
          memcpy(&vectorx[j * n], &pop.x[i * n], n * sizeof(double));
          j++;
        }

      objf(n, j, vectorx, lb, ub, vectorfx);
      stats.objfunctions += j;

      for (j = 0, i = 0; i < opt.s;
           i++) /* we could avoid a second cycle if we saved the indices */
        if (pop.active[i] &&
            feasible_p(n, &pop.x[i * n], lincons, A, b, lb, ub)) {
          pop.fx[i] = vectorfx[j];
          j++;
        }

      free(vectorx);
      free(vectorfx);
    } else {
      for (i = 0; i < opt.s; i++) {
        if (pop.active[i]) {
          if (feasible_p(n, &pop.x[i * n], lincons, A, b, lb, ub)) {
            objf(n, 1, &pop.x[i * n], lb, ub, &pop.fx[i]);
            stats.objfunctions++;
          } else {
            pop.fx[i] = +Inf;
          }
        }
      }
    }

    process = 0;
    for (i = 0; i < opt.s; i++) {
      if (pop.active[i]) {
        if (pop.fy[i] > pop.fx[i]) { /* progress obtained */
          pop.fy[i] = pop.fx[i];     /* Search step */
          memcpy(&pop.y[i * n], &pop.x[i * n], n * sizeof(double));

          /* check if a new leader is available or if a progress was
          obtained on the leader */
          if (pop.fy[gbest] > pop.fy[i] || gbest == i) {
            gbest = i;           /* global best indice */
            success = 1;         /* success for leader obtained */
            last_success = NULL; /* reset successful direction on poll step */
          }
        }
      }
    }

    if (!success) { /* no success for the gbest particle in one generation, so
                       performe a poll step */
      if (pop.delta >= opt.tol) {
        /* performe a poll step, update y and delta */
        pollstep(n, lincons, gbest, objf, lb, ub, A, b, &last_success);

        stats.pollsteps++;
        iterunsuc = 0;
      } else {
        iterunsuc++;
        // printf("Consecutive unsuccesseful iterations %d\n", iterunsuc);
      }
    } else { /* success for the gbest particle */
      iterunsuc = 0;
      // printf("Success in Swarm iteration\n");
      /* increase delta */
      if (pop.delta < opt.delta) {
        pop.delta *= opt.idelta;
        //  printf("Increasing delta in search step\n");
      }
      /* allow at least one more iteration */
      if (pop.delta < opt.tol) pop.delta = 2 * opt.tol;
    }

    /* inertia factor is a linear interpolation from iweight to fweight */
    weight =
        opt.iweight -
        (opt.iweight - opt.fweight) * ((double)(iter)) / ((double)opt.maxiter);

    // printf("Before step computation\n");
    // check_feasible_pop(n, opt.s, lincons, &pop, lb, ub, A, b);
    // printf("Before step computation\n");

    for (i = 0; i < opt.s; i++) { /* for each particle */

      if (pop.active[i]) { /* active particle */
        /* update velocity */
        for (j = 0; j < n; j++) {
          pop.v[i * n + j] =
#if SYS_RANDOM == 1
              projection(weight * pop.v[i * n + j] +
                             opt.mu * (rand() / (RAND_MAX + 1.0)) *
                                 (pop.y[i * n + j] - pop.x[i * n + j]) +
                             opt.nu * (rand() / (RAND_MAX + 1.0)) *
                                 (pop.y[gbest * n + j] - pop.x[i * n + j]),
                         -maxv[j], maxv[j]);
#else
              projection(weight * pop.v[i * n + j] +
                             opt.mu * (randflt(&rand_seed)) *
                                 (pop.y[i * n + j] - pop.x[i * n + j]) +
                             opt.nu * (randflt(&rand_seed)) *
                                 (pop.y[gbest * n + j] - pop.x[i * n + j]),
                         -maxv[j], maxv[j]);
#endif

          AlphaMax[j] = 1.0; /* a step no longer than 1 */
        }

        for (j = 0; j < n; j++) {
          if (pop.v[i * n + j] < 0.0) {
            if (AlphaMax[j] > (lb[j] - pop.x[i * n + j]) / pop.v[i * n + j])
              AlphaMax[j] = (lb[j] - pop.x[i * n + j]) / pop.v[i * n + j];
          }
          if (pop.v[i * n + j] > 0.0) {
            if (AlphaMax[j] > (ub[j] - pop.x[i * n + j]) / pop.v[i * n + j])
              AlphaMax[j] = (ub[j] - pop.x[i * n + j]) / pop.v[i * n + j];
          }
        }

        for (j = 0; j < n; j++)
          if (AlphaMax[j] < 0.0) {
            // printf("Também não deveria de acontecer\n");
            AlphaMax[j] = 0.0;
          }

#ifdef LINEAR
        /* account for linear constraints */
        AlphaMaxLinCons = 1.0;

        for (j = 0; j < n; j++) velocity[j] = AlphaMax[j] * pop.v[i * n + j];

        if (lincons) {
          /* -Ax */
          dgemv_(&Normal, &lincons, &n, &none, A, &lincons, &pop.x[i * n],
                 &oneI, &zero, tmpm1, &oneI);

          /* Av */
          dgemv_(&Normal, &lincons, &n, &one, A, &lincons, velocity, &oneI,
                 &zero, tmpm2, &oneI);

          for (j = 0; j < lincons; j++) {
            if (tmpm2[j] > 0 &&
                ((b[j] + tmpm1[j]) / tmpm2[j]) < AlphaMaxLinCons)
              AlphaMaxLinCons = ((b[j] + tmpm1[j]) / tmpm2[j]);
          }
        }
        // printf("AlphaMax=%.5f\n",AlphaMaxLinCons);
        if (AlphaMaxLinCons > 0.0) {
          for (j = 0; j < n; j++) {
            /* update particle and check bound limits */
            pop.x[i * n + j] = projection(
                pop.x[i * n + j] + AlphaMaxLinCons * velocity[j], lb[j], ub[j]);
          }

        } /* AlphaMax is zero and particle is not updated */

#else
        /* update particle and check bound limits */
        for (j = 0; j < n; j++) {
          pop.x[i * n + j] = projection(
              pop.x[i * n + j] + AlphaMax[j] * pop.v[i * n + j], lb[j], ub[j]);
        }
#endif
      }
    }

    // printf("After step computation\n");
    // check_feasible_pop(n, opt.s, lincons, &pop, lb, ub, A, b);
    // printf("After step computation\n");

    /* check for all norm velocities to zero */

    /* first for gbest */
    normtmp = 0.0;
    for (j = 0; j < n; j++) normtmp += pow(pop.v[gbest * n + j], 2.0);
    maxnormv = sqrt(normtmp);

    /* remove particle close to gbest and compute maximum velocity */
    actives = 0;
    for (i = 0; i < opt.s; i++) {        /* for each particle */
      if (pop.active[i] && i != gbest) { /* active particle and not the gbest */
        normtmp = 0.0;
        normtmp2 = 0.0;
        for (j = 0; j < n; j++) {
          normtmp += pow(pop.y[i * n + j] - pop.y[gbest * n + j], 2.0);
          normtmp2 += pow(pop.v[i * n + j], 2.0);
        }
        normtmp = sqrt(normtmp);
        normtmp2 = sqrt(normtmp2);
        if (normtmp < opt.delta &&
            normtmp2 <
                opt.delta) {  //(fabs((double)(iter-(opt.maxiter)/100.0)))){
          pop.active[i]--;    /* remove particle */
          // printf("Particle %d inactive iter=%d\n", i, iter);
        } else { /* particle not removed, so account for maxnormv */
          if (maxnormv < normtmp2) maxnormv = normtmp2;
        }
      }
      if (pop.active[i]) actives++; /* account in actives */
    }

//    printf("Maximum velocity norm: %f\n", maxnormv);

// printf("%d;%.20f\n",stats.objfunctions,pop.fy[gbest]);

#ifdef LINEAR
//  check_feasible_pop(n, opt.s, lincons, &pop, lb, ub, A, b);
#endif
  }

  if (opt.IPrint >= 0) {
    if (opt.IPrint > 0 && iter != 0 && iter % opt.IPrint != 0)
      opt.outfcn(n, opt.s, iter, gbest, &pop);

    if (iter >= opt.maxiter)
      printf("\n\nStopping due to maximum number of iterations reached\n\n");

    if (stats.objfunctions >= opt.maxf)
      printf(
          "\n\nStopping due to maximum number of function evaluations "
          "reached\n\n");

    // print_pop(n, gbest, opt.s, &pop);
    // print_best(n, gbest, opt.s, &pop, iter);

    // matlab_write_pop(n, gbest, opt.s, &pop, iter);

    printf("maxnormv=%.20f\n", maxnormv);
    printf("delta=%.20f\n", pop.delta);
    printf("%d iterations\n", iter);
    printf("%d function evaluations\n", stats.objfunctions);
    printf("%d poll steps performed\n", stats.pollsteps);
    printf("%d poll steps performed with success\n", stats.sucpollsteps);
    printf("%d & %d & %d & %d & %.4f\n", iter, stats.objfunctions,
           stats.pollsteps, stats.sucpollsteps, pop.fy[gbest]);

    /* some printf can be done here */
  };

  clean_pattern();

  time(&tt);
  buff = ctime(&tt);
  if (opt.IPrint >= 0) printf("Final time: %s\n", buff);

  free(AlphaMax);
#ifdef LINEAR
  free(velocity);
  free(tmpm1);
  free(tmpm2);
#endif

  /* returning the best of the population */
  *sol = pswarm_malloc(n * sizeof(double));
  if (!(*sol)) return EXIT_MEM;

  memcpy(*sol, &pop.y[gbest * n], n * sizeof(double));

  if (f) *f = pop.fy[gbest];

  /* free allocated memory */
  free(pop.x);
  free(pop.y);
  free(pop.v);
  free(pop.fx);
  free(pop.fy);
  free(pop.active);
  free(maxv);

  return EXIT_OK;
}

int init_pop(int n, int s, int lincons, struct swarm *pop, double *lb,
             double *ub, double *A, double *b, int ninitials,
             double *initials) {
  int i, j;
  double normtmp = 10.0; /* should never be used by default */
  int feas, accepted;
#ifdef LINEAR
  double one = 1.0;
  int oneI = 1;
  double zero = 0.0;
  char Normal = 'N';
  char Upper = 'U';
  double maxub, minlb;
  double tmpptnorm;
  double randscale;
  int msg, info;
  int reallincons;
  int tmplincons;
  int recovered = 0;

  double *Ai;
  double *Aext;
  double *Aextext;
  double *bext;
  double *center_ini;
  double *center;
  double *E2;
  double *E;
  double *tmppt;
  double *E2tmppt;
#endif

  /*   The user can provide an initial guess to include in the initial swarm
  A reset in the population can also proposed a fixed number of point
  to be in the next swarm */

  /* Do simple check in the simple bound limits */

  if (ninitials > s) {
    printf("Populations should be increased to %d particles\n", ninitials);
    ninitials = s;
  }

  accepted = 0;

  if (ninitials && initials) {
    for (i = 0; i < ninitials; i++) {
      // printf("Init %d: %f\n", i, objfun(n,&initials[i*n]));
      for (j = 0; j < n; j++)
        pop->x[accepted * n + j] =
            projection(initials[i * n + j], lb[j], ub[j]);
      feas = 1;

#ifdef LINEAR
      if (lincons > 0) {
        /* check for linear feasibility */
        Ai = pswarm_malloc(lincons * sizeof(double));
        dgemv_(&Normal, &lincons, &n, &one, A, &lincons, &initials[i * n],
               &oneI, &zero, Ai, &oneI);

        for (j = 0; j < lincons; j++)
          if (Ai[j] > b[j]) feas = 0;
        free(Ai);
      }
#endif
      if (feas) {
        pop->fy[accepted] = +Inf * 10; /* in first iteration, y will be set */
        pop->active[accepted] = 1;     /* chances to be near the gbest */
        accepted++;
      }
    }
  }

  if (accepted) {
    /* compute standard deviation of first particle */
    normtmp = 0.0;
    for (j = 0; j < n; j++) normtmp += pow(pop->x[j], 2.0);
    if (normtmp < 10) normtmp = opt.blim;
  }

  if (accepted < ninitials)
    printf("Only %d particles were linear feasible after projection\n",
           accepted);

#ifdef LINEAR

  if (lincons > 0) {
    if (opt.IPrint >= 0) printf("Generating ellipsoid\n");

    /* compute ellipsoide with maximum volume */
    Aext = pswarm_malloc((lincons + 2 * n) * n * sizeof(double));
    bext = pswarm_malloc((lincons + 2 * n) * sizeof(double));
    center_ini = pswarm_malloc(n * sizeof(double));
    center = pswarm_malloc(n * sizeof(double));

    memset(Aext, 0, (lincons + 2 * n) * n * sizeof(double));
    memset(bext, 0, (lincons + 2 * n) * sizeof(double));
    memcpy(bext, b, lincons * sizeof(double));

    maxub = -Inf;
    minlb = Inf;

    reallincons = lincons;
    for (i = 0; i < n; i++) {
      if (ub[i] < Inf) {
        if (ub[i] > maxub) maxub = ub[i];
        reallincons++;
      }
      if (lb[i] > -Inf) {
        if (lb[i] < minlb) minlb = lb[i];
        reallincons++;
      }
    }

    for (i = 0; i < lincons; i++) {
      for (j = 0; j < n; j++) {
        Aext[i + j * reallincons] = A[i + j * lincons];
      }
    }

    tmplincons = 0;

    for (i = 0; i < n; i++) {
      if (ub[i] < Inf) {
        Aext[lincons + tmplincons + i * reallincons] =
            1.0; /* remaining were set to zero */
        bext[lincons + tmplincons] = ub[i];
        tmplincons++;
      }
      if (lb[i] > -Inf) {
        Aext[lincons + tmplincons + i * reallincons] =
            -1.0; /* remaining were set to zero */
        bext[lincons + tmplincons] = -lb[i];
        tmplincons++;
      }
    }

    /* solve for initial guess */
    msg = mve_presolve(reallincons, n, Aext, bext, 80, 1e-8, center_ini);
    if (msg) {
      printf("MVE error %d\n", msg);
      printf("Trying to recover by adding fictitious bounds\n");

      /* add fictitious bounds */

      if (reallincons >= lincons + 2 * n) {
        printf(
            "MVE cannot recover. All bounds already considered!\nAborting\n");

        free(Aext);
        free(bext);
        free(center_ini);
        free(center);

        return 1;
        ;
      }

      Aextext = pswarm_malloc((lincons + 2 * n) * n * sizeof(double));

      for (i = 0; i < reallincons; i++) {
        for (j = 0; j < n; j++) {
          Aextext[i + j * (lincons + 2 * n)] = Aext[i + j * reallincons];
        }
      }

      tmplincons = 0;

      for (i = 0; i < n; i++) {
        if (ub[i] >= Inf) {
          Aextext[reallincons + tmplincons + i * (lincons + 2 * n)] =
              1.0; /* remaining were set to zero */
          if (lb[i] > -Inf) {
            bext[reallincons + tmplincons] = 100;
            if (bext[reallincons + tmplincons] < lb[i] + 3 * fabs(lb[i]))
              bext[reallincons + tmplincons] = lb[i] + 3 * fabs(lb[i]);
            tmplincons++;
          } else { /* both limits not defined */
            bext[reallincons + tmplincons] = 100;
            if (bext[reallincons + tmplincons] < 10 * maxub)
              bext[reallincons + tmplincons] = 10 * maxub;
            tmplincons++;
          }
        }
        if (lb[i] <= -Inf) {
          Aextext[reallincons + tmplincons + i * (lincons + 2 * n)] =
              -1.0; /* remaining were set to zero */
          if (ub[i] < Inf) {
            bext[reallincons + tmplincons] = -100;
            if (bext[reallincons + tmplincons] > lb[i] + 3 * fabs(lb[i]))
              bext[reallincons + tmplincons] = ub[i] - 3 * fabs(ub[i]);
            bext[reallincons + tmplincons] = -bext[reallincons + tmplincons];
            tmplincons++;
          } else { /* both limits not defined */
            /* this code is not reached!! */
            bext[reallincons + tmplincons] = -100;
            if (bext[reallincons + tmplincons] > -10 * minlb)
              bext[reallincons + tmplincons] = -10 * minlb;
            bext[reallincons + tmplincons] = -bext[reallincons + tmplincons];
            tmplincons++;
          }
        }
      }

      tmplincons = lincons + 2 * n;

      msg = mve_presolve(tmplincons, n, Aextext, bext, 80, 1e-8, center_ini);
      if (msg) {
        printf("MVE error %d\n", msg);
        printf("Cannot recover!\nAborting\n");

        free(Aext);
        free(bext);
        free(center_ini);
        free(center);

        return 1;
      }
      recovered = 1;
    }

    /* print center */
    /*printf("Center= ");
    for(j=0;j<n;j++)
            printf("%.3f ",center_ini[j]);*/

    /* compute ellipsoid */
    /* try first with real constraints */
    E2 = pswarm_malloc(n * n * sizeof(double));
    msg = mve_solver(reallincons, n, Aext, bext, center_ini, 80, 1e-6, center,
                     E2);
    if (msg) {
      printf("MVEsolver error %d\n", msg);

      printf("Trying to recover by adding fictitious bounds\n");

      /* add fictitious bounds */

      if (reallincons >= lincons + 2 * n) {
        printf(
            "MVE cannot recover. All bounds already considered!\nAborting\n");

        free(Aext);
        free(bext);
        free(center_ini);
        free(center);

        return 1;
      }

      Aextext = pswarm_malloc((lincons + 2 * n) * n * sizeof(double));

      for (i = 0; i < reallincons; i++) {
        for (j = 0; j < n; j++) {
          Aextext[i + j * (lincons + 2 * n)] = Aext[i + j * reallincons];
        }
      }

      tmplincons = 0;

      for (i = 0; i < n; i++) {
        if (ub[i] >= Inf) {
          Aextext[reallincons + tmplincons + i * (lincons + 2 * n)] =
              1.0; /* remaining were set to zero */
          if (lb[i] > -Inf) {
            bext[reallincons + tmplincons] = 100;
            if (bext[reallincons + tmplincons] < lb[i] + 3 * fabs(lb[i]))
              bext[reallincons + tmplincons] = lb[i] + 3 * fabs(lb[i]);
            tmplincons++;
          } else { /* both limits not defined */
            bext[reallincons + tmplincons] = 100;
            if (bext[reallincons + tmplincons] < 10 * maxub)
              bext[reallincons + tmplincons] = 10 * maxub;
            tmplincons++;
          }
        }
        if (lb[i] <= -Inf) {
          Aextext[reallincons + tmplincons + i * (lincons + 2 * n)] =
              -1.0; /* remaining were set to zero */
          if (ub[i] < Inf) {
            bext[reallincons + tmplincons] = -100;
            if (bext[reallincons + tmplincons] > lb[i] + 3 * fabs(lb[i]))
              bext[reallincons + tmplincons] = ub[i] - 3 * fabs(ub[i]);
            bext[reallincons + tmplincons] = -bext[reallincons + tmplincons];
            tmplincons++;
          } else { /* both limits not defined */
            /* this code is not reached!! */
            bext[reallincons + tmplincons] = -100;
            if (bext[reallincons + tmplincons] > -10 * minlb)
              bext[reallincons + tmplincons] = -10 * minlb;
            bext[reallincons + tmplincons] = -bext[reallincons + tmplincons];
            tmplincons++;
          }
        }
      }

      tmplincons = lincons + 2 * n;

      msg = mve_solver(tmplincons, n, Aextext, bext, center_ini, 80, 1e-6,
                       center, E2);
      if (msg) {
        printf("MVE cannot recover.\nAborting\n");

        free(Aext);
        free(bext);
        free(Aextext);
        free(center_ini);
        free(center);

        return 1;
        ;
      }
    }

    dpotrf_(&Upper, &n, E2, &n, &info);

    E = pswarm_malloc(n * n * sizeof(double));

    memset(E, 0, n * n * sizeof(double));
    for (i = 0; i < n; i++)
      for (j = i; j < n; j++) E[i * n + j] = E2[j * n + i];
  }

#endif /* LINEAR */

#ifdef LINEAR
  if (lincons > 0) {
    tmppt = pswarm_malloc(n * sizeof(double));
    E2tmppt = pswarm_malloc(n * sizeof(double));
  }
#endif

  for (i = accepted; i < s; i++) { /* for all remaining particle   */
#ifdef LINEAR

    if (lincons > 0) {
      for (j = 0; j < n; j++) {
#if SYS_RANDOM == 1
        tmppt[j] = ((rand() / (RAND_MAX + 1.0)) - 0.5) * 2;
#else
        tmppt[j] = ((randflt(&rand_seed)) - 0.5) * 2;
#endif
      }
      tmpptnorm = dnrm2_(&n, tmppt, &oneI);
      for (j = 0; j < n; j++) tmppt[j] /= tmpptnorm;

      dgemv_(&Normal, &n, &n, &one, E, &n, tmppt, &oneI, &zero, E2tmppt, &oneI);

#if SYS_RANDOM == 1
      randscale = pow((rand() / (RAND_MAX + 1.0)), 1.0 / n);
#else
      randscale = pow((randflt(&rand_seed)), 1.0 / n);
#endif

      for (j = 0; j < n; j++) {
        pop->x[i * n + j] = center[j] + randscale * E2tmppt[j];
      }
    } else {
#endif

      for (j = 0; j < n; j++) { /* for all dimensions */
        if (lb[j] > -Inf && ub[j] < Inf) {
#if SYS_RANDOM == 1
          pop->x[i * n + j] =
              (rand() / (RAND_MAX + 1.0)) * (ub[j] - lb[j]) + lb[j];
#else
        pop->x[i * n + j] = (randflt(&rand_seed)) * (ub[j] - lb[j]) + lb[j];
#endif
        } else {
          if (lb[j] <= -Inf && ub[j] >= Inf) { /* both limits infinite */
#if SYS_RANDOM == 1
            pop->x[i * n + j] =
                2 * (rand() / (RAND_MAX + 1.0) - 0.5) * normtmp + initials[j];
#else
          pop->x[i * n + j] =
              2 * (randflt(&rand_seed) - 0.5) * normtmp + initials[j];
#endif
          } else {
            if (ninitials <= 0) {
              printf(
                  "No finite bounds on all variables and no initial guess "
                  "provided\n");
              printf("Unable to obtain initial population\n");
              return 1;
            }
            if (lb[j] <= -Inf) { /* lower infinite and upper finite */
#if SYS_RANDOM == 1
              pop->x[i * n + j] = 2 * (rand() / (RAND_MAX + 1.0) - 0.5) *
                                      (ub[j] - initials[j]) +
                                  initials[j];
#else
            pop->x[i * n + j] =
                2 * (randflt(&rand_seed) - 0.5) * (ub[j] - initials[j]) +
                initials[j];
#endif
            } else { /* upper infinite and lower finite */
#if SYS_RANDOM == 1
              pop->x[i * n + j] = 2 * (rand() / (RAND_MAX + 1.0) - 0.5) *
                                      (initials[j] - lb[j]) +
                                  initials[j];
#else
            pop->x[i * n + j] =
                2 * (randflt(&rand_seed) - 0.5) * (initials[j] - lb[j]) +
                initials[j];
#endif
            }
          }
        }
      }

#ifdef LINEAR
    }
#endif

    pop->fy[i] = +Inf * 10; /* in first iteration, y will be set */
    pop->active[i] = 1;     /* chances to be near the gbest */
  }

#ifdef LINEAR
  if (lincons > 0) {
    free(Aext);
    free(bext);
    free(center_ini);
    free(center);
    free(tmppt);
    free(E2tmppt);
    free(E2);
    free(E);
  }
#endif

  pop->delta = opt.delta;
  /* initialize  velocity */
  memset(pop->v, 0, s * n * sizeof(double));

  // print_pop(n, 0, s, pop);
  return 0;
}

double projection(double xi, double lbi, double ubi) {
  if (xi < lbi) {
    return lbi;
  }
  if (xi > ubi) {
    return ubi;
  }
  return xi;
}

/* print the best of each particle in a population */
void print_pop(int n, int gbest, int s, struct swarm *pop) {
  int i, j, inactive;

  printf("Printing the best so far for each particle\n");

  inactive = 0;
  for (i = 0; i < s; i++) { /* for each particle */
    if (pop->active[i]) {   /* active particle */
      printf("x(%d)=[", i);
      for (j = 0; j < n - 1; j++) /* for each coordinate */
        printf("%.4f,", pop->x[i * n + j]);
      printf("%.4f];\n", pop->x[i * n + n - 1]);
      //      printf("y(%d)=[", i);
      //      for(j=0;j<n-1;j++) /* for each coordinate */
      //          printf("%.4f,", pop->y[i*n+j]);
      //      printf("%.4f];\n", pop->y[i*n+n-1]);
      //      printf("v(%d)=[", i);
      //      for(j=0;j<n-1;j++) /* for each coordinate */
      //          printf("%.4f,", pop->v[i*n+j]);
      //      printf("%.4f];\n", pop->v[i*n+n-1]);
      //      printf("f(%d)=%.20f\n", i, pop->fy[i]);
    } else {
      inactive++;
    }
  }

  printf("%d inactive particles\n", inactive);

  printf("\n The very best\n");
  printf("p(%d)=[", gbest);
  for (j = 0; j < n - 1; j++) /* for each coordinate */
    printf("%.10f,", pop->y[gbest * n + j]);
  printf("%.10f];\n", pop->y[gbest * n + n - 1]);
  printf("f(%d)=%.10f\n", gbest, pop->fy[gbest]);
}

/* write the best of each particle in a population */
void matlab_write_pop(int n, int gbest, int s, struct swarm *pop, int iter) {
  int i;
  FILE *mfile;

  if (n != 2) return;

  if (iter == 1)
    mfile = fopen("pop.m", "w");
  else
    mfile = fopen("pop.m", "a");

  if (!mfile) return;

  fprintf(mfile, "xa1=[");
  for (i = 0; i < s; i++) { /* for each particle */
    if (pop->active[i]) fprintf(mfile, "%.2f,", pop->y[i * n]);
  }
  fprintf(mfile, "];");

  fprintf(mfile, "xa2=[");
  for (i = 0; i < s; i++) { /* for each particle */
    if (pop->active[i]) fprintf(mfile, "%.2f,", pop->y[i * n + 1]);
  }
  fprintf(mfile, "];");
  fprintf(mfile, "hold off;\nir2;\nhold on;\nplot(xa1,xa2,'or');\n");
  // fprintf(mfile,"plot(xa1,xa2,'ob');\n");

  //  fprintf(mfile,"xb1=[");
  //  for(i=0;i<s;i++){ /* for each particle */
  //    if(!pop->active[i])
  //        fprintf(mfile,"%.2f,", pop->y[i*n]);
  //  }
  //  fprintf(mfile,"];");

  //  fprintf(mfile,"xb2=[");
  //  for(i=0;i<s;i++){ /* for each particle */
  //    if(!pop->active[i])
  //        fprintf(mfile,"%.2f,", pop->y[i*n+1]);
  //  }
  //  fprintf(mfile,"];");
  //  fprintf(mfile,"hold off;\nir2;\nhold
  //  on;\nplot(xa1,xa2,'or');\nplot(xb1,xb2,'ob');\n");
  //  fprintf(mfile,"plot(xb1,xb2,'or');\n");
  fprintf(mfile,
          "title('iter=%d, best fx=%.4f, pollsteps=%d, suc=%d, delta=%.8f "
          "nfx=%d');\npause;\n",
          iter, pop->fy[gbest], stats.pollsteps, stats.sucpollsteps, pop->delta,
          stats.objfunctions);

  fclose(mfile);
}

/* write the best swarm particle */
void print_best(int n, int gbest, int s, struct swarm *pop, int iter) {
  int i, nactive;
  FILE *file;

  if (iter == 1)
    file = fopen("results.txt", "w");
  else
    file = fopen("results.txt", "a");

  nactive = 0;
  for (i = 0; i < s; i++)
    if (pop->active[i]) nactive++;

  if (!file) return;

  fprintf(file, "x=[");
  for (i = 0; i < n; i++) { /* for each dimension */
    fprintf(file, "%.8f,", pop->y[gbest * n + i]);
  }
  fprintf(file, "]  fx=%lf\n", pop->fy[gbest]);

  fprintf(file, "Nobj=%d  Npoll=%d  Nsucpoll=%d Active=%d\n",
          stats.objfunctions, stats.pollsteps, stats.sucpollsteps, nactive);

  fclose(file);
}

/* Output function */
double outfcn(int n, int s, int iter, int gbest, struct swarm *pop) {
  if (iter == 0) {
    printf("\n  Iter     Leader     Objective  ");
    printf("\n  -------------------------------\n");
  }

  printf("    %4d   %4d   %4.6e\n", iter, gbest, pop->fy[gbest]);

  return 1.0;
}

/* write array to stdout */
void print_array(int n, double *x) {
  int i;

  printf("=[");
  for (i = 0; i < n - 1; i++) { /* for each dimension */
    printf("%.4f,", x[i]);
  }
  printf("%.4f]\n", x[i]);
}

#if SYS_RANDOM != 1

#define LONG_INT long

/*
The next two functions, myrand and randflt, were copied from
user.c in ASA.
*/

#define MULT ((LONG_INT)25173)
#define MOD ((LONG_INT)65536)
#define INCR ((LONG_INT)13849)
#define FMOD ((double)65536.0)

/***********************************************************************
* double myrand - returns random number between 0 and 1
*   This routine returns the random number generator between 0 and 1
***********************************************************************/

static double myrand(LONG_INT *rand_seed) {
#if TRUE /* (change to FALSE for alternative RNG) */
  *rand_seed = (LONG_INT)((MULT * (*rand_seed) + INCR) % MOD);
  return ((double)(*rand_seed) / FMOD);
#else
/* See "Random Number Generators: Good Ones Are Hard To Find,"
Park & Miller, CACM 31 (10) (October 1988) pp. 1192-1201.
***********************************************************
THIS IMPLEMENTATION REQUIRES AT LEAST 32 BIT INTEGERS
*********************************************************** */
#define _A_MULTIPLIER 16807L
#define _M_MODULUS 2147483647L /* (2**31)-1 */
#define _Q_QUOTIENT 127773L    /* 2147483647 / 16807 */
#define _R_REMAINDER 2836L     /* 2147483647 % 16807 */
  long lo;
  long hi;
  long test;

  hi = *rand_seed / _Q_QUOTIENT;
  lo = *rand_seed % _Q_QUOTIENT;
  test = _A_MULTIPLIER * lo - _R_REMAINDER * hi;
  if (test > 0) {
    *rand_seed = test;
  } else {
    *rand_seed = test + _M_MODULUS;
  }
  return ((double)*rand_seed / _M_MODULUS);
#endif                         /* alternative RNG */
}

/***********************************************************************
* double randflt
***********************************************************************/

static double randflt(LONG_INT *rand_seed) {
  return (resettable_randflt(rand_seed, 0));
}

/***********************************************************************
* double resettable_randflt
***********************************************************************/
static double resettable_randflt(LONG_INT *rand_seed, int reset)

/* shuffles random numbers in random_array[SHUFFLE] array */

{
  /* This RNG is a modified algorithm of that presented in
  * %A K. Binder
  * %A D. Stauffer
  * %T A simple introduction to Monte Carlo simulations and some
  *    specialized topics
  * %B Applications of the Monte Carlo Method in statistical physics
  * %E K. Binder
  * %I Springer-Verlag
  * %C Berlin
  * %D 1985
  * %P 1-36
  * where it is stated that such algorithms have been found to be
  * quite satisfactory in many statistical physics applications. */

  double rranf;
  unsigned kranf;
  int n;
  static int randflt_initial_flag = 0;
  LONG_INT initial_seed;
  static double random_array[SHUFFLE]; /* random variables */

  if (*rand_seed < 0) *rand_seed = -*rand_seed;

  if ((randflt_initial_flag == 0) || reset) {
    initial_seed = *rand_seed;

    for (n = 0; n < SHUFFLE; ++n) random_array[n] = myrand(&initial_seed);

    randflt_initial_flag = 1;

    for (n = 0; n < 1000; ++n) /* warm up random generator */
      rranf = randflt(&initial_seed);

    rranf = randflt(rand_seed);

    return (rranf);
  }

  kranf = (unsigned)(myrand(rand_seed) * SHUFFLE) % SHUFFLE;
  rranf = *(random_array + kranf);
  *(random_array + kranf) = myrand(rand_seed);

  return (rranf);
}
#endif

#ifdef LINEAR
void check_feasible_pop(int n, int s, int lincons, struct swarm *pop,
                        double *lb, double *ub, double *A, double *b) {
  int feasible;
  double how_many;
  int pi, i;
  char Normal = 'N';
  double none = -1.0;
  int oneI = 1;
  double zero = 0.0;
#ifdef linux
  double tmpm[lincons];
#else
  double *tmpm;

  tmpm = pswarm_malloc(lincons * sizeof(double));
#endif

  for (pi = 0; pi < s; pi++) {
    feasible = 1;
    /* check for linear feasibility */
    dgemv_(&Normal, &lincons, &n, &none, A, &lincons, &pop->x[pi * n], &oneI,
           &zero, tmpm, &oneI); /* tmpm = -Ax */
    for (i = 0; i < lincons; i++)
      if (b[i] + tmpm[i] < -1e-5) {
        feasible = 0;
        how_many = b[i] + tmpm[i];
      }

    if (!feasible) {
      printf(
          "Particle %d is not linear feasible %.15f\nThis should not "
          "happen!!!\n",
          pi, how_many);
    }

    feasible = 1;
    for (i = 0; i < n; i++)
      if (pop->x[pi * n + i] > ub[i] || pop->x[pi * n + i] < lb[i])
        feasible = 0;

    if (!feasible) {
      printf(
          "Particle %d is not simple bound feasible\nThis should not "
          "happen!!!\n",
          pi);
    }
  }

#ifndef linux
  free(tmpm);
#endif
}
#endif

/* check is particle is feasible */
int feasible_p(int n, double *x, int lincons, double *A, double *b, double *lb,
               double *ub) {
  int j;
#ifdef LINEAR
  char Normal = 'N';
  double none = -1.0;
  int oneI = 1;
  double zero = 0.0;
#ifdef linux
  double tmpm1[lincons];
#else
  double *tmpm1;

  tmpm1 = pswarm_malloc(lincons * sizeof(double));
#endif

  /* check for linear feasibility */
  if (lincons && A && b) {
    dgemv_(&Normal, &lincons, &n, &none, A, &lincons, x, &oneI, &zero, tmpm1,
           &oneI); /* tmpm = -Ap */
    for (j = 0; j < lincons; j++)
      if (b[j] + tmpm1[j] < 0) {
#ifndef linux
        free(tmpm1);
#endif
        return 0;
      }
  }
#ifndef linux
  free(tmpm1);
#endif

#endif /* LINEAR */
  /* check for bound feasibility */
  for (j = 0; j < n; j++)
    if (x[j] > ub[j] || x[j] < lb[j]) return 0;

  return 1;
}

void *pswarm_malloc(size_t size) {
  void *pointer;

  pointer = malloc(size);

  if (!pointer) {
    printf("Error allocating memory\nAborting\n");
    exit(1);
  }

  return pointer;
}
