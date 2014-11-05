
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

#include "pswarm.h"

extern struct Stats stats;

/**********************************************/
/*                                            */
/* User defined global variables              */
/*                                            */
/**********************************************/

#define DIM 6         /* a linear problem with 6 varibles */
#define CONS 5        /* a linear problem with 5 linear constraints */
double obj_coef[DIM]; /* objective function coefficients */

int n_dims() {
  return DIM;
}

int n_cons() {
  return CONS;
}

/**********************************************/
/*                                            */
/* End of user defined global variables       */
/*                                            */
/**********************************************/

/*******************************************************/
/*                                                     */
/* User define objective function                      */
/*                                                     */
/*******************************************************/

void objfun(int n, int m, double *x, double *lb, double *ub, double *fx) {
  int i, j;

  for (j = 0; j < m; j++) {
    fx[j] = 0.0;
    for (i = 0; i < n; i++) fx[j] += obj_coef[i] * x[j * n + i];
  }
}

/*******************************************************/
/*                                                     */
/* End of user define objective function               */
/*                                                     */
/*******************************************************/

/*******************************************************/
/*                                                     */
/* User define problem data                            */
/*                                                     */
/*******************************************************/

void set_problem_dimension(int *n, int *lincons) {
  *n = DIM;        /* problem dimension -- number of variables */
  *lincons = CONS; /* number of linear constraints */
}

void set_problem(double *x, double *lb, double *ub, double *A, double *b) {
  /* initial guess */
  x[0] = 10;
  x[1] = 11;
  x[2] = 12;
  x[3] = 13;
  x[4] = 14;
  x[5] = 15;

  /* lower bounds */
  lb[0] = 0.0;
  lb[1] = 0.0;
  lb[2] = 0.0;
  lb[3] = 0.0;
  lb[4] = 0.0;
  lb[5] = 0.0;

  /* upper bounds */
  ub[0] = 50;
  ub[1] = 60;
  ub[2] = 85;
  ub[3] = 70;
  ub[4] = 40;
  ub[5] = 24;

  /* A - linear constraints defined in a Frontran way (transposed) */

  /* b - columnwise vector of independent terms */

  /*
  A=[-1  0  0   0.9    0   0;
      1  0  0  -1.15   0   0;
      2  2  1   1     -1  -1;
      2  4  1   2      1   0;
      3  6  1   5      0   1]
  */

  memset(A, 0, DIM * CONS * sizeof(double)); /* clear memory. We only need to
                                                set the nonzero elements */
  A[0] = -1;
  A[1] = 1;
  A[2] = 2;
  A[3] = 2;
  A[4] = 3;

  A[7] = 2;
  A[8] = 4;
  A[9] = 6;

  A[12] = 1;
  A[13] = 1;
  A[14] = 1;

  A[15] = 0.9;
  A[16] = -1.15;
  A[17] = 1;
  A[18] = 2;
  A[19] = 5;

  A[22] = -1;
  A[23] = 1;

  A[27] = -1;
  A[29] = 1;

  /*
    b=[0; 0; 160; 200; 80];
  */

  b[0] = b[1] = 0;
  b[2] = 160.0;
  b[3] = 200.0;
  b[4] = 80.0;
}

/*******************************************************/
/*                                                     */
/* End of user define problem data                     */
/*                                                     */
/*******************************************************/

/*******************************************************/
/*                                                     */
/* User initialize problem data                        */
/*                                                     */
/*******************************************************/

void user_init() {
  /* use this function if some global data must be initialized */

  obj_coef[0] = -10.0;
  obj_coef[1] = -15.0;
  obj_coef[2] = -22.0;
  obj_coef[3] = -17.0;
  obj_coef[4] = 0.0;
  obj_coef[5] = 0.0;
}

/*******************************************************/
/*                                                     */
/* End of user initialize problem data                 */
/*                                                     */
/*******************************************************/
