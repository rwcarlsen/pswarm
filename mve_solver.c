/***************************************************************************
*        MVE C implementation
*        See: ----
*        For details
***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Parameters:
A - m*n coefficients matriz
b - independent term
x0 - initial guess
maxiter - maximum allowed iterations
tol - tolerance (stopping criteria)
x - Ellipsoide center
E2 - Ellipsoide radious

*/

#ifndef max
#define max(a, b) (((a) > (b)) ? (a) : (b))
#endif

extern double dnrm2_();
extern void dgemv_();
extern void daxpy_();
extern void dgemm_();
extern void dgetri_();
extern void dgesv_();
extern void dgetrf_();

extern void *pswarm_malloc(size_t size);

int mve_solver(int m, int n, double *A, double *b, double *x0, double maxiter,
               double tol, double *x, double *E2) {
  double one = 1.0;
  double none = -1.0;
  double zero = 0.0;
  char Normal = 'N';
  char Transpose = 'T';
  double bnrm;

  int oneI = 1;

  double minmu = 1e-8;
  double tau0 = 0.75;
  double tau;
  double pp;

  int i, j, iter, info;

  double astep;
  int lmywork = max(n, m);  //*n*m*m;
  double res;
  double t, t2;
  double gap;
  double rmu;

  double *bmAx0;
  double *tmpnm;
  double *tmpmm;
  double *tmpnn;
  double *R1;
  double *R2;
  double *R3;
  double *dy;
  double *dyDy;
  double *dz;
  double *R23;
  double *R3Dy;
  double *yz;
  double *yh;
  double *y2h;
  double *G;
  double *Y;
  double *YA;
  double *T;
  double *Q;
  double *YQ;
  double *YQQY;
  double *mywork;
  int *myworkI;
  double *bmAx;
  double *nA;
  double *Adx;
  double *nb;
  double *h;
  double *y;
  double *z;

  bnrm = dnrm2_(&m, b, &oneI);

  if (maxiter <= 0) maxiter = 50;

  if (tol <= 0) tol = 1.e-4;

  bmAx0 = pswarm_malloc(sizeof(double) * m);
  memcpy(bmAx0, b, sizeof(double) * m);

  dgemv_(&Normal, &m, &n, &none, A, &m, x0, &oneI, &one, bmAx0,
         &oneI); /* bmAx0=b-Ax0 */

  for (i = 0; i < m; i++) {
    if (bmAx0[i] <= 0.0) {
      printf("x0 is not interior\n\nAborting\n\n");
      free(bmAx0);

      return (2);
    }
  }

  lmywork = max(n, m);
  tmpnm = pswarm_malloc(sizeof(double) * n * m);
  R1 = pswarm_malloc(sizeof(double) * n);
  R2 = pswarm_malloc(sizeof(double) * m);
  R3 = pswarm_malloc(sizeof(double) * m);
  dy = pswarm_malloc(sizeof(double) * m);
  dyDy = pswarm_malloc(sizeof(double) * m);
  dz = pswarm_malloc(sizeof(double) * m);
  R23 = pswarm_malloc(sizeof(double) * m);
  R3Dy = pswarm_malloc(sizeof(double) * m);
  yz = pswarm_malloc(sizeof(double) * m);
  yh = pswarm_malloc(sizeof(double) * m);
  y2h = pswarm_malloc(sizeof(double) * m);
  G = pswarm_malloc(sizeof(double) * m * m);
  Y = pswarm_malloc(sizeof(double) * m * m);
  YA = pswarm_malloc(sizeof(double) * m * n);
  T = pswarm_malloc(sizeof(double) * m * n);
  Q = pswarm_malloc(sizeof(double) * m * m);
  YQ = pswarm_malloc(sizeof(double) * m * m);
  YQQY = pswarm_malloc(sizeof(double) * m * m);
  mywork = pswarm_malloc(sizeof(double) * lmywork);
  myworkI = pswarm_malloc(sizeof(int) * lmywork);
  tmpmm = pswarm_malloc(sizeof(double) * m * m);
  tmpnn = pswarm_malloc(sizeof(double) * n * n);
  bmAx = pswarm_malloc(sizeof(double) * m);
  nA = pswarm_malloc(sizeof(double) * m * n);
  Adx = pswarm_malloc(sizeof(double) * m);
  nb = pswarm_malloc(sizeof(double) * m);
  h = pswarm_malloc(sizeof(double) * m);
  y = pswarm_malloc(sizeof(double) * m);
  z = pswarm_malloc(sizeof(double) * m);

  for (i = 0; i < m; i++) {   /* for each row */
    for (j = 0; j < n; j++) { /* for each column */
      nA[i + j * m] = A[i + j * m] / bmAx0[i];
    }
    nb[i] = bmAx[i] = 1.0;
    y[i] = 1.0;
  }

  for (j = 0; j < n; j++) x[j] = 0.0;

  res = 1;

  for (iter = 0; iter < maxiter; iter++) { /* Main loop */

    if (iter) {
      astep = -astep;
      daxpy_(&m, &astep, Adx, &oneI, bmAx, &oneI); /* atencao astep = -astep */
    }

    memset(Y, 0, sizeof(double) * m * m);

    for (i = 0; i < m; i++) Y[i * m + i] = y[i];

    dgemm_(&Transpose, &Normal, &n, &m, &m, &one, nA, &m, Y, &m, &zero, tmpnm,
           &n);

    dgemm_(&Normal, &Normal, &n, &n, &m, &one, tmpnm, &n, nA, &m, &zero, E2,
           &n);

    dgetrf_(&n, &n, E2, &n, myworkI, &info);
    dgetri_(&n, E2, &n, myworkI, mywork, &lmywork,
            &info); /* we should check for info==0 */

    dgemm_(&Normal, &Normal, &m, &n, &n, &one, nA, &m, E2, &n, &zero, tmpnm,
           &m); /* tmpnm matriz m por n */

    dgemm_(&Normal, &Transpose, &m, &m, &n, &one, tmpnm, &m, nA, &m, &zero, Q,
           &m);

    for (i = 0; i < m; i++) h[i] = sqrt(Q[i * m + i]);

    if (!iter) {
      /* min */
      t = bmAx[0] / h[0];
      for (i = 1; i < m; i++) {
        if (t > bmAx[i] / h[i]) t = bmAx[i] / h[i];
      }

      for (i = 0; i < m; i++) {
        y[i] /= t * t;
        h[i] *= t;
      }

      for (i = 0; i < m; i++) {
        z[i] = 1e-1;
        if (z[i] < bmAx[i] - h[i]) z[i] = bmAx[i] - h[i];
      }

      t2 = t * t; /* t=t^2 */
      dgemm_(&Normal, &Normal, &m, &m, &m, &zero, Y, &m, Y, &m, &t2, Q,
             &m); /* Q=t^2*Q , Q=Y*/

      t2 = 1 / t2; /*invt=1/(t^2) */
      dgemm_(&Normal, &Normal, &m, &m, &m, &zero, Y, &m, Y, &m, &t2, Y,
             &m); /* Y=Y/t^2 */
    }

    gap = 0.0;
    for (i = 0; i < m; i++) {
      yz[i] = y[i] * z[i];
      yh[i] = y[i] * h[i];
      gap += yz[i];
    }

    gap /= m; /* sum(yz)/m */

    rmu = 0.5;
    if (rmu > gap) rmu = gap;

    rmu *= gap;

    if (rmu < minmu) rmu = minmu;

    dgemv_(&Transpose, &m, &n, &none, nA, &m, yh, &oneI, &zero, R1,
           &oneI); /* -yh*A' */

    memcpy(R2, bmAx, m * sizeof(double));   /* R2=bmAx */
    daxpy_(&m, &none, h, &oneI, R2, &oneI); /* R2=R2-h*/
    daxpy_(&m, &none, z, &oneI, R2, &oneI); /* R2=R2-z*/

    for (i = 0; i < m; i++) R3[i] = rmu - yz[i];

    res = fabs(R1[0]);

    for (i = 1; i < n; i++)
      if (res < fabs(R1[i])) res = fabs(R1[i]);

    for (i = 0; i < m; i++)
      if (res < fabs(R2[i])) res = fabs(R2[i]);

    for (i = 0; i < m; i++)
      if (res < fabs(R3[i])) res = fabs(R3[i]);

    if (res < tol * (1 + bnrm) && rmu <= minmu) {
      /* converged */
      for (i = 0; i < n; i++) x[i] += x0[i];

      return 0;
    }

    dgemm_(&Normal, &Normal, &m, &m, &m, &one, Y, &m, Q, &m, &zero, YQ, &m);

    for (i = 0; i < m; i++)
      for (j = 0; j < m; j++)
        YQQY[i * m + j] = YQ[i * m + j] * YQ[i + j * m]; /* YQQY = YQ .* YQ' */

    for (i = 0; i < m; i++) y2h[i] = 2 * yh[i];

    dgemm_(&Normal, &Normal, &m, &n, &m, &one, Y, &m, nA, &m, &zero, YA, &m);

    memcpy(G, YQQY, m * m * sizeof(double));
    for (i = 0; i < m; i++) {
      if (y2h[i] * z[i] > 1e-12) {
        G[i * m + i] += y2h[i] * z[i];
      } else {
        G[i * m + i] += 1e-12;
      }
    }

    memset(tmpmm, 0, m * m * sizeof(double));
    for (i = 0; i < m; i++) tmpmm[i * m + i] = h[i] + z[i];

    memset(T, 0, m * n * sizeof(double));
    dgemm_(&Normal, &Normal, &m, &n, &m, &one, tmpmm, &m, YA, &m, &zero, T, &m);

    memcpy(tmpmm, G, m * m * sizeof(double));
    dgesv_(&m, &n, tmpmm, &m, myworkI, T, &m,
           &info); /* we should check for info */

    memset(tmpmm, 0, m * m * sizeof(double));
    for (i = 0; i < m; i++) tmpmm[i * m + i] = y2h[i];

    dgemm_(&Normal, &Normal, &m, &n, &m, &one, tmpmm, &m, T, &m, &none, YA,
           &m); /* ATP = YA' */

    for (i = 0; i < m; i++) {
      R3Dy[i] = R3[i] / y[i];
      R23[i] = R2[i] - R3Dy[i];
    }

    dgemm_(&Transpose, &Normal, &n, &n, &m, &one, YA, &m, nA, &m, &zero, tmpnn,
           &n); /* ATP*A = tmpnn */

    dgemv_(&Transpose, &m, &n, &one, YA, &m, R23, &oneI, &one, R1, &oneI);

    dgesv_(&n, &oneI, tmpnn, &n, myworkI, R1, &n, &info);
        /* R1 = dx */ /* we should check for info */

    dgemv_(&Normal, &m, &n, &one, nA, &m, R1, &oneI, &zero, Adx, &oneI);

    for (i = 0; i < m; i++) {
      dyDy[i] = y2h[i] * (Adx[i] - R23[i]);
    }

    dgesv_(&m, &oneI, G, &m, myworkI, dyDy, &m,
           &info); /* we should check for info */

    for (i = 0; i < m; i++) {
      dy[i] = y[i] * dyDy[i];
      dz[i] = R3Dy[i] - z[i] * dyDy[i];
    }

    astep = 1;

    for (i = 1; i < n; i++) {
      pp = -0.5;
      if (pp > -Adx[i] / bmAx[i]) pp = -Adx[i] / bmAx[i];
      if (astep > (-1.0 / pp)) astep = -1.0 / pp;
    }

    for (i = 1; i < m; i++) {
      pp = -0.5;
      if (pp > dyDy[i]) pp = dyDy[i];
      if (astep > (-1.0 / pp)) astep = -1.0 / pp;
    }

    for (i = 1; i < m; i++) {
      pp = -0.5;
      if (pp > dz[i] / z[i]) pp = dz[i] / z[i];
      if (astep > (-1.0 / pp)) astep = -1.0 / pp;
    }

    tau = tau0;
    if (tau < 1 - res) tau = 1 - res;

    astep *= tau;

    daxpy_(&n, &astep, R1, &oneI, x, &oneI);
    daxpy_(&m, &astep, dy, &oneI, y, &oneI);
    daxpy_(&m, &astep, dz, &oneI, z, &oneI);
  }

  for (i = 0; i < n; i++) x[i] += x0[i];

  printf("Maximum number of iterations in MVE reached\n");

  free(bmAx0);
  free(tmpnm);
  free(R1);
  free(R2);
  free(R3);
  free(dy);
  free(dyDy);
  free(dz);
  free(R23);
  free(R3Dy);
  free(yz);
  free(yh);
  free(y2h);
  free(G);
  free(Y);
  free(YA);
  free(T);
  free(Q);
  free(YQ);
  free(YQQY);
  free(mywork);
  free(myworkI);
  free(tmpmm);
  free(tmpnn);
  free(bmAx);
  free(nA);
  free(Adx);
  free(nb);
  free(h);
  free(y);
  free(z);

  return 1;
}

/*
int main(int argc, char **argv)
{
#define  m 6
#define  n 3

double A[m*n];
double b[m];
double x0[n];
double x[n];
double E2[n*n];


memset(A, 0, m*n*sizeof(double));
memset(b, 0, m*sizeof(double));

// A declared in a Fortran way (transpose)

A[0+0*m]=1.0;
A[1+1*m]=1.0;
A[2+2*m]=1.0;

A[3+0*m]=-1.0;
A[4+1*m]=-1.0;
A[5+2*m]=-1.0;

b[0]=3.0;
b[1]=2.0;
b[2]=2.0;
b[3]=0.0;
b[4]=0.0;
b[5]=0.0;


x0[0]=1.0;
x0[1]=1.0;
x0[2]=1.0;

mve_solver(m, n, A, b, x0, 10, 1e-2, x, E2);

return 0;
}
*/

/*

DO 70 I = 1, N
DO 80 J = 1, LDA
write(*,*)'A(',I,',',J,')',A(I,J)
80      continue
70   continue
*/
