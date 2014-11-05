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
   maxiter - maximum allowed iterations
   tol - tolerance (stopping criteria)
   x - interior point
*/

extern double dnrm2_();
extern double ddot_();
extern void dgemv_();
extern void daxpy_();
extern void dgemm_();
extern void dgetri_();
extern void dgesv_();
extern void dpotrs_();

extern void *pswarm_malloc(size_t size);

void calcstep(int m, int n, double *A, double *B, double *s, double *y,
              double *r1, double *r2, double r3, double *r4, double *dx,
              double *ds, double *dt, double *dy) {
  char Transpose = 'T';
  char Normal = 'N';
  int n1 = n + 1;
  int oneI = 1;
  double none = -1.0;
  double one = 1.0;
  int info;
  int i;

  int *myworkI;
  double *dxdt;
  double *tmp;
  double *tmpB;

  tmp = pswarm_malloc(m * sizeof(double));
  dxdt = pswarm_malloc(n1 * sizeof(double));

  memset(dxdt, 0, n1 * sizeof(double));

  dxdt[n] = 0.0;
  for (i = 0; i < m; i++) {
    tmp[i] = (r1[i] * y[i] - r4[i]) / s[i];
    dxdt[n] += tmp[i];
  }

  memcpy(dxdt, r2, n * sizeof(double));
  dgemv_(&Transpose, &m, &n, &one, A, &m, tmp, &oneI, &one, dxdt, &oneI);

  /*  dpotrs_(&Upper, &n1, &oneI, B, &n1, dxdt, &n1, &info); */

  free(tmp);

  tmpB = pswarm_malloc(n1 * n1 * sizeof(double));
  myworkI = pswarm_malloc(n1 * sizeof(int));

  memcpy(tmpB, B, n1 * n1 * sizeof(double));

  dgesv_(&n1, &oneI, tmpB, &n1, myworkI, dxdt, &n1, &info);

  memcpy(dx, dxdt, n * sizeof(double));
  *dt = dxdt[n];

  memcpy(ds, r1, m * sizeof(double));
  dgemv_(&Normal, &m, &n, &none, A, &m, dx, &oneI, &one, ds, &oneI);

  for (i = 0; i < m; i++) {
    ds[i] -= (*dt);
    dy[i] = (r4[i] - y[i] * ds[i]) / s[i];
  }

  free(myworkI);
  free(dxdt);
  free(tmpB);
}

int mve_presolve(int m, int n, double *A, double *b, int maxiter, double tol,
                 double *x) {
  double t;
  double dt, dtc;
  double tau0 = 0.995;
  double sigma0 = 0.2;
  double tau;
  double sigma;
  double mu;
  double bnrm;
  double gap;
  double rgap;
  double total_err;
  double prif;
  double drif;
  double alphap;
  double alphad;
  double ratio;

  char Normal = 'N';
  char Transpose = 'T';

  double r3;

  double one = 1.0;
  double none = -1.0;
  int oneI = 1;
  double zero = 0.0;

  int iter, i, j;
  int n1 = n + 1;

  double *y;
  double *tmpmm;
  double *tmpnn;
  double *tmpm1;
  double *tmpm2;
  double *B;
  double *AtD;
  double *AtDe;
  double *e_m;
  double *s;
  double *d;
  double *dx;
  double *dxc;
  double *ds;
  double *dsc;
  double *dy;
  double *dyc;
  double *r1;
  double *r2;
  double *r23;
  double *r4;

  y = pswarm_malloc(m * sizeof(double));
  tmpmm = pswarm_malloc(m * m * sizeof(double));
  tmpnn = pswarm_malloc(n * n * sizeof(double));
  tmpm1 = pswarm_malloc(m * sizeof(double));
  tmpm2 = pswarm_malloc(m * sizeof(double));
  B = pswarm_malloc(n1 * n1 * sizeof(double));
  AtD = pswarm_malloc(n * m * sizeof(double));
  AtDe = pswarm_malloc(n * sizeof(double));
  e_m = pswarm_malloc(m * sizeof(double));
  s = pswarm_malloc(m * sizeof(double));
  d = pswarm_malloc(m * sizeof(double));
  dx = pswarm_malloc(n * sizeof(double));
  dxc = pswarm_malloc(n * sizeof(double));
  ds = pswarm_malloc(m * sizeof(double));
  dsc = pswarm_malloc(m * sizeof(double));
  dy = pswarm_malloc(m * sizeof(double));
  dyc = pswarm_malloc(m * sizeof(double));
  r1 = pswarm_malloc(m * sizeof(double));
  r2 = pswarm_malloc(n * sizeof(double));
  r23 = pswarm_malloc(n1 * sizeof(double));
  r4 = pswarm_malloc(m * sizeof(double));

  /* initialize */
  memset(x, 0, n * sizeof(double));
  memset(dx, 0, n * sizeof(double));
  memset(dxc, 0, n * sizeof(double));

  memset(dy, 0, m * sizeof(double));
  memset(dyc, 0, m * sizeof(double));
  memset(ds, 0, m * sizeof(double));
  memset(dsc, 0, m * sizeof(double));

  dt = dtc = 0.0;

  bnrm = dnrm2_(&m, b, &oneI);

  t = b[0];
  for (i = 0; i < m; i++) {
    if (t > b[i]) t = b[i];
    y[i] = 1.0 / m;
  }

  t -= 1.0;

  for (i = 0; i < m; i++) {
    e_m[i] = 1.0;
    s[i] = b[i] - t;
  }

  for (iter = 0; iter < maxiter; iter++) {
    memcpy(r1, s, m * sizeof(double));
    dgemv_(&Normal, &m, &n, &one, A, &m, x, &oneI, &one, r1,
           &oneI); /* r1= s + Ax */

    for (i = 0; i < m; i++) { /* r1 = r1 +t;*/
      r1[i] += t;
      r1[i] = -r1[i];
    }

    daxpy_(&m, &one, b, &oneI, r1, &oneI); /* r1 = b + r1*/

    dgemv_(&Transpose, &m, &n, &none, A, &m, y, &oneI, &zero, r2, &oneI);

    r3 = 1.0;
    for (i = 0; i < m; i++) r3 -= y[i];

    gap = 0.0;
    for (i = 0; i < m; i++) {
      r4[i] = -s[i] * y[i];
      gap -= r4[i];
    }

    prif = dnrm2_(&m, r1, &oneI) / (1 + bnrm);

    memcpy(r23, r2, n * sizeof(double));
    r23[n] = r3;
    drif = dnrm2_(&n1, r23, &oneI);

    rgap = ddot_(&m, b, &oneI, y, &oneI);
    rgap = fabs(rgap - t) / (1 + fabs(t));

    total_err = rgap;
    if (prif > total_err) total_err = prif;
    if (drif > total_err) total_err = drif;

    if (total_err < tol) {
      free(y);
      free(tmpmm);
      free(tmpnn);
      free(tmpm1);
      free(tmpm2);
      free(B);
      free(AtD);
      free(AtDe);
      free(e_m);
      free(s);
      free(d);
      free(dx);
      free(dxc);
      free(ds);
      free(dsc);
      free(dy);
      free(dyc);
      free(r1);
      free(r2);
      free(r23);
      free(r4);

      return 0; /* sucess */
    }

    if (dt > (1e+3) * bnrm || t > (1e+6) * bnrm) {
      free(y);
      free(tmpmm);
      free(tmpnn);
      free(tmpm1);
      free(tmpm2);
      free(B);
      free(AtD);
      free(AtDe);
      free(e_m);
      free(s);
      free(d);
      free(dx);
      free(dxc);
      free(ds);
      free(dsc);
      free(dy);
      free(dyc);
      free(r1);
      free(r2);
      free(r23);
      free(r4);

      return 1; /* unbounded ? */
    }

    for (i = 0; i < m; i++) {
      d[i] = 5.0e+15;
      if (d[i] > y[i] / s[i]) d[i] = y[i] / s[i];
    }

    memset(tmpmm, 0, m * m * sizeof(double));
    for (i = 0; i < m; i++) tmpmm[i * m + i] = d[i];

    dgemm_(&Transpose, &Normal, &n, &m, &m, &one, A, &m, tmpmm, &m, &zero, AtD,
           &n);

    dgemv_(&Normal, &n, &m, &one, AtD, &n, e_m, &oneI, &zero, AtDe, &oneI);

    /* construct B matrix */

    dgemm_(&Normal, &Normal, &n, &n, &m, &one, AtD, &n, A, &m, &zero, tmpnn,
           &n);

    for (i = 0; i < n; i++) {
      B[i * (n + 1) + n] = B[i + n * (n + 1)] = AtDe[i];
      for (j = 0; j < n; j++) {
        B[i + j * (n + 1)] = tmpnn[i + j * n];
      }
    }

    B[n + n * (n + 1)] = d[0];
    for (i = 1; i < m; i++) B[n + n * (n + 1)] += d[i];

    for (i = 0; i < n + 1; i++) B[i * (n + 1) + i] += 1e-14;

    /*    dpotrf_(&Upper, &n1, B, &n1, &info); No Cholesky decomposition */

    calcstep(m, n, A, B, s, y, r1, r2, r3, r4, dx, ds, &dt, dy);

    alphap = -1.0;
    alphad = -1.0;

    for (i = 0; i < m; i++) {
      if (alphap > ds[i] / s[i]) alphap = ds[i] / s[i];
      if (alphad > dy[i] / y[i]) alphad = dy[i] / y[i];
    }

    alphad = -1.0 / alphad;
    alphap = -1.0 / alphap;

    memcpy(tmpm1, s, m * sizeof(double));
    daxpy_(&m, &alphap, ds, &oneI, tmpm1, &oneI);

    memcpy(tmpm2, y, m * sizeof(double));
    daxpy_(&m, &alphad, dy, &oneI, tmpm2, &oneI);

    ratio = ddot_(&m, tmpm1, &oneI, &tmpm2, &oneI) / gap;

    sigma = sigma0;
    if (sigma > ratio * ratio) sigma = ratio * ratio;
    mu = sigma * gap / m;

    memset(r1, 0, m * sizeof(double));
    memset(r2, 0, n * sizeof(double));
    r3 = 0.0;
    for (i = 0; i < m; i++) r4[i] = mu - ds[i] * dy[i];

    calcstep(m, n, A, B, s, y, r1, r2, r3, r4, dxc, dsc, &dtc, dyc);

    daxpy_(&n, &one, dxc, &oneI, dx, &oneI);
    daxpy_(&m, &one, dsc, &oneI, ds, &oneI);
    daxpy_(&m, &one, dyc, &oneI, dy, &oneI);
    dt += dtc;

    alphap = -0.5;
    alphad = -0.5;

    for (i = 0; i < m; i++) {
      if (alphap > ds[i] / s[i]) alphap = ds[i] / s[i];
      if (alphad > dy[i] / y[i]) alphad = dy[i] / y[i];
    }

    alphap = -1.0 / alphap;
    alphad = -1.0 / alphad;

    tau = tau0;

    if (tau < 1 - gap / m) tau = 1 - gap / m;
    if (tau * alphap < 1) {
      alphap = tau * alphap;
    } else {
      alphap = 1.0;
    }

    if (tau * alphad < 1) {
      alphad = tau * alphad;
    } else {
      alphad = 1.0;
    }

    daxpy_(&n, &alphap, dx, &oneI, x, &oneI);
    daxpy_(&m, &alphap, ds, &oneI, s, &oneI);
    t += alphap * dt;
    daxpy_(&m, &alphad, dy, &oneI, y, &oneI);
  }

  free(y);
  free(tmpmm);
  free(tmpnn);
  free(tmpm1);
  free(tmpm2);
  free(B);
  free(AtD);
  free(AtDe);
  free(e_m);
  free(s);
  free(d);
  free(dx);
  free(dxc);
  free(ds);
  free(dsc);
  free(dy);
  free(dyc);
  free(r1);
  free(r2);
  free(r23);
  free(r4);

  if (t < 1e-16) return 2; /* no volume */

  return 0;
}

/*
int main(int argc, char **argv)
{
  int m=6, n=3;
  int msg;

#ifdef linux
  double A[m*n];
  double b[m];
  double x[n];
#else
  double *A;
  double *b;
  double *x;

  A=malloc(m*n*sizeof(double));
  b=malloc(m*sizeof(double));
  x=malloc(n*sizeof(double));
#endif


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

  msg=mve_presolve(m, n, A, b, 10, 1e-2, x);

#ifndef linux
  free(A);
  free(b);
  free(x);
#endif

  return 0;
}
*/
