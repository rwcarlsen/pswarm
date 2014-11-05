#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef linux
#include <unistd.h>
#else
#include <io.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <math.h>

typedef struct point {
  double *key;
  double *value;
  struct point *next;
} point;

static struct point *cache_list = NULL;

int check_cache(int n, double *key, int m, double *value) {
  struct point *tmp;

  //  printf("Check cache\n");
  //  printf("x=[%lf,%lf,%lf,%lf,%lf,%lf]-> [%lf,%lf,%lf,%lf]\n",key[0],
  //	 key[1], key[2], key[3], key[4], key[5],
  //	 value[0], value[1], value[2], value[3]);

  if (!cache_list) { /* list not initialized */
    return 0;
  }

  tmp = cache_list;
  while (tmp) {
    if (!memcmp(key, tmp->key, n * sizeof(double))) {
      memcpy(value, tmp->value, m * sizeof(double));
      return 1;
    };

    tmp = tmp->next;
  }

  return 0;
}

void insert_cache(int n, double *key, int m, double *value) {
  struct point *tmp;

  //  printf("Insert cache\n");
  //  printf("x=[%lf,%lf,%lf,%lf,%lf,%lf]-> [%lf,%lf,%lf,%lf]\n",key[0],
  //	 key[1], key[2], key[3], key[4], key[5],
  //	 value[0], value[1], value[2], value[3]);

  if (!cache_list) { /* initialize list */
    cache_list = malloc(sizeof(struct point));
    cache_list->key = malloc(n * sizeof(double));
    cache_list->value = malloc(m * sizeof(double));
    cache_list->next = NULL;
    memcpy((void *)cache_list->key, (void *)key, n * sizeof(double));
    memcpy((void *)cache_list->value, (void *)value, m * sizeof(double));
    return;
  }

  tmp = malloc(sizeof(struct point));
  tmp->key = malloc(n * sizeof(double));
  tmp->value = malloc(m * sizeof(double));
  memcpy(tmp->key, key, n * sizeof(double));
  memcpy(tmp->value, value, m * sizeof(double));

  tmp->next = cache_list;
  cache_list = tmp;

  return;
}

void save_cache_file(int n, int m) {
  struct point *tmp;
  int fid;

#ifdef linux
  if ((fid = open("cache_file", O_RDWR | O_TRUNC | O_CREAT,
                  S_IRUSR | S_IWUSR)) == -1) {
#else
  if ((fid = _open("cache_file", _O_RDWR | _O_TRUNC | _O_CREAT,
                   _S_IREAD | _S_IWRITE)) == -1) {
#endif
    perror("Open cache file error\n");
    return;
  }

  tmp = cache_list;
  while (tmp) {
    write(fid, tmp->key, n * sizeof(double));
    write(fid, tmp->value, m * sizeof(double));
    tmp = tmp->next;
  }

  close(fid);

  return;
}

void load_cache_file(int n, int m) {
  struct point *aux;
  int fid;
#ifdef linux
  double key[n], value[m];
#else
  double *key, *value;

  key = malloc(n * sizeof(double));
  value = malloc(m * sizeof(double));
#endif

  cache_list = NULL;

  if ((fid = open("cache_file", O_RDONLY)) == -1) {
    perror("Open cache file error\n");
    return;
  }

  while (read(fid, key, n * sizeof(double)) == n * sizeof(double)) {
    read(fid, value, m * sizeof(double));
    aux = malloc(sizeof(struct point));
    aux->key = malloc(n * sizeof(double));
    aux->value = malloc(m * sizeof(double));
    memcpy(aux->key, key, n * sizeof(double));
    memcpy(aux->value, value, m * sizeof(double));
    aux->next = cache_list;
    cache_list = aux;
  }

  close(fid);

#ifndef linux
  free(key);
  free(value);
#endif

  return;
}

/* int main() */
/* { */
/*   double x[3], y[2]; */

/*   x[0]=1.3;x[1]=2.4;x[2]=4.0; */
/*   y[0]=2.0;y[1]=4.9; */

/*   load_cache_file(3, 2); */

/*   if(!check_cache(3, x, 2, y)){ */
/*     printf("Point was not on cache\n"); */
/*     insert_cache(3, x, 2, y); */
/*   } */

/*   save_cache_file(3, 2); */

/*   return 0; */
/* } */

int read_cesam_cache(int n, double *x, double *age, double *teff, double *lum,
                     double *r) {
  double y[4];
#ifdef linux
  double xtmp[n];
#else
  double *xtmp;

  xtmp = malloc(n * sizeof(double));
#endif

  xtmp[0] = (double)((long int)(x[0] * 1000 + 0.5) / 1000.0);
  xtmp[1] = (double)((long int)(x[1] + 0.5));
  xtmp[2] = (double)((long int)(x[2] * 10000 + 0.5) / 10000.0);
  xtmp[3] = (double)((long int)(x[3] * 10000 + 0.5) / 10000.0);
  xtmp[4] = (double)((long int)(x[4] * 100 + 0.5) / 100.0);
  xtmp[5] = (double)((long int)(x[5] * 100 + 0.5) / 100.0);

  if (check_cache(n, xtmp, 4, y)) {
    *age = y[0];
    *teff = y[1];
    *lum = y[2];
    *r = y[3];
#ifdef linux
    free(xtmp);
#endif
    return 1;
  }

#ifdef linux
  free(xtmp);
#endif

  return 0;
}

void write_cesam_cache(int n, double *x, double *age, double *teff, double *lum,
                       double *r) {
  double y[4];
#ifdef linux
  double xtmp[n];
#else
  double *xtmp;

  xtmp = malloc(n * sizeof(double));
#endif

  xtmp[0] = (double)((long int)(x[0] * 1000 + 0.5) / 1000.0);
  xtmp[1] = (double)((long int)(x[1] + 0.5));
  xtmp[2] = (double)((long int)(x[2] * 10000 + 0.5) / 10000.0);
  xtmp[3] = (double)((long int)(x[3] * 10000 + 0.5) / 10000.0);
  xtmp[4] = (double)((long int)(x[4] * 100 + 0.5) / 100.0);
  xtmp[5] = (double)((long int)(x[5] * 100 + 0.5) / 100.0);

  y[0] = *age;
  y[1] = *teff;
  y[2] = *lum;
  y[3] = *r;

  insert_cache(n, xtmp, 4, y);

#ifdef linux
  free(xtmp);
#endif

  return;
}
