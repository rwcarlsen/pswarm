#ifndef Pattern_included

#define Pattern_included


/* A vector in a linked list of vectors */
struct poll_vector {
	double *vector;    /* vector */
	struct poll_vector *next;
};




void pollstep(int n, int lincons, int pi, void (*objf)(), double *lb, double *ub,
	      double *A, double *b, struct poll_vector **last_sucess);
void init_D(int);
void print_D(int);
void print_TC(int);
void print_poll_vector(int, double *);
void init_pattern(int);
void clean_pattern();
void clean_D();
void clean_TC();
void insert_TC(int n, double *dir);
double try_poll_point(int, double *, double (*objf)(), double *, double *);

#endif
