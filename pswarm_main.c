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


#ifdef MPE
/* for MPE */
int ComputeID_begin, ComputeID_end, SendID_begin, SendID_end, RecvID_begin,
  RecvID_end;
#endif


typedef struct { char *msg; } Exit_code;
extern struct Stats stats;


static Exit_code exit_codes[] = {
  {"Normal exit"},                /* 0 */
  {"Abnormal exit"},              /* 1 */
  {"Failed to allocate memory"},  /* 2 */
  {"Unable to initalize population"},  /* 3 */
};

extern void *pswarm_malloc(size_t size);


void objfun(int, int, double *, double *, double *, double *);

#ifdef MPI
void MPI_objfun_deamon(int n, double *, double *, int);
#endif


extern int PSwarm(int, void (*)(),double *, double *, int, double *, double *, double **, double *, double *);
extern struct Options opt;

extern void save_cache_file(int n, int m);
extern void load_cache_file(int n, int m);
extern int read_cesam_cache(int n, double *x, double *age, 
                double *teff, double *lum, double *r);
extern void write_cesam_cache(int n, double *x, double *age,
                  double *teff, double *lum, double *r);

extern void set_problem(double *x, double *lb, double *ub, double *A, double *b);
extern void set_problem_dimension(int *n, int *lincons);

#ifdef MPI

extern void user_init_MPI(int);

#else

extern void user_init();

#endif

static jmp_buf Jb;

void catchfpe(int n)
{
  printf("\nFloating point error.\n");
  fflush(stdout);
  longjmp(Jb,1);
}


/************************************************
                    Main      
*************************************************/
int main(int argc, char **argv)
{
  int exit_code, i;
  double *sol=NULL;
  double *f=NULL;
  double *lb, *ub;
  double *A, *b; /* linear constraints defined in a Fortran way */
  int lincons; /* number of linear constraints */

#ifdef MPI
  int stop=0, MPI_myrank, MPI_numprocs;
  /*char hostname[256]; */
#endif

  int n;
  double *X0;
  
  if (!setjmp(Jb)){
    signal(SIGFPE, catchfpe);

    /* lower and upper bounds on variables */
    lb=pswarm_malloc(n*sizeof(double));    


    ub=pswarm_malloc(n*sizeof(double));
    
    if(!lb || !ub){
      printf("Unable to allocate memory for variable bounds\n");
      exit(1);
    }

    X0=pswarm_malloc(n*sizeof(double));
    if(!X0){
      printf("Unable to allocate memory for initial guess\n");
      exit(1);
    }

    if(lincons){
      A = pswarm_malloc(lincons*n*sizeof(double));
      b = pswarm_malloc(lincons*sizeof(double));
    
      memset(A, 0, (lincons)*n*sizeof(double));
      memset(b, 0, (lincons)*sizeof(double));
    } else {
      A=b=NULL;
    }

    set_problem(X0, lb, ub, A, b);


    /* check for finite bound on the variables */
/*    for(i=0;i<n;i++){
      if(lb[i]<=-Inf || ub[i]>=Inf){
    printf("Variable %i has no finite bounds\nPlease provide finite bounds for all variables\n",i);
    exit(1);
      }
    }*/

    

#ifdef MPI
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_numprocs);

    user_init_MPI(MPI_myrank);

#else 

    user_init();

#endif

#ifdef MPE
    MPE_Init_log();

    ComputeID_begin = MPE_Log_get_event_number();
    ComputeID_end   = MPE_Log_get_event_number();
    SendID_begin    = MPE_Log_get_event_number();
    SendID_end      = MPE_Log_get_event_number();
    RecvID_begin    = MPE_Log_get_event_number();
    RecvID_end      = MPE_Log_get_event_number();
#endif

  

#ifdef MPI
    if(!MPI_myrank){ /* I am the algorithm */
#else
    if(1){
#endif

#ifdef MPE      
      MPE_Describe_state(ComputeID_begin, ComputeID_end, "Compute", "red");
      MPE_Describe_state(SendID_begin, SendID_end, "Send", "blue");
      MPE_Describe_state(RecvID_begin, RecvID_end, "Recv", "green");
      
      MPE_Start_log();
#endif

      //      load_cache_file(n, 4);
    
      exit_code=PSwarm(n, &objfun, lb, ub, lincons, A, b, &sol, f, X0);

#ifdef MPI
      /* Goodbye to objective function processes */
      stop=0;
      for(i=1;i<MPI_numprocs;i++)
    MPI_Send(&stop, 1, MPI_INT, i, 99, MPI_COMM_WORLD);
#endif
      
      //      save_cache_file(n, 4);
      
	  if(opt.IPrint>=0)
        printf("\n%s\n", exit_codes[exit_code].msg);


  
    } else { /* I am an objective function process */

#ifdef MPI
      stats.objfunctions=0;
  
      MPI_objfun_deamon(n, lb, ub, MPI_myrank);

      printf("Deamon %d computed %d objectives\n", MPI_myrank, stats.objfunctions);
#endif

      exit_code=0; /* never executed if not MPI */
    }

#ifdef MPE
    MPE_Finish_log("PPswarm");
#endif

#ifdef MPI
    MPI_Finalize();
#endif

  }

  if(lb)
	  free(lb);
  if(ub)
	  free(ub);
  if(A)
	  free(A);
  if(b)
	  free(b);

  return exit_code;
}

#ifdef MPI

void MPI_objfun_deamon(int n, double *lb, double *ub, int MPI_myrank)
{
  int action=1, particle;
  double fx, x[n];
  MPI_Status status;


  do {
#ifdef MPE
    MPE_Log_event(RecvID_begin, 0, NULL);
#endif
    MPI_Recv(&action, 1, MPI_INT, 0, 99, MPI_COMM_WORLD, &status);
    if(action){ /* Objective function request */
      if(action==2)
    MPI_Recv(&particle, 1, MPI_INT, 0, 99, MPI_COMM_WORLD, &status);
      MPI_Recv(&x, n, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD, &status);

#ifdef MPE
      MPE_Log_event(RecvID_end, 0, NULL);
 
      MPE_Log_event(ComputeID_begin, 0, NULL);
#endif

      objfun(n, 1, x, lb, ub, &fx);
      stats.objfunctions++;



#ifdef MPE
      MPE_Log_event(ComputeID_end, 0, NULL);

      MPE_Log_event(SendID_begin, 0, NULL);
#endif
      if(action==2)
    MPI_Send(&particle, 1, MPI_INT, 0, 99, MPI_COMM_WORLD);
      MPI_Send(&fx, 1, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD);
#ifdef MPE
      MPE_Log_event(SendID_end, 0, NULL);
#endif
    } else {
#ifdef MPE
      MPE_Log_event(RecvID_end, 0, NULL); /* received instruction to stop */
#endif
    }
  } while(action);


  return;
}

#endif

