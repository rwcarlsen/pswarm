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

#ifndef AMPL

extern void set_problem(double *x, double *lb, double *ub, double *A, double *b);
extern void set_problem_dimension(int *n, int *lincons);

#ifdef MPI

extern void user_init_MPI(int);

#else

extern void user_init();

#endif

#endif


#ifdef AMPL

static double objsign;
static fint NERROR = -1;
#define asl cur_ASL

char pswarm_version[]="PSwarm v1.4";

/* This struct member names must be in alphabetic order,
   for binary search */

keyword keywds[] = {
  KW("cognitial"   , pswarm_opt_d, (Char*)&opt.mu,
     "Cognitial parameter"),
  KW("ddelta"   , pswarm_opt_d, (Char*)&opt.ddelta,
     "Deacresing delta factor (<1)"),
  KW("delta"   , pswarm_opt_d, (Char*)&opt.delta,
     "Initial delta"),
  KW("fweight"   , pswarm_opt_d, (Char*)&opt.fweight,
     "Final weight (inercial parameter)"),
  KW("idelta"   , pswarm_opt_d, (Char*)&opt.idelta,
     "Increase delta factor (>1)"),
  KW("iprint"   , pswarm_opt_i, (Char*)&opt.IPrint,
     "Print for every iprint iterations (<0 no print, =0 print final; >0 print)"),
  KW("iweight"   , pswarm_opt_d, (Char*)&opt.iweight,
     "Initial weight (inercial parameter)"),
  KW("maxf", pswarm_opt_i, (Char*)&opt.maxf,
     "Maximum number of function evaluations times problem dimension"),
  KW("maxit", pswarm_opt_i, (Char*)&opt.maxiter,
     "Maximum number of iterations times problem dimension"),
  KW("size"        , pswarm_opt_i, (Char*)&opt.s,
     "Swarm size"),
  KW("social"        , pswarm_opt_d, (Char*)&opt.nu,
     "Social parameter"),
  KW("tol"        , pswarm_opt_d, (Char*)&opt.tol,
     "Stopping tolerance parameter"),
  KW("vectorized"    , pswarm_opt_i, (Char*)&opt.vectorized,
     "Vectorized call to the objective function"),
};


struct Option_Info Oinfo = { "pswarm", "PSwarm", "pswarm_options",
                 keywds, nkeywds, 1, pswarm_version, 0, NULL};




/**********************************************
Set options. String type
**********************************************/
char *pswarm_opt_s(Option_Info *oi, keyword *kw, char *value)
{
  char *s;


  /* never echo options */
  oi->option_echo &= ~ASL_OI_echo; 


  if(!strncmp("method", kw->name, 6)){
    s=value;
    while(*s!=' ' && *s!=0)
      s++;
    if(s<=value)
      return value;

    if(!strncmp("disc_hett", value, 9)){
      *(int *)kw->info=0; /*DISC_METHOD;*/
      printf("Discretization method selected Hettich version\n");
      return s;
    }

    /* unknown method */
    return value;
  }


  /* not implemented option */
  return value;
}


/**********************************************
Set options. Integer type
**********************************************/
char *pswarm_opt_i(Option_Info *oi, keyword *kw, char *value)
{
  long optval;
  char *s;

  /* never echo options */
  oi->option_echo &= ~ASL_OI_echo;

  optval=strtol(value, &s, 10);
  if(s > value){
    /* existing integer number */
    *(int *)kw->info=(int)optval;
    printf("\nDefault option %s=%d changed\n", kw->name, *(int *)kw->info);
    return s;
  }

return value;
}


/**********************************************
Set options. Double type
**********************************************/
char *pswarm_opt_d(Option_Info *oi, keyword *kw, char *value)
{
  double optval;
  char *s;

  /* never echo options */
  oi->option_echo &= ~ASL_OI_echo;  

  optval=strtod(value, &s);

  if(s > value){
    /* existing double number */
    *(double *)kw->info=optval;
    printf("\nDefault option %s=%.6f changed\n", kw->name, *(double *)kw->info);
    return s;
  }

return value;
}

#endif /* AMPL */


static jmp_buf Jb;

void catchfpe(int n)
{
#ifdef AMPL
  report_where(asl);
#endif /* AMPL */
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

#ifdef AMPL
  char *stub;
  double *lbt, *ubt, *tmp;
  ASL *asl;
  FILE *nl;
  fint m, n, no, nz, mxr, mxc;
  cgrad *cg;
  int tmplincons;
#else
  int n;
  double *X0;
#endif
  

#ifdef AMPL
  asl = ASL_alloc(ASL_read_fg);
  
  stub = getstops(argv, &Oinfo);
  
  nl=jac_dim_ASL(asl, stub, &m, &n, &no, &nz, &mxr, &mxc, (fint)strlen(stub));
  if (!nl){
    printf("Can't read problem\n");
    exit(1);
  }

  want_deriv = 0;
  want_derivs = 0; /* no derivs */
  want_xpi0=1;    /* initial guess, if available */
  
  
  fg_read_ASL(asl, nl, 0);
  
  if(n_obj<1){
    printf("At least one objective is requested\n");
    exit(1);
  }

  if(n_obj>1){
    printf("Current implementation only supports one objective function\n");
    printf("Considering first objective\n");
  }

  if(nlc){
    printf("\nIgnoring %d nonlinear constraint%c\n",
       nlc, nlc==1? ' ':'s');
  }
  
  dense_j();
  
  objsign = objtype[0] ? -1. : 1.; // default is minimization
  /* objsign =  1 minimization problem */
  /* objsign = -1 maximization problem */

  /* Deal with linear constraints */
  if(n_con-nlc>0){
    /* Ax<=b constraint type */
    /* account for how many */
    lincons=0;
    for(i=nlc;i<n_con;i++){
      if(LUrhs[2*i]>-Inf)
		  lincons++;
      if(LUrhs[2*i+1]<Inf)
		  lincons++;
    }

    A = pswarm_malloc(lincons*n_var*sizeof(double));
    b = pswarm_malloc(lincons*sizeof(double));

    memset(A, 0, (lincons)*n_var*sizeof(double));
    memset(b, 0, (lincons)*sizeof(double));


    tmplincons=0;
    for(i=nlc;i<n_con;i++){
      if(LUrhs[2*i]>-Inf){
    b[tmplincons]=-LUrhs[2*i];
    for(cg=Cgrad[i];cg;cg=cg->next)
      A[tmplincons+lincons*cg->varno]=-cg->coef;    
    tmplincons++;
      }
      if(LUrhs[2*i+1]<Inf){
    b[tmplincons]=LUrhs[2*i+1];
    for(cg=Cgrad[i];cg;cg=cg->next)
      A[tmplincons+lincons*cg->varno]=cg->coef;
    tmplincons++;
      }

    }

  } else {
    lincons=0;
    A=NULL;
    b=NULL;
  }
  
  //opt.s = 5*n;

  //printf("**** Pop size = %d ****\n", opt.s);
#endif /* AMPL */
  
  if (!setjmp(Jb)){
    signal(SIGFPE, catchfpe);

#ifndef AMPL
    set_problem_dimension(&n,&lincons);
#endif
    
    /* lower and upper bounds on variables */

#ifdef AMPL
    lbt=
#endif
      lb=pswarm_malloc(n*sizeof(double));    


#ifdef AMPL
    ubt=
#endif
      ub=pswarm_malloc(n*sizeof(double));
    
    if(!lb || !ub){
      printf("Unable to allocate memory for variable bounds\n");
      exit(1);
    }

#ifdef AMPL
    tmp=LUv;
    for(i=0;i<n;i++){
      *lbt++=*tmp++;
      *ubt++=*tmp++;
    }
#else
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

#endif /*AMPL*/

    /* check for finite bound on the variables */
/*    for(i=0;i<n;i++){
      if(lb[i]<=-Inf || ub[i]>=Inf){
    printf("Variable %i has no finite bounds\nPlease provide finite bounds for all variables\n",i);
    exit(1);
      }
    }*/

    

#ifndef AMPL    
#ifdef MPI
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_numprocs);

    user_init_MPI(MPI_myrank);

#else 

    user_init();

#endif
#endif /* AMPL */



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
      
#ifdef AMPL
      write_sol("Exiting PSwarm", sol, NULL, &Oinfo);
#endif
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


#ifdef AMPL
void objfun(int n, int m, double *x, double *lb, double *ub, double *fx)
{
  int j;

  for(j=0;j<m;j++){
	  fx[j]=objval(0, &x[j*n], &NERROR);
  }
  

}

#endif /* AMPL */

