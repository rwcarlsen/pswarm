#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<memory.h>

#include "pattern.h"
#include "pswarm.h"

extern struct swarm pop;
extern struct Options opt;
extern struct Stats stats;

struct poll_vector *PVectors=NULL; /* lninked list of poll vectors */
struct poll_vector *D=NULL; /* linked list of poll vectors maximal basis */
struct poll_vector *last_D=NULL;
struct poll_vector *TC=NULL; /* extra poll vectors for tangent cone */

#ifdef LINEAR
extern void dgemv_();
extern void dgemm_();
extern void dorgqr_();
extern void dgeqrf_();
extern void dgesv_();
extern double dnrm2_();

void tangent_cone(int n, int lincons, double *A, double *b, double*x, double *lb, double *ub);
#endif

extern int feasible_p(int n, double *x, int lincons, double *A, double *b, double *lb, double *ub);

extern void *pswarm_malloc(size_t size);


/* performe a poll step for y[pi] position */
void pollstep(int n, int lincons, int pi, void (*objf)(), double *lb, double *ub, double *A, double *b, struct poll_vector **last_success)
{
  int i,j;
  double *poll_point, fx, minfx;
  struct poll_vector *tmp, *minvector;
  double *vectorx, *vectorfx;

#ifdef LINEAR
  int feasible;

  double none=-1.0;
  double zero=0.0;
  char Normal='N';
  int oneI=1;

  double *tmpm;
#endif

#ifdef LINEAR /* get tangent cone for active constraints */

    tmpm=pswarm_malloc(lincons*sizeof(double));

	feasible=1;
	/* check for linear feasibility */
	if(lincons){
		dgemv_(&Normal, &lincons, &n, &none, A, &lincons, &pop.y[pi*n], &oneI, &zero, tmpm, &oneI); /* tmpm = -Ax */
		for(i=0;i<lincons;i++)
			if(b[i]+tmpm[i]<-opt.tol){
				printf("Linear constraint %d is %f and should be greater than zero\n", i, b[i]+tmpm[i]);
				feasible=0;
			}

			if(!feasible){
				printf("Leader particle %d is not linear feasible\nThis should not happen!!!\n", pi);
			}
	}


    free(tmpm);
    
	for(i=0;i<n;i++)
		if(pop.y[pi*n+i]>ub[i] || pop.y[pi*n+i]<lb[i]){
			printf("Leader Particle is not bound feasible\nThis should not happen!!!\n");
			feasible=0;
		}

		if(lincons) /* TC is already a NULL pointer */
			tangent_cone(n, lincons, A, b, &pop.y[pi*n], lb, ub);

#endif



  poll_point=pswarm_malloc(n*sizeof(double));

  if(!poll_point){
    printf("Unable to alocate memory in pattern.c\n");
    exit(1);
  }




  /* performe a poll step for each poll vector in D or TC */

  /* Use D if TC is NULL */
  if(TC){
    PVectors=TC;
  } else {
    PVectors=D;
  }

/*   if(TC){ */
/*     printf("*********************************************\n"); */
/*     print_TC(n); */
/*     printf("*********************************************\n"); */
/*   } */


if(opt.vectorized){

	/* account for the number of trial points */
	tmp=PVectors;
	j=0;
	while(tmp){
		j++;
		tmp=tmp->next;
	}

	//printf("I have %d possible points to evaluate\n",j);

	/* allocate memory, we may be allocating more memory than necessary */
	vectorx=(double *)pswarm_malloc(j*n*sizeof(double));
	vectorfx=(double *)pswarm_malloc(j*sizeof(double));

	/* compute trial points */
	tmp=PVectors;
	j=0;
	while(tmp){
		for(i=0;i<n;i++)
			vectorx[j*n+i]=pop.y[pi*n+i]+pop.delta*tmp->vector[i];

		/* is it linear feasible ? */
		if(feasible_p(n, &vectorx[j*n], lincons, A, b, lb, ub))
			j++; /* yes */

		tmp=tmp->next;
	}

	//printf("And %d are feasible\n", j);

	if(j>0){
		objf(n, j, vectorx, lb, ub, vectorfx);
		stats.objfunctions+=j;
	}

	minvector=NULL;
	minfx=Inf;
	j=0;
	tmp=PVectors;
	while(tmp){
		for(i=0;i<n;i++)
			poll_point[i]=pop.y[pi*n+i]+pop.delta*tmp->vector[i];

		if(feasible_p(n, poll_point, lincons, A, b, lb, ub)){
			if(minfx>vectorfx[j]){
				minfx=vectorfx[j];
				minvector=tmp; /* don't break. We are in a non oportunistic version */
			}
			j++;
		}
		tmp=tmp->next;
	}

	free(vectorx);
	free(vectorfx);

} else {
	tmp=PVectors;
	minvector=NULL;
	minfx=Inf;
	while(tmp){
		for(i=0;i<n;i++)
			poll_point[i]=pop.y[pi*n+i]+pop.delta*tmp->vector[i];

		if(feasible_p(n, poll_point, lincons, A, b, lb, ub)){
			objf(n, 1, poll_point, lb, ub, &fx);
			stats.objfunctions++;
			if(minfx>fx){
				minfx=fx;
				minvector=tmp; /* if oportunistic then break if fx < leader fx */
				if(pop.fy[pi]>minfx)
					break;
			}
		}
		tmp=tmp->next;
	}
}
  
  if(pop.fy[pi]>minfx){ /* a successful poll point */
    stats.sucpollsteps++;
    for(i=0;i<n;i++)
      pop.y[pi*n+i]=pop.y[pi*n+i]+pop.delta*minvector->vector[i];
    pop.fy[pi]=minfx;
    if(*last_success==minvector){ /* last success vector twice, increase delta */
      pop.delta*=opt.idelta;
      //printf("Increasing delta in poll step %f\n", pop.delta);
    } else { /* last success vector is different */
      *last_success=minvector;
    }
  } else {
    pop.delta*=opt.ddelta;
    //printf("Decreasing delta\n");
    *last_success=NULL;
  }
  
  
  /* free variables */
  free(poll_point);
  
  PVectors=NULL;
  
}


void init_pattern(int n)
{
  
  init_D(n);
  //    print_D(n);

  TC=NULL;
  PVectors=NULL;
  
}

void clean_pattern()
{
  
  clean_D();
  clean_TC();
  PVectors=NULL;
  
  
}


void clean_D()
{
  struct poll_vector *tmp1, *tmp2;
  
  tmp1=D;
  while(tmp1){
    tmp2=tmp1->next;
    free(tmp1->vector);
    free(tmp1);
    tmp1=tmp2;
  }
  
  D=NULL;
}


void clean_TC()
{
  struct poll_vector *tmp1, *tmp2;
  
  tmp1=TC;
  while(tmp1){
    tmp2=tmp1->next;
    if(tmp1->vector)
      free(tmp1->vector);
    free(tmp1);
    tmp1=tmp2;
  }
  
  TC=NULL;
}


#ifdef LINEAR /* get tangent cone for active constraints */

void tangent_cone(int n, int lincons, double *A, double *b, double*x, double *lb, double *ub)
{
  double ActiveThreshold;
  double ActiveThresholdLimit;

  double none=-1.0;
  double zero=0.0;
  double one=1.0;
  int oneI=1;
  char Normal='N';
  double *tmpnm;
  double *tmpmn;
  double *Active;
  double *R;
  double tmpnrm;
  char Transpose='T';
  int mylwork=2*n, info;
  int i, j, activetmp, NActive;
  

  int *myworkI;
  double *tmpm;
  double *tau;
  double *N;
  double *tmpdir;
  double *mywork;
  
  myworkI=pswarm_malloc(n*sizeof(int));
  tmpm=pswarm_malloc(lincons*sizeof(double));
  tau=pswarm_malloc(n*sizeof(double));
  N=pswarm_malloc(n*n*sizeof(double));
  tmpdir=pswarm_malloc(n*sizeof(double));
  mywork=pswarm_malloc(2*n*sizeof(double));

  if(TC) /* clean previous tangent cone directions */
    clean_TC();

  ActiveThreshold=10*pop.delta;
  if(ActiveThreshold>opt.EpsilonActive)
    ActiveThreshold=opt.EpsilonActive;

  ActiveThresholdLimit=pow(ActiveThreshold, 2.0);
  if(ActiveThresholdLimit>0.1)
    ActiveThresholdLimit=0.1;

  while(ActiveThreshold>=ActiveThresholdLimit){

    dgemv_(&Normal, &lincons, &n, &none, A, &lincons, x, &oneI, &zero, tmpm, &oneI); /* tmpm = -Ax */

    NActive=0;
    for(i=0;i<lincons;i++)
      if(b[i]+tmpm[i]<=ActiveThreshold*pop.scale) /* i is epsilon active */
    NActive++;

    for(i=0;i<n;i++){
      if(x[i]-lb[i]<=ActiveThreshold*pop.scale)
    NActive++;
      if(ub[i]-x[i]<=ActiveThreshold*pop.scale)
    NActive++;
    }

    if(NActive>0 && NActive<=n){ /* we have at least a squared matrix (may be singular) */
      Active=pswarm_malloc(NActive*n*sizeof(double));
      tmpnm=pswarm_malloc(NActive*n*sizeof(double));
      tmpmn=pswarm_malloc(NActive*n*sizeof(double));
      R=pswarm_malloc(NActive*NActive*sizeof(double));
      
      if(!Active || !tmpnm || !tmpmn || !R){
    printf("Unable to allocate memory for the Active linear constraints\n");
    exit(1); /* we shouldn't be so drastic */
      }
      
      activetmp=0;
      memset(Active, 0, NActive*n*sizeof(double)); /* zero out Active */
      /* form Active matrix - transposed */
      for(i=0;i<lincons;i++){
    if(b[i]+tmpm[i]<=ActiveThreshold*pop.scale){ /* i is epsilon active */
      for(j=0;j<n;j++)
        Active[activetmp*n+j]=A[i+j*lincons]; /* transpose */
      activetmp++;
    }
      }

      for(i=0;i<n;i++){
    if(x[i]-lb[i]<=ActiveThreshold*pop.scale){ /* - identity */
      Active[activetmp*n+i]=-1.0;
      activetmp++;
    }
    if(ub[i]-x[i]<=ActiveThreshold*pop.scale){ /* identity */
      Active[activetmp*n+i]=1.0;
      activetmp++;
    }
      }

      memcpy(tmpnm, Active, n*activetmp*sizeof(double));
      /* perform a QR factorization */
      dgeqrf_(&n, &NActive, tmpnm, &n, tau, mywork, &mylwork, &info);

      if(!info){
    memset(R, 0, NActive*NActive*sizeof(double));
    for(i=0;i<NActive;i++){
      for(j=i;j<NActive;j++)
        R[j*NActive+i]=tmpnm[j*n+i];
    }
    
    /* min(NActive,n)=NActive */
    
    dorgqr_(&n, &NActive, &NActive, tmpnm, &n, tau, mywork, &mylwork, &info); /* tmpnm = Q(n,m) */

    if(!info){
    
      /* Transpose tmpnm to obtain Q'=tmpmn*/
      for(i=0;i<NActive;i++)
        for(j=0;j<n;j++){
          tmpmn[j*NActive+i]=tmpnm[i*n+j];
        }     
      
      /* Solve linear system R*B'=Q' */
      
      for(i=0;i<n;i++)
        myworkI[i]=i+1;
      
      dgesv_(&NActive, &n, R, &NActive, myworkI, tmpmn, &NActive, &info); /* tmpmn = B' */

      if(!info){ /* solution computed */
      
        /* Build N */
        
        memset(N, 0, n*n*sizeof(double));
        for(i=0;i<n;i++) /* N = Identity */
          N[i*n+i]=1.0;
        dgemm_(&Transpose, &Transpose, &n, &n, &NActive, &none, tmpmn, &NActive, Active, &n, &one, N, &n); 
        
        /* generate all directions -- [N -N B -B]*/
        
        for(i=0;i<NActive;i++){ /* insert B and -B */
          for(j=0;j<n;j++)
        tmpdir[j]=tmpmn[i+j*NActive];
          tmpnrm=dnrm2_(&n, tmpdir, &oneI);
          if(tmpnrm>=opt.tol){
        insert_TC(n, tmpdir);
          }
        }
        
        
        for(i=0;i<n;i++){ /* insert N and -N */
          for(j=0;j<n;j++)
        tmpdir[j]=N[i+j*n];
          tmpnrm=dnrm2_(&n, tmpdir, &oneI);
          if(tmpnrm>=opt.tol){
        insert_TC(n, tmpdir);
          }
        }

      } /* Closes linear system solution */

    } /* Closes Q computation */
      
      } /* Closes QR factorization */
      
      /* generate structure with directions */
      
      free(Active);
      free(tmpnm);
      free(tmpmn);
      free(R);
      break;
    } else {
      if(NActive<=0)
    break;
      ActiveThreshold/=2;
    }
    
  }


  free(myworkI);
  free(tmpm);
  free(tau);
  free(N);
  free(tmpdir);
  free(mywork);

}


void insert_TC(int n, double *dir)
{
  struct poll_vector *tmp1, *tmp2;
  int i;

  if(!dir)
    return;

  /* should we check for repeated directions ? */

  tmp1=pswarm_malloc(sizeof(struct poll_vector));
  tmp2=pswarm_malloc(sizeof(struct poll_vector));
  if(!tmp1 || !tmp2){
    printf("Unable to allocate memory for vector in tangent cone\nAborting!!\n");
    exit(1);
  }

  tmp1->vector=pswarm_malloc(n*sizeof(double));
  tmp2->vector=pswarm_malloc(n*sizeof(double));
  if(!tmp1->vector || !tmp2->vector){
    printf("Unable to allocate memory for vector in tangent cone\nAborting!!\n");
    exit(1);
  }

  memcpy(tmp1->vector, dir, n*sizeof(double));
  for(i=0;i<n;i++)
    tmp2->vector[i]=-dir[i];

  tmp1->next=tmp2;
  if(TC){
    tmp2->next=TC;
  } else {
    tmp2->next=NULL;
  }

  TC=tmp1; /* insert in top of list */
    
}

#endif

void init_D(int n)
{
  int i;
  struct poll_vector *tmp1, *tmp2;
  
  
  if(D) /* D already initialized */
    return;
  
  switch(opt.pollbasis){
  default:
    printf("\n Poll basis order not defined\n");
    printf("\n Using I -I order\n");
  case N2: /* I -I */
    D=pswarm_malloc(sizeof(struct poll_vector));
    D->next=NULL;
    D->vector=pswarm_malloc(n*sizeof(double));
    memset(D->vector, 0, n*sizeof(double));
    D->vector[0]=+1.0;
    tmp2=D;
    for(i=1;i<2*n;i++){
      tmp1=pswarm_malloc(sizeof(struct poll_vector));
      tmp2->next=tmp1;
      tmp1->vector=pswarm_malloc(n*sizeof(double));
      memset(tmp1->vector, 0, n*sizeof(double)); /* clear memory */
      if(i<n)
    tmp1->vector[i]=+1.0;
      else
    tmp1->vector[i-n]=-1.0;
      tmp1->next=NULL ;
      tmp2=tmp1;
    }
    /* add e -e to basis */
    /*     tmp1=malloc(sizeof(struct poll_vector)); */
    /*     tmp2->next=tmp1; */
    /*     tmp1->vector=malloc(n*sizeof(double)); */
    /*     for(i=0;i<n;i++) */
    /*       tmp1->vector[i]=+1.0; */
    /*     tmp1->next=NULL; */
    /*     tmp2=tmp1; */
    /*     tmp1=malloc(sizeof(struct poll_vector)); */
    /*     tmp2->next=tmp1; */
    /*     tmp1->vector=malloc(n*sizeof(double)); */
    /*     for(i=0;i<n;i++) */
    /*       tmp1->vector[i]=-1.0; */
    /*     tmp1->next=NULL; */
    /*     tmp2=tmp1; */
    last_D=tmp2;
  }
}


void print_D(n)
{
  struct poll_vector *tmp;
  
  tmp=D;
  while(tmp){
    print_poll_vector(n, tmp->vector);
    tmp=tmp->next;
  }
  
}

void print_TC(n)
{
  struct poll_vector *tmp;
  
  tmp=TC;
  while(tmp){
    print_poll_vector(n, tmp->vector);
    tmp=tmp->next;
  }
  
}


void print_poll_vector(int n, double *vector)
{
  int i;
  
  
  if(!vector)
    return;
  
  printf("D=(%.2f, ", vector[0]);
  for(i=1;i<n;i++)
    printf("%.2f ", vector[i]);
  printf(")\n");
  
}
