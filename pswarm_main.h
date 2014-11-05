#ifdef MPI
#include "mpi.h"
#endif
#ifdef MPE
#include "mpe.h"
#endif
#ifdef AMPL
#include "nlp.h"
#include "asl.h"
#include "getstub.h"
#else
#include <stdio.h>
#endif

#include <signal.h>



#ifdef AMPL
char *pswarm_opt_s(Option_Info *, keyword *, char *);
char *pswarm_opt_i(Option_Info *, keyword *, char *);
char *pswarm_opt_d(Option_Info *, keyword *, char *);
#endif

