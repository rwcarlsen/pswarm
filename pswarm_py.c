#include "pswarm_py.h"



static jmp_buf Jb;

static PyMethodDef pswarm_methods[]={
  {"pswarm", pswarm_py, METH_VARARGS,
   "Call PSwarm solver."},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initpswarm_py(void) {

  Py_InitModule("pswarm_py", pswarm_methods);
  import_array();
}


void catchfpe(int n)
{
  printf("\nFloating point error.\n");
  fflush(stdout);
  longjmp(Jb,1);
}

void PrintRealVector(char *txt, int n, double *x)
{
  int i;

  if(txt==NULL || x==NULL)
    return;

  printf("%s=[", txt);

  for(i=0;i<n-1;i++)
    printf("%f,",x[i]);

  printf("%f]\n",x[i]);

}

void getRealOption(PyObject *Options, char *opt_name, double *option)
{
  PyObject *py_opt=NULL;
    
  py_opt = PyDict_GetItemString(Options, opt_name);
  if(!py_opt)
    return;
  
  if(!PyFloat_Check(py_opt)){
    printf("%s option must be a float\n", opt_name);
    return;
  }
  
  *option = PyFloat_AsDouble(py_opt);
}

void getIntOption(PyObject *Options, char *opt_name, int *option)
{
  PyObject *py_opt=NULL;

  py_opt = PyDict_GetItemString(Options, opt_name);
  if(!py_opt)
    return;
  
  if(!PyInt_Check(py_opt)){
    printf("%s option must be integer\n", opt_name);
    return;
  }
  
  *option = PyInt_AsLong(py_opt);
}


int getPyRealVector(char *txt, int n, PyObject *obj, double *x) {

  int i;
  PyArrayObject *array = NULL;
  char msn[256];

  if(x==NULL || obj==NULL)
    return 1;
  
  Py_INCREF(obj);

  array = (PyArrayObject *) PyArray_ContiguousFromAny(obj, PyArray_DOUBLE, 0,0);

  if(array == NULL){
    PyErr_SetString(PyExc_ValueError, "Null array");    
    goto clean_exit;
  }

  if(PyArray_DIM(array,0) != n){
    sprintf(msn, "Array '%s' is of wrong size. Expected %d and got %d.",
        txt, n, (int)PyArray_DIM(array,0));
    PyErr_SetString(PyExc_ValueError, msn);
    goto clean_exit;
  }

  for (i = 0; i < n; i++)
    x[i] = ((double*)PyArray_DATA(array))[i];

  Py_XDECREF(obj);
  Py_XDECREF(array);
  return 0;

 clean_exit:
  
  Py_XDECREF(obj);
  Py_XDECREF(array);
  return 1;
}


int getPyRealMatrix(char *txt, int n, int m, PyObject *obj, double *A) {

  int i,j;
  char msn[256];
  PyObject *py_row=NULL;
#ifdef linux
  double row[n];
#else
  double *row;
  
  row=malloc(n*sizeof(double));
#endif

  if(A==NULL || obj==NULL)
    return 1;
  
  Py_INCREF(obj);

  for(i=0;i<m;i++){
    sprintf(msn, "%s[%d]", txt, i);
    py_row=PyList_GetItem(obj,i);
      if(getPyRealVector(msn, n, py_row, row)){
    goto clean_exit;
      }
      for(j=0;j<n;j++)
    A[i+j*m]=row[j];
  }

#ifndef linux
  free(row);
#endif
  Py_XDECREF(obj);
  return 0;

 clean_exit:
#ifndef linux
  free(row);
#endif
  Py_XDECREF(obj);
  return 1;
}


void py_objfun(int n, int m, double *x, double *lb, double *ub, double *fx)
{
  PyObject *result = NULL, *py_x = NULL;
  int i,j,k;
  npy_intp dim[2];
  char msn[256];

  if(x==NULL || m==0)
    return;

  /* pswarm controls bound feasibility, but just in case... */
  for(j=0;j<m;j++){
	for(i=0;i<n;i++){
		if(x[j*n+i]<lb[i] || x[j*n+i]>ub[i]){
			PySys_WriteStdout("Error computing objective function for unfeasible bound point\nReturning all as infinity.");
			for(k=0;k<m;k++)
				fx[k]=1e20;
			return;
		}
	}
  }

  /* Build a list for calling the objective function */
  dim[0]=m;
  dim[1]=n;
  py_x=PyArray_SimpleNewFromData(2, dim, PyArray_DOUBLE, x);
  if(py_x==NULL){
	  PySys_WriteStdout("Error making objective argument for objective function\n");
	  for(j=0;j<m;j++)
		  fx[j]=1e20;
	  return;
  }

  if((result=PyEval_CallFunction(py_objf, "(O)", py_x))==NULL){
	  PySys_WriteStdout("Error calling Python objective function\n");
	  for(j=0;j<m;j++)
		  fx[j]=1e20;
	  Py_DECREF(py_x);
	  return;
  }

  if(PyArray_DIM(result,0) != m){
    sprintf(msn, "Objective function returned bad function vector. Expected %d and got %d.",
        m, (int)PyArray_DIM(result,0));
    PyErr_SetString(PyExc_ValueError, msn);
	for(j=0;j<m;j++)
		fx[j]=1e20;
	Py_DECREF(result);
	Py_DECREF(py_x);
	return;
  }

  for (j = 0; j < m; j++)
    fx[j] = ((double*)PyArray_DATA(result))[j];

  Py_XDECREF(py_x);
  Py_XDECREF(result);
    
  return;
}


/* returning a negative value causes pswarm to stop */
double py_outfcn(int n, int s, int iter, int gbest, struct swarm *pop)
{
  PyObject *py_iter = NULL, *py_leader=NULL, *py_fx=NULL, *py_x = NULL, *result=NULL;
  npy_intp dim[1];
  double fx;


  fx=1.0; /* do not exit PSwarm */

  if(py_outf){
	  /* Build a list for calling the output function */
	  dim[0]=1;
      py_iter=PyArray_SimpleNewFromData(1, dim, PyArray_INT, &iter);
      if(py_iter==NULL)
        goto clean_exit;

      py_leader=PyArray_SimpleNewFromData(1, dim, PyArray_INT, &gbest);
      if(py_leader==NULL)
		  goto clean_exit;

      py_fx=PyArray_SimpleNewFromData(1, dim, PyArray_DOUBLE, &(pop->fy[gbest]));
      if(py_fx==NULL)
		  goto clean_exit;

      dim[0]=n;
      py_x=PyArray_SimpleNewFromData(1, dim, PyArray_DOUBLE, &(pop->y[gbest*n]));
      if(py_x==NULL)
		  goto clean_exit;

      if((result=PyEval_CallFunction(py_outf, "(OOOO)", py_iter, py_leader, py_fx, py_x))==NULL){
	    PySys_WriteStdout("Error calling outputfcn\n");
		goto clean_exit;
	  }
	  
	  if(!PyFloat_Check(result) && !PyInt_Check(result)){
		goto clean_exit;
	  }
	  
	  if(PyFloat_Check(result)){
		  fx=PyFloat_AsDouble(result);
	  } else {
		  fx=(double)PyInt_AsLong(result);
	  }

clean_exit:
	  Py_XDECREF(py_iter);
      Py_XDECREF(py_leader);
      Py_XDECREF(py_fx);
      Py_XDECREF(py_x);
	  Py_XDECREF(result);

	  return fx;
  
  } else {
	  /* print to stdout */
	  	if(iter==0){
			PySys_WriteStdout("\n  Iter     Leader     Objective  ");
			PySys_WriteStdout("\n  -------------------------------\n");
		}

		PySys_WriteStdout("    %4d   %4d   %4.6e\n", iter, gbest, pop->fy[gbest]);

		return 1.0;
  }
  

  return 1.0;  /* never reached */
}


static PyObject *pswarm_py(PyObject *self, PyObject *args)
{
  double *lb=NULL, *ub=NULL;
  double *A=NULL, *b=NULL, *x0=NULL;
  int n, i, lincons, exit_code=0;
  double *sol=NULL;
  double f;

  PyObject *py_n, *py_lb, *py_ub, *py_x0, *py_A, *py_b, *py_res, *py_f,*py_ret;
  PyArrayObject *py_x;
  npy_intp dim[1];
  
  if(!saved_options){
    memcpy(&opt_backup, &opt, sizeof(struct Options));
    saved_options++;
    // PySys_WriteStdout("Saving options\n");
  } else {
    memcpy(&opt, &opt_backup, sizeof(struct Options));
    // PySys_WriteStdout("Recovering saved options\n");
  }


  /* Process input arguments */

  if(!PyArg_ParseTuple(args,"O!O!", &PyDict_Type, &Problem, &PyDict_Type, &Options))
    return NULL; /* return error */

  Py_INCREF(Problem);
  Py_INCREF(Options);

  //PyObject_Print(Problem,stdout,0);

  /* Process Problem definition */
  py_objf = PyDict_GetItemString(Problem,"objf");
  if(!py_objf || !PyFunction_Check(py_objf)){
    PyErr_SetString(PyExc_ValueError,"objf must be defined as a function object");
    goto clean_exit;
  } else {
    Py_INCREF(py_objf);
  }
  
  py_n = PyDict_GetItemString(Problem,"Variables");
  if(!py_n || !PyInt_Check(py_n)){
    PyErr_SetString(PyExc_ValueError,"The number of Variables must be provided as an integer");
    goto clean_exit;
  } else {
    Py_INCREF(py_n);
  }
  
  n = PyInt_AsLong(py_n);
  if(n<=0){
    PyErr_SetString(PyExc_ValueError,"The number of Variables must be positive");
    Py_DECREF(py_n);
    goto clean_exit;
  }
  Py_DECREF(py_n);


  if((lb=(double *) malloc(n*sizeof(double)))==NULL){
    PyErr_SetString(PyExc_ValueError,"Unable to allocate memory for lb");
    goto clean_exit;
  }

  py_lb = PyDict_GetItemString(Problem,"lb");
  if(!py_lb){
    for(i=0;i<n;i++)
      lb[i]=-1e20;
  } else {
    Py_INCREF(py_lb);

    if(getPyRealVector("lb", n, py_lb, lb)){
      Py_DECREF(py_lb);
      goto clean_exit;
    } else {
      Py_DECREF(py_lb);
    }
  }


  if((ub=(double *) malloc(n*sizeof(double)))==NULL){
    PyErr_SetString(PyExc_ValueError,"Unable to allocate memory for ub");
    goto clean_exit;
  }

  py_ub = PyDict_GetItemString(Problem,"ub");
  if(!py_ub){
    for(i=0;i<n;i++)
      ub[i]=1e20;
  } else {
    Py_INCREF(py_ub);
    
    if(getPyRealVector("ub", n, py_ub, ub)){
      Py_DECREF(py_ub);
      goto clean_exit;
    } else {
      Py_DECREF(py_ub);
    }
  }


  py_x0 = PyDict_GetItemString(Problem,"x0");
  if(!py_x0){
    x0=NULL;
  } else {
    Py_INCREF(py_x0);
    
    if((x0=(double *) malloc(n*sizeof(double)))==NULL){
      PyErr_SetString(PyExc_ValueError,"Unable to allocate memory for x0");
      Py_DECREF(py_x0);
      goto clean_exit;
    }

    if(getPyRealVector("x0", n, py_x0, x0)){
      Py_DECREF(py_x0);
      goto clean_exit;
    } else {
      Py_DECREF(py_x0);
    }
  }


  py_A = PyDict_GetItemString(Problem,"A");
  if(!py_A){
    lincons=0;
    A=NULL;
    b=NULL;
  } else {
    Py_INCREF(py_A);

    if(!PyList_Check(py_A)){
      PyErr_SetString(PyExc_ValueError,"A must be defined as a list of arrays (rows)");
      Py_DECREF(py_A);
      goto clean_exit;
    }

    lincons=(int)PyList_Size(py_A);

    if((A=(double *) malloc(lincons*n*sizeof(double)))==NULL){
      PyErr_SetString(PyExc_ValueError,"Unable to allocate memory for A");
      Py_DECREF(py_A);
      goto clean_exit;
    }

    if(getPyRealMatrix("A", n, lincons, py_A, A)){
      Py_DECREF(py_A);
      goto clean_exit;
    } else {
      Py_DECREF(py_A);
    }
    

    py_b = PyDict_GetItemString(Problem,"b");
    if(!py_b){
      PyErr_SetString(PyExc_ValueError,"b must be defined as an array");
      goto clean_exit;
    } else {
      Py_INCREF(py_b);
      
      if((b=(double *) malloc(lincons*sizeof(double)))==NULL){
    PyErr_SetString(PyExc_ValueError,"Unable to allocate memory for b");
    Py_DECREF(py_b);
    goto clean_exit;
      }

      if(getPyRealVector("b", lincons, py_b, b)){
    Py_DECREF(py_b);
    goto clean_exit;
      } else {
    Py_DECREF(py_b);
      }
    }

    
  }

  py_outf = PyDict_GetItemString(Options,"outputfcn");
  if(py_outf && !PyFunction_Check(py_outf)){
		  PyErr_SetString(PyExc_ValueError,"outputfcn must be defined as a function object");
		  goto clean_exit;
  };


  if(py_outf)
	  Py_INCREF(py_outf);
  opt.outfcn=&py_outfcn;
		  

  //PrintRealVector("b", lincons, b);

  //PyObject_Print(Options,stdout,0);
  getRealOption(Options, "cognitial", &(opt.mu));
  getRealOption(Options, "fweight",   &(opt.fweight));
  getRealOption(Options, "iweight",   &(opt.iweight));
  getIntOption (Options, "maxf",      &(opt.maxf));
  getIntOption (Options, "maxit",     &(opt.maxiter));
  getIntOption (Options, "size",      &(opt.s));
  getIntOption (Options, "iprint",    &(opt.IPrint));
  getRealOption(Options, "social",    &(opt.nu));
  getRealOption(Options, "tol",       &(opt.tol));
  getRealOption(Options, "delta",     &(opt.delta));
  getRealOption(Options, "ddelta",    &(opt.ddelta));
  getRealOption(Options, "idelta",    &(opt.idelta));
  getIntOption (Options, "vectorized",&(opt.vectorized));



  if (!setjmp(Jb)){
    signal(SIGFPE, catchfpe);

    exit_code=PSwarm(n, &py_objfun, lb, ub, lincons, A, b, &sol, &f, x0);
    
/*	if(opt.IPrint>=0){
      if(exit_code){
        PySys_WriteStdout("Abnormal exit\n");
	  } else {
        PySys_WriteStdout("Normal exit\n");
	  }
	}*/
  
  }

  /* Build a list with the solution */
  py_res=PyDict_New();

  Py_INCREF(py_res);

    
  py_ret=PyInt_FromLong(exit_code);
  if(py_ret==NULL){
    PyErr_SetString(PyExc_ValueError,"Unable to create the return value object");
    goto clean_exit;
  }
  if(PyDict_SetItemString(py_res, "ret", py_ret)){
    PyErr_SetString(PyExc_ValueError,"Unable to insert the return value object into dict");
    goto clean_exit;
  }
  Py_INCREF(py_ret);


  
  py_f=PyFloat_FromDouble(f);
  if(py_f==NULL){
    PyErr_SetString(PyExc_ValueError,"Unable to create the solution objective value object");
    goto clean_exit;
  }
  if(PyDict_SetItemString(py_res, "f", py_f)){
    PyErr_SetString(PyExc_ValueError,"Unable to insert the solution objective value object into dict");
    goto clean_exit;
  }
  Py_INCREF(py_f);


  dim[0]=n;
  py_x=(PyArrayObject *)PyArray_SimpleNewFromData(1, dim, PyArray_DOUBLE, sol);
  if(py_x==NULL){
    PyErr_SetString(PyExc_ValueError,"Unable to return the solution value");
    goto clean_exit;
  }
  if(PyDict_SetItemString(py_res, "x", (PyObject *)py_x)){
    PyErr_SetString(PyExc_ValueError,"Unable to insert the solution value object into dict");
    goto clean_exit;
  }

  Py_INCREF((PyObject *)py_x);

  py_x->flags |= NPY_OWNDATA; /* to be dealocated with py_x */


  /*
  PrintRealVector("sol=", n, sol);
  PyObject_Print(py_x,stdout,0);
  printf ("\n");
  PyObject_Print(py_res,stdout,0);
  printf ("\n");
  */

  if(lb)
    free(lb);
  if(ub)
    free(ub);
  if(A)
    free(A);
  if(b)
    free(b);
  if(x0)
    free(x0);

  /* We may not dealocate sol, as it may be used in an object */


  Py_XDECREF(py_objf);
  if(py_outf)
	  Py_XDECREF(py_outf);
  Py_XDECREF(Problem);
  Py_XDECREF(Options);


  Py_INCREF(py_res);
  return py_res;

 clean_exit:
  if(lb)
    free(lb);
  if(ub)
    free(ub);
  if(A)
    free(A);
  if(b)
    free(b);
  if(x0)
    free(x0);

  Py_XDECREF(py_objf);
  Py_XDECREF(Problem);
  Py_XDECREF(Options);
  return NULL;
  
}
