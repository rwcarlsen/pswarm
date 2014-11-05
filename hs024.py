from math import sqrt;
from numpy import zeros;

def py_objf(x):
	"""The objective function."""
	fx = zeros(len(x))
	for i in range(len(x)):
		fx[i]= (pow(x[i][0]-3.0,2.0)-9.0)*pow(x[i][1],3.0)/(2.07*sqrt(3.0));
	return fx;

def py_outf(it, leader, fx, x):
	"""The output function. Return a negative number to stop PSwarm"""
	if(it==0):
		print "  Iter     Leader     Objective"
		print "  ------------------------------"
	print '    %4d   %4d   %4.6e' % (it[0], leader[0],  fx[0])
	return 1.0; # do not stop



Problem = {
	'Variables':  2,
	'objf': py_objf,
	'lb': [0.0, 0.0],
#    'ub': [1e20,1e20], # same as not defined
	'A': [[-1.0/sqrt(3), 1], [-1.0, sqrt(3)], [1.0, sqrt(3)]],
	'b': [0,0,6],
	'x0': [1, 0.5]
	}

# Optimal solution
# Problem$x0 <- array(c(3, 1.73205))

# Options list already defined

Options = {
# These values are the defaults
	'maxf': 2000,
	'maxit': 2000,
	'social': 0.5,
	'cognitial': 0.5,
	'fweight': 0.4,
	'iweight': 0.9,
	'size': 42,
	'iprint': 10,
	'tol': 1E-5,
#    'delta': 2, # delta is problem specific (computed)
	'ddelta': 0.5,
	'idelta': 2.0,
	'outputfcn': py_outf,
        'vectorized': 1
	}
