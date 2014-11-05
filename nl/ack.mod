############################################################################
#
# Ackley problem
# as described in
#   M.M.  Ali,  C.  Khompatraporn  and  Z.B.  Zabinsky,   "A numerical
#   evaluation of several stochastic algorithms on selected continuous
#   global optimization test problems", Journal of Global Optimization
#   (2005) 31:635-672
#
# Coded for AMPL by: aivaz@dps.uminho.pt
# Minho University, Portugal
# 14 September 2005
############################################################################

# parameters
param n := 10;
param pi := 4*atan(1);

# variables
var x{1..n};

# objective function, minimization problem
minimize fx:
    -20*exp(-0.02*sqrt((1/n)*(sum {i in 1..n} (x[i]^2))))
    -exp((1/n)*sum {i in 1..n} (cos(2*pi*x[i])))
    +20+exp(1);

subject to bounds {i in 1..n}:
    -30 <= x[i] <= 30;

option solver pswarm;
solve;

display fx;
display x;

printf "Unknown number of local minima and a global minima at the origin with fx=0\n";
