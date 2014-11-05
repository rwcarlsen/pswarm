##################################################################
#
#  Run Pswarm with a test problem
#
##################################################################
#
# 02-07-2008 aivaz@dps.uminho.pt
#
# See README file for instructions
#
##################################################################

# Read problem definition, functions and options
from hs024 import *

# Load the solver
from pswarm_py import pswarm


# Call the solver
result = pswarm(Problem,Options)

print "Returning value" # zero means successful
print result['ret']

print "Best particle "
print result['x']

print "Objective value "
print result['f']
