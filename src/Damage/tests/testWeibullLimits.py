#-------------------------------------------------------------------------------
# Test out our derivation of how the expected limits of a given Weibull
# realization should go.
#-------------------------------------------------------------------------------
from Spheral import *
from math import *
import random
import matplotlib.pyplot as plt

# Use amsmath with latex in matplotlib
plt.rc('text', usetex=True)
plt.rc('text.latex',
       preamble=r'\usepackage{amsmath}')

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(V = 100.0,          # volume
            k = 1e10,           # k Weibull constant (coefficient)
            m = 9.0,            # m Weibull constant (exponent)
            Nflaws = 5000,      # number of flaws
            Ntrials = 5000,     # number of trials
            seed = 49489212)    # random number seed

mInv = 1.0/m
chi = (Nflaws/(k*V))**mInv

#-------------------------------------------------------------------------------
# Generate an example flaw distribution
#-------------------------------------------------------------------------------
random.seed(seed)
minflaws0, maxflaws0 = [2.0]*Ntrials, [-1.0]*Ntrials
for itrial in xrange(Ntrials):
    for i in xrange(Nflaws):
        fi = chi*random.random()**mInv
        minflaws0[itrial] = min(minflaws0[itrial], fi)
        maxflaws0[itrial] = max(maxflaws0[itrial], fi)
    #print "    %5i/%5i: Min/max flaws: [%g : %g]" % (itrial, Ntrials, minflaws0[itrial], maxflaws0[itrial])
print "After creating flaw distributions the hard way:"
print "  Range of minflaws: [%g : %g]" % (min(minflaws0), max(minflaws0))
print "  Range of maxflaws: [%g : %g]" % (min(maxflaws0), max(maxflaws0))

#-------------------------------------------------------------------------------
# Generate distributions of the expected min/max flaws for comparison
#-------------------------------------------------------------------------------
minflaws = [chi * (1.0 - random.random()**(1.0/Nflaws))**mInv for i in xrange(Ntrials)]
maxflaws = [chi * random.random()**(1.0/(m*Nflaws)) for i in xrange(Ntrials)]
print "Using math to estimate min/max flaw ranges:"
print "  Range of minflaws: [%g : %g]" % (min(minflaws), max(minflaws))
print "  Range of maxflaws: [%g : %g]" % (min(maxflaws), max(maxflaws))

#-------------------------------------------------------------------------------
# Compute the ratio of the effective number of flaws to our target number
#-------------------------------------------------------------------------------
Neffratio = [(1.0 + k*V*(epsmax**m - epsmin**m))/Nflaws for epsmin, epsmax in zip(minflaws, maxflaws)]
print "  Range of Neff/Nflaws: [%g : %g]" % (min(Neffratio), max(Neffratio))

#-------------------------------------------------------------------------------
# Now plot the results.
#-------------------------------------------------------------------------------
fig1 = plt.figure(1)
fig2 = plt.figure(2)
fig3 = plt.figure(3)

plt.figure(1)
plt.hist(minflaws0,
         bins = max(10, Ntrials/100),
         color = "black",
         label = "Randomly realized actual distribution")
plt.hist(minflaws,
         bins = max(10, Ntrials/100),
         color = "red",
         alpha = 0.5,
         label = "Predicted distribution")
#plt.axvline(x=minflaw, color="red", linewidth=3)
plt.xlabel(r"Min flaw activation energy ($\varepsilon_{\min}$)")
plt.ylabel(r"N($\varepsilon_{\min}$)")
plt.legend(loc="upper left")

plt.figure(2)
plt.hist(maxflaws0,
         bins = max(10, Ntrials/100),
         color = "black",
         label = "Randomly realized actual distribution")
plt.hist(maxflaws,
         bins = max(10, Ntrials/100),
         color = "red",
         alpha = 0.5,
         label = "Predicted distribution")
#plt.axvline(x=maxflaw, color="red", linewidth=3)
plt.yscale("log")
plt.xlabel(r"Max flaw activation energy ($\varepsilon_{\max}$)")
plt.ylabel(r"N($\varepsilon_{\max}$)")
plt.legend(loc="upper left")

plt.figure(3)
plt.hist(Neffratio,
         bins = max(10, Ntrials/100),
         color = "black",
         label = "Neff/Nflaws")
plt.xlabel(r"$N_{\text{eff}}/N_{\text{flaws}}$")
plt.ylabel(r"$N(N_{\text{eff}}/N_{\text{flaws}})$")
#plt.legend(loc="upper left")
