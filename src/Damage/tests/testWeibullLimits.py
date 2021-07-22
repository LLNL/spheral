#-------------------------------------------------------------------------------
# Test out our derivation of how the expected limits of a given Weibull
# realization should go.
#-------------------------------------------------------------------------------
from Spheral import *
from SpheralTestUtilities import *
from math import *
import random
import matplotlib.pyplot as plt

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(V = 100.0,          # volume
            k = 1e10,           # k Weibull constant (coefficient)
            m = 9.0,            # m Weibull constant (exponent)
            Nflaws = 10000,     # number of flaws
            Ntrials = 10000,    # number of trials
            seed = 49489212)    # random number seed

#-------------------------------------------------------------------------------
# Generate an example flaw distribution
#-------------------------------------------------------------------------------
random.seed(seed)
minflaw, maxflaw = 2.0, -1.0
for i in xrange(Nflaws):
    fi = (Nflaws/(k*V)*random.random())**(1.0/m)
    minflaw = min(minflaw, fi)
    maxflaw = max(maxflaw, fi)
print "Min/max flaws: [%g : %g]" % (minflaw, maxflaw)

#-------------------------------------------------------------------------------
# Generate distributions of the expected min/max flaws for comparison
#-------------------------------------------------------------------------------
minflaws = [(Nflaws/(k*V))**(1.0/m) * (1.0 - random.random()**(1.0/Nflaws))**(1.0/m) for i in xrange(Ntrials)]
print "Range of minflaws: [%g : %g]" % (min(minflaws), max(minflaws))

maxflaws = [(Nflaws/(k*V))**(1.0/m) * random.random()**(1.0/(m*Nflaws)) for i in xrange(Ntrials)]
print "Range of maxflaws: [%g : %g]" % (min(maxflaws), max(maxflaws))

#-------------------------------------------------------------------------------
# Now plot the results.
#-------------------------------------------------------------------------------
fig1 = plt.figure(1)
fig2 = plt.figure(2)

plt.figure(1)
plt.hist(minflaws,
         bins = max(10, Ntrials/100))
plt.axvline(x=minflaw, color="red", linewidth=3)
plt.xlabel(r"Min flaw activation energy ($\varepsilon_{\min}$)")
plt.ylabel(r"N($\varepsilon_{\min}$)")

plt.figure(2)
plt.hist(maxflaws,
         bins = max(10, Ntrials/100))
plt.axvline(x=maxflaw, color="red", linewidth=3)
plt.xlabel(r"Max flaw activation energy ($\varepsilon_{\max}$)")
plt.ylabel(r"N($\varepsilon_{\max}$)")
