#-------------------------------------------------------------------------------
# Experimenting with algorithms for seeding the flaw distributions to follow
# the Weibull distribution.
#-------------------------------------------------------------------------------
from math import *
import Gnuplot
import random

# This is all the psyco magic.
import psyco
psyco.full()

gen = random.Random()

#-------------------------------------------------------------------------------
# The algorithm outlined in Benz & Asphaug.
#-------------------------------------------------------------------------------
def benzasphaug(n, k, m, vol):
    result = [[] for i in xrange(n)]
    ienergy = 1
    while (min([len(x) for x in result]) == 0):
        i = gen.randint(0, n - 1)
        result[i].append((ienergy/(k*vol))**(1.0/m))
        ienergy += 1
    return result

#-------------------------------------------------------------------------------
# Use the Weibull variate random algorithm.
#-------------------------------------------------------------------------------
def wei(n, k, m, vol):
    result = [[] for i in xrange(n)]
    beta = (vol/n*k)**(-1.0/m)
    while (min([len(x) for x in result]) == 0):
        i = gen.randint(0, n - 1)
        result[i].append((gen.weibullvariate(beta, m)))
    return result

#-------------------------------------------------------------------------------
# Subdividing algorithm.
#-------------------------------------------------------------------------------
def log2(x):
    return log(x)/log(2.0)

def subdivide(n, k, m, vol):
    result = [[] for i in xrange(n)]
    nenergy = int(log2(n))
    assert nenergy > 0
    while (min([len(x) for x in result]) == 0):
        for ienergy in xrange(nenergy):
            eps = (ienergy/(k*vol))**(1.0/m)
            for j in xrange(ienergy):
                i = gen.randint(0, n - 1)
                result[i].append(eps)
    return result

#-------------------------------------------------------------------------------
# Similar to Benz & Asphaug, but one flaw per node.
#-------------------------------------------------------------------------------
def benzasphaug2(n, k, m, vol):
    result = [[] for i in xrange(n)]
    pool = range(n)
    gen.shuffle(pool)
    ienergy = 1
    for i in pool:
        result[i].append((ienergy/(k*vol))**(1.0/m))
        ienergy += 1
    assert sum([len(x) for x in result]) == n
    return result

#-------------------------------------------------------------------------------
# Generate the cumulative number distribution.
#-------------------------------------------------------------------------------
def cumulativeDistribution(population):
    xresult = []
    for x in population:
        xresult.extend(x)
    assert len(xresult) == sum([len(x) for x in population])
    xresult.sort()
    nresult = [0]*len(xresult)
    nresult[0] = 1
    for i in xrange(1, len(xresult)):
        nresult[i] = nresult[i - 1] + 1
    return xresult, nresult

#-------------------------------------------------------------------------------
# Plot the various algorithms.
#-------------------------------------------------------------------------------
n = 1000
k = 1.12e4 # 27.0
m = 2.63
vol = pi*10.0**2 * 0.3
p = Gnuplot.Gnuplot()

# The Benz & Asphaug method.
xba, yba = cumulativeDistribution(benzasphaug(n, k, m, vol))
dba = Gnuplot.Data(xba, yba,
                   with = "points",
                   title = "Benz \& Asphaug",
                   inline = True)
p.plot(dba)

# Benz & Asphaug with single value per node.
xba2, yba2 = cumulativeDistribution(benzasphaug2(n, k, m, vol))
dba2 = Gnuplot.Data(xba2, yba2,
                    with = "points",
                    title = "Benz \& Asphaug (single flaw per node)",
                    inline = True)
p.replot(dba2)

# The python weibull variate distribution.
xwei, ywei = cumulativeDistribution(wei(n, k, m, vol))
dwei = Gnuplot.Data(xwei, ywei,
                    with = "points",
                    title = "Weibull variate",
                    inline = True)
p.replot(dwei)

# Subdivision
xsub, ysub = cumulativeDistribution(subdivide(n, k, m, vol))
dsub = Gnuplot.Data(xsub, ysub,
                    with = "points",
                    title = "Subdivide",
                    inline = True)
p.replot(dsub)

# The analytic expectation.
xans = xba
yans = [vol*k*x**m for x in xans]
dans = Gnuplot.Data(xans, yans,
                    with = "lines",
                    title = "Analytic expecation",
                    inline = True)
p.replot(dans)

p("set logscale xy")
p.refresh()

