#-------------------------------------------------------------------------------
# Plot the fragment mass distribution from three resolutions of the TensileRod
# simulation.
#-------------------------------------------------------------------------------
import pickle
import pylab
import Gnuplot

crapPile = (("dumps-TensileRod-2d-96x24/fragment-properties-t=500.0microsec", "96 x 24"),
            ("dumps-TensileRod-2d-192x48/fragment-properties-t=500.0microsec", "192 x 48"),
            ("dumps-TensileRod-2d-384x96/fragment-properties-t=500.0microsec", "384 x 96"))
            
useGnuplot = True

#-------------------------------------------------------------------------------
# Create the cumulative distribution from a discrete set.
#-------------------------------------------------------------------------------
def cumulativeDistribution(fragMasses):
    total = 0.0
    massDistribution = []
    for m in fragMasses:
        total += m
        massDistribution.append(total)
    assert len(massDistribution) == len(fragMasses)
    return massDistribution

def barwidth(x):
    assert len(x) >= 2
    result = [x[1] - x[0]]
    for i in xrange(1, len(x) - 1):
        result.append(x[i + 1] - x[i])
    result.append(x[-1] - x[-2])
    return result

#-------------------------------------------------------------------------------
# Prepare the plot.
#-------------------------------------------------------------------------------
if useGnuplot:
    p = Gnuplot.Gnuplot()
    p.title("Cumulative fragment mass distribution function")
    p.xlabel('m_{frag}')
    p.ylabel('\Sum m_{frag}')
    p("set logscale x")
    p("set key top left")
else:
    pylab.figure(0)
    pylab.title("Cumulative fragment mass distribution function")
    pylab.xlabel('m_{frag}')
    pylab.ylabel('\Sum m_{frag}')
    pylab.semilogx()

#-------------------------------------------------------------------------------
# Load the data 
#-------------------------------------------------------------------------------
for filename, label in crapPile:

    # Load the fragment properties dictionary.
    f = open(filename, "r")
    fragProps = pickle.load(f)
    f.close()

    # Compute total mass, and normalize the fragment masses by the total.
    mfrags = [fragProps[i]["mass"] for i in fragProps]
    mfrags.sort()
    totalMass = sum(mfrags)
    assert totalMass > 0.0
    for i in xrange(len(mfrags)):
        mfrags[i] /= totalMass
    print "Total mass of fragments: ", totalMass

    # Compute the cumulative mass distribution.
    mdist = cumulativeDistribution(mfrags)

    # Plot the sucker.
    if useGnuplot:
        d = Gnuplot.Data(mfrags, mdist,
                         with = "hist",
                         title = label,
                         inline = True)
        p.replot(d)
    else:
        pylab.plot(mfrags, mdist,
                   label = label)

if not useGnuplot:
    pylab.legend(loc = "best")
    pylab.show()
