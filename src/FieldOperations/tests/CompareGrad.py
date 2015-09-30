from math import *
from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *

################################################################################
n1 = 100
rho1 = 1.0
m1 = rho1/n1
v1 = Vector1d(1.0)
x0, x1 = 0.0, 1.0

gamma = 2.0
mu = 1.0

intercept = 0.0
multiplier = 1.0

neighborSearchType = 3 # GatherScatter
numGridLevels = 10
topGridCellSize = 2.0
origin = Vector1d(0.0)

################################################################################
# Linear node distribution
def linearSeed(n,
               distort = 0.1,
               xmin = 0.0,
               xmax = 1.0,
               nPerh = 2.01):
    assert n > 0

    x = []
    H = []
    dx0 = (xmax - xmin)/((1.0 + distort)*n)
    for i in xrange(n):
        dx = dx0*(1.0 + distort*(i + 0.5)/n)
        x.append(xmin + (i + 0.5)*dx)
        h = nPerh*dx
        H.append(1.0/h)
    return x, H

################################################################################
# Apply a sinusoidal distortion to a node distribution
def sinSeed(n,
            distort = 0.1,
            xmin = 0.0,
            xmax = 1.0,
            nPerh = 2.01):
    assert n > 0

    x = []
    H = []
    dx0 = (xmax - xmin)/n
    for i in xrange(n):
        dx = dx0*(1.0 + distort*sin(2.0*pi*(i + 0.5)/n))
        x.append(xmin + (i + 0.5)*dx)
        h = nPerh*dx
        H.append(1.0/h)
    return x, H

################################################################################
# Compute the given moment comparison between two distributions.
def LNnorm(Lnum, array0, array1):
    assert len(array0) == len(array1)
    assert len(array0) > 0

    result = 0.0
    for i in xrange(len(array0)):
        result += pow(abs(array1[i] - array0[i]), Lnum)
    result = (result/len(array0))**(1.0/Lnum)
    return result

################################################################################
def linearFunc(x):
    return intercept + multiplier*x

def linearGrad(x):
    return multiplier

################################################################################
def quadraticFunc(x):
    return intercept + multiplier*x**2

def quadraticGrad(x):
    return 2*multiplier*x

################################################################################
def cubicFunc(x):
    return intercept + multiplier*x**3

def cubicGrad(x):
    return 3*multiplier*x**2

################################################################################
def xComponent(list):
    return map(lambda x: x.x, list)

################################################################################
eos = GammaLawGasMKS1d(gamma, mu)

W = BSplineKernel1d()
#W = W4SplineKernel1d()
#W = GaussianKernel1d()
#W = SuperGaussianKernel1d()
#W = PiGaussianKernel1d(1.0)
kernelExtent = W.kernelExtent

# Set the table kernel for the FieldList divergence.
WT = TableKernel1d()
WT.setTableData(W, 100)

import random
generator = random.Random()
ranMag = 0.1

x, H = linearSeed(n1,
                  distort = 0.5,
                  xmin = x0,
                  xmax = x1)

# Create an appropriate set of boundary conditions.
xPlane0, xPlane1 = Plane1d((0.0), (1.0)), Plane1d((1.0), (-1.0))
xbc0 = ReflectingBoundary1d(xPlane0)
xbc1 = ReflectingBoundary1d(xPlane1)

################################################################################
# Set up an SPH node list and data base.
sphNodes = SphNodeList1d(eos, n1)

for i in xrange(n1):
    sphNodes.positions[i] = x[i]
    sphNodes.Hfield[i].xx = H[i]
    #sphNodes.positions[i] = (i + 0.5 + ranMag*generator.uniform(-1.0,1.0))*dx1

sphNodes.mass[:] = [m1]*sphNodes.numNodes
sphNodes.velocity[:] = [v1]*sphNodes.numNodes

sphNeighbor = NestedGridNeighbor1d(sphNodes,
                                   neighborSearchType,
                                   numGridLevels,
                                   topGridCellSize,
                                   origin,
                                   kernelExtent)
sphNodes.neighbor = sphNeighbor

sphNodes.massDensity[:] = [rho1]*sphNodes.numNodes

sphDb = DataBase1d()
sphDb.appendNodeList(sphNodes)

sphDb.updateFluidMassDensity()
sphNodes.updateWeight()

sphPressure = sphDb.fluidPressure

sphMass = sphDb.fluidMass
sphMassDensity = sphDb.fluidMassDensity
sphPosition = sphDb.fluidPosition
sphWeight = sphDb.fluidWeight
sphHfield = sphDb.fluidHfield

xbc0.setGhostNodes(sphDb)
xbc1.setGhostNodes(sphDb)

for bc in [xbc0, xbc1]:
    bc.applyFieldListGhostBoundary(sphMass)
    bc.applyFieldListGhostBoundary(sphMassDensity)
    bc.applyFieldListGhostBoundary(sphWeight)
    bc.applyFieldListGhostBoundary(sphHfield)
    bc.applyFieldListGhostBoundary(sphPressure)

sphDb.updateFluidMassDensity()
sphNodes.updateWeight()

for bc in [xbc0, xbc1]:
    bc.applyFieldListGhostBoundary(sphMassDensity)
    bc.applyFieldListGhostBoundary(sphWeight)

sphNodes.neighbor.updateNodes()

################################################################################
# Set up an MASH node list and data base.
mashNodes = MashNodeList1d(eos, n1)

for i in xrange(n1):
    mashNodes.positions[i] = sphNodes.positions[i]
    mashNodes.Hfield[i].xx = sphNodes.Hfield[i].xx

mashNodes.mass[:] = [m1]*mashNodes.numNodes
mashNodes.weight[:] = [m1]*mashNodes.numNodes
mashNodes.velocity[:] = [v1]*mashNodes.numNodes

mashNeighbor = NestedGridNeighbor1d(mashNodes,
                                   neighborSearchType,
                                   numGridLevels,
                                   topGridCellSize,
                                   origin,
                                   kernelExtent)
mashNodes.neighbor = mashNeighbor

mashNodes.massDensity[:] = [rho1]*mashNodes.numNodes

mashDb = DataBase1d()
mashDb.appendNodeList(mashNodes)

xbc0.setGhostNodes(mashDb)
xbc1.setGhostNodes(mashDb)

mashPressure = mashDb.fluidPressure

mashMass = mashDb.fluidMass
mashMassDensity = mashDb.fluidMassDensity
mashPosition = mashDb.fluidPosition
mashWeight = mashDb.fluidWeight
mashHfield = mashDb.fluidHfield

for bc in [xbc0, xbc1]:
    bc.applyFieldListGhostBoundary(mashMass)
    bc.applyFieldListGhostBoundary(mashMassDensity)
    bc.applyFieldListGhostBoundary(mashWeight)
    bc.applyFieldListGhostBoundary(mashHfield)
    bc.applyFieldListGhostBoundary(mashPressure)

mashNodes.neighbor.updateNodes()

# Set up the linear correction for the MASH nodes.
mashNodes.linearCorrection = 1
bcVector = vector_of_Boundary1dPtr(2)
bcVector[0] = xbc0
bcVector[1] = xbc1
mashNodes.initialize(mashDb, WT, bcVector)

################################################################################
title('Linear field test')

# Set up a linear pressure field, including ghost nodes.
for i in xrange(sphNodes.numNodes):
    sphPressure[0][i] = linearFunc(sphNodes.positions[i].x)
    mashPressure[0][i] = linearFunc(mashNodes.positions[i].x)

# Generate the SPH smoothed estimate of the pressure.
sphSmooth = sphPressure.smoothFields(sphPosition,
                                     sphWeight,
                                     sphMass,
                                     sphMassDensity,
                                     sphHfield,
                                     WT)

# Generate the MASH smoothed estimate of the pressure.
mashSmooth = mashPressure.smoothFieldsMash(mashPosition,
                                           mashMass,
                                           mashHfield,
                                           WT)

mashSmooth2 = mashPressure.smoothFieldsMash2(mashPosition,
                                             mashMass,
                                             mashMassDensity,
                                             mashHfield,
                                             WT)

# Plot the SPH and MASH smoothed fields against the analytic solution.
linearPlot = plotFieldList(sphPressure,
                           lineTitle='Analytic solution',
                           winTitle='Linear field')
plotFieldList(sphSmooth,
              plot=linearPlot,
              plotStyle='points',
              lineTitle='SPH smoothed estimate')
plotFieldList(mashSmooth,
              plot=linearPlot,
              lineTitle='MASH smoothed estimate')
plotFieldList(mashSmooth2,
              plot=linearPlot,
              lineTitle='MASH corrected estimate')

# Calculate the norms of the smoothed vs actual answer.
print 'Interpolate in the Field:'
print 'SPH estimates:  L1 norm = %g, L2 norm = %g' % \
      (LNnorm(1, sphPressure[0][:n1], sphSmooth[0][:n1]),
       LNnorm(2, sphPressure[0][:n1], sphSmooth[0][:n1]))
print 'MASH estimates: L1 norm = %g, L2 norm = %g' % \
      (LNnorm(1, mashPressure[0][:n1], mashSmooth[0][:n1]),
       LNnorm(2, mashPressure[0][:n1], mashSmooth[0][:n1]))
print 'MASH estimates: L1 norm = %g, L2 norm = %g' % \
      (LNnorm(1, mashPressure[0][:n1], mashSmooth2[0][:n1]),
       LNnorm(2, mashPressure[0][:n1], mashSmooth2[0][:n1]))
print '\n'

# Evaluate the SPH and MASH gradients.
sphGrad = sphPressure.gradient(sphPosition, sphWeight, sphMass, sphMassDensity,
                               sphHfield, WT)
mashGrad = mashPressure.gradientMash(mashPosition, mashMass, mashHfield, WT)
mashGrad2 = mashPressure.gradientMash2(mashPosition,
                                       mashMass,
                                       mashMassDensity,
                                       mashHfield,
                                       WT)

# Calculate the analytic gradient.
xans = array([0.0]*n1)
linearGradAnswer = array([0.0]*n1)
for i in xrange(n1):
    xans[i] = sphNodes.positions[i].x
    linearGradAnswer[i] = linearGrad(xans[i])
ansData = Gnuplot.Data(xans, linearGradAnswer, linearGradAnswer, with='lines')

# Plot gradient values.
linearGradPlot = plotFieldList(sphGrad,
                               yFunction='%s.x',
                               plotStyle='points',
                               userXRange=[0, 1],
                               userYRange=[-1, 2],
                               lineTitle='SPH gradient',
                               winTitle='Linear gradient')
plotFieldList(mashGrad,
              yFunction='%s.x',
              plot=linearGradPlot,
              userXRange=[0, 1],
              userYRange=[0.5, 1.5],
              lineTitle='MASH gradient')
plotFieldList(mashGrad2,
              yFunction='%s.x',
              plot=linearGradPlot,
              userXRange=[0, 1],
              userYRange=[0.5, 1.5],
              lineTitle='MASH corrected gradient')
linearGradPlot.replot(ansData)
linearGradPlot.refresh()

# Calculate the norms of the smoothed vs actual answer for the gradient.
print 'Calculate the gradient.'
print 'SPH estimates:  L1 norm = %g, L2 norm = %g' % \
      (LNnorm(1, linearGradAnswer[:n1], xComponent(sphGrad[0][:n1])),
       LNnorm(2, linearGradAnswer[:n1], xComponent(sphGrad[0][:n1])))
print 'MASH estimates: L1 norm = %g, L2 norm = %g' % \
      (LNnorm(1, linearGradAnswer[:n1], xComponent(mashGrad[0][:n1])),
       LNnorm(2, linearGradAnswer[:n1], xComponent(mashGrad[0][:n1])))
print 'MASH estimates: L1 norm = %g, L2 norm = %g' % \
      (LNnorm(1, linearGradAnswer[:n1], xComponent(mashGrad2[0][:n1])),
       LNnorm(2, linearGradAnswer[:n1], xComponent(mashGrad2[0][:n1])))
print '\n'

# ################################################################################
# title('Quadratic field test')

# # Set up a linear pressure field, including ghost nodes.
# for i in xrange(nodes1.numNodes):
#     pressure[0][i] = quadraticFunc(nodes1.positions[i].x)

# # Generate the SPH smoothed estimate of the pressure.
# SPHsmooth = pressure.smoothFields(fluidPosition,
#                                   fluidWeight,
#                                   fluidHfield,
#                                   WT)

# # Generate the MASH smoothed estimate of the pressure.
# MASHsmooth = pressure.smoothFieldsMash(fluidPosition,
#                                             fluidMass,
#                                             fluidHfield,
#                                             WT)

# MASHsmooth2 = pressure.smoothFieldsMash2(fluidPosition,
#                                          fluidMass,
#                                          fluidMassDensity,
#                                          fluidHfield,
#                                          WT)

# # Plot the SPH and MASH smoothed fields against the analytic solution.
# quadraticPlot = plotFieldList(pressure,
#                            lineTitle='Analytic solution',
#                            winTitle='Quadratic field')
# plotFieldList(SPHsmooth,
#               plot=quadraticPlot,
#               plotStyle='points',
#               lineTitle='SPH smoothed estimate')
# plotFieldList(MASHsmooth,
#               plot=quadraticPlot,
#               lineTitle='MASH smoothed estimate')
# plotFieldList(MASHsmooth2,
#               plot=quadraticPlot,
#               lineTitle='MASH corrected estimate')

# # Calculate the norms of the smoothed vs actual answer.
# print 'Interpolate in the Field:'
# print 'SPH estimates:  L1 norm = %g, L2 norm = %g' % \
#       (LNnorm(1, pressure[0][:n1], SPHsmooth[0][:n1]),
#        LNnorm(2, pressure[0][:n1], SPHsmooth[0][:n1]))
# print 'MASH estimates: L1 norm = %g, L2 norm = %g' % \
#       (LNnorm(1, pressure[0][:n1], MASHsmooth[0][:n1]),
#        LNnorm(2, pressure[0][:n1], MASHsmooth[0][:n1]))
# print 'MASH estimates: L1 norm = %g, L2 norm = %g' % \
#       (LNnorm(1, pressure[0][:n1], MASHsmooth2[0][:n1]),
#        LNnorm(2, pressure[0][:n1], MASHsmooth2[0][:n1]))
# print '\n'

# # Evaluate the SPH and MASH gradients.
# SPHgrad = pressure.gradient(fluidPosition, fluidWeight, fluidHfield, WT)
# MASHgrad = pressure.gradientMash(fluidPosition, fluidMass, fluidHfield, WT)

# # Calculate the analytic gradient.
# xans = array([0.0]*n1)
# quadraticGradAnswer = array([0.0]*n1)
# for i in xrange(n1):
#     xans[i] = nodes1.positions[i].x
#     quadraticGradAnswer[i] = quadraticGrad(xans[i])
# ansData = Gnuplot.Data(xans, quadraticGradAnswer, quadraticGradAnswer, with='lines')

# # Plot gradient values.
# quadraticGradPlot = plotFieldList(SPHgrad,
#                                   yFunction='%s.x',
#                                   plotStyle='points',
#                                   userXRange=[0, 1],
#                                   userYRange=[-1, 3],
#                                   lineTitle='SPH gradient',
#                                   winTitle='Quadratic gradient')
                                  
# plotFieldList(MASHgrad,
#               yFunction='%s.x',
#               plot=quadraticGradPlot,
#               userXRange=[0, 1],
#               userYRange=[-1, 3],
#               lineTitle='MASH gradient')
# quadraticGradPlot.replot(ansData)
# quadraticGradPlot.refresh()

# # Calculate the norms of the smoothed vs actual answer for the gradient.
# print 'Calculate the gradient.'
# print 'SPH estimates:  L1 norm = %g, L2 norm = %g' % \
#       (LNnorm(1, quadraticGradAnswer[:n1], xComponent(SPHgrad[0][:n1])),
#        LNnorm(2, quadraticGradAnswer[:n1], xComponent(SPHgrad[0][:n1])))
# print 'MASH estimates: L1 norm = %g, L2 norm = %g' % \
#       (LNnorm(1, quadraticGradAnswer[:n1], xComponent(MASHgrad[0][:n1])),
#        LNnorm(2, quadraticGradAnswer[:n1], xComponent(MASHgrad[0][:n1])))
# print '\n'

# ################################################################################
# title('Cubic field test')

# # Set up a linear pressure field, including ghost nodes.
# for i in xrange(nodes1.numNodes):
#     pressure[0][i] = cubicFunc(nodes1.positions[i].x)

# # Generate the SPH smoothed estimate of the pressure.
# SPHsmooth = pressure.smoothFields(fluidPosition,
#                                   fluidWeight,
#                                   fluidHfield,
#                                   WT)

# # Generate the MASH smoothed estimate of the pressure.
# MASHsmooth = pressure.smoothFieldsMash(fluidPosition,
#                                             fluidMass,
#                                             fluidHfield,
#                                             WT)

# MASHsmooth2 = pressure.smoothFieldsMash2(fluidPosition,
#                                          fluidMass,
#                                          fluidMassDensity,
#                                          fluidHfield,
#                                          WT)

# # Plot the SPH and MASH smoothed fields against the analytic solution.
# cubicPlot = plotFieldList(pressure,
#                            lineTitle='Analytic solution',
#                            winTitle='Cubic field')
# plotFieldList(SPHsmooth,
#               plot=cubicPlot,
#               plotStyle='points',
#               lineTitle='SPH smoothed estimate')
# plotFieldList(MASHsmooth,
#               plot=cubicPlot,
#               lineTitle='MASH smoothed estimate')
# plotFieldList(MASHsmooth2,
#               plot=cubicPlot,
#               lineTitle='MASH corrected estimate')

# # Calculate the norms of the smoothed vs actual answer.
# print 'Interpolate in the Field:'
# print 'SPH estimates:  L1 norm = %g, L2 norm = %g' % \
#       (LNnorm(1, pressure[0][:n1], SPHsmooth[0][:n1]),
#        LNnorm(2, pressure[0][:n1], SPHsmooth[0][:n1]))
# print 'MASH estimates: L1 norm = %g, L2 norm = %g' % \
#       (LNnorm(1, pressure[0][:n1], MASHsmooth[0][:n1]),
#        LNnorm(2, pressure[0][:n1], MASHsmooth[0][:n1]))
# print 'MASH estimates: L1 norm = %g, L2 norm = %g' % \
#       (LNnorm(1, pressure[0][:n1], MASHsmooth2[0][:n1]),
#        LNnorm(2, pressure[0][:n1], MASHsmooth2[0][:n1]))
# print '\n'

# # Evaluate the SPH and MASH gradients.
# SPHgrad = pressure.gradient(fluidPosition, fluidWeight, fluidHfield, WT)
# MASHgrad = pressure.gradientMash(fluidPosition, fluidMass, fluidHfield, WT)

# # Calculate the analytic gradient.
# xans = array([0.0]*n1)
# cubicGradAnswer = array([0.0]*n1)
# for i in xrange(n1):
#     xans[i] = nodes1.positions[i].x
#     cubicGradAnswer[i] = cubicGrad(xans[i])
# ansData = Gnuplot.Data(xans, cubicGradAnswer, cubicGradAnswer, with='lines')

# # Plot gradient values.
# cubicGradPlot = plotFieldList(SPHgrad,
#                                   yFunction='%s.x',
#                                   plotStyle='points',
#                                   userXRange=[0, 1],
#                                   userYRange=[-1, 4],
#                                   lineTitle='SPH gradient',
#                                   winTitle='Cubic gradient')
                                  
# plotFieldList(MASHgrad,
#               yFunction='%s.x',
#               plot=cubicGradPlot,
#               userXRange=[0, 1],
#               userYRange=[-1, 4],
#               lineTitle='MASH gradient')
# cubicGradPlot.replot(ansData)
# cubicGradPlot.refresh()

# # Calculate the norms of the smoothed vs actual answer for the gradient.
# print 'Calculate the gradient.'
# print 'SPH estimates:  L1 norm = %g, L2 norm = %g' % \
#       (LNnorm(1, cubicGradAnswer[:n1], xComponent(SPHgrad[0][:n1])),
#        LNnorm(2, cubicGradAnswer[:n1], xComponent(SPHgrad[0][:n1])))
# print 'MASH estimates: L1 norm = %g, L2 norm = %g' % \
#       (LNnorm(1, cubicGradAnswer[:n1], xComponent(MASHgrad[0][:n1])),
#        LNnorm(2, cubicGradAnswer[:n1], xComponent(MASHgrad[0][:n1])))
# print '\n'

