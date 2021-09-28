from SolidSpheral3d import *
from SpheralGnuPlotUtilities import *

# We'll work in CGuS units.
units = PhysicalConstants(0.01,     # Unit length in meters
                          0.001,    # Unit mass in kg
                          1.0e-6)   # Unit time in sec

#-------------------------------------------------------------------------------
# Build an ANEOS SiO2 like thing.
#-------------------------------------------------------------------------------
izetl = [-1] #, -2, -3]
initializeANEOS("ANEOS.SIO2", "ANEOS.barf", izetl)
#initializeANEOS("ANEOS.INPUT", "ANEOS.barf", izetl)
etaMin, etaMax = 0.2, 5.0
Tmin, Tmax = 1.0e0, 1.0e8
rhoMin, rhoMax = 0.1, 20.0

eosSiO2 =       ANEOS(1,                 # Material number (offset sequentially from ANEOS.INPUT)
                      200,               # num rho vals
                      400,               # num T vals
                      rhoMin,            # minimum density (kg/m^3)
                      rhoMax,            # maximum density (kg/m^3)
                      10.0*Tmin,              # minimum temperature (K)
                      Tmax,              # maximum temperature (K)
                      units)
# eosForsterite = ANEOS(2,                 # Material number (offset sequentially from ANEOS.INPUT)
#                       100,               # num rho vals
#                       100,               # num T vals
#                       rhoMin,            # minimum density (kg/m^3)
#                       rhoMax,            # maximum density (kg/m^3)
#                       Tmin,              # minimum temperature (K)
#                       Tmax,              # maximum temperature (K)
#                       units)
# eosWater =      ANEOS(3,                 # Material number (offset sequentially from ANEOS.INPUT)
#                       100,               # num rho vals
#                       100,               # num T vals
#                       rhoMin,            # minimum density (kg/m^3)
#                       rhoMax,            # maximum density (kg/m^3)
#                       Tmin,              # minimum temperature (K)
#                       Tmax,              # maximum temperature (K)
#                       units)

eosTillotsonBasalt = TillotsonEquationOfState(materialName = "basalt",
                                              etamin = etaMin,
                                              etamax = etaMax,
                                              units = units)

#-------------------------------------------------------------------------------
# Plot the pressure, entropy, & sound speed as a function of (rho, eps)
#-------------------------------------------------------------------------------
n = 100
drho = (rhoMax - rhoMin)/n
rho = [rhoMin + i*drho for i in xrange(n + 1)]

plots = []
gdata = []
def plotIt(data, xlabel, ylabel, title):
    gdata.append(Gnuplot.Data(data))
    plots.append(generateNewGnuPlot())
    plots[-1].xlabel(xlabel)
    plots[-1].ylabel(ylabel)
    plots[-1].splot(gdata[-1], title=title)

for eos, label in ((eosSiO2, "SiO2"),
                   # (eosForsterite, "Forsterite"),
                   # (eosWater, "water"),
                   (eosTillotsonBasalt, "Tillotson")):
    epsMin = eos.specificThermalEnergy(rhoMin, Tmin)
    epsMax = eos.specificThermalEnergy(rhoMax, Tmax)
    deps = (epsMax - epsMin)/n
    dT = (Tmax - Tmin)/n
    eps = [epsMin + i*deps for i in xrange(n + 1)]
    T = [Tmin + i*dT for i in xrange(n + 1)]

    PA, csA, sA, gA, epsT, Teps, epsTratio = [], [], [], [], [], [], []
    for rhoi in rho:
        for epsi in eps:
            PA.append((rhoi, epsi, eos.pressure(rhoi,epsi)))
            csA.append((rhoi, epsi, eos.soundSpeed(rhoi,epsi)))
            sA.append((rhoi, epsi, eos.entropy(rhoi,epsi)))
            gA.append((rhoi, epsi, eos.gamma(rhoi,epsi)))
            Teps.append((rhoi, epsi, log10(eos.temperature(rhoi,epsi))))
        for Ti in T:
            epsi = eos.specificThermalEnergy(rhoi,Ti)
            Tii = eos.temperature(rhoi,epsi)
            epsT.append((rhoi, Ti, epsi))
            epsTratio.append((rhoi, Ti, Tii/Ti))

    print "Pressure range for %s    : [%g, %g]" % (label, min([x[2] for x in PA]), max([x[2] for x in PA]))
    print "Sound speed range for %s : [%g, %g]" % (label, min([x[2] for x in csA]), max([x[2] for x in csA]))
    print "Entropy range for %s     : [%g, %g]" % (label, min([x[2] for x in sA]), max([x[2] for x in sA]))
    print "Gamma range for %s       : [%g, %g]" % (label, min([x[2] for x in gA]), max([x[2] for x in gA]))
    print "eps lookup range for %s  : [%g, %g]" % (label, min([x[2] for x in epsT]), max([x[2] for x in epsT]))
    print "T lookup range for %s    : [%g, %g]" % (label, 10.0**min([x[2] for x in Teps]), 10.0**max([x[2] for x in Teps]))
    print "T(rho,eps)/T range for %s: [%g, %g]" % (label, min([x[2] for x in epsTratio]), max([x[2] for x in epsTratio]))

    plotIt(PA,
           "rho (g/cm^3)",
           "eps (Mb cm^2/g)",
           "Pressure %s" % label)

    plotIt(csA,
           "rho (g/cm^3)",
           "eps (Mb cm^2/g)",
           "sound speed %s" % label)

    plotIt(sA,
           "rho (g/cm^3)",
           "eps (Mb cm^2/g)",
           "entropy %s" % label)

    plotIt(gA,
           "rho (g/cm^3)",
           "eps (Mb cm^2/g)",
           "gamma %s" % label)

    plotIt(epsT,
           "rho (g/cm^3)",
           "T",
           "eps(rho,T) %s" % label)

    plotIt(Teps,
           "rho (g/cm^3)",
           "eps (Mb cm^2/gm)",
           "log10[T(rho,eps)] %s" % label)

    plotIt(epsTratio,
           "rho (g/cm^3)",
           "T",
           "T(rho,eps(rho,T))/T %s" % label)
