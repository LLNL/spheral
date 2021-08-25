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
Tmin, Tmax = 1.0e3, 1.0e8
rhoMin, rhoMax = 0.1, 20.0

eosSiO2 =       ANEOS(1,                 # Material number (offset sequentially from ANEOS.INPUT)
                      100,               # num rho vals
                      100,               # num T vals
                      rhoMin,            # minimum density (kg/m^3)
                      rhoMax,            # maximum density (kg/m^3)
                      Tmin,              # minimum temperature (K)
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
n = 50
drho = (rhoMax - rhoMin)/n
rho = [rhoMin + i*drho for i in xrange(n + 1)]

plots = []
for eos, label in ((eosSiO2, "SiO2"),
                   # (eosForsterite, "Forsterite"),
                   # (eosWater, "water"),
                   (eosTillotsonBasalt, "Tillotson")):
    epsMin = eos.specificThermalEnergy(rhoMin, Tmin)
    epsMax = eos.specificThermalEnergy(rhoMax, Tmax)
    deps = (epsMax - epsMin)/n
    eps = [epsMin + i*deps for i in xrange(n + 1)]

    PA, csA, sA, gA = [], [], [], []
    for rhoi in rho:
        for epsi in eps:
            PA.append((rhoi, epsi, eos.pressure(rhoi,epsi)))
            csA.append((rhoi, epsi, eos.soundSpeed(rhoi,epsi)))
            sA.append((rhoi, epsi, eos.entropy(rhoi,epsi)))
            gA.append((rhoi, epsi, eos.gamma(rhoi,epsi)))

    print "Pressure range for %s    : [%g, %g]" % (label, min([x[2] for x in PA]), max([x[2] for x in PA]))
    print "Sound speed range for %s : [%g, %g]" % (label, min([x[2] for x in PA]), max([x[2] for x in csA]))
    print "Entropy range for %s     : [%g, %g]" % (label, min([x[2] for x in PA]), max([x[2] for x in sA]))
    print "Gamma range for %s       : [%g, %g]" % (label, min([x[2] for x in PA]), max([x[2] for x in gA]))

    plots.append(generateNewGnuPlot())
    plots[-1].xlabel("rho (g/cm^3)")
    plots[-1].ylabel("eps (Mb cm^2/g)")
    PAdata = Gnuplot.Data(PA)
    plots[-1].splot(PAdata, title="Pressure %s" % label)

    plots.append(generateNewGnuPlot())
    plots[-1].xlabel("rho (g/cm^3)")
    plots[-1].ylabel("eps (Mb cm^2/g)")
    csAdata = Gnuplot.Data(csA)
    plots[-1].splot(csAdata, title="sound speed %s" % label)
    plots.append(plots[-1])

    plots.append(generateNewGnuPlot())
    plots[-1].xlabel("rho (g/cm^3)")
    plots[-1].ylabel("eps (Mb cm^2/g)")
    sAdata = Gnuplot.Data(sA)
    plots[-1].splot(sAdata, title="entropy %s" % label)
    plots.append(plots[-1])

    plots.append(generateNewGnuPlot())
    plots[-1].xlabel("rho (g/cm^3)")
    plots[-1].ylabel("eps (Mb cm^2/g)")
    gAdata = Gnuplot.Data(gA)
    plots[-1].splot(gAdata, title="gamma %s" % label)
    plots.append(plots[-1])
