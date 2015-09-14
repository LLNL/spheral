#-------------------------------------------------------------------------------
# Compare the Osborne and Gruneisen EOS surfaces for Beryllium.
#-------------------------------------------------------------------------------
from math import *
import Gnuplot
from SolidSpheral3d import *
from SpheralGnuPlotUtilities import *

#-------------------------------------------------------------------------------
# Build the EOS's we're going to consider.
#-------------------------------------------------------------------------------
units = PhysicalConstants(0.01,     # Unit length in meters
                          0.001,    # Unit mass in kg
                          1.0e-6)   # Unit time in seconds
etamin, etamax = 0.01, 100.0
eosG = GruneisenEquationOfState(1.85, etamin, etamax, 
                                8.00000000e-01,
                                1.12400000e+00,
                                0.00000000e+00,
                                0.00000000e+00,
                                1.11000000e+00,
                                1.60000000e-01,
                                9.015,
                                units)
eosO = OsborneEquationOfState( 1.85, etamin, etamax,     # Parameters from Howell & Ball 2002
                               0.951168,   # a1
                               0.345301,   # a2pos
                              -0.345301,   # asneg
                               0.926914,   # b0
                               2.948420,   # b1
                               0.507979,   # b2pos
                               0.507979,   # b2neg
                               0.564362,   # c0
                               0.620422,   # c1
                               0.0,        # c2pos
                               0.0,        # c2neg
                               0.8,
                               9.015,
                               units)

#-------------------------------------------------------------------------------
# Plot the pressure and sound speed for each EOS.
#-------------------------------------------------------------------------------
n = 50
rhoMin, rhoMax = 0.05, 50.0
drho = (rhoMax - rhoMin)/n
rho = [rhoMin + i*drho for i in xrange(n + 1)]

epsMin, epsMax = 1.0e-5, 1e10
deps = (log(epsMax) - log(epsMin))/n
eps = [exp(log(epsMin) + i*deps) for i in xrange(n + 1)]

# Make a fake NodeList so we can call the EOS with Fields.
nodes = makeVoidNodeList("nodes", numInternal=1)

for (label, eos) in (("Gruneisen", eosG),
                     ("Osborne", eosO)):
    rhof = ScalarField("rho", nodes)
    epsf = ScalarField("eps", nodes)
    Pf = ScalarField("P", nodes)
    csf = ScalarField("cs", nodes)
    P, cs = [], []
    for rhoi in rho:
        for epsi in eps:
            rhof[0] = rhoi
            epsf[0] = epsi
            eos.setPressure(Pf, rhof, epsf)
            eos.setSoundSpeed(csf, rhof, epsf)
            P.append((rhoi, epsi, Pf[0]))
            cs.append((rhoi, epsi, csf[0]))
    exec("plot%s_P  = generateNewGnuPlot(); Pplot = plot%s_P" % (label, label))
    exec("plot%s_cs = generateNewGnuPlot(); csPlot = plot%s_cs"  % (label, label))
    Pplot("set logscale y; set logscale z")
    csPlot("set logscale y; set logscale z")
    Pdata, csData = Gnuplot.Data(P, inline=True), Gnuplot.Data(cs, inline=True)
    Pplot.splot(Pdata, title="%s pressure" % label)
    Pplot.xlabel("rho")
    Pplot.ylabel("eps")
    csPlot.splot(csData, title="%s sound speed" % label)
    csPlot.xlabel("rho")
    csPlot.ylabel("eps")
