#-------------------------------------------------------------------------------
# Plot out pressure and sound speed surfaces for various materials in our
# little Tillotson library of values.
#-------------------------------------------------------------------------------
from math import *
import Gnuplot
from SolidSpheral3d import *

#-------------------------------------------------------------------------------
# Build the EOS's we're going to consider.
#-------------------------------------------------------------------------------
mats = ["Granite", "Pumice", "Nylon", "Glass", "copper"]
units = PhysicalConstants(0.01,     # Unit length in meters
                          0.001,    # Unit mass in kg
                          1.0e-6)   # Unit time in seconds
etamin, etamax = 0.01, 100.0
EOSes = [TillotsonEquationOfState(mat, etamin, etamax, units) for mat in mats]

#-------------------------------------------------------------------------------
# Plot the pressure and sound speed for each EOS.
#-------------------------------------------------------------------------------
n = 50
rhoMin, rhoMax = 0.2, 50.0
#rhoMin, rhoMax = 200.0, 5e4
drho = (rhoMax - rhoMin)/n
rho = [rhoMin + i*drho for i in xrange(n + 1)]

#epsMin, epsMax = 1.0, 1e15
epsMin, epsMax = 1.0e-2, 1e2
deps = (log(epsMax) - log(epsMin))/n
eps = [exp(log(epsMin) + i*deps) for i in xrange(n + 1)]

# Make a fake NodeList so we can call the EOS with Fields.
nodes = makeVoidNodeList("nodes", numInternal=1)

plots = []
for matLabel, eos in zip(mats, EOSes):
    rhof = ScalarField("rho", nodes)
    epsf = ScalarField("eps", nodes)
    Pf = ScalarField("P", nodes)
    csf = ScalarField("cs", nodes)
    gamf = ScalarField("gamma", nodes)
    P, cs, gam = [], [], []
    for rhoi in rho:
        for epsi in eps:
            rhof[0] = rhoi
            epsf[0] = epsi
            eos.setPressure(Pf, rhof, epsf)
            eos.setSoundSpeed(csf, rhof, epsf)
            eos.setGammaField(gamf, rhof, epsf)
            P.append((rhoi, epsi, Pf[0]))
            cs.append((rhoi, epsi, csf[0]))
            gam.append((rhoi, epsi, gamf[0]))
    plots.append(Gnuplot.Gnuplot())
    plots.append(Gnuplot.Gnuplot())
    plots.append(Gnuplot.Gnuplot())
    Pplot = plots[-3]
    csPlot = plots[-2]
    gamPlot = plots[-1]
    Pplot("set logscale y; set logscale z")
    csPlot("set logscale y; set logscale z")
    gamPlot("set logscale y; set logscale z")
    Pdata, csData, gamData = Gnuplot.Data(P, inline=True), Gnuplot.Data(cs, inline=True), Gnuplot.Data(gam, inline=True)
    Pplot.splot(Pdata, title="%s pressure" % matLabel)
    Pplot.xlabel("rho")
    Pplot.ylabel("eps")
    csPlot.splot(csData, title="%s sound speed" % matLabel)
    csPlot.xlabel("rho")
    csPlot.ylabel("eps")
    gamPlot.splot(gamData, title="%s gam" % matLabel)
    gamPlot.xlabel("rho")
    gamPlot.ylabel("eps")
