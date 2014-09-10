from Spheral3d import *
import Gnuplot
from SpheralTestUtilities import *


rho0 = 100.0
Pmin = 1e-6
Pmax = 1e35
Tmin = 100.0

units = PhysicalConstants(0.01,
                          0.001,
                          1.0e-6)
nodes = SphNodeList3d(1)

eos = HelmholtzEquationOfState(nodes, # !!!!
                               units,
                               Pmin,
                               Pmax,
                               Tmin)
                               
                               

#-------------------------------------------------------------------------------
# Plot the pressure as a 
#-------------------------------------------------------------------------------
n = 50.0
rhoMin, rhoMax = 100.0, 1.0e9
drho = (1.0/n) * log1p(rhoMax/rhoMin)
rho = [rhoMin * exp(drho*i) for i in xrange(n + 1)]

eMin, eMax = 1.0e10, 1.0e20
de = (1.0/n) * log1p(eMax,eMin)
e = [eMin * exp(de*j) for j in xrange(n+1)]

P, cs = [], []
for rhoi in rho:
    for ei in e:
        nRho = ScalarField3d(nodes, rhoi)
        ne = ScalarField3d(nodes, ei)
        Pr = ScalarField3d(nodes)
        soundSpeed = ScalarField3d(nodes)
        eos.setPressure(Pr,nRho,ne)
        eos.setSoundSpeed(soundSpeed,nRho,ne)
        P.append((rhoi, ei, Pr[0]))
        cs.append((rhoi, ei, soundSpeed[0]))

Pplot = Gnuplot.Gnuplot()
Pplot.xlabel("rho/rho0")
Pplot.ylabel("eps (J/kg)")
Pdata = Gnuplot.Data(P)
Pplot.splot(Pdata, title="Pressure")

csplot = Gnuplot.Gnuplot()
csplot.xlabel("rho/rho0")
csplot.ylabel("eps (J/kg)")
csdata = Gnuplot.Data(cs)
csplot.splot(csdata, title="sound speed")
