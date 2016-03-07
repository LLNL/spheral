from Spheral1d import *
import Gnuplot
from SpheralTestUtilities import *


rho0 = 100.0
Pmin = 1e-6
Pmax = 1e35
Tmin = 100.0

units = PhysicalConstants(0.01,
                          0.001,
                          1.0e-6)


eos = HelmholtzEquationOfState(units,
                               Pmin,
                               Pmax,
                               Tmin)

#eos = GammaLawGas(4.0/3.0,
#                  13.6,
#                  units)

hmin, hmax, nPerh = 1,1,1

nodes1 = makeFluidNodeList("nodes1", eos,
                           hmin = hmin,
                           hmax = hmax,
                           nPerh = nPerh)

nodes1.numInternalNodes = 1

n = 20
rhoMin, rhoMax = 100.0, 1.0e9
drho = (1.0/n) * log1p(rhoMax/rhoMin)
rho = [rhoMin * exp(drho*i) for i in xrange(n + 1)]

eMin, eMax = 1.0e17, 1.0e20
de = (1.0/n) * log1p(eMax/eMin)
e = [eMin * exp(de*j) for j in xrange(n+1)]



P, cs, gam = [], [], []
for rhoi in rho:
    for ei in e:
        nRho = ScalarField("testRho", nodes1, rhoi)
        ne = ScalarField("testU", nodes1, ei)
        Pr = ScalarField("testP", nodes1)
        ng = ScalarField("testGamma", nodes1)
        soundSpeed = ScalarField("testCs", nodes1)
        eos.setPressure(Pr,nRho,ne)
        eos.setSoundSpeed(soundSpeed,nRho,ne)
        eos.setGammaField(ng,nRho,ne)
        P.append((rhoi, ei, Pr[0]))
        cs.append((rhoi, ei, soundSpeed[0]))
        gam.append((rhoi,ei,ng[0]))
        

Pplot = Gnuplot.Gnuplot()
Pplot("set term x11")
Pplot("set logscale xy")
Pplot.xlabel("rho/rho0")
Pplot.ylabel("eps (J/kg)")
Pdata = Gnuplot.Data(P)
Pplot.splot(Pdata, title="Pressure")

csplot = Gnuplot.Gnuplot()
csplot("set term x11")
csplot("set logscale xy")
csplot.xlabel("rho/rho0")
csplot.ylabel("eps (J/kg)")
csdata = Gnuplot.Data(cs)
csplot.splot(csdata, title="sound speed")

gamplot = Gnuplot.Gnuplot()
gamplot("set term x11")
gamplot("set logscale xy")
gamplot.xlabel("rho/rho0")
gamplot.ylabel("eps (J/kg)")
gamdata = Gnuplot.Data(gam)
gamplot.splot(gamdata, title="gamma")
