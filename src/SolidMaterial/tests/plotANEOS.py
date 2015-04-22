from SolidSpheral3d import *
import Gnuplot

# We'll work in CGS units.
units = PhysicalConstants(0.01,  # Unit length in meters
                          0.001, # Unit mass in kg
                          1.0)   # Unit time in sec

#-------------------------------------------------------------------------------
# Build a Gruneisen for SiO2 like materials.
#-------------------------------------------------------------------------------
mconv = 1.0
lconv = 1.0
tconv = 1.0

rho0 = 2.65               # g/cc
C0 = 0.36839e6            # cm/s
S1 = 1.8954               # dimensionless
S2 = 0.0                  # dimensionless
S3 = 0.0                  # dimensionless
gamma0 = 0.9              # dimensionless
b = 1.0                   # dimensionless
etaMin = 0.2
etaMax = 5.0

eosSiO2 = GruneisenEquationOfState(rho0,            # ref density (g/cc)
                                   etaMin,          # etamin             
                                   etaMax,          # etamax
                                   C0,
                                   S1,
                                   S2,
                                   S3,
                                   gamma0,
                                   b,
                                   60.0843,         # atomic weight
                                   units)

#-------------------------------------------------------------------------------
# Build an ANEOS SiO2 like thing.
#-------------------------------------------------------------------------------
izetl = vector_of_int(1, -1)
initializeANEOS("ANEOS.INPUT", "ANEOS.barf", izetl)
rhoMin, rhoMax = 0.9*etaMin*rho0, 1.1*etaMax*rho0
Tmin, Tmax = 1.0, 1.0e8
eosANEOS = ANEOS(0,                 # Material number
                 1000,              # num rho vals
                 1000,              # num T vals
                 rhoMin,            # minimum density (kg/m^3)
                 rhoMax,            # maximum density (kg/m^3)
                 Tmin,              # minimum temperature (K)
                 Tmax,              # maximum temperature (K)
                 units)
eps0ANEOS = eosANEOS.specificThermalEnergy(rho0, 1.0)  # Specific energy at 1K, reference density
print "eps0ANEOS = ", eps0ANEOS

#-------------------------------------------------------------------------------
# Plot the pressure as a function of (rho, eps)
#-------------------------------------------------------------------------------
n = 50
drho = (rhoMax - rhoMin)/n
rho = [rhoMin + i*drho for i in xrange(n + 1)]

epsMin = eosANEOS.specificThermalEnergy(rho0, 0.1*Tmin)
epsMax = eosANEOS.specificThermalEnergy(rho0, 1.1*Tmax)
deps = (epsMax - epsMin)/n
eps = [epsMin + i*deps for i in xrange(n + 1)]

# Write the (rho, eps, P, cs) set to a file.
f = open("SiOS_ANEOS.txt", "w")
f.write("""
# ANEOS vs. Gruneisen EOS dump for SiO2 like material (all units CGS).
#
# ANEOS eps(1K) = %g
#
""" % (eps0ANEOS))
f.write((6*'"%20s "' + "\n") % ("rho (g/cm^3)", "eps (erg/g)", 
                                "P Grun (dyne)", "cs Grun (cm/sec)", 
                                "P ANEOS (dyne)", "cs ANEOS (cm/sec)"))

PG, csG, PA, csA = [], [], [], []
for rhoi in rho:
    for epsi in eps:
        PG.append((rhoi, epsi, eosSiO2.pressure(rhoi, epsi - epsMin)))
        csG.append((rhoi, epsi, eosSiO2.soundSpeed(rhoi, epsi - epsMin)))
        PA.append((rhoi, epsi, eosANEOS.pressure(rhoi, epsi)))
        csA.append((rhoi, epsi, eosANEOS.soundSpeed(rhoi, epsi)))
        f.write((6*"%20g " + "\n") % (rhoi, epsi, PG[-1][-1], csG[-1][-1], PA[-1][-1], csA[-1][-1]))
f.close()

PGplot = Gnuplot.Gnuplot()
PGplot.xlabel("rho (g/cm^3)")
PGplot.ylabel("eps (erg/g)")
PGdata = Gnuplot.Data(PG)
PGplot.splot(PGdata, title="Pressure (Gruneisen)")


csGplot = Gnuplot.Gnuplot()
csGplot.xlabel("rho (g/cm^3)")
csGplot.ylabel("eps (erg/g)")
csGdata = Gnuplot.Data(csG)
csGplot.splot(csGdata, title="sound speed (Gruneisen)")

PAplot = Gnuplot.Gnuplot()
PAplot.xlabel("rho (g/cm^3)")
PAplot.ylabel("eps (erg/g)")
PAdata = Gnuplot.Data(PA)
PAplot.splot(PAdata, title="Pressure (ANEOS)")


csAplot = Gnuplot.Gnuplot()
csAplot.xlabel("rho (g/cm^3)")
csAplot.ylabel("eps (erg/g)")
csAdata = Gnuplot.Data(csA)
csAplot.splot(csAdata, title="sound speed (ANEOS)")
