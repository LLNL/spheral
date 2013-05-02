from SolidSpheral3d import *
import Gnuplot

# We'll work in MKS units.
units = PhysicalConstants(1.0,   # Unit length in meters
                          1.0,   # Unit mass in kg
                          1.0)   # Unit time in sec

#-------------------------------------------------------------------------------
# Build a Tillotson for Granite like materials.
#-------------------------------------------------------------------------------
mconv = 1.0
lconv = 1.0
tconv = 1.0
MPa = mconv/(lconv*tconv**2) * 1e6  # Convert MPa
GPa = mconv/(lconv*tconv**2) * 1e9  # Convert GPa
MJperKg = (mconv/tconv)**2 * 1e6    # Convert MJ/Kg
GJperKg = (mconv/tconv)**2 * 1e9    # Convert GJ/Kg

rho0Granite = 2.65e3
etaMinGranite = 0.1
etaMaxGranite = 10.0

eosGranite = TillotsonEquationOfState(rho0Granite,     # ref density (kg/m^3)
                                      etaMinGranite,   # etamin             
                                      etaMaxGranite,   # etamax             
                                      0.5,             # a      (dimensionless)
                                      1.5,             # b      (dimensionless)
                                      60.0 * GPa,      # A      (Pa)
                                      40.0 * GPa,      # B      (Pa)
                                      5.0,             # alpha  (dimensionless)
                                      5.0,             # beta   (dimensionless)
                                      10.0  * MJperKg, # eps0   (J/Kg) -- energies for Dolomite
                                      250.0 * MJperKg, # epsLiq (J/Kg) -- energies for Dolomite
                                      1.4   * GJperKg, # epsVap (J/Kg) -- energies for Dolomite
                                      55.350,          # atomic weight -- complete punt
                                      units)

#-------------------------------------------------------------------------------
# Build an ANEOS SiO2 like thing.
#-------------------------------------------------------------------------------
izetl = vector_of_int(1, -1)
initializeANEOS("ANEOS.INPUT", "ANEOS.barf", izetl)
rhoMin, rhoMax = 0.9*etaMinGranite*rho0Granite, 1.1*etaMaxGranite*rho0Granite
eosANEOS = ANEOSEquationOfState(0,                 # Material number
                                1000,              # num rho vals
                                1000,              # num T vals
                                rhoMin,            # minimum density (kg/m^3)
                                rhoMax,            # maximum density (kg/m^3)
                                1.0,               # minimum temperature (K)
                                1e8,               # maximum temperature (K)
                                units)

#-------------------------------------------------------------------------------
# Plot the pressure as a function of (rho, eps)
#-------------------------------------------------------------------------------
n = 50
drho = (rhoMax - rhoMin)/n
rho = [rhoMin + i*drho for i in xrange(n + 1)]

epsMin, epsMax = 0.0, 1.1*eosGranite.epsVapor
deps = (epsMax - epsMin)/n
eps = [epsMin + i*deps for i in xrange(n + 1)]

# Write the (rho, eps, P, cs) set to a file.
f = open("Granite_TillotsonEOS.txt", "w")
f.write("""
# Tillotson EOS dump for a granite like material (all units MKS).
# rho0 = %g kg/m^3
# a = %g (dimensionless)
# b = %g (dimensionless)
# A = %g (Pa)
# B = %g (Pa)
# alpha = %g (dimensionless)
# beta = %g (dimensionless)
# eps0 = %g (J/kg)
# epsLiq = %g (J/kg)
# epsVap = %g (J/kg)
#
""" % (eosGranite.referenceDensity,
       eosGranite.a,
       eosGranite.b,
       eosGranite.A,
       eosGranite.B,
       eosGranite.alpha,
       eosGranite.beta,
       eosGranite.eps0,
       eosGranite.epsLiquid,
       eosGranite.epsVapor))
f.write((6*'"%20s "' + "\n") % ("rho (kg/m^3)", "eps (J/kg)", 
                                "P Till (Pa)", "cs Till (m/sec)", 
                                "P ANEOS (Pa)", "cs ANEOS (m/sec)"))

PT, csT, PA, csA = [], [], [], []
for rhoi in rho:
    for epsi in eps:
        PT.append((rhoi/rho0Granite, epsi/eosGranite.eps0, eosGranite.pressure(rhoi, epsi)))
        csT.append((rhoi/rho0Granite, epsi/eosGranite.eps0, eosGranite.soundSpeed(rhoi, epsi)))
        PA.append((rhoi/rho0Granite, epsi/eosGranite.eps0, eosANEOS.pressure(rhoi, epsi)))
        csA.append((rhoi/rho0Granite, epsi/eosGranite.eps0, eosANEOS.soundSpeed(rhoi, epsi)))
        f.write((6*"%20g " + "\n") % (rhoi, epsi, PT[-1][-1], csT[-1][-1], PA[-1][-1], csA[-1][-1]))
f.close()

PTplot = Gnuplot.Gnuplot()
PTplot.xlabel("rho/rho0")
PTplot.ylabel("eps (J/kg)")
PTdata = Gnuplot.Data(PT)
PTplot.splot(PTdata, title="Pressure (Tillotson)")


csTplot = Gnuplot.Gnuplot()
csTplot.xlabel("rho/rho0")
csTplot.ylabel("eps (J/kg)")
csTdata = Gnuplot.Data(csT)
csTplot.splot(csTdata, title="sound speed (Tillotson)")

PAplot = Gnuplot.Gnuplot()
PAplot.xlabel("rho/rho0")
PAplot.ylabel("eps (J/kg)")
PAdata = Gnuplot.Data(PA)
PAplot.splot(PAdata, title="Pressure (ANEOS)")


csAplot = Gnuplot.Gnuplot()
csAplot.xlabel("rho/rho0")
csAplot.ylabel("eps (J/kg)")
csAdata = Gnuplot.Data(csA)
csAplot.splot(csAdata, title="sound speed (ANEOS)")
