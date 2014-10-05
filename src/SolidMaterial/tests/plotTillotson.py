from SolidSpheral3d import *
import Gnuplot

#-------------------------------------------------------------------------------
# Granite (solid) material parameters.
# Extracted from information in
# (a) Saito, Kaiho, Abe, Katayama, & Takayama, Int. J. of Imp. Eng., 35,
#     1770-1777, 2008
# (b) Ivanov, Badukov, Yakovlev, Gerasimov, Dikov, Pope, & Ocampo, Geological
#     Society of America, Special paper 307, The Cretacious-Tertiary Event and
#     other catastrophies in Earth History, 1996
# (c) www.matweb.com
# (d) http://www.efunda.com/materials/common_matl/common_matl.cfm?matlphase=solid&matlprop=mechanical
#
# Look up:
#   Young's modulus           : 20-60 GPa  : using 20.0 GPa
#   Poison's ratio            : 0.2-0.3    : using 0.25
#   Ultimate tensile strength : 7-25 MPa   : using 7.0
#
# Derive:
#   Shear modulus   : 8.0 GPa
#
# Moduli related by E = 2*G*(1 + nu) = 3*K*(1 - 2*nu)
#-------------------------------------------------------------------------------
mconv = 1.0
lconv = 1.0
tconv = 1.0
MPa = mconv/(lconv*tconv**2) * 1e6  # Convert MPa
GPa = mconv/(lconv*tconv**2) * 1e9  # Convert GPa
MJperKg = (mconv/tconv)**2 * 1e6    # Convert MJ/Kg
GJperKg = (mconv/tconv)**2 * 1e9    # Convert GJ/Kg

rho0Granite = 2.63 * 1e3
etaMinGranite = 0.5
etaMaxGranite = 1.5

units = PhysicalConstants(1.0,  # unit length in meters
                          1.0,  # unit mass in kg
                          1.0)  # unit time in sec
eosGranite = TillotsonEquationOfState(rho0Granite,     # ref density (kg/m^3)
                                      etaMinGranite,   # etamin             
                                      etaMaxGranite,   # etamax             
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
# Plot the pressure as a 
#-------------------------------------------------------------------------------
n = 50
rhoMin, rhoMax = 0.9*etaMinGranite*rho0Granite, 1.1*etaMaxGranite*rho0Granite
#rhoMin, rhoMax = rho0Granite, 1.1*etaMaxGranite*rho0Granite
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
f.write((4*"%20s " + "\n") % ("rho (kg/m^3)", "eps (J/kg)", "P (Pa)", "cs (m/sec)"))

P, cs = [], []
for rhoi in rho:
    for epsi in eps:
        P.append((rhoi/rho0Granite, epsi/eosGranite.eps0, eosGranite.pressure(rhoi, epsi)))
        cs.append((rhoi/rho0Granite, epsi/eosGranite.eps0, eosGranite.soundSpeed(rhoi, epsi)))
        f.write((4*"%20g " + "\n") % (rhoi, epsi, P[-1][-1], cs[-1][-1]))
f.close()

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
