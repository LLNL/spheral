from Spheral3d import *
from SpheralMatplotlib import plotSurface
import matplotlib.pyplot as plt
import numpy as np

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
epsMin, epsMax = 1e-5, 1.1*eosGranite.epsVapor
epsOff = min(2.0*epsMin, 0.0)
rho = np.geomspace(rhoMin, rhoMax, num = n)
eps = np.geomspace(epsMin - epsOff, epsMax, num = n) + epsOff

rho_grid, eps_grid = np.meshgrid(rho, eps)

shape = rho_grid.shape
P_grid, cs_grid = np.zeros(shape), np.zeros(shape)
s_grid, g_grid = np.zeros(shape), np.zeros(shape)
dPdR_grid, dPdU_grid = np.zeros(shape), np.zeros(shape)

# Write the (rho, eps, P, cs) set to a file.
with open("Granite_TillotsonEOS.txt", "w") as f:
    f.write(f"""
# Tillotson EOS dump for a granite like material (all units MKS).
# rho0 = {eosGranite.referenceDensity} kg/m^3
# a = {eosGranite.a} (dimensionless)
# b = {eosGranite.b} (dimensionless)
# A = {eosGranite.A} (Pa)
# B = {eosGranite.B} (Pa)
# alpha = {eosGranite.alpha} (dimensionless)
# beta = {eosGranite.beta} (dimensionless)
# eps0 = {eosGranite.eps0} (J/kg)
# epsLiq = {eosGranite.epsLiquid} (J/kg)
# epsVap = {eosGranite.epsVapor} (J/kg)
#
""")
    f.write((4*"%20s " + "\n") % ("rho (kg/m^3)", "eps (J/kg)", "P (Pa)", "cs (m/sec)"))

    P, cs = [], []
    for j in range(n):
        for i in range(n):
            rhoi = rho_grid[j][i]
            epsi = rho_grid[j][i]
            Pi, dPdUi, dPdRi = eosGranite.pressureAndDerivs(rhoi, epsi)
            P_grid[j][i] = Pi
            cs_grid[j][i] = eosGranite.soundSpeed(rhoi, epsi)
            s_grid[j][i] = eosGranite.entropy(rhoi, epsi)
            g_grid[j][i] = eosGranite.gamma(rhoi, epsi)
            dPdR_grid[j][i] = dPdRi
            dPdU_grid[j][i] = dPdUi
            f.write((4*"%20g " + "\n") % (rhoi, epsi, P_grid[j][i], cs_grid[j][i]))

plotSurface(np.log10(rho_grid), np.log10(eps_grid - epsOff), P_grid,
             xlabel = r"$\log(\rho)$ (kg/m$^3$)",
             ylabel = r"$\log(\varepsilon)$ (J/kg)",
             zlabel = r"$\log(P)$ (Pa)",
             title = "Pressure")

plotSurface(np.log10(rho_grid), np.log10(eps_grid - epsOff), dPdR_grid,
             xlabel = r"$\log(\rho)$ (kg/m$^3$)",
             ylabel = r"$\log(\varepsilon)$ (J/kg)",
             zlabel = r"$\partial P/\partial \rho$ (Pa m$^3$/kg)",
             title = r"$\partial P/\partial \rho$")

plotSurface(np.log10(rho_grid), np.log10(eps_grid - epsOff), dPdU_grid,
             xlabel = r"$\log(\rho)$ (kg/m$^3$)",
             ylabel = r"$\log(\varepsilon)$ (J/kg)",
             zlabel = r"$\partial P/\partial \varepsilon$ (Pa kg/J)",
             title = r"$\partial P/\partial \varepsilon$")

plotSurface(np.log10(rho_grid), np.log10(eps_grid - epsOff), cs_grid,
             xlabel = r"$\log(\rho)$ (kg/m$^3$)",
             ylabel = r"$\log(\varepsilon)$ (J/kg)",
             zlabel = r"$c_s$ (m/sec)",
             title = "Sound speed")

plotSurface(np.log10(rho_grid), np.log10(eps_grid - epsOff), s_grid,
             xlabel = r"$\log(\rho)$ (kg/m$^3$)",
             ylabel = r"$\log(\varepsilon)$ (J/kg)",
             zlabel = r"$s$",
             title = "entropy")

plotSurface(np.log10(rho_grid), np.log10(eps_grid - epsOff), g_grid,
             xlabel = r"$\log(\rho)$ (kg/m$^3$)",
             ylabel = r"$\log(\varepsilon)$ (J/kg)",
             zlabel = r"$\gamma$",
             title = "gamma")
