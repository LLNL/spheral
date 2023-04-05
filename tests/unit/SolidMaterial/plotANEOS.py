from Spheral3d import *
from SpheralMatplotlib import plotSurface
import matplotlib.pyplot as plt
import numpy as np

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
                      10.0*Tmin,         # minimum temperature (K)
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
rho = np.geomspace(rhoMin, rhoMax, num = n)

stuff = []

for eos, label in ((eosSiO2, "SiO2"),
                   # (eosForsterite, "Forsterite"),
                   # (eosWater, "water"),
                   (eosTillotsonBasalt, "Tillotson")
                   ):
    rho0 = eos.referenceDensity
    epsMin = eos.specificThermalEnergy(rho0, Tmin)
    epsMax = eos.specificThermalEnergy(rho0, Tmax)
    epsOff = min(2.0*epsMin, 0.0)
    eps = np.geomspace(epsMin - epsOff, epsMax, num = n) + epsOff
    T = np.geomspace(Tmin, Tmax, num = n)

    rho_grid, eps_grid = np.meshgrid(rho, eps)
    rho_grid, T_grid = np.meshgrid(rho, T)
    shape = rho_grid.shape
    PA_grid, csA_grid = np.zeros(shape), np.zeros(shape)
    sA_grid, gA_grid = np.zeros(shape), np.zeros(shape)
    Teps_grid, epsT_grid, epsTratio_grid = np.zeros(shape), np.zeros(shape), np.zeros(shape)

    for j in range(n):
        for i in range(n):
            rhoi = rho_grid[j][i]
            epsi = rho_grid[j][i]
            PA_grid[j][i] = eos.pressure(rhoi,epsi)
            csA_grid[j][i] = eos.soundSpeed(rhoi,epsi)
            sA_grid[j][i] = eos.entropy(rhoi,epsi)
            gA_grid[j][i] = eos.gamma(rhoi,epsi)
            Teps_grid[j][i] = eos.temperature(rhoi,epsi)

            Ti = T_grid[j][i]
            epsi = eos.specificThermalEnergy(rhoi, Ti)
            Tii = eos.temperature(rhoi, epsi)
            epsT_grid[j][i] = epsi
            epsTratio_grid[j][i] = Tii/Ti

    print("Pressure range for %s    : [%g, %g]" % (label, np.min(PA_grid), np.max(PA_grid)))
    print("Sound speed range for %s : [%g, %g]" % (label, np.min(csA_grid), np.max(csA_grid)))
    print("Entropy range for %s     : [%g, %g]" % (label, np.min(sA_grid), np.max(sA_grid)))
    print("Gamma range for %s       : [%g, %g]" % (label, np.min(gA_grid), np.max(gA_grid)))
    print("eps lookup range for %s  : [%g, %g]" % (label, np.min(epsT_grid), np.max(epsT_grid)))
    print("T lookup range for %s    : [%g, %g]" % (label, np.min(Teps_grid), np.max(Teps_grid)))
    print("T(rho,eps)/T range for %s: [%g, %g]" % (label, np.min(epsTratio_grid), np.max(epsTratio_grid)))

    stuff.append(plotSurface(np.log10(rho_grid), np.log10(eps_grid - epsOff), np.log10(PA_grid),
                             xlabel = "$\log(\\rho)$ (g/cm$^3$)",
                             ylabel = "$\log(\\varepsilon)$ (Mb cm$^2$/g)",
                             zlabel = "$\log(P)$ (dyne)",
                             title = "Pressure %s" % label))

    stuff.append(plotSurface(np.log10(rho_grid), np.log10(eps_grid - epsOff), np.log10(csA_grid),
                             xlabel = "$\log(\\rho)$ (g/cm$^3$)",
                             ylabel = "$\log(\\varepsilon)$ (Mb cm$^2$/g)",
                             zlabel = "$\log(c_s)$ (cm/sec)",
                             title = "Sound speed %s" % label))

    stuff.append(plotSurface(np.log10(rho_grid), np.log10(eps_grid - epsOff), np.log10(sA_grid),
                             xlabel = "$\log(\\rho)$ (g/cm$^3$)",
                             ylabel = "$\log(\\varepsilon)$ (Mb cm$^2$/g)",
                             zlabel = "$\log(s)$",
                             title = "entropy %s" % label))

    stuff.append(plotSurface(np.log10(rho_grid), np.log10(eps_grid - epsOff), np.log10(gA_grid),
                             xlabel = "$\log(\\rho)$ (g/cm$^3$)",
                             ylabel = "$\log(\\varepsilon)$ (Mb cm$^2$/g)",
                             zlabel = "$\log(\\gamma)$",
                             title = "gamma %s" % label))

    stuff.append(plotSurface(np.log10(rho_grid), np.log10(eps_grid - epsOff), np.log10(epsT_grid),
                                             xlabel = "$\log(\\rho)$ (g/cm$^3$)",
                                             ylabel = "$\log(\\varepsilon)$ (Mb cm$^2$/g)",
                                             zlabel = "$\log(\\varepsilon(\\rho, T))$",
                                             title = "eps(rho,T) %s" % label))

    stuff.append(plotSurface(np.log10(rho_grid), np.log10(eps_grid - epsOff), np.log10(Teps_grid),
                             xlabel = "$\log(\\rho)$ (g/cm$^3$)",
                             ylabel = "$\log(\\varepsilon)$ (Mb cm$^2$/g)",
                             zlabel = "$\log(T(\\rho, \\varepsilon))$",
                             title = "log10[T(rho,eps)] %s" % label))

    stuff.append(plotSurface(np.log10(rho_grid), np.log10(eps_grid - epsOff), np.log10(epsTratio_grid),
                             xlabel = "$\log(\\rho)$ (g/cm$^3$)",
                             ylabel = "$\log(\\varepsilon)$ (Mb cm$^2$/g)",
                             zlabel = "$\log(T(\\rho, \\varepsilon)/T)$",
                             title = "T(rho,eps(rho,T))/T %s" % label))
