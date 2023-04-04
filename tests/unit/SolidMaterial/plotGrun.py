from math import *
from Spheral1d import *
from SpheralTestUtilities import *
import matplotlib.pyplot as plt
from SpheralMatplotlib import plotSurface
import numpy as np

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
izetl = [-1]
initializeANEOS("ANEOS.INPUT", "ANEOS.barf", izetl)
rhoMin, rhoMax = 0.9*etaMin*rho0, 1.1*etaMax*rho0
Tmin, Tmax = 1.0, 1.0e8
eosANEOS = ANEOS(1,                 # Material number
                 100,               # num rho vals
                 100,               # num T vals
                 rhoMin,            # minimum density (kg/m^3)
                 rhoMax,            # maximum density (kg/m^3)
                 Tmin,              # minimum temperature (K)
                 Tmax,              # maximum temperature (K)
                 units)
eps0ANEOS = eosANEOS.specificThermalEnergy(rho0, 1.0)  # Specific energy at 1K, reference density
print("eps0ANEOS = ", eps0ANEOS)

#-------------------------------------------------------------------------------
# Plot the pressure as a function of (rho, eps)
#-------------------------------------------------------------------------------
n = 50
rho = np.geomspace(rhoMin, rhoMax, num = n)

epsMin = eosANEOS.specificThermalEnergy(rho0, 0.1*Tmin)
epsMax = eosANEOS.specificThermalEnergy(rho0, 1.1*Tmax)
eps = np.geomspace(epsMin, epsMax, num = n)

rho_grid, eps_grid = np.meshgrid(rho, eps)
shape = rho_grid.shape
PG_grid, csG_grid = np.zeros(shape), np.zeros(shape)
PA_grid, csA_grid = np.zeros(shape), np.zeros(shape)

# Write the (rho, eps, P, cs) set to a file.
with open("SiOS_ANEOS.txt", "w") as f:
    f.write("""
    # ANEOS vs. Gruneisen EOS dump for SiO2 like material (all units CGS).
    #
    # ANEOS eps(1K) = %g
    #
    """ % (eps0ANEOS))
    f.write((6*'"%20s "' + "\n") % ("rho (g/cm^3)", "eps (erg/g)", 
                                    "P Grun (dyne)", "cs Grun (cm/sec)", 
                                    "P ANEOS (dyne)", "cs ANEOS (cm/sec)"))

    for j in range(n):
        for i in range(n):
            PG_grid[j][i] = eosSiO2.pressure(rho_grid[j][i], eps_grid[j][i] - epsMin)
            csG_grid[j][i] = eosSiO2.soundSpeed(rho_grid[j][i], eps_grid[j][i])
            PA_grid[j][i] = eosANEOS.pressure(rho_grid[j][i], eps_grid[j][i])
            csA_grid[j][i] = eosANEOS.soundSpeed(rho_grid[j][i], eps_grid[j][i])
            f.write((6*"%20g " + "\n") % (rho_grid[j][i], eps_grid[j][i],
                                          PG_grid[j][i], csG_grid[j][i],
                                          PA_grid[j][i], csA_grid[j][i]))

PGplot, PGax, PGsurf = plotSurface(rho_grid, eps_grid, PG_grid,
                                   xlabel = "$\\rho$ (g/cm$^3$)",
                                   ylabel = "$\\varepsilon$ (erg/g)",
                                   zlabel = "$P$ (dynes)",
                                   title = "Pressure (Gruneisen)")

csGplot, csGax, csGsurf = plotSurface(rho_grid, eps_grid, csG_grid,
                                      xlabel = "$\\rho$ (g/cm$^3$)",
                                      ylabel = "$\\varepsilon$ (erg/g)",
                                      zlabel = "$c_s$ (cm/sec)",
                                      title = "Sound speed (Gruneisen)")

PAplot, PAax, PAsurf = plotSurface(rho_grid, eps_grid, PA_grid,
                                   xlabel = "$\\rho$ (g/cm$^3$)",
                                   ylabel = "$\\varepsilon$ (erg/g)",
                                   zlabel = "$P$ (dynes)",
                                   title = "Pressure (ANEOS)")

csAplot, csAax, csAsurf = plotSurface(rho_grid, eps_grid, csA_grid,
                                      xlabel = "$\\rho$ (g/cm$^3$)",
                                      ylabel = "$\\varepsilon$ (erg/g)",
                                      zlabel = "$c_s$ (cm/sec)",
                                      title = "Sound speed (ANEOS)")

plt.show()
