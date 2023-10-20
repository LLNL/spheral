#-------------------------------------------------------------------------------
# Compute a Hugoniot curve for the given equation of state from a starting
# (density, energy) point.
# Arguments:
#     eos: the equation of state to profile
#    rho0: initial mass density
#    eps0: initial specific thermal energy
#  epsmax: maximum specific thermal energy to probe
# The result is returned as a tuple of arrays:
#-------------------------------------------------------------------------------
import numpy as np
from scipy.optimize import fsolve

def computeHugoniot(eos, rho0, eps0, upiston, n = 101):
    assert upiston >= 0.0
    assert n > 1
    dupiston = upiston/(n - 1)

    # How should we query the single pressure response from the EOS?  Some
    # Spheral equations of state support this, but others require the Field
    # interface
    if hasattr(eos, "pressure"):
        Pfunc = eos.pressure
    else:
        raise RuntimeError("computeHugoniot does not work yet with EOS's that can't compute a single pressure")

    # Starting point
    P0 = Pfunc(rho0, eps0)

    # The Ranking-Hugoniot conservation relations across a shock front (in the frame of the shock)
    class RKjumpRelations:
        def __init__(self, upiston):
            self.upiston = upiston
        def __call__(self, args):
            u0, u1, rho1, eps1 = args
            P1 = Pfunc(rho1, eps1)
            return np.array([np.abs(u1 - u0) - self.upiston,
                             rho1*u1 - rho0*u0,
                             P1 + rho1*u1*u1 - P0 - rho0*u0*u0,
                             eps1 + P1/np.maximum(1e-10, rho1) + 0.5*np.square(u1) - eps0 - P0/np.maximum(1.0e-10, rho0) - 0.5*u0*u0])

    u0_vals, u1_vals, rho_vals, eps_vals = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)
    for i in range(n):
        up = i*dupiston
        last_guess = (-1.5*up, -0.5*up, 4.0*rho0, 0.5*up*up)
        current_solution = fsolve(RKjumpRelations(up), last_guess)
        u0_vals[i], u1_vals[i], rho_vals[i], eps_vals[i] = current_solution
        #last_guess = current_solution
        print(up, " --> ", current_solution, " ===> ", RKjumpRelations(up)(current_solution))

    return u0_vals, u1_vals, rho_vals, eps_vals, np.array([Pfunc(rhoi, epsi) for rhoi, epsi in zip(rho_vals, eps_vals)])

if __name__ == "__main__":
    from Spheral1d import *
    #eos = GammaLawGas(5.0/3.0, 1.0, CGS())
    #eos = TillotsonEquationOfState("aluminum melosh89", etamin = 0.1, etamax = 10.0, units=CGuS())
    eos = ANEOS("dunite", constants=CGuS())
    u0, u1, rho, eps, P = computeHugoniot(eos, 0.5*eos.referenceDensity, 0.0, 4.0)
    from matplotlib import pyplot as plt
    plt.plot(rho, P)
    
