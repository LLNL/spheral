#-------------------------------------------------------------------------------
# Compute Hugoniot shock jump values for the given equation of state based on an
# initial (density, energy, upiston) state.  upiston can be a single value or
# a numpy array.
# Arguments:
#       eos: the equation of state to profile
#      rho0: initial mass density
#      eps0: initial specific thermal energy
#   upiston: piston velocity (single value or numpy array)
# 
# Returns:
#        us: shock velocity
#      rhos: post-shock density
#      epss: post-shock specific thermal energy
#        Ps: post-shock pressure
# 
# If upiston is a numpy array then the returned values are also arrays corresponding
# to each piston value.
#-------------------------------------------------------------------------------
import mpi
import numpy as np
from numpy import linalg as LA
import scipy.integrate, scipy.optimize, scipy.stats

#-------------------------------------------------------------------------------
# A wrapper for iterating the solution of a system combining scipy minimize and
# fsolve
#-------------------------------------------------------------------------------
def solve(func, initial_guess, bounds,
          acceptance_tol = 1.0e-8,
          maxiter = 100,
          verbose = False):
    solution = None
    if mpi.rank == 0:
        nparams = len(initial_guess)
        iter = 0
        lasterr = 2.0*acceptance_tol
        next_guess = initial_guess
        while iter < maxiter and lasterr > acceptance_tol:
            iter += 1
            opt_guess = scipy.optimize.minimize(func.norm,
                                                x0 = next_guess,
                                                bounds = bounds)
            solution = scipy.optimize.fsolve(func, opt_guess.x, full_output = True)
            lasterr = np.max(np.abs(func(solution[0])))
            if verbose:
                print("--> solve (iteration, err) ", iter, lasterr, func(solution[0]))
            if lasterr > acceptance_tol:
                next_guess = solution[0] + 0.1*scipy.stats.gennorm.rvs(0.1, size=nparams)
                if bounds:
                    for i, bd in enumerate(bounds):
                        if bd[0]:
                            next_guess[i] = max(next_guess[i], bd[0])
                        if bd[1]:
                            next_guess[i] = min(next_guess[i], bd[1])
    solution = mpi.bcast(solution, 0)
    return solution[0]

#-------------------------------------------------------------------------------
# Hugoniot solution
#-------------------------------------------------------------------------------
def computeHugoniot(eos, rho0, eps0, upiston):
    """Compute Hugoniot shock jump values for the given equation of state based on an
initial (density, energy, upiston) state.  upiston can be a single value or
a numpy array.
Arguments:
      eos: the equation of state to profile
     rho0: initial mass density
     eps0: initial specific thermal energy
  upiston: piston velocity (single value or numpy array)

Returns:
       us: shock velocity
     rhos: post-shock density
     epss: post-shock specific thermal energy
       Ps: post-shock pressure

If upiston is a numpy array then the returned values are also arrays corresponding
to each piston value.
    """

    if type(upiston) is float:
        upiston = np.array([upiston])
    assert np.min(upiston) >= 0.0
    n = len(upiston)

    # How should we query the single pressure response from the EOS?  Some
    # Spheral equations of state support this, but others require the Field
    # interface
    if hasattr(eos, "pressure"):
        Pfunc = eos.pressure
    else:
        raise RuntimeError("computeHugoniot does not work yet with EOS's that can't compute a single pressure")

    # Initial pressure
    P0 = Pfunc(rho0, eps0)

    # Couple of useful functions
    def sgn(x):
        return -1.0 if x < 0.0 else 1.0

    def safeInv(x, fuzz = 1e-50):
        return sgn(x)/max(fuzz, abs(x))

    # Functor to help us solve the Rankine-Hugoniot jump relations including a porosity
    class RankineHugoniotJumpRelations:
        def __init__(self, upiston, u0, rho0, eps0, P0):
            self.upiston = upiston
            self.u0 = u0
            self.rho0 = rho0
            self.eps0 = eps0
            self.P0 = P0
            return

        def __call__(self, args):
            us, eps1 = args
            m1 = self.rho0*(us - self.u0)            # mass/time
            rho1 = m1*safeInv(us - self.upiston)
            P1 = Pfunc(rho1, eps1)
            return np.array([m1*(self.upiston - self.u0) - (P1 - self.P0),                                               # Conservation of momentum
                             m1*(eps1 - self.eps0 + 0.5*(self.upiston - self.u0)**2) - P1*(self.upiston - self.u0)])     # Conservation of energy

        def norm(self, args):
            return LA.norm(self(args))
        
    # Prepare to return the arrays of values.  We return these in the same frame as the piston velocity was given, so presumably lab
    # s => post shock region
    US, RHOS, EPSS, PS = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)
    for i, up in enumerate(upiston):
        us, epss = solve(func = RankineHugoniotJumpRelations(up, 0.0, rho0, eps0, P0),
                         initial_guess = (1.5*up,                     # us
                                          0.5*up**2),                 # epss
                         bounds = [(up, np.inf),       # us
                                   (0.0, np.inf)],     # epss
                         verbose = False)
        rhos = rho0*us*safeInv(us - up)
        # print("  Shock conditions:  us = ", us, "\n",
        #       "                   rhos = ", rhos, "\n",
        #       "                   epss = ", epss, "\n",
        #       "                     up = ", up)
        US[i] = us
        RHOS[i] = rhos
        EPSS[i] = epss
        PS[i] = Pfunc(rhos, epss)

    result = US, RHOS, EPSS, PS
    if n == 1:
        return [x[0] for x in result]
    else:
        return result

#-------------------------------------------------------------------------------
# If run directly, provide some example Hugoniot curves
#-------------------------------------------------------------------------------
if __name__ == "__main__":
    from Spheral1d import *
    from matplotlib import pyplot as plt
    from SpheralMatplotlib import *
    def plotIt(x, y, ylabel, xlabel=r"$x$", ylog=False, title=None, style="k-"):
        fig = newFigure()
        if ylog:
            fig.semilogy(x, y, style)
        else:
            fig.plot(x, y, style)
        fig.set_xlabel(xlabel)
        fig.set_ylabel(ylabel)
        if title:
            fig.set_title(title)
        else:
            fig.set_title(ylabel)
        return

    #eos = GammaLawGas(5.0/3.0, 1.0, CGS())
    eos = TillotsonEquationOfState("aluminum melosh89", etamin = 0.1, etamax = 10.0, units=CGuS())
    #eos = ANEOS("dunite", constants=CGuS())
    rho0 = eos.referenceDensity
    eps0 = 0.0
    up = np.linspace(1e-4, 0.1, 50)
    us, rhos, epss, Ps = computeHugoniot(eos, rho0, eps0, up)
    plotIt(rhos, Ps, xlabel=r"$\rho_s$ (g/cc)", ylabel=r"$P_s$ (Mbar)", title=r"$\rho$-$P$ shock Hugoniot", style="ro-")
    plotIt(us, Ps, xlabel=r"$v_s$ (cm/$\mu$sec)", ylabel=r"$P_s$ (Mbar)", title=r"$v_s$-$P$ shock Hugoniot", style="ro-")
    plotIt(up, us, xlabel=r"$v_p$ (cm/$\mu$sec)", ylabel=r"$v_s$ (cm/$\mu$sec)", title=r"$v_p$-$v_s$ shock Hugoniot", style="ro-")
