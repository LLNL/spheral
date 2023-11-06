#-------------------------------------------------------------------------------
# Compute the analytic solution to the piston driven shock compaction of
# porosity.  Based on the 1D compaction wave problem (5.1) in 
# Jutzi, M., Benz, W., & Michel, P. (2008). Numerical simulations of impacts
# involving porous bodies.I. Implementing sub-resolution porosity in a 3D SPH
# hydrocode. Icarus, 198(1), 242–255.
#
# We are using the solution method described more fully in the section on
# Hugoniots and Isentropes from
# Herrmann, W. (1969). Constitutive Equation for the Dynamic Compaction of 
# Ductile Porous Materials. Journal of Applied Physics, 40(6), 2490–2499.
#-------------------------------------------------------------------------------
import numpy as np
from numpy import linalg as LA
import scipy.integrate, scipy.optimize

#-------------------------------------------------------------------------------
# The P-alpha crush curve from Jutzi 2008
#-------------------------------------------------------------------------------
class PalphaCrushCurve:

    def __init__(self,
                 rhoS0,
                 P0,
                 alpha0,   # alpha(P0)
                 alphat,
                 Pe,
                 Pt,
                 Ps,
                 cS0,
                 c0,
                 n1,
                 n2):
        self.rhoS0 = rhoS0
        self.P0 = P0
        self.alpha0 = alpha0
        self.alphat = alphat
        self.Pe = Pe
        self.Pt = Pt
        self.Ps = Ps
        self.cS0 = cS0
        self.c0 = c0
        self.K0 = cS0*cS0*rhoS0
        self.n1 = n1
        self.n2 = n2

        # Find alphae = alpha(Pe)
        self.alphae = alpha0   # Starting point
        last_alphae = 0.0
        iter = 0
        while abs(self.alphae - last_alphae) > 1.0e-15 and iter < 1000:
            iter += 1
            last_alphae = self.alphae
            print(" --> ", iter, self.alphae, self.c0, self.cS0, self.K0)
            self.alphae = scipy.integrate.solve_ivp(self.Dalpha_elasticDP,
                                                    t_span = [self.P0, self.Pe],
                                                    y0 = [self.alpha0],
                                                    t_eval = [self.Pe]).y[0][0]
        print("alphae: ", self.alphae, last_alphae, iter)

        if self.alphat is None:
            self.alphat = self.alphae    # Reduces to Eq 8 in Jutzi 2008

        assert self.rhoS0 > 0.0
        assert 1.0 <= self.alphat <= self.alphae 
        assert self.Pe <= self.Pt <= self.Ps
        assert self.c0 <= self.cS0
        return

    def h(self, alpha):
        assert self.alphae > 1.0 and self.c0 < self.cS0, "alphae={}, c0={}, cS0={}".format(self.alphae, self.c0, self.cS0)
        return 1.0 + (alpha - 1.0)*(self.c0 - self.cS0)/(self.cS0*(self.alphae - 1.0))

    def Dalpha_elasticDP(self, P, alpha):
        return alpha*alpha/self.K0*(1.0 - 1.0/self.h(alpha)**2)

    def ce(self, alpha):
        return self.h(alpha)*self.cS0

    def alphaElastic(self, P):
        "Just the elastic curve response"
        stuff = scipy.integrate.solve_ivp(self.Dalpha_elasticDP,
                                          t_span = [self.P0, P],
                                          y0 = [self.alpha0],
                                          t_eval = [P])
        return stuff.y[0][0]

    def __call__(self, P):
        if P <= self.P0:
            return self.alpha0

        if (P >= self.Ps):

            # Solid limit, all porosity squeezed out
            return 1.0
        
        elif P < self.Pe:
            assert self.c0 < self.cS0

            # Elastic regime -- integrate elastic Dalpha/DP equation
            stuff = scipy.integrate.solve_ivp(self.Dalpha_elasticDP,
                                              t_span = [self.P0, P],
                                              y0 = [self.alpha0],
                                              t_eval = [P])
            return stuff.y[0][0]

        else:

            # Plastic limit
            if P < self.Pt:
                return (self.alphae - self.alphat) * ((self.Pt - P)/(self.Pt - self.Pe))**self.n1 + (self.alphat - 1.0) * ((self.Ps - P)/(self.Ps - self.Pe))**self.n2 + 1.0
            else:
                return (self.alphat - 1.0) * ((self.Ps - P)/(self.Ps - self.Pe))**self.n2 + 1.0
                
#-------------------------------------------------------------------------------
# Hugoniot solution with porosity
#-------------------------------------------------------------------------------
def computeHugoniotWithPorosity(eos, rho0, eps0, upiston, crushCurve, n = 101):
    assert upiston >= 0.0
    assert n > 0
    dupiston = upiston/max(1, n - 1)

    # How should we query the single pressure response from the EOS?  Some
    # Spheral equations of state support this, but others require the Field
    # interface
    if hasattr(eos, "pressure"):
        Pfunc = eos.pressure
    else:
        raise RuntimeError("computeHugoniotWithPorosity does not work yet with EOS's that can't compute a single pressure")

    def relativeDiff(a, b):
        return abs(a - b)/(abs(a) + abs(b))

    # Extract a few constants from the crush curve
    alpha0 = crushCurve.alpha0
    P0 = Pfunc(alpha0*rho0, eps0)/alpha0
    Pe = crushCurve.Pe
    alphae = crushCurve.alphae
    ce = crushCurve.ce(alpha0)    # Elastic wave speed

    # print("Initial conditions: u0: ", 0.0, "\n",
    #       "                  rho0: ", rho0, "\n",
    #       "                  eps0: ", eps0, "\n",
    #       "                    P0: ", P0, "\n",
    #       "                alpha0: ", alpha0, "\n")

    # Functor to help us solve the Rankine-Hugoniot jump relations including a porosity
    class RankineHugoniotJumpRelations:
        def __init__(self, upiston, u0, rho0, eps0, P0, alpha0, alphaPfunc):
            self.upiston = upiston
            self.u0 = u0
            self.rho0 = rho0
            self.eps0 = eps0
            self.P0 = P0
            self.alpha0 = alpha0
            self.alphaPfunc = alphaPfunc
            return

        def __call__(self, args):
            us, rho1, eps1, alpha1 = args
            P1 = Pfunc(alpha1*rho1, eps1)/alpha1
            alpha1_new = self.alphaPfunc(P1)
            m1 = rho1*(us - self.upiston)       # mass/time
            return np.array([self.rho0*(us - self.u0) - m1,                                                              # Conservation of mass
                             m1*(self.upiston - self.u0) - (P1 - self.P0),                                               # Conservation of momentum
                             m1*(eps1 - self.eps0 + 0.5*(self.upiston**2 - self.u0**2)) - (P1 - self.P0)*self.upiston,   # Conservation of energy
                             alpha1_new - alpha1])                                                                       # Convergence of alpha(P)

        def norm(self, args):
            return LA.norm(self(args))
        
    # Solve for elastic wave properties (wave moves at speed ce).
    # We need to know what is the maximum piston speed that just generates an elastic wave,
    # but does not create a shock. This is denoted upe.
    #   Known: Pe, ce, alphae
    # Unknown: ue, rhoe, epse, upe

    # Solve for upe
    def elasticLimitConditions(upi):

        # Solve for the jump conditions at this piston speed
        RKfunc = RankineHugoniotJumpRelations(upi, 0.0, rho0, eps0, P0, alpha0, crushCurve.alphaElastic)
        wild_guess = (1.5*upi, 2.0*rho0, 0.5*upi**2, alpha0)
        # opt_guess = scipy.optimize.minimize(RKfunc.norm,
        #                                     x0 = wild_guess,
        #                                     bounds = [(0.0, np.inf),   # us (ce)
        #                                               (rho0, np.inf),  # rho1
        #                                               (eps0, np.inf),  # eps1
        #                                               (1.0, alpha0)])  # alpha1
        solution = scipy.optimize.fsolve(RKfunc, wild_guess, xtol = 1.0e-10, full_output = True)
        us, rho1, eps1, alpha1 = solution[0]
        P1 = Pfunc(alpha1*rho1, eps1)/alpha1
        return (us/ce - 1.0)**2 + (P1/Pe - 1.0)**2 + (alpha1/alphae - 1.0)**2

    stuff = scipy.optimize.minimize_scalar(elasticLimitConditions, bracket = (0.0, upiston), tol = 1e-20)
    upe = stuff.x
    # print("Found elastic limit piston velocity upe: ", upe, "\n", stuff)

    # With upe recover the full elastic limit conditions
    RKfunc = RankineHugoniotJumpRelations(upe, 0.0, rho0, eps0, P0, alpha0, crushCurve.alphaElastic)
    wild_guess = (ce, 2.0*rho0, 0.5*upe**2, alphae)
    opt_guess = scipy.optimize.minimize(RKfunc.norm,
                                        x0 = wild_guess,
                                        bounds = [(0.99*ce, 1.01*ce),           # ce
                                                  (rho0, np.inf),     # rhoe
                                                  (eps0, np.inf),     # epse
                                                  (0.99*alphae, 1.01*alphae)])  # alphae
    solution = scipy.optimize.fsolve(RKfunc, opt_guess.x, xtol = 1.0e-10, full_output = True)
    u1, rhoe, epse, alpha1 = solution[0]
    # print("Elastic conditions:  ce = ", ce, u1, abs(u1 - ce)/ce, "\n",
    #       "                   rhoe = ", rhoe, "\n",
    #       "                   epse = ", epse, "\n",
    #       "                     Pe = ", Pe, Pfunc(alphae*rhoe, epse)/alphae, abs(Pfunc(alphae*rhoe, epse)/alphae - Pe)/Pe, "\n",
    #       "                 alphae = ", alphae, alpha1, abs(alpha1 - alphae)/alphae)
    # print("Elastic opt_guess was ", opt_guess)

    # Prepare to return the arrays of values.  We return these in the same frame as the piston velocity was given, so presumably lab
    # e => elastic region
    # s => post shock region
    US, RHOS, EPSS, PS, ALPHAS = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)
    UE, RHOE, EPSE, PE, ALPHAE = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)
    for i in range(n):
        up = i*dupiston if n > 1 else upiston

        if up <= upe:

            # There is no shock, just the elastic wave
            RKfunc = RankineHugoniotJumpRelations(up, 0.0, rho0, eps0, P0, alpha0, crushCurve)
            wild_guess = (1.5*up, 2.0*rho0, 0.5*up**2, 1.0 + 0.5*(alpha0 - 1.0))
            opt_guess = scipy.optimize.minimize(RKfunc.norm,
                                                x0 = wild_guess,
                                                bounds = [(0.0, ce),        # u
                                                          (rho0, rhoe),     # rho
                                                          (eps0, epse),     # eps
                                                          (1.0, alphae)])   # alpha
            solution = scipy.optimize.fsolve(RKfunc, opt_guess.x, xtol = 1.0e-10, full_output = True)
            u1, rho1, eps1, alpha1 = solution[0]
            UE[i] = u1
            RHOE[i] = rho1
            EPSE[i] = eps1
            PE[i] = Pfunc(alpha1*rho1, eps1)/alpha1
            ALPHAE[i] = alpha1

        else:

            # Now we have the full elastic wave properties, so solve for the shock jump conditions encountering the
            # post-elastic wave initial conditions.
            RKfunc = RankineHugoniotJumpRelations(up, upe, rhoe, epse, Pe, alphae, crushCurve)
            wild_guess = (1.5*up, 2.0*rhoe, 0.5*up**2, 1.0)
            opt_guess = scipy.optimize.minimize(RKfunc.norm,
                                                x0 = wild_guess,
                                                bounds = [(0.0, np.inf),      # us
                                                          (rhoe, np.inf),     # rhos
                                                          (epse, np.inf),     # epss
                                                          (1.0, alphae)])     # alphas
            solution = scipy.optimize.fsolve(RKfunc, opt_guess.x, xtol = 1.0e-10, full_output = True)
            us, rhos, epss, alphas = solution[0]
            # print("  Shock conditions:  us = ", us, "\n",
            #       "                   rhos = ", rhos, "\n",
            #       "                   epss = ", epss, "\n",
            #       "                 alphas = ", alphas)

    # If the shock speed exceeds ce, then the shock front overtakes the elastic wave and there is no
    # elastic region
    if us >= ce:
        RKfunc = RankineHugoniotJumpRelations(upiston, 0.0, rho0, eps0, P0, alpha0, crushCurve)
        wild_guess = (1.5*upiston, 2.0*rhoe, 0.5*upiston**2, 1.0)
        opt_guess = scipy.optimize.minimize(RKfunc.norm,
                                            x0 = wild_guess,
                                            bounds = [(0.0, np.inf),      # us
                                                      (rho0, np.inf),     # rhos
                                                      (eps0, np.inf),     # epss
                                                      (1.0, alpha0)])     # alphas
        solution = scipy.optimize.fsolve(RKfunc, opt_guess.x, xtol = 1.0e-10, full_output = True)
        us, rhos, epss, alphas = solution[0]
        # print("  Shock conditions:  us = ", us, "\n",
        #       "                   rhos = ", rhos, "\n",
        #       "                   epss = ", epss, "\n",
        #       "                 alphas = ", alphas)
        
        US[i] = us
        RHOS[i] = rhos
        EPSS[i] = epss
        PS[i] = Pfunc(alphas*rhos, epss)/alphas
        ALPHAS[i] = alphas

    else:

        # There is both an elastic wave and shocked region
        US[i] = us
        RHOS[i] = rhos
        EPSS[i] = epss
        PS[i] = Pfunc(alphas*rhos, epss)/alphas
        ALPHAS[i] = alphas
        UE[i] = upe
        RHOE[i] = rhoe
        EPSE[i] = epse
        PE[i] = Pe
        ALPHAE[i] = alphae

    result = US, RHOS, EPSS, PS, ALPHAS, UE, RHOE, EPSE, PE, ALPHAE
    if len(result[0] == 1):
        return [x[0] for x in result]
    else:
        return result
 
#-------------------------------------------------------------------------------
# Compute the solution profiles for the test
#-------------------------------------------------------------------------------
class PlanarCompactionSolution:

    def __init__(self,
                 eos,
                 vpiston,
                 eps0,
                 alpha0,
                 alphat,
                 Pe,
                 Pt,
                 Ps,
                 n1,
                 n2,
                 cS0 = None,
                 c0 = None,
                 h0 = None,
                 nPoints = 101,
                 pistonFrame = False):
        self.eos = eos
        self.vpiston = abs(vpiston)
        self.eps0 = eps0
        self.alpha0 = alpha0
        rhoS0 = eos.referenceDensity
        self.rho0 = rhoS0/alpha0
        PS0 = eos.pressure(rhoS0, eps0)
        self.P0 = PS0/alpha0
        if cS0 is None:
            cS0 = eos.soundSpeed(rhosS0, eps0)
        if c0 is None:
            c0 = cS0
        if h0 is None:
            self.h0 = 2.0/nPoints
        else:
            self.h0 = h0
        self.nPoints = nPoints
        self.pistonFrame = pistonFrame
        self.crushCurve = PalphaCrushCurve(rhoS0, self.P0, alpha0, alphat, Pe, Pt, Ps, cS0, c0, n1, n2)
        return

    # Find the shock and elastic wave properties as a function of time
    def waveProperties(self, t):

        # Compute the shock jump conditions
        us, rhos, epss, Ps, alphas, ue, rhoe, epse, Pe, alphae = computeHugoniotWithPorosity(self.eos, self.rho0, self.eps0, abs(self.vpiston), self.crushCurve, n=1)
        ce = self.crushCurve.ce(self.alpha0)
        xs = 1.0 - us*t    # position of the shock
        xe = 1.0 - ce*t    # position of the elastic wave

        # Conditions behind shock
        v1 = -self.vpiston
        h1 = self.h0 * self.rho0/rhos

        # Conditions behind the elastic wave
        v2 = -ue
        h2 = self.h0 * self.rho0/rhoe
        
        return us, rhos, epss, Ps, alphas, ue, rhoe, epse, Pe, alphae, xs, xe, v1, h1, v2, h2

    # Return the solution profiles as x, v, eps, rho, P, h
    def solution(self, t,
                 x = None):

        # Get the current regional properties
        us, rhos, epss, Ps, alphas, ue, rhoe, epse, Pe, alphae, xs, xe, v1, h1, v2, h2 = self.waveProperties(t)

        # Did the user specify the x positions?
        xpiston = 1.0 - self.vpiston*t
        if x is None:
            x = np.linspace(-1.0, xpiston, self.nPoints)

        # Generate the solution profiles
        n = len(x)
        v, eps, rho, P, h = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)
        for i in range(n):
            if x[i] >= xs:
                v[i], eps[i], rho[i], P[i], h[i] = v1, epss, rhos, Ps, h1
            elif x[i] >= xe:
                v[i], eps[i], rho[i], P[i], h[i] = v2, epse, rhoe, Pe, h2
            else:
                v[i], eps[i], rho[i], P[i], h[i] = 0.0, self.eps0, self.rho0, self.P0, self.h0

        # In the piston frame we actually look like the fluid is flowing into a stationary wall
        if self.pistonFrame:
            x += self.vpiston*t
            v += self.vpiston

        return x, v, eps, rho, P, h

    # Return the alpha solution
    def alpha_solution(self, t,
                       x = None):
        
        # Get the current regional properties
        us, rhos, epss, Ps, alphas, ue, rhoe, epse, Pe, alphae, xs, xe, v1, h1, v2, h2 = self.waveProperties(t)

        # Did the user specify the x positions?
        if x is None:
            xpiston = 1.0 - self.vpiston*t
            x = np.linspace(-1.0, xpiston, self.nPoints)

        # Generate the solution profiles
        n = len(x)
        alpha = np.zeros(n)
        for i in range(n):
            if x[i] >= xs:
                alpha[i] = alphas
            elif x[i] >= xe:
                alpha[i] = alphae
            else:
                alpha[i] = self.alpha0

        # In the piston frame we actually look like the fluid is flowing into a stationary wall
        if self.pistonFrame:
            x += self.vpiston*t

        return x, alpha

    # Return the sound speed solution
    def soundSpeed_solution(self, t,
                            x = None):
        
        # Get the current regional properties
        us, rhos, epss, Ps, alphas, ue, rhoe, epse, Pe, alphae, xs, xe, v1, h1, v2, h2 = self.waveProperties(t)
        def _soundSpeed(rhoi, epsi, alphai):
            cS0 = self.eos.soundSpeed(alphai*rhoi, epsi)
            return cS0 + (alphai - 1.0)/(self.alpha0 - 1.0)*(self.crushCurve.c0 - cS0)
        c_s = _soundSpeed(rhos, epss, alphas)
        c_e = _soundSpeed(rhoe, epse, alphae)
        c_0 = self.crushCurve.c0

        # Did the user specify the x positions?
        if x is None:
            xpiston = 1.0 - self.vpiston*t
            x = np.linspace(-1.0, xpiston, self.nPoints)

        # Generate the solution profiles
        n = len(x)
        cs = np.zeros(n)
        for i in range(n):
            if x[i] >= xs:
                cs[i] = c_s
            elif x[i] >= xe:
                cs[i] = c_e
            else:
                cs[i] = c_0

        # In the piston frame we actually look like the fluid is flowing into a stationary wall
        if self.pistonFrame:
            x += self.vpiston*t

        return x, cs

#-------------------------------------------------------------------------------
# Main if run directly
#-------------------------------------------------------------------------------
if __name__ == "__main__":
    from Spheral1d import *
    units = CGuS()
    cgs = CGS()
    mCGSconv = cgs.unitMassKg/units.unitMassKg
    lCGSconv = cgs.unitLengthMeters/units.unitLengthMeters
    tCGSconv = cgs.unitTimeSec/units.unitTimeSec
    PCGSconv = mCGSconv/(tCGSconv*tCGSconv*lCGSconv)
    vCGSconv = lCGSconv/tCGSconv

    eos = TillotsonEquationOfState("aluminum melosh89", etamin = 0.1, etamax = 10.0, units=units)
    #eos = ANEOS("dunite", constants=CGuS())

    alpha0 = 1.275
    rhoS0 = eos.referenceDensity
    eps0 = 0.0
    PS0 = eos.pressure(rhoS0, eps0)
    Pe = 8e8 * PCGSconv
    Ps = 7e9 * PCGSconv
    cS0 = 5.35e5 * vCGSconv
    c0 = 4.11e5 *vCGSconv
    n1 = 0.0
    n2 = 2.0
    alphat = None # (alpha0 - 1.0)*((Ps - Pe)/(Ps - PS0))**2 + 1.0

    alpha_curve = PalphaCrushCurve(rhoS0, PS0/alpha0, alpha0, alphat, Pe, Pe, Ps, cS0, c0, n1, n2)

    P = np.linspace(PS0, 1.5*Ps, 1000)
    alpha = np.array([alpha_curve(x) for x in P])
    from matplotlib import pyplot as plt
    from SpheralMatplotlib import *
    fig = newFigure()
    fig.plot(P, alpha, "k-")
    fig.set_xlabel(r"$P$")
    fig.set_ylabel(r"$\alpha$")
    fig.set_title(r"$\alpha(P)$ crush curve")

    vpiston = -45.8e-3
    #print("Hugoniot solution: ", computeHugoniotWithPorosity(eos, rhoS0/alpha0, eps0, abs(vpiston), alpha_curve, n=1))

    solution = PlanarCompactionSolution(eos,
                                        vpiston,
                                        eps0,
                                        alpha0,
                                        alphat,
                                        Pe,
                                        Pe,
                                        Ps,
                                        n1,
                                        n2,
                                        cS0,
                                        c0)
    x, v, eps, rho, P, h = solution.solution(t=3.5)
    xx, alpha = solution.alpha_solution(t=3.5)
    def plotIt(x, y, label, ylog=False):
        fig = newFigure()
        if ylog:
            fig.semilogy(x, y, "k-")
        else:
            fig.plot(x, y, "k-")
        fig.set_xlabel(r"$x$")
        fig.set_ylabel(label)
        fig.set_title(label)
    plotIt(x, v, r"$v$")
    plotIt(x, eps, r"$\varepsilon$")
    plotIt(x, rho, r"$\rho$")
    plotIt(x, P, r"$P$", True)
    plotIt(x, alpha, r"$\alpha$")

    
