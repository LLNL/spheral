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
import scipy.integrate, scipy.optimize

#-------------------------------------------------------------------------------
# The P-alpha crush curve from Jutzi 2008
#-------------------------------------------------------------------------------
class PalphaCrushCurve:

    def __init__(self,
                 rhoS0,
                 P0,
                 alpha0,
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

        # Solve for alphae
        self.alphae = alpha0   # Starting point
        last_alphae = 0.0
        iter = 0
        while abs(self.alphae - last_alphae) > 1.0e-6 and iter < 1000:
            iter += 1
            last_alphae = self.alphae
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
    Pe = crushCurve.Pe
    alphae = crushCurve.alphae
    ce = crushCurve.ce(alphae)    # Elastic wave speed

    # The Ranking-Hugoniot conservation relations across a shock front (in the frame of the shock)
    class RKjumpRelations:
        def __init__(self, upiston):
            self.upiston = upiston
        def __call__(self, args):
            u0, u1, rho1, eps1, alpha1 = args

            # Look up the current pressure and distension
            P1 = Pfunc(alpha1*rho1, eps1)/alpha1
            alpha1_new = crushCurve(P1)

            # The exact Hugoniot jump conditions depend on what regime of crushing porosity we're in
            if P1 <= Pe:

                # Elastic regime
                P0 = Pfunc(alpha0*rho0, eps0)/alpha0
                rho00 = rho0

            elif abs(u0) < ce:

                # Plastic regime with an elastic wave moving in front of the shock (shock speed slower than elastic wave speed).
                # In this case the shock is encountering the elastically compressed state rather than virgin state.
                P0 = Pfunc(alpha0*rho0, eps0)/alphae
                rho00 = alpha0*rho0/alphae

            else:

                # The shock speed exceeds the elastic wave speed, so there is no elastic region and we're back to
                # the shock encountering virgin initial state material
                P0 = Pfunc(alpha0*rho0, eps0)/alpha0
                rho00 = rho0

            return np.array([u1 - u0 - self.upiston,
                             rho1*u1 - rho00*u0,
                             P1 + rho1*u1*u1 - P0 - rho00*u0*u0,
                             eps1 + P1/np.maximum(1e-10, rho1) + 0.5*np.square(u1) - eps0 - P0/np.maximum(1.0e-10, rho00) - 0.5*u0*u0,
                             alpha1 - alpha1_new])

    # Prepare to return the arrays of values.  We return these in the same frame as the piston velocity was given, so presumably lab
    # e => elastic region
    # s => post shock region
    US, RHOS, EPSS, PS, ALPHAS = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)
    UE, RHOE, EPSE, PE, ALPHAE = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)
    for i in range(n):
        up = i*dupiston if n > 1 else upiston
        last_guess = (-1.5*up, -0.5*up, 2.0*rho0, 0.5*up*up, alpha0)
        current_solution = scipy.optimize.fsolve(RKjumpRelations(up), last_guess, full_output = True)
        #print("current_solution: ", current_solution)

        # The shock speed vs. the elastic wave speed tells us if we have an elastic compaction region
        # ahead of the shock or not.
        u0, u1, rho1, eps1, alpha1 = current_solution[0]
        print("              up: ", up, u1 - u0)
        assert relativeDiff(u1 - u0, up) < 1.0e-5, current_solution
        US[i] = -u0
        RHOS[i] = rho1
        EPSS[i] = eps1
        PS[i] = Pfunc(alpha1*rho1, eps1)/alpha1
        ALPHAS[i] = alpha1
        if US[i] < ce:
            UE[i] = ce
            RHOE[i] = alpha0*rho0/alphae
            EPSE[i] = eps0
            PE[i] = Pfunc(alpha0*rho0, eps0)/alphae
            ALPHAE[i] = alphae
        else:
            UE[i] = 0.0
            RHOE[i] = rho0
            EPSE[i] = eps0
            PE[i] = Pfunc(alpha0*rho0, eps0)/alpha0
            ALPHAE[i] = alpha0

        # u0_vals[i], u1_vals[i], rho_vals[i], eps_vals[i], alpha_vals[i] = current_solution
        # #last_guess = current_solution
        # #print(up, " --> ", current_solution, " ===> ", RKjumpRelations(up)(current_solution))

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
                 nPoints = 101):
        self.eos = eos
        self.vpiston = abs(vpiston)
        self.eps0 = eps0
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
        self.crushCurve = PalphaCrushCurve(rhoS0, self.P0, alpha0, alphat, Pe, Pt, Ps, cS0, c0, n1, n2)
        return

    # Return the solution profiles as x, v, eps, rho, P, h
    def solution(self, t,
                 x = None):

        # The current piston position
        xpiston = 1.0 - self.vpiston*t

        # Did the user specify the x positions?
        if x is None:
            x = np.linspace(-1.0, xpiston, self.nPoints)

        # Compute the shock jump conditions
        us, rhos, epss, Ps, alphas, ue, rhoe, epse, Pe, alphae = computeHugoniotWithPorosity(self.eos, self.rho0, self.eps0, abs(self.vpiston), self.crushCurve, n=1)
        xs = 1.0 - us*t    # position of the shock
        xe = 1.0 - ue*t    # position of the elastic wave

        # Conditions behind shock
        v1 = -self.vpiston
        h1 = self.h0 * self.rho0/rhos

        # Conditions behind the elastic wave
        v2 = -ue
        h2 = self.h0 * self.rho0/rhoe

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

        return x, v, eps, rho, P, h

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
    print("Hugoniot solution: ", computeHugoniotWithPorosity(eos, rhoS0/alpha0, eps0, abs(vpiston), alpha_curve, n=1))

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
    def plotIt(x, y, label):
        fig = newFigure()
        fig.plot(x, y, "k-")
        fig.set_xlabel(r"$x$")
        fig.set_ylabel(label)
        fig.set_title(label)
    plotIt(x, v, r"$v$")
    plotIt(x, eps, r"$\varepsilon$")
    plotIt(x, rho, r"$\rho$")
    plotIt(x, P, r"$P$")

    
