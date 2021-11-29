#!/usr/bin/env python
#-------------------------------------------------------------------------------
# RiemannSolution
#
# Adapted from code I got from Chad Noble, which in turn was based on code from
# Toro as described in the following comments.
#
# Exact Riemann solver for the Euler equations in one dimension
# Translated from the Fortran code er1pex.f and er1pex.ini
# by Dr. E.F. Toro downloaded from
# http://www.numeritek.com/numerica_software.html#freesample
#-------------------------------------------------------------------------------
from math import *
import numpy as np
import argparse

# Several standard tests listed as (x0, x1, xdiaph, gamma_gas, out_time, dl, vl, pl, dr, vr, pr)
# The bounds are chosen so that a good comparison can be made in the range x \in [0,1]
Riemann_packaged_problems = {
    "sod"                 : ( 0.0, 1.0, 0.5, 1.4, 0.20, 1.0, 0.0, 1.0, 0.125, 0.0, 0.1),                         # TEST 1 (Modified Sod)
    "123"                 : ( 0.0, 1.0, 0.5, 1.4, 0.15, 1.0, -2.0, 0.4, 1.0, 2.0, 0.4),                          # TEST 2 (123 problem)
    "leftwc"              : ( 0.0, 1.0, 0.5, 1.4, 0.012, 1.0, 0.0, 1000.0, 1.0, 0.0, 0.01),                      # TEST 3 (Left Woodward & Colella)
    "2shock_collision"    : (-1.0, 2.0, 0.4, 1.4, 0.035, 5.99924, 19.5975, 460.894, 5.99242, -6.19633, 46.0950), # TEST 4 (Collision of 2 shocks)
    "stationary_contact"  : (-0.5, 1.5, 0.8, 1.4, 0.012, 1.0, -19.59745, 1000.0, 1.0, -19.59745, 0.01),          # TEST 5 (Stationary contact)
    "slow_shock"          : (-2.0, 8.0, 0.5, 1.4, 1.0, 3.857143, -0.810631, 10.33333, 1.0, -3.44, 1.0),          # TEST 6 (Slow shock)
    "shock_contact_shock" : (-1.0, 2.0, 0.5, 1.4, 0.3, 1.0, 0.5, 1.0, 1.25, -0.5, 1.0),                          # TEST 7 (Shock-Contact-Shock)
    "leblanc"             : ( 0.0, 1.0, 0.3, 1.4, 0.5, 1.0, 0.0, 2.0e-1/3.0, 0.01, 0.0, 2.0e-10/3.0),            # TEST 8 (LeBlanc)
}

#-------------------------------------------------------------------------------
# The main object.
#-------------------------------------------------------------------------------
class RiemannSolution:

    def __init__(self, 
                 problem = "Sod",  # ("", "Sod", "123", "Stationary_contact", "Slow_shock", "Slow_contact_shock", "LeBlanc")
                 n = 1000,         # number of points in evaluating exact solution
                 x0 = None,        # box min coordinate
                 x1 = None,        # box max coordinate
                 xdiaph = None,    # position of diaphragm xdiaph \in [x0, x1]
                 gamma_gas = None, # ratio of specific heats
                 out_time = None,  # default time of solution
                 dl = None,        # density (left state)
                 vl = None,        # velocity (left state)
                 pl = None,        # pressure (left state)
                 hl = None,        # smoothing scale (left state)
                 dr = None,        # density (right state)
                 vr = None,        # velocity (right state)
                 pr = None,        # pressure (right state)
                 hr = None):       # smoothing scale (right state)

        assert problem or (x0 and x1 and out_time and xdiaph and gamma_gas and dl and vl and pl and dr and vr and pr)

        # Get the ICs.
        if problem:
            assert problem.lower() in Riemann_packaged_problems
            _x0, _x1, _xdiaph, _gamma_gas, _out_time, _dl, _vl, _pl, _dr, _vr, _pr = Riemann_packaged_problems[problem.lower()]
            if x0 is None:
                x0 = _x0
            if x1 is None:
                x1 = _x1
            if xdiaph is None:
                xdiaph = _xdiaph
            if gamma_gas is None:
                gamma_gas = _gamma_gas
            if out_time is None:
                out_time = _out_time
            if dl is None:
                dl = _dl
            if vl is None:
                vl = _vl
            if pl is None:
                pl = _pl
            if dr is None:
                dr = _dr
            if vr is None:
                vr = _vr
            if pr is None:
                pr = _pr

        # Store the variables
        self.n = n
        self.x0 = x0
        self.x1 = x1
        self.xdiaph = xdiaph
        self.gamma_gas = gamma_gas
        self.out_time = out_time
        self.dl = dl
        self.vl = vl
        self.pl = pl
        self.hl = hl
        self.dr = dr
        self.vr = vr
        self.pr = pr
        self.hr = hr
        return

    #---------------------------------------------------------------------------
    # Compute the solution.
    #---------------------------------------------------------------------------
    def solution(self, 
                 time = None,
                 x = None):

        n = self.n
        x0 = self.x0
        x1 = self.x1
        xdiaph = self.xdiaph
        gamma_gas = self.gamma_gas
        out_time = self.out_time
        dl = self.dl
        vl = self.vl
        pl = self.pl
        hl = self.hl
        dr = self.dr
        vr = self.vr
        pr = self.pr
        hr = self.hr

        # Solution time
        if not time is None:
            out_time = time
        else:
            out_time = self.out_time

        # Sampling positions
        if x is None:
            assert n > 0
            assert x1 > x0
            x = np.linspace(x0, x1, n)
        else:
            n = len(x)

        # Did we get the initial (left, right) h?
        if hl is None:
            hl = x[1] - x[0]
        if hr is None:
            hr = x[-1] - x[-2]
        assert hl > 0 and hr > 0
        
        # compute gamma related constants
        g1 = (gamma_gas - 1.0)/(2.0*gamma_gas)
        g2 = (gamma_gas + 1.0)/(2.0*gamma_gas)
        g3 = 2.0*gamma_gas/(gamma_gas - 1.0)
        g4 = 2.0/(gamma_gas - 1.0)
        g5 = 2.0/(gamma_gas + 1.0)
        g6 = (gamma_gas - 1.0)/(gamma_gas + 1.0)
        g7 = (gamma_gas - 1.0)/2.0
        g8 = gamma_gas - 1.0
        
        # compute sound speeds
        cl = sqrt(gamma_gas*pl/dl)
        cr = sqrt(gamma_gas*pr/dr)

        #---------------------------------------------------------------------------
        # purpose: to provide a guessed value for pressure
        #          pm in the Star Region. The choice is made
        #          according to adaptive Riemann solver using
        #          the PVRS, TRRS and TSRS approximate
        #          Riemann solvers. See Sect. 9.5 of Chapt. 9 of Ref. 1
        #---------------------------------------------------------------------------
        def guessp():
            quser = 2.0
        
            # compute guess pressure from PVRS Riemann solver
            cup = 0.25*(dl + dr)*(cl + cr)
            ppv = 0.5*(pl + pr) + 0.5*(vl - vr)*cup
            ppv = max(0.0, ppv)
            pmin = min(pl,  pr)
            pmax = max(pl,  pr)
            qmax = pmax/pmin
        
            if (qmax <= quser and (pmin <= ppv and ppv <= pmax)):
                pm = ppv     # select PVRS Riemann solver
            else:
                if (ppv < pmin):
                    # select Two-Rarefaction Riemann solver
                    pq = pow(pl/pr, g1)
                    vm = (pq*vl/cl + vr/cr + g4*(pq - 1.0))/(pq/cl + 1.0/cr)
                    ptl = 1.0 + g7*(vl - vm)/cl
                    ptr = 1.0 + g7*(vm - vr)/cr
                    pm = 0.5*(pow(pl*ptl, g3) + pow(pr*ptr, g3))
                else:
                    # select Two-Shock Riemann solver with PVRS as estimate
                    gel = sqrt((g5/dl)/(g6*pl + ppv))
                    ger = sqrt((g5/dr)/(g6*pr + ppv))
                    pm = (gel*pl + ger*pr - (vr - vl))/(gel + ger)
            return pm
        
        #---------------------------------------------------------------------------
        # purpose: to evaluate the pressure functions
        #          fl and fr in exact Riemann solver
        #          and their first derivatives
        #---------------------------------------------------------------------------
        def prefun(p, dk, pk, ck):
            if (p <= pk):
                # rarefaction wave
                pratio = p/pk
                f = g4*ck*(pow(pratio, g1) - 1.0)
                fd = (1.0/(dk*ck))*pow(pratio, -g2)
            else:
                #  shock wave
                ak = g5/dk
                bk = g6*pk
                qrt = sqrt(ak/(bk + p))
                f = (p - pk)*qrt
                fd = (1.0 - 0.5*(p - pk)/(bk + p))*qrt
            return f, fd
        
        #---------------------------------------------------------------------------
        # purpose: to compute the solution for pressure and
        #          velocity in the Star Region
        #---------------------------------------------------------------------------
        def starpu(pscale):
        
            nriter = 20
            tolpre = 1.0e-6
        
            # guessed value pstart is computed
            pstart = guessp()
            pold = pstart
            udiff = vr - vl
        
            print ("----------------------------------------\n"
                   "   Iteration number     Change\n"
                   "----------------------------------------")
        
            i = 1
            change = 10.0*tolpre
            while i <= nriter and change > tolpre:
                fl, fld = prefun(pold, dl, pl, cl)
                fr, frd = prefun(pold, dr, pr, cr)
                p = pold - (fl + fr + udiff)/(fld + frd)
                change = 2.0*abs((p - pold)/(p + pold))
                print '\t', i, "\t\t", change
                if (p < 0.0):
                    p = tolpre
                pold = p
                i += 1
            if (i > nriter):
                print "divergence in Newton-Raphson iteration"
        
            # compute velocity in star region
            u = 0.5*(vl + vr + fr - fl)
            print "----------------------------------------\n" \
                  "     Pressure           Velocity\n"         \
                  "----------------------------------------\n" \
                  "     ", p/pscale, "\t\t", u, '\n'           \
                  "----------------------------------------"
            return p, u
        
        #---------------------------------------------------------------------------
        # purpose: to sample the solution throughout the wave
        #          pattern. Pressure pm and velocity vm in the
        #          star region are known. Sampling is performed
        #          in terms of the 'speed' s = x/t. Sampled
        #          values are d, v, p
        #---------------------------------------------------------------------------
        def sample(pm, vm, s):
        
            if (s <= vm):
                # sampling point lies to the left of the contact discontinuity
                if (pm <= pl):
                    # left rarefaction
                    shl = vl - cl
                    if (s <= shl):
                        # sampled point is left data state
                        d = dl
                        v = vl
                        p = pl
                        h = hl
                    else:
                        cml = cl*pow(pm/pl, g1)
                        stl = vm - cml
                        if (s > stl):
                            # sampled point is star left state
                            d = dl*pow(pm/pl, 1.0/gamma_gas)
                            v = vm
                            p = pm
                            h = hl*dl/d
                        else:
                            # sampled point is inside left fan
                            v = g5*(cl + g7*vl + s)
                            c = g5*(cl + g7*(vl - s))
                            d = dl*pow(c/cl, g4)
                            p = pl*pow(c/cl, g3)
                            h = hl*dl/d
                else:
                    # left shock
                    pml = pm/pl
                    sl = vl - cl*sqrt(g2*pml + g1)
                    if (s <= sl):
                        # sampled point is left data state
                        d = dl
                        v = vl
                        p = pl
                        h = hl
                    else:
                        # sampled point is star left state
                        d = dl*(pml + g6)/(pml*g6 + 1.0)
                        v = vm
                        p = pm
                        h = hl*dl/d
            else:
                # sampling point lies to the right of the contact discontinuity
                if (pm > pr):
                    # right shock
                    pmr = pm/pr
                    sr  = vr + cr*sqrt(g2*pmr + g1)
                    if (s >= sr):
                        # sampled point is right data state
                        d = dr
                        v = vr
                        p = pr
                        h = hr
                    else:
                        # sampled point is star right state
                        d = dr*(pmr + g6)/(pmr*g6 + 1.0)
                        v = vm
                        p = pm
                        h = hr*dr/d
                else:
                    # right rarefaction
                    shr = vr + cr
                    if (s >= shr):
                        # sampled point is right data state
                        d = dr
                        v = vr
                        p = pr
                        h = hr
                    else:
                        cmr = cr*pow(pm/pr, g1)
                        str = vm + cmr
                        if (s <= str):
                            # sampled point is star right state
                            d = dr*pow(pm/pr, 1.0/gamma_gas)
                            v = vm
                            p = pm
                            h = hr*dr/d
                        else:
                            # sampled point is inside left fan
                            v = g5*(-cr + g7*vr + s)
                            c = g5*(cr - g7*(vr - s))
                            d = dr*pow(c/cr, g4)
                            p = pr*pow(c/cr, g3)
                            h = hr*dr/d
            return d, v, p, h
        
        
        # the pressure positivity condition is tested for
        if (g4*(cl+cr) <= (vr-vl)):
            raise RunTimeError, ("the initial data is such that vacuum is generated"
                                 "\nstopping program")
        
        # exact solution for pressure and velocity in star region is found
        pm, vm = starpu(1.0)
        
        # complete solution at time out_time is found
        d = np.empty(n)
        v = np.empty(n)
        p = np.empty(n)
        eps = np.empty(n)
        A = np.empty(n)
        h = np.empty(n)
        for i in xrange(n):
            s = (x[i] - xdiaph)/max(1e-10, out_time)
            ds, vs, ps, hs = sample(pm, vm, s)
            d[i] = ds
            v[i] = vs
            p[i] = ps
            eps[i] = ps/(g8*ds)
            A[i] = ps/pow(ds, gamma_gas)
            h[i] = hs
        
        return x, v, eps, d, p, A, h

#-------------------------------------------------------------------------------
# Provide a way to call this script as a standalone executable.
#-------------------------------------------------------------------------------
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description = "Compute the Riemann solution, with optional output to a file or plotted to the screen.")
    ap.add_argument("--problem", 
                    default = "Sod",
                    help = """
Use one of the canned Riemann initial conditions: (Sod, 123, Stationary_contact, Slow_shock, Slow_contact_shock, LeBlanc).
If specified as the empty string "" (or None), the full state must be specified explicitly.""")
    ap.add_argument("--n",
                    default = 1000,
                    type = int,
                    help = "Number of points to generate in the solution.")
    ap.add_argument("--x0", 
                    default = None,
                    type = float,
                    help = "Minimum spatial coordinate in the tube.")
    ap.add_argument("--x1", 
                    default = None,
                    type = float,
                    help = "Maximum spatial coordinate in the tube.")
    ap.add_argument("--xdiaph", 
                    default = None,
                    type = float,
                    help = "Coordinate of the diaphragm.")
    ap.add_argument("--gamma_gas",
                    default = None,
                    type = float,
                    help = "Ratio of specific heats.")
    ap.add_argument("--out_time", 
                    default = None,
                    type = float,
                    help = "Solution time.")
    ap.add_argument("--dl",
                    default = None,
                    type = float,
                    help = "Initial density for left state.")
    ap.add_argument("--vl",
                    default = None,
                    type = float,
                    help = "Initial velocity for left state.")
    ap.add_argument("--pl",
                    default = None,
                    type = float,
                    help = "Initial pressure for left state.")
    ap.add_argument("--hl",
                    default = None,
                    type = float,
                    help = "Initial smoothing scale for left state.")
    ap.add_argument("--dr",
                    default = None,
                    type = float,
                    help = "Initial density for right state.")
    ap.add_argument("--vr",
                    default = None,
                    type = float,
                    help = "Initial velocity for right state.")
    ap.add_argument("--pr",
                    default = None,
                    type = float,
                    help = "Initial pressure for right state.")
    ap.add_argument("--hr",
                    default = None,
                    type = float,
                    help = "Initial smoothing scale for right state.")
    ap.add_argument("--file",
                    default = None,
                    help = "Write profiles to given file.")
    ap.add_argument("--noheader",
                    action = "store_true",
                    help = "Do not write a header at the top of the output file.")
    ap.add_argument("--plot",
                    action = "store_true",
                    help = "Plot the profiles to the screen.")
    ap.add_argument("--plotsize",
                    default = 10,
                    type = float,
                    help = "Set the size of the figure (in inches) when plotting.")
    args = ap.parse_args()
    globals().update(vars(args))

    # Compute the solution.
    answer = RiemannSolution(problem = problem,
                             n = n,
                             x0 = x0,
                             x1 = x1,
                             xdiaph = xdiaph,
                             gamma_gas = gamma_gas,
                             out_time = out_time,
                             dl = dl,
                             vl = vl,
                             pl = pl,
                             hl = hl,
                             dr = dr,
                             vr = vr,
                             pr = pr,
                             hr = hr)
    x, v, eps, d, p, A, h = answer.solution(time = out_time)

    # Write the output to a text file.
    if file:
        with open(file, "w") as f:
            # Write a header
            if not noheader:
                f.write(
"""# Output from RiemannSolution using the arguments:
#          problem = %(problem)s
#                n = %(n)s
#               x0 = %(x0)s
#               x1 = %(x0)s
#           xdiaph = %(x0)s
#        gamma_gas = %(x0)s
#         out_time = %(out_time)s
#   dl, vl, pl, hl = %(dl)s, %(vl)s, %(pl)s, %(hl)s
#       dr, vr, pr = %(dr)s, %(vr)s, %(pr)s, %(hr)s
#
#            x                       rho                      vel                       P                        eps                       A                       h
""" % {"problem" : problem,
       "n"       : n,
       "x0"      : x0,
       "x1"      : x1,
       "xdiaph"  : xdiaph,
       "gamma_gas" : gamma_gas,
       "out_time"  : out_time,
       "dl"        : dl,
       "vl"        : vl,
       "pl"        : pl,
       "hl"        : hl,
       "dr"        : dr,
       "vr"        : vr,
       "pr"        : pr,
       "hr"        : hr})
            for xi, di, vi, pi, epsi, Ai, hi in zip(x, d, v, p, eps, A, h):
                f.write((7*"%20.17e     ") % (xi, di, vi, pi, epsi, Ai, hi) + "\n")

    # Plot the results to the screen via matplotlib (if available)
    if plot:
        try:
            import matplotlib.pyplot as plt
            fig = plt.figure(figsize=(plotsize, 2.0/3.0*plotsize))
            axes = []
            for i, (q, label) in enumerate([(d, "Density"),
                                            (v, "Velocity"),
                                            (p, "Pressure"),
                                            (eps, "Specific Thermal Energy"),
                                            (A, "Entropy"),
                                            (h, "Smoothing scale")]):
                axes.append(fig.add_subplot(2, 3, i + 1))
                plt.plot(x, q, linewidth=3)
                plt.title(label)
                qmin = min(q)
                qmax = max(q)
                qdiff = qmax - qmin
                axes[i].set_ylim(qmin - 0.1*qdiff, qmax + 0.1*qdiff)
            plt.show()

        except:
            print "ERROR: unable to import matplotlib for graphics."
            pass
