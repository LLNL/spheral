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

double                // density, velocity, pressure, speed of sound
    dl, ul, pl, cl,   // in left region
    dr, ur, pr, cr;   // in right region

def RiemannSolution(problem = "Sod",  # ("", "Sod", "Einfeldt", "Stationary_contact", "Slow_shock", "Slow_contact_shock", "LeBlanc")
                    n = 1000,         # number of points in evaluating exact solution
                    x0 = None,        # box min coordinate
                    x1 = None,        # box max coordinate
                    xdiaph = None,    # position of diaphragm xdiaph \in [x0, x1]
                    gamma_gas = None, # ratio of specific heats
                    out_time = None,  # output time
                    dl = None,        # density (left state)
                    vl = None,        # velocity (left state)
                    pl = None,        # pressure (left state)
                    dr = None,        # density (right state)
                    vr = None,        # velocity (right state)
                    pr = None)        # pressure (right state)

    assert problem or (x0 and x1 and out_time and xdiaph and gamma_gas and dl and vl and pl and dr and vr and pr)

    # Check for one of the prepackaged test cases.
    twothirds = 2.0/3.0
    leblanc1 = twothirds*1e-1
    leblanc2 = twothirds*1e-10
    packaged_problems = {"sod"                 : (0.0, 1.0, 0.3, 1.4, 0.20, 1.0, 0.75, 1.0, 0.125, 0.0, 0.1),                        # TEST 1 (Modified Sod)
                         "123"                 : (0.0, 1.0, 0.5, 1.4, 0.15, 1.0, -2.0, 0.4, 1.0, 2.0, 0.4),                          # TEST 2 (123 problem)
                         "leftwc"              : (0.0, 1.0, 0.5, 1.4, 0.012, 1.0, 0.0, 1000.0, 1.0, 0.0, 0.01),                      # TEST 3 (Left Woodward & Colella)
                         "2shock_collision"    : (0.0, 1.0, 0.4, 1.4, 0.035, 5.99924, 19.5975, 460.894, 5.99242, -6.19633, 46.0950), # TEST 4 (Collision of 2 shocks)
                         "stationary_contact"  : (0.0, 1.0, 0.8, 1.4, 0.012, 1.0, -19.59745, 1000.0, 1.0, -19.59745, 0.01),          # TEST 5 (Stationary contact)
                         "slow_shock"          : (0.0, 1.0, 0.5, 1.4, 1.0, 3.857143, -0.810631, 10.33333, 1.0, -3.44, 1.0),          # TEST 6 (Slow shock)
                         "shock_contact_shock" : (0.0, 1.0, 0.5, 1.4, 0.3, 1.0, 0.5, 1.0, 1.25, -0.5, 1.0),                          # TEST 7 (Shock-Contact-Shock)
                         "leblanc"             : (0.0, 1.0, 0.3, 1.4, 0.5, 1.0, 0.0, leblanc1, 0.01, 0.0, leblanc2),                 # TEST 8 (LeBlanc)
                         }
    if problem:
        assert problem.lower() in packaged_problems
        x0, x1, xdiaph, gamma_gas, out_time, dl, vl, pl, dr, vr, pr = packaged_problems[problem.lower()]
    else:
        assert dl

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

    boxsize = x1 - x0

    #---------------------------------------------------------------------------
    # pvrpose: to provide a guessed value for pressure
    #          pm in the Star Region. The choice is made
    #          according to adaptive Riemann solver using
    #          the PVRS, TRRS and TSRS approximate
    #          Riemann solvers. See Sect. 9.5 of Chapt. 9 of Ref. 1
    #---------------------------------------------------------------------------
    def guessp:
        quser = 2.0

        # compute guess pressure from PVRS Riemann solver
        cup = 0.25*(dl + dr)*(cl + cr)
        ppv = 0.5*(pl + pr) + 0.5*(vl - vr)*cup
        ppv = max(0.0, ppv)
        pmin = min(pl,  pr)
        pmax = max(pl,  pr)
        qmax = pmax/pmin

        if (qmax <= quser && (pmin <= ppv && ppv <= pmax)):
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
            i += 1
            fl, fld = prefun(pold, dl, pl, cl)
            fr, frd = prefun(pold, dr, pr, cr)
            p = pold - (fl + fr + udiff)/(fld + frd)
            change = 2.0*abs((p - pold)/(p + pold))
            print '\t', i, "\t\t" << change
            if (p < 0.0):
                p = tolpre
            pold = p
        if (i > nriter):
            print "divergence in Newton-Raphson iteration"

        # compute velocity in star region
        u = 0.5*(vl + vr + fr - fl)
        print ("----------------------------------------\n"
               "     Pressure           Velocity\n"
               "----------------------------------------\n"
               "     ", p/pscale, "\t\t", u, '\n'
               "----------------------------------------")
        return p, u

    #---------------------------------------------------------------------------
    # purpose: to sample the solution throughout the wave
    #          pattern. Pressure pm and velocity vm in the
    #          star region are known. Sampling is performed
    #          in terms of the 'speed' s = x/t. Sampled
    #          values are d, u, p
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
                    u = vl
                    p = pl
                else:
                    cml = cl*pow(pm/pl, g1)
                    stl = vm - cml
                if (s > stl):
                    # sampled point is star left state
                    d = dl*pow(pm/pl, 1.0/gamma_gas)
                    u = vm
                    p = pm
                else:
                    # sampled point is inside left fan
                    u = g5*(cl + g7*vl + s)
                    c = g5*(cl + g7*(vl - s))
                    d = dl*pow(c/cl, g4)
                    p = pl*pow(c/cl, g3)
            else:
                # left shock
                pml = pm/pl
                sl = vl - cl*sqrt(g2*pml + g1)
                if (s <= sl):
                    # sampled point is left data state
                    d = dl
                    u = vl
                    p = pl
                else:
                    # sampled point is star left state
                    d = dl*(pml + g6)/(pml*g6 + 1.0)
                    u = vm
                    p = pm
        else:
            # sampling point lies to the right of the contact discontinuity
            if (pm > pr):
                # right shock
                pmr = pm/pr
                sr  = vr + cr*sqrt(g2*pmr + g1)
                if (s >= sr):
                    # sampled point is right data state
                    d = dr
                    u = vr
                    p = pr
                else:
                    # sampled point is star right state
                    d = dr*(pmr + g6)/(pmr*g6 + 1.0)
                    u = vm
                    p = pm
            else:
                # right rarefaction
                shr = vr + cr
                if (s >= shr):
                    # sampled point is right data state
                    d = dr
                    u = vr
                    p = pr
                else:
                    cmr = cr*pow(pm/pr, g1)
                    str = vm + cmr
                    if (s <= str):
                        # sampled point is star right state
                        d = dr*pow(pm/pr, 1.0/gamma_gas)
                        u = vm
                        p = pm
                    else:
                        # sampled point is inside left fan
                        u = g5*(-cr + g7*vr + s)
                        c = g5*(cr - g7*(vr - s))
                        d = dr*pow(c/cr, g4)
                        p = pr*pow(c/cr, g3)
        return d, u, p


    # the pressure positivity condition is tested for
    if (g4*(cl+cr) <= (vr-vl)):
        raise RunTimeError, ("the initial data is such that vacuum is generated"
                             "\nstopping program")

    # exact solution for pressure and velocity in star region is found
    pm, vm = starpu(1.0)
    dx = boxsize/n

    # complete solution at time out_time is found
    with open(filename) as f:
        for i in xrange(n):
            xpos = x0 + (i - 0.5)*dx
            s  = (xpos - xdiaph)/out_time

            # solution at point (x,t) = (xpos-xdiaph, out_time) is found
            ds, us, ps = sample(pm, vm, s)

            # exact solution profiles are written to data file
            f.write((5*"%20.17e     ") % (xpos, ds, us, ps, ps/ds/g8) + "\n")
    return
