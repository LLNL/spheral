#!/usr/bin/env python3

import subprocess

# Helper function for executing commands stolen from uberenv
def sexe(cmd,ret_output=False,echo=True):
    """ Helper for executing shell commands. """
    if echo:
        print("[exe: {0}]".format(cmd))

    # If we want to return the output as string a print to stdout
    # in real-time we need to let subprocess print as normal to
    # PIPE and STDOUT. We then need to read it back ourselves and
    # append to an ouput string of our own making. There is no way
    # to do this with subprocess currently.
    if ret_output:
        p = subprocess.Popen(cmd,
                             shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             encoding='utf8')
        out = "";
        while True:
            realtime_output = p.stdout.readline()

            if realtime_output == '' and p.poll() is not None:
                break

            if realtime_output:
                print(realtime_output.strip(), flush=True)
                out += realtime_output

        if echo:
            print(out)
        return out
    else:
        # If we do not need to return the output as a string, run()
        # will suffice.
        p = subprocess.run(cmd, shell=True,
                           check=True, text=True)

def num_3d_cyl_nodes(rmin, rmax, zmin, zmax, thetamin, thetamax, nr0, nz0, Ntot):
    """
    This routine estimates values for nr and nz for a 3D cylindrical distribution
    to closely match a total desired number of SPH nodes

    Parameters:
    rmin, rmax: Inner/outer radial location
    zmin, zmax: Lower/upper axial locations
    thetamin, thetamax: Lower/upper polar angles
    nr0, nz0: Initial guesses for number of radial and axial SPH nodes
    Ntot: Desired total number of nodes
    Returns:
    nr, nz: Optimal number of radial and axial SPH nodes
    """
    import numpy as np
    from scipy.optimize import basinhopping

    class takestep:
        def __init__(self):
            self.rng = np.random.default_rng()
        def __call__(self, x):
            x += int(self.rng.uniform(0, 100)/10)
            return x
    is_2d = False
    if (zmax < zmin):
        ig = [nr0]
        zlen = 0.
        is_2d = True
    else:
        ig = [nr0, nz0]
        zlen = zmax - zmin
        if nr0*nz0 > Ntot:
            raise Exception("Initial guesses are too high")
    rlen = rmax - rmin
    Dtheta = thetamax - thetamin
    mintheta = int(Dtheta/(0.5*np.pi) + 0.5)
    def calc_ntot(x):
        if (len(x) > 1):
            nr = int(x[0])
            nz = int(x[1])
            dz = zlen / nz
            dr = rlen / nr
            maxd = max(dr, dz)
        else:
            nr = int(x)
            nz = 1
            dr = rlen / nr
            maxd = dr
        new_tot = 0
        for i in range(nr):
            ri = rmin + (i + 0.5)*dr
            n_th = max(mintheta, int(ri*Dtheta/maxd))
            new_tot += n_th
        return abs(new_tot*nz - Ntot)
    result = basinhopping(calc_ntot, ig, take_step=takestep())
    if is_2d:
        return result.x, 0
    else:
        nr, nz = result.x[0], result.x[1]
        return nr, nz

def num_2d_cyl_nodes(rmin, rmax, thetamax, nr0, Ntot):
    nr, nz = num_3d_cyl_nodes(rmin, rmax, 0., -1., 0., thetamax, nr0, 1, Ntot)
    return nr
