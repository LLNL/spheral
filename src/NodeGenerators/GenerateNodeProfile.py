from math import *
import numpy as np

from NodeGeneratorBase import *

from Spheral import Vector1d, Tensor1d, SymTensor1d, \
     rotationMatrix1d, testPointInBox1d
from SpheralTestUtilities import fuzzyEqual

#-------------------------------------------------------------------------------
# Class to generate 1-D node positions for a fixed node mass to fit the given
# density profile in a range (xmin, xmax).
#-------------------------------------------------------------------------------
class GenerateNodeProfile1d(NodeGeneratorBase):

    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def __init__(self,
                 mi,                 # target mass per point
                 rho,                # density profile
                 xmin,
                 xmax,
                 nNodePerh = 2.01,
                 numbins = 10000):

        assert mi > 0.0
        assert xmin < xmax
        assert nNodePerh > 0.0

        # If the user provided a constant for rho, then use the constantRho
        # class to provide this value.
        if type(rho) == type(1.0):
            self.rhofunc = ConstantRho(rho)
        else:
            self.rhofunc = rho

        # Build the evenly sampled cumulative mass as a function of position.
        ok = False
        while not ok:
            dx = (xmax - xmin)/numbins
            mcum = np.cumsum(np.array([0.0] + [0.5*dx*(self.rhofunc(xmin + i*dx) + self.rhofunc(xmin + (i + 1)*dx)) for i in xrange(numbins)]))
            if mcum[-1]/mi > 0.5*numbins:
                numbins = int(2*mcum[-1]/mi)
                print "Warning, boosting numbins to %i to increase mass resolution for interpolation" % numbins
            else:
                ok = True

        # Now go through and bisect for positions to get the mass per point we want.
        xi = xmin
        self.x = []
        self.rho = []
        mtarget = -0.5*mi
        while xi < xmax:
            mtarget += mi
            if mtarget <= mcum[-1]:
                i = np.searchsorted(mcum, mtarget) - 1
                assert mtarget >= mcum[i] and mtarget <= mcum[i+1]
                xi = xmin + (i + (mtarget - mcum[i])/(mcum[i+1] - mcum[i]))*dx
                assert (xi >= xmin + i*dx) and (xi <= xmin + (i+1)*dx)
                self.x.append(xi)
                self.rho.append(self.rhofunc(xi))
            else:
                xi = xmax
        n = len(self.x)
        print "Generated %i 1D points." % n
        self.m = [mi]*n

        # Figure out the H.
        self.H = []
        for i in xrange(n):
            if i == 0:
                dxavg = self.x[i+1] - self.x[i]
            elif i == n-1:
                dxavg = self.x[i] - self.x[i-1]
            else:
                dxavg = 0.5*(self.x[i+1] - self.x[i-1])
            self.H.append(SymTensor1d(1.0/(nNodePerh*dxavg)))

        # Have the base class break up the serial node distribution
        # for parallel cases.
        NodeGeneratorBase.__init__(self, True,
                                   self.x, self.m, self.rho, self.H)

        return

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localPosition(self, i):
        assert i >= 0 and i < len(self.x)
        return Vector1d(self.x[i])

    #---------------------------------------------------------------------------
    # Get the mass for the given node index.
    #---------------------------------------------------------------------------
    def localMass(self, i):
        assert i >= 0 and i < len(self.m)
        return self.m[i]

    #---------------------------------------------------------------------------
    # Get the mass density for the given node index.
    #---------------------------------------------------------------------------
    def localMassDensity(self, i):
        assert i >= 0 and i < len(self.x)
        return self.rho[i]

    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        assert i >= 0 and i < len(self.H)
        return self.H[i]

