from math import *

from NodeGeneratorBase import *

from Spheral import Vector2d, Tensor2d, SymTensor2d, CylindricalBoundary, \
     rotationMatrix2d, testPointInBox2d
from SpheralTestUtilities import fuzzyEqual

#-------------------------------------------------------------------------------
# Specialized node generator to generate nodes ratioing from a 1D slab surface.
#-------------------------------------------------------------------------------
class GenerateRatioSlab2d(NodeGeneratorBase):
    
    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def __init__(self,
                 dxSurface, xratio,
                 dySurface, yratio,
                 rho,
                 xmin,
                 xmax,
                 nNodePerh = 2.01,
                 SPH = False,
                 flipx = False,
                 flipy = False):
        
        assert dxSurface > 0.0
        assert xratio > 0.0
        assert dySurface > 0.0
        assert yratio > 0.0
        assert xmin[0] < xmax[0]
        assert xmin[1] < xmax[1]
        assert nNodePerh > 0.0
        
        # If the user provided a constant for rho, then use the constantRho
        # class to provide this value.
        if type(rho) == type(1.0):
            self.rho = ConstantRho(rho)
        else:
            self.rho = rho

        self.x, self.y, self.m, self.H = [], [], [], []

        # Decide the actual ratios we're going to use to arrive at an integer number of radial bins.
        def adjustRatio(drStart, drRatio, rmin, rmax):
            if abs(drRatio - 1.0) > 1e-4:
                neff = max(1, int(log(1.0 - (rmax - rmin)*(1.0 - drRatio)/drStart)/log(drRatio) + 0.5))
                drStart = (rmax - rmin)*(1.0 - drRatio)/(1.0 - drRatio**neff)
            else:
                neff = max(1, int((rmax - rmin)/drStart + 0.5))
                drStart = (rmax - rmin)/neff
            return drStart, neff
        dxSurface, nxeff = adjustRatio(dxSurface, xratio, xmin[0], xmax[0])
        dySurface, nyeff = adjustRatio(dySurface, yratio, xmin[1], xmax[1])
        print "Adjusting initial spacing to (%g, %g) in order to create integer numbers of bins (%i, %i) to edges." % (dxSurface, dySurface, nxeff, nyeff)

        # Work our way in from the y surface.
        y1 = xmax[1]
        dy = dySurface
        while y1 > xmin[1]:
            y0 = max(xmin[1], y1 - dy)
            yi = 0.5*(y0 + y1)
            hy = nNodePerh*dy

            if flipy:
                yi = xmin[1] + xmax[1] - yi

            # Work our way in from the x surface.
            x1 = xmax[0]
            dx = dxSurface
            while x1 > xmin[0]:
                x0 = max(xmin[0], x1 - dx)
                xi = 0.5*(x0 + x1)
                hx = nNodePerh*dx

                if flipx:
                    xi = xmin[0] + xmax[0] - xi

                self.x.append(xi)
                self.y.append(yi)
                self.m.append(self.rho(Vector2d(xi, yi))*dx*dy)
                self.H.append(SymTensor2d(1.0/hx, 0.0, 0.0, 1.0/hy))
                
                x1 = x0
                dx *= xratio

            y1 = y0
            dy *= yratio

        # # Make sure the total mass is what we intend it to be, by applying
        # # a multiplier to the particle masses.
        # M0 = 0.0
        # for m in self.m:
        #     M0 += m
        # assert M0 > 0.0
        # M1 = (xmax[0] - xmin[0]) * (xmax[1] - xmin[1]) * rho
        # massCorrection = M1/M0
        # for i in xrange(len(self.m)):
        #     self.m[i] *= massCorrection
        # print "Applied a mass correction of %f to ensure total mass is %f." % (massCorrection, M1)

        # Have the base class break up the serial node distribution
        # for parallel cases.
        NodeGeneratorBase.__init__(self, True,
                                   self.x, self.y, self.m, self.H)
        return


    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localPosition(self, i):
        assert i >= 0 and i < len(self.x)
        assert len(self.x) == len(self.y)
        return Vector2d(self.x[i], self.y[i])
    
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
        return self.rho(self.localPosition(i))
    
    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        assert i >= 0 and i < len(self.H)
        return self.H[i]
    
