from math import *

from NodeGeneratorBase import *
from fitspline import *

from Spheral import Vector2d
from Spheral import Tensor2d
from Spheral import SymTensor2d

#-------------------------------------------------------------------------------
# Generate a node distribution by interpolating between two spline fitted
# curves.
#
#   p02  p12  p22
#   p01  p11  p21
#   p00  p10  p20
#-------------------------------------------------------------------------------
class PolarSplineNodeGenerator2d(self,
                                 radcurve1,
                                 radcurve2,
                                 thetamin,
                                 thetamax,
                                 origin,
                                 ncurve1,
                                 nbetween,
                                 rho,
                                 distributionType = "constantPerShell",
                                 nNodePerh = 2.01,
                                 SPH = False):

    # Pre-conditions.
    assert len(radcurve1) > 1
    assert len(radcurve1) == len(radcurve2)
    assert thetamin < thetamax
    assert ncurve1 > 1
    assert nbetween > 1
    assert rho > 0.0
    assert distributionType in ("constantPerShell",
                                "constantLinearStep")

    self.radcurve1 = radcurve1
    self.radcurve2 = radcurve2
    self.thetamin = thetamin
    self.thetamax = thetamax
    self.origin = origin
    self.ncurve1 = ncurve1
    self.nbetween = nbetween
    self.rho = rho
    self.distributionType = distributionType
    self.nNodePerh = nNodePerh
    self.SPH = SPH

    # The length of curve 1.
    nc = len(radcurve1)
    dtheta = (thetamax - thetamin)/(nc - 1)*pi
    thetas = [thetamin + i*dtheta for i in range(nc)]
    assert min([x >= thetamin and x <= thetamax for x in thetas])
    coefs = fitspline(radcurve1, thetas)
    L0 = self.computeLength(coefs, thetamin, thetamax)
    dL0 = L0/ncurve1

    # Build up the node positions by interpolating between the spline curves.
    self.x = []
    self.y = []
    self.m = []
    self.rho = []
    drad = [(radcurve2[i] - radcurve1[i])/(nbetween - 1) for i in range(nc)]
    for ishell in range(nbetween):

        # Compute the spline curves for this shell.
        rsp0 = [radcurve1[i] + drad[i]*ishell         for i in range(nc)]
        rsp1 = [radcurve1[i] + drad[i]*(ishell + 0.5) for i in range(nc)]
        rsp2 = [radcurve1[i] + drad[i]*(ishell + 1)   for i in range(nc)]
        assert min([rsp0[i] >= radcurve1[i] and rsp0[i] <= radcurve2[i] for i in range(nc)])
        assert min([rsp1[i] >= radcurve1[i] and rsp1[i] <= radcurve2[i] for i in range(nc)])
        assert min([rsp2[i] >= radcurve1[i] and rsp2[i] <= radcurve2[i] for i in range(nc)])
        coefs0 = fitspline(rsp0, thetas)
        coefs1 = fitspline(rsp1, thetas)
        coefs2 = fitspline(rsp2, thetas)

        # The length of the curves for this shell, and the step size for
        # points along the curves.
        L0 = self.computeLength(coefs0, thetamin, thetamax)
        L1 = self.computeLength(coefs1, thetamin, thetamax)
        L2 = self.computeLength(coefs2, thetamin, thetamax)
        if distributionType == "constantPerShell":
            nshell = ncurve1
        elif distributionType == "constantLinearStep":
            nshell = int(L0/dL0 + 0.5)
        dL0 = L0/nshell
        dL1 = L1/nshell
        dL2 = L2/nshell

        # The starting points on these curves.
        p20 = origin + Vector2d(cos(thetamin), sin(thetamin))*evaluniformspline(coefs0, thetamin)
        p21 = origin + Vector2d(cos(thetamin), sin(thetamin))*evaluniformspline(coefs1, thetamin)
        p22 = origin + Vector2d(cos(thetamin), sin(thetamin))*evaluniformspline(coefs2, thetamin)

        # Integrate along the splines, filling in the values for this shell.
        for i in range(nshell):
            p00 = p20
            p01 = p21
            p02 = p22

            # Find the points defining the area.
            p10 = self.integrateToLength(coefs0, p00, 0.5*dL0)
            p20 = self.integrateToLength(coefs0, p00,     dL0)
            p11 = self.integrateToLength(coefs1, p01, 0.5*dL1)
            p21 = self.integrateToLength(coefs1, p01,     dL1)
            p12 = self.integrateToLength(coefs2, p02, 0.5*dL2)
            p22 = self.integrateToLength(coefs2, p02,     dL2)

            # Approximate the area, assuming line segments.
            area = (self.quadArea(p00, p10, p11, p01) +
                    self.quadArea(p10, p20, p21, p11) +
                    self.quadArea(p01, p11, p12, p02) +
                    self.quadArea(p11, p21, p22, p12))

            # Append the node values.
            self.x.append(p11.x)
            self.y.append(p11.y)
            self.m.append(area * rho)
            self.rho.append(rho)
