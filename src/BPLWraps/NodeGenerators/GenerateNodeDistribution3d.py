from math import *

from NodeGeneratorBase import *

from Spheral import Vector3d
from Spheral import Tensor3d
from Spheral import SymTensor3d

from loadmpi import loadmpi
mpi, procID, nProcs = loadmpi()

#-------------------------------------------------------------------------------
# Class to generate 3-D node positions.
#-------------------------------------------------------------------------------
class GenerateNodeDistribution3d(NodeGeneratorBase):

    #-------------------------------------------------------------------------------
    # Constructor
    #-------------------------------------------------------------------------------
    def __init__(self, n1, n2, n3, rho,
                 distributionType = 'lattice',
                 xmin = (0.0, 0.0, 0.0),
                 xmax = (1.0, 1.0, 1.0),
                 rmin = None,
                 rmax = None,
                 thetamin = None,
                 thetamax = None,
                 phimin = None,
                 phimax = None,
                 zmin = None,
                 zmax = None,
                 nNodePerh = 2.01,
                 SPH = False,
                 rejecter = None):

        # Check that the input parameters are good.
        assert n1 > 0
        assert n2 > 0
        assert n3 > 0 or distributionType == "cylindrical"
        assert rho > 0.0
        assert ((distributionType == 'lattice' and
                 xmin is not None and
                 xmax is not None) or
                (distributionType == 'line' and 
                 xmin is not None and 
                 xmax is not None) or
                (distributionType == 'cylindrical' and
                 rmin is not None and
                 rmax is not None and
                 thetamin is not None and
                 thetamax is not None and
                 zmin is not None and
                 zmax is not None) or
                (distributionType == 'constantNTheta' and
                 rmin is not None and
                 rmax is not None and
                 thetamin is not None and
                 thetamax is not None and
                 zmin is not None and
                 zmax is not None)                 
                )

        self.n1 = n1
        self.n2 = n2
        self.n3 = n3
        self.rho = rho
        self.distributionType = distributionType
        self.xmin = xmin
        self.xmax = xmax
        self.rmin = rmin
        self.rmax = rmax
        self.thetamin = thetamin
        self.thetamax = thetamax
        self.phimin = phimin
        self.phimax = phimax
        self.zmin = zmin
        self.zmax = zmax

        assert nNodePerh > 0.0

        self.nNodePerh = nNodePerh

        # Generate the internal sets of positions, masses, and H
        if distributionType == 'optimal':
            self.x, self.y, self.z, self.m, self.H = \
                    self.optimalSphericalDistribution(self.nRadial,
                                                      self.nTheta,
                                                      self.nPhi,
                                                      self.rho,
                                                      self.rmin,
                                                      self.rmax,
                                                      self.thetamin,
                                                      self.thetamax,
                                                      self.phimin,
                                                      self.phimax,
                                                      self.nNodePerh)
        elif distributionType == 'constantDAngle':
            self.x, self.y, self.z, self.m, self.H = \
                    self.constantDAngleSphericalDistribution(self.nRadial,
                                                             self.nTheta,
                                                             self.nPhi,
                                                             self.rho,
                                                             self.rmin,
                                                             self.rmax,
                                                             self.thetamin,
                                                             self.thetamax,
                                                             self.phimin,
                                                             self.phimax,
                                                             self.nNodePerh)
        elif distributionType == 'constantNShell':
            self.x, self.y, self.z, self.m, self.H = \
                    self.constantNShellSphericalDistribution(self.nRadial,
                                                             self.nTheta,
                                                             self.nPhi,
                                                             self.rho,
                                                             self.rmin,
                                                             self.rmax,
                                                             self.thetamin,
                                                             self.thetamax,
                                                             self.phimin,
                                                             self.phimax,
                                                             self.nNodePerh)

        elif distributionType == 'lattice':
            self.x, self.y, self.z, self.m, self.H = \
                    self.latticeDistribution(self.n1,      # nx
                                             self.n2,      # ny
                                             self.n3,      # nz
                                             self.rho,
                                             self.xmin,    # (xmin, ymin, zmin)
                                             self.xmax,    # (xmax, ymax, zmax)
                                             self.rmin,
                                             self.rmax,
                                             self.nNodePerh)

        elif distributionType == 'cylindrical':
            self.x, self.y, self.z, self.m, self.H = \
                    self.cylindrical(self.n1,      # nr
                                     self.n2,      # nz
                                     self.rho,
                                     self.rmin,
                                     self.rmax,
                                     self.thetamin,
                                     self.thetamax,
                                     self.zmin,
                                     self.zmax,
                                     self.nNodePerh)

        elif distributionType == 'constantNTheta':
            self.x, self.y, self.z, self.m, self.H = \
                    self.constantNTheta(self.n1,      # nr
                                        self.n2,      # ntheta
                                        self.n3,      # nz
                                        self.rho,
                                        self.rmin,
                                        self.rmax,
                                        self.thetamin,
                                        self.thetamax,
                                        self.zmin,
                                        self.zmax,
                                        self.nNodePerh)

        elif distributionType == 'line':
            self.x, self.y, self.z, self.m, self.H = \
                    self.lineDistribution(self.n1,      # nx
                                          self.n2,      # ny
                                          self.n3,      # nz
                                          self.rho,
                                          self.xmin,    # (xmin, ymin, zmin)
                                          self.xmax,    # (xmax, ymax, zmax)
                                          self.rmin,
                                          self.rmax,
                                          self.nNodePerh)


        # If the user provided a "rejecter", give it a pass
        # at the nodes.
        if rejecter:
            self.x, self.y, self.z, self.m, self.H = rejecter(self.x,
                                                              self.y,
                                                              self.z,
                                                              self.m,
                                                              self.H)

        # Have the base class break up the serial node distribution
        # for parallel cases.
        NodeGeneratorBase.__init__(self, True,
                                   self.x, self.y, self.z, self.m, self.H)

        # If SPH has been specified, make sure the H tensors are round.
        if SPH:
            self.makeHround()

        return

    #-------------------------------------------------------------------------------
    # Get the position for the given node index.
    #-------------------------------------------------------------------------------
    def localPosition(self, i):
        assert i >= 0 and i < len(self.x)
        assert len(self.x) == len(self.y)
        return Vector3d(self.x[i], self.y[i], self.z[i])

    #-------------------------------------------------------------------------------
    # Get the mass for the given node index.
    #-------------------------------------------------------------------------------
    def localMass(self, i):
        assert i >= 0 and i < len(self.m)
        return self.m[i]

    #-------------------------------------------------------------------------------
    # Get the mass density for the given node index.
    #-------------------------------------------------------------------------------
    def localMassDensity(self, i):
        assert i >= 0 and i < len(self.x)
        return self.rho

    #-------------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #-------------------------------------------------------------------------------
    def localHtensor(self, i):
        assert i >= 0 and i < len(self.H)
        return self.H[i]

    #-------------------------------------------------------------------------------
    # Force the H tensors to be round.
    #-------------------------------------------------------------------------------
    def makeHround(self):
        for i in xrange(len(self.H)):
            h0 = (self.H[i].Determinant())**(1.0/3.0)
            self.H[i] = SymTensor3d(h0, 0.0, 0.0,
                                    0.0, h0, 0.0,
                                    0.0, 0.0, h0)
        return

    #-------------------------------------------------------------------------------
    # This version tries to combine the optimal configuration of constant spacing
    # for the innermost node distribution, and constant nodes per ring for the outer
    # regions.
    #-------------------------------------------------------------------------------
    def optimalSphericalDistribution(self, nRadial, nTheta, rho,
                                     rmin = 0.0,
                                     rmax = 1.0,
                                     nNodePerh = 2.01,
                                     theta = pi/2.0):

        # Determine the radius at which we want to shift from constant spacing
        # to constant nodes per radius.
        dr = (rmax - rmin)/nRadial
        dTheta = theta/nTheta

        r0 = rmin + 2**(log(nTheta*theta/(2.0*pi))/log(2))*dr/theta
        #r0 = rmin + dr/dTheta

        print 'optimalSphericalDistribution: Cutoff radius is ', r0
        nRadial0 = max(0, min(nRadial, int((r0 - rmin)/dr)))
        nRadial1 = nRadial - nRadial0
        r0 = rmin + nRadial0*dr
        print 'Shifted to ', r0
        print nRadial0, nRadial1

        if nRadial0 and not nRadial1:
            # Only use constant spacing for the nodes.
            x, y, m, H = self.constantDThetaSphericalDistribution(nRadial0, rho,
                                                                  rmin,
                                                                  rmax,
                                                                  nNodePerh,
                                                                  theta)

        elif not nRadial0 and nRadial1:
            # Only use constant nodes per radial bin
            x, y, m, H = self.constantNThetaSphericalDistribution(nRadial1,
                                                                  nTheta,
                                                                  rho,
                                                                  rmin,
                                                                  rmax,
                                                                  nNodePerh,
                                                                  theta)

        else:
            # Combine the two schemes.
            x0, y0, m0, H0 = self.constantDThetaSphericalDistribution(nRadial0, rho,
                                                                      rmin,
                                                                      r0,
                                                                      nNodePerh,
                                                                      theta)
            x1, y1, m1, H1 = self.constantNThetaSphericalDistribution(nRadial1, nTheta, rho,
                                                                      r0,
                                                                      rmax,
                                                                      nNodePerh,
                                                                      theta)
            x, y, m, H = x0 + x1, y0 + y1, m0 + m1, H0 + H1

        return x, y, m, H

    #-------------------------------------------------------------------------------
    # Volume of a spherical section.
    #-------------------------------------------------------------------------------
    def sphericalSectionVolume(self,
                               rmin, rmax,
                               thetamin, thetamax,
                               phimin, phimax):
        return ((thetamax - thetamin)*(cos(phimin) - cos(phimax))*
                (rmax**3 - rmin**3)/3.0)

    #-------------------------------------------------------------------------------
    # Seed positions/masses for spherical symmetry.
    # This version tries to seed nodes in spherical shells with roughly constant
    # mass per radial shell.
    #-------------------------------------------------------------------------------
    def constantMassShells(self, nNodes, nRadialShells, rho,
                           rmin, rmax,
                           thetamin, thetamax,
                           phimin, phimax,
                           nNodePerh):

        assert nNodes > 0
        assert nRadialShells > 0 and nRadialShells < nNodes
        assert rho > 0.0
        assert rmin < rmax
        assert thetamin < thetamax
        assert phimin < phimax

        assert thetamin >= 0.0 and thetamin <= 2.0*pi
        assert thetamax >= 0.0 and thetamax <= 2.0*pi

        assert phimin >= 0.0 and phimin <= 2.0*pi
        assert phimax >= 0.0 and phimax <= 2.0*pi

        # Return values.
        x = []
        y = []
        z = []
        m = []
        H = []
        
        # Nominal mass per node.
        totalVolume = sphericalSectionVolume(rmin, rmax,
                                             thetamin, thetamax,
                                             phimin, phimax)
        totalMass = rho*totalVolume
        nominalmNode = totalMass/nNodes
        assert nominalmNode > 0.0
        nominalmShell = toalMass/nRadialShells
        assert nominalmShell > 0.0
        nominalnShell = int(nominalmShell/nominalmNode)
        assert nominalnShell > 0
        nominalVShell = totalVolume/nRadialShells
        assert nominalVShell > 0.0

        # The fixed H tensor value.
        averageSpacing = (totalVolume/nNodes)**(1.0/3.0)
        h0 = 1.0/(nNodePerh*averageSpacing)
        nominalH = SymTensor3d(h0, 0.0, 0.0,
                               0.0, h0, 0.0,
                               0.0, 0.0, h0)

        # Loop over each shell.
        rinner = 0.0
        router = 0.0
        for ishell in xrange(nRadial):

            # Radii of this shell.
            rinner = router
            deltar = (3.0*nominalVShell/((thetamax - thetamin)*
                                         (cos(phimin) - cos(phimax))) +
                      rinner**3)**(1.0/3.0)
            router = rinner + deltar
            rshell = rinner + 0.5*deltar
            assert router <= 1.00001*rmax
            assert rshell < rmax

            # How many nodes in this shell?
            if ishell < nRadial - 1:
                n = nshell
            else:
                n = nNodes - nshell*(nRadial - 1)
            assert n > 0

            # Mass per node.
            mnode = nominalmShell/n

            # Call the routine to generate the node positions for this
            # shell.
            shellx, shelly, shellz = self.shellPositions(rshell, n,
                                                         thetamin, thetamax,
                                                         phimin, phimax)

            # Append this shells results.
            x += shellx
            y += shelly
            z += shellz
            m += [mnode]*n
            H += [nominalH]*n

        # That's it, return the result.
        return x, y, z, m, H

    #-------------------------------------------------------------------------------
    # Helper method.  For a given shell, seed the requested number of nodes.
    # This method functions by plunking points down initially randomly in the given
    # spherical section, and then iterating over each point and having it push all
    # other points inside it's optimal area radially away.  Theoretically this
    # should converge on a pretty uniform distribution.
    #-------------------------------------------------------------------------------
    def shellPositions(self, r, n,
                       thetamin, thetamax,
                       phimin, phimax,
                       maxIterations = 100,
                       tolerance = 0.01,
                       seed = None):

        assert r > 0.0
        assert n > 0
        assert thetamin < thetamax
        assert phimin < phimax

        assert thetamin >= 0.0 and thetamin <= 2.0*pi
        assert thetamax >= 0.0 and thetamax <= 2.0*pi

        assert phimin >= 0.0 and phimin <= 2.0*pi
        assert phimax >= 0.0 and phimax <= 2.0*pi

        # The total area of this spherical section.
        totalArea = (thetamax - thetamin)*(cos(phimin) - cos(phimax))*r**2

        # The nominal area for each point.
        nominalArea = totalArea/n
        meanSpacing = sqrt(nominalArea/pi)

        # Seed the requested number of points randomly in the given spherical
        # section.  For now we just store the (theta, phi) coordinatess of each
        # point.
        theta = []
        phi = []
        import random
        g = random.Random(seed)
        for i in xrange(n):
            theta.append(g.uniform(thetamin, thetamax))
            phi.append(g.uniform(phimin, phimax))

        # Iterate over looping over the points, and push them apart so that each
        # gets the expected circular area.  There are two stopping criteria:
        # 1.  maxIterations: only up up to maxIterations iterations through the loop.
        # 2.  tolerance: stop if the maximum displacement of a node during
        #     an iteration is <= tolerance*meanSpacing.
        iter = 0
        tolDisplacement = tolerance*meanSpacing
        maxDisplacement = 2.0*tolDisplacement
        while iter < maxIterations and maxDisplacement > tolDisplacement:
            print 'Iteration %i for shell @ r=%f' % (iter, r)
            maxDisplacement = 0.0
            iter += 1

            # Loop over each point in the shell.
            for i in xrange(n):
                xi = r*theta[i]
                yi = r*phi[i]

                # Loop over every other point, and push it away if it's
                # within the mean seperation from this point.
                for j in xrange(n):
                    dx = r*theta[j] - xi
                    dy = r*phi[j] - yi
                    sep = sqrt(dx*dx + dy*dy)
                    if sep < meanSpacing:
                        xx = xi + meanSpacing*dx/sep
                        yy = yi + meanSpacing*dy/sep
                        theta[j] = xx/r
                        phi[j] = yy/r
                        maxDisplacement = max(maxDisplacement,
                                              meanSpacing - sep)

        # Convert the final (r, theta, phi) distribution to the (x,y,z)
        # return values.
        x = []
        y = []
        z = []
        for i in xrange(n):
            x.append(r*sin(phi[i])*cos(theta[i]))
            y.append(r*sin(phi[i])*sin(theta[i]))
            z.append(r*cos(phi[i]))

        return x, y, z

    #-------------------------------------------------------------------------------
    # Seed positions/masses for circular symmetry.
    # This version seeds a constant number of nodes per radial bin.
    #-------------------------------------------------------------------------------
    def constantNThetaSphericalDistribution(self, nRadial, nTheta, rho,
                                              rmin = 0.0,
                                              rmax = 1.0,
                                              nNodePerh = 2.01,
                                              theta = pi/2.0):

        from Spheral import Tensor3d
        from Spheral import SymTensor3d

        dr = (rmax - rmin)/nRadial
        dTheta = theta/nTheta
        hr = 1.0/(nNodePerh*dr)

        x = []
        y = []
        m = []
        H = []

        for i in xrange(0, nRadial):
            rInner = rmin + i*dr
            rOuter = rmin + (i + 1)*dr
            ri = rmin + (i + 0.5)*dr

            mRing = (rOuter**2 - rInner**2) * theta/2.0 * rho
            mi = mRing/nTheta

            hTheta = 1.0/(nNodePerh*ri*dTheta)

            for j in xrange(nTheta):
                thetai = (j + 0.5)*dTheta
                x.append(ri*cos(thetai))
                y.append(ri*sin(thetai))
                m.append(mi)

                Hi = SymTensor3d(hr, 0.0, 0.0, hTheta)
                rot = Tensor3d(cos(thetai), sin(thetai), -sin(thetai), cos(thetai))
                rotInv = rot.Transpose()
                H.append(((rotInv.dotsym(Hi)).dot(rot)).Symmetric())

        return x, y, m, H

    #-------------------------------------------------------------------------------
    # Seed positions/masses on a lattice
    #-------------------------------------------------------------------------------
    def latticeDistribution(self, nx, ny, nz, rho,
                            xmin,
                            xmax,
                            rmin,
                            rmax,
                            nNodePerh = 2.01):

        assert nx > 0
        assert ny > 0
        assert nz > 0
        assert rho > 0

        dx = (xmax[0] - xmin[0])/nx
        dy = (xmax[1] - xmin[1])/ny
        dz = (xmax[2] - xmin[2])/nz

        n = nx*ny*nz
        volume = ((xmax[0] - xmin[0])*
                  (xmax[1] - xmin[1])*
                  (xmax[2] - xmin[2]))
        m0 = rho*volume/n

        hx = 1.0/(nNodePerh*dx)
        hy = 1.0/(nNodePerh*dy)
        hz = 1.0/(nNodePerh*dz)
        H0 = SymTensor3d(hx, 0.0, 0.0,
                         0.0, hy, 0.0,
                         0.0, 0.0, hz)

        x = []
        y = []
        z = []
        m = []
        H = []

        for k in xrange(nz):
            for j in xrange(ny):
                for i in xrange(nx):
                    xx = xmin[0] + (i + 0.5)*dx
                    yy = xmin[1] + (j + 0.5)*dy
                    zz = xmin[2] + (k + 0.5)*dz
                    r = sqrt(xx*xx + yy*yy + zz*zz)
                    if ((r >= rmin or rmin is None) and
                        (r <= rmax or rmax is None)):
                        x.append(xx)
                        y.append(yy)
                        z.append(zz)
                        m.append(m0)
                        H.append(H0)

        return x, y, z, m, H

    #-------------------------------------------------------------------------------
    # Seed positions/masses on a line
    #-------------------------------------------------------------------------------
    def lineDistribution(self, nx, ny, nz, rho,
                         xmin,
                         xmax,
                         rmin,
                         rmax,
                         nNodePerh = 2.01):

        # Make sure at least two of the dimensions are 1 (and that the third
        # is greater than 1).
        assert (((nx > 1) and (ny == 1) and (nz == 1)) or \
                ((nx == 1) and (ny > 1) and (nz == 1)) or \
                ((nx == 1) and (ny == 1) and (nz > 1)))
        assert rho > 0

        dx = (xmax[0] - xmin[0])/nx
        dy = (xmax[1] - xmin[1])/ny
        dz = (xmax[2] - xmin[2])/nz
        ns = max([nx, ny, nz])
        ds = min([dx, dy, dz])

        volume = ns * ds
        m0 = rho*volume/ns

        # Compute the H tensor.
        hs = 1.0/(nNodePerh*ds)
        H0 = SymTensor3d(hs, 0.0, 0.0,
                         0.0, hs, 0.0,
                         0.0, 0.0, hs)

        x = []
        y = []
        z = []
        m = []
        H = []

        for k in xrange(nz):
            for j in xrange(ny):
                for i in xrange(nx):
                    xx = xmin[0] + (i + 0.5)*dx
                    yy = xmin[1] + (j + 0.5)*dy
                    zz = xmin[2] + (k + 0.5)*dz
                    r = sqrt(xx*xx + yy*yy + zz*zz)
                    if ((r >= rmin or rmin is None) and
                        (r <= rmax or rmax is None)):
                        x.append(xx)
                        y.append(yy)
                        z.append(zz)
                        m.append(m0)
                        H.append(H0)

        return x, y, z, m, H

    #-------------------------------------------------------------------------------
    # Seed positions/masses on a cylindrical section about the z-axis.
    #-------------------------------------------------------------------------------
    def cylindrical(self, nr, nz, rho,
                    rmin,
                    rmax,
                    thetamin,
                    thetamax,
                    zmin,
                    zmax,
                    nNodePerh = 2.01):

        assert nr > 0
        assert nz > 0
        assert rho > 0
        assert rmin < rmax
        assert zmin < zmax

        dr = (rmax - rmin)/nr
        dz = (zmax - zmin)/nz
        Dtheta = thetamax - thetamin

        hr = 1.0/(nNodePerh*dr)
        hz = 1.0/(nNodePerh*dz)

        minNtheta = int(Dtheta/(0.5*pi) + 0.5)

        x = []
        y = []
        z = []
        m = []
        H = []

        for iz in xrange(nz):
            zi = zmin + (iz + 0.5)*dz
            for ir in xrange(nr):
                rInner = rmin + ir*dr
                rOuter = rInner + dr
                ri = 0.5*(rInner + rOuter)
                circumference = ri*Dtheta
                nominalNtheta = int(circumference/max(dr, dz))
                ntheta = max(minNtheta, nominalNtheta)
                dtheta = Dtheta/ntheta
                if ir > 0:
                    mring = 0.5*(rOuter**2 - rInner**2)*Dtheta*dz*rho
                else:
                    mring = 0.5*(rOuter**2)*Dtheta*dz*rho
                mi = mring/ntheta
                htheta = 1.0/(nNodePerh*dtheta*ri)
                for itheta in xrange(ntheta):
                    thetai = thetamin + (itheta + 0.5)*dtheta
                    x.append(ri*cos(thetai))
                    y.append(ri*sin(thetai))
                    z.append(zi)
                    m.append(mi)
                    Hi = SymTensor3d(hr, 0.0, 0.0,
                                     0.0, htheta, 0.0,
                                     0.0, 0.0, hz)
                    R = Tensor3d(cos(thetai), sin(thetai), 0.0,
                                 -sin(thetai), cos(thetai), 0.0,
                                 0.0, 0.0, 1.0)
                    Hi.rotationalTransform(R)
                    H.append(Hi)

        return x, y, z, m, H

    #-------------------------------------------------------------------------------
    # Seed positions/masses on a cylindrical section about the z-axis.
    # In this case we maintain a constant number of nodes per angle.
    #-------------------------------------------------------------------------------
    def constantNTheta(self, nr, ntheta, nz, rho,
                       rmin,
                       rmax,
                       thetamin,
                       thetamax,
                       zmin,
                       zmax,
                       nNodePerh = 2.01):

        assert nr > 0
        assert ntheta > 0
        assert nz > 0
        assert rho > 0
        assert rmin < rmax
        assert zmin < zmax

        dr = (rmax - rmin)/nr
        Dtheta = thetamax - thetamin
        dtheta = Dtheta/ntheta
        dz = (zmax - zmin)/nz

        hr = 1.0/(nNodePerh*dr)
        hz = 1.0/(nNodePerh*dz)

        x = []
        y = []
        z = []
        m = []
        H = []

        for iz in xrange(nz):
            zi = zmin + (iz + 0.5)*dz
            for ir in xrange(nr):
                rInner = rmin + ir*dr
                rOuter = rInner + dr
                ri = 0.5*(rInner + rOuter)
                circumference = ri*Dtheta
                mring = 0.5*(rOuter**2 - rInner**2)*Dtheta*dz*rho
                mi = mring/ntheta
                htheta = 1.0/(nNodePerh*dtheta*ri)
                for itheta in xrange(ntheta):
                    thetai = thetamin + (itheta + 0.5)*dtheta
                    x.append(ri*cos(thetai))
                    y.append(ri*sin(thetai))
                    z.append(zi)
                    m.append(mi)
                    Hi = SymTensor3d(hr, 0.0, 0.0,
                                     0.0, htheta, 0.0,
                                     0.0, 0.0, hz)
                    R = Tensor3d(cos(thetai), sin(thetai), 0.0,
                                 -sin(thetai), cos(thetai), 0.0,
                                 0.0, 0.0, 1.0)
                    Hi.rotationalTransform(R)
                    H.append(Hi)

        return x, y, z, m, H

#--------------------------------------------------------------------------------
# Create cylindrical distributions similar to the RZ generator.
#--------------------------------------------------------------------------------
from GenerateNodeDistribution2d import GenerateNodeDistributionRZ
from NodeGenerators import generateCylDistributionFromRZ
class GenerateCylindricalNodeDistribution3d(GenerateNodeDistributionRZ):

    def __init__(self, nRadial, nTheta, rho,
                 distributionType = 'optimal',
                 xmin = None,
                 xmax = None,
                 rmin = None,
                 rmax = None,
                 nNodePerh = 2.01,
                 theta = pi/2.0,
                 phi = 2.0*pi,
                 SPH = False):
        GenerateNodeDistributionRZ.__init__(self,
                                            nRadial,
                                            nTheta,
                                            rho,
                                            distributionType,
                                            xmin,
                                            xmax,
                                            rmin,
                                            rmax,
                                            nNodePerh,
                                            theta,
                                            SPH)
        from Spheral import Vector3d, CylindricalBoundary

        # The base class already split the nodes up between processors, but
        # we want to handle that ourselves.  Distribute the full set of RZ
        # nodes to every process, then redecompose them below.
        self.x = mpi.allreduce(self.x[:], mpi.SUM)
        self.y = mpi.allreduce(self.y[:], mpi.SUM)
        self.m = mpi.allreduce(self.m[:], mpi.SUM)
        self.z = [0.0]*len(self.x)
        self.globalIDs = [0]*len(self.x)

        # Convert the 2-D H tensors to 3-D, and correct the masses.
        n = len(self.x)
        H2dlist = mpi.allreduce(self.H[:], mpi.SUM)
        assert len(H2dlist) == n
        self.H = []
        for i in xrange(n):
            xi = self.x[i]
            yi = self.y[i]
            H2d = H2dlist[i]

            hxy0 = 1.0/(H2d.eigenValues().minElement())
            dphi = CylindricalBoundary.angularSpacing(yi, hxy0, nNodePerh, 2.0)

            hz = dphi*yi*nNodePerh
            self.H.append(SymTensor3d(H2d.xx, H2d.xy, 0.0,
                                      H2d.yx, H2d.yy, 0.0,
                                      0.0,    0.0,    1.0/hz))
            if SPH:
                h0 = self.H[-1].Determinant()**(1.0/3.0)
                self.H[-1] = SymTensor3d.one * h0

            # Convert the mass to the full hoop mass, which will then be used in
            # generateCylDistributionFromRZ to compute the actual nodal masses.
            mi = self.m[i] # * xi/(xi*xi + 1.0e-60) # Undo what was done in the base class
            circ = 2.0*pi*yi
            mhoop = mi*circ
            self.m[i] = mhoop
            #nhoop = 2.0*pi*dphi/(dphi*dphi + 1.0e-60)
            #self.m[i] = mhoop*nhoop/(nhoop*nhoop + 1.0e-60)

        assert len(self.m) == n
        assert len(self.H) == n

        # Duplicate the nodes from the xy-plane, creating rings of nodes about
        # the x-axis.  We use a C++ helper method for the sake of speed.
        kernelExtent = 2.0
        extras = []
        generateCylDistributionFromRZ(self.x, self.y, self.z,
                                      self.m, self.H, self.globalIDs,
                                      extras,
                                      nNodePerh, kernelExtent, phi,
                                      procID, nProcs)
        return

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localPosition(self, i):
        from Spheral import Vector3d
        return Vector3d(self.x[i], self.y[i], self.z[i])

    #-------------------------------------------------------------------------------
    # Get the mass for the given node index.
    #-------------------------------------------------------------------------------
    def localMass(self, i):
        return self.m[i]

    #-------------------------------------------------------------------------------
    # Get the mass density for the given node index.
    #-------------------------------------------------------------------------------
    def localMassDensity(self, i):
        return self.rho(self.localPosition(i))

    #-------------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #-------------------------------------------------------------------------------
    def localHtensor(self, i):
        return self.H[i]
