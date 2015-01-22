from math import *

from NodeGeneratorBase import *

from Spheral import Vector3d
from Spheral import Tensor3d
from Spheral import SymTensor3d

from Spheral import vector_of_int, vector_of_double, vector_of_SymTensor3d, vector_of_vector_of_double

import mpi
procID = mpi.rank
nProcs = mpi.procs

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
                 origin = (0.0, 0.0, 0.0),  # For rmin, rmax
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
        assert isinstance(rho,(float,int)) and  (rho>0.0)
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
                 zmax is not None) or 
                (distributionType == 'hcp' and
                 xmin is not None and
                 xmax is not None)
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
        self.origin = origin
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
            serialInitialization = True
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
            serialInitialization = True
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
            serialInitialization = True
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
            serialInitialization = False
            self.x, self.y, self.z, self.m, self.H = \
                self.latticeDistribution(self.n1,      # nx
                                         self.n2,      # ny
                                         self.n3,      # nz
                                         self.rho,
                                         self.xmin,    # (xmin, ymin, zmin)
                                         self.xmax,    # (xmax, ymax, zmax)
                                         self.rmin,
                                         self.rmax,
                                         self.origin,
                                         self.nNodePerh)

        elif distributionType == 'cylindrical':
            serialInitialization = True
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
            serialInitialization = True
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
            serialInitialization = True
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

        elif distributionType == 'hcp':
            serialInitialization = False
            self.x, self.y, self.z, self.m, self.H = \
                self.hcpDistribution(self.n1,      # nx
                                     self.n2,      # ny
                                     self.n3,      # nz
                                     self.rho,
                                     self.xmin,    # (xmin, ymin, zmin)
                                     self.xmax,    # (xmax, ymax, zmax)
                                     self.rmin,
                                     self.rmax,
                                     self.origin,
                                     self.nNodePerh)

        # If the user provided a "rejecter", give it a pass
        # at the nodes.
        if rejecter:
            self.x, self.y, self.z, self.m, self.H = rejecter(self.x,
                                                              self.y,
                                                              self.z,
                                                              self.m,
                                                              self.H)

        # Make rho a list
        self.rho = [self.rho] * len(self.m)

        # Initialize the base class.  If "serialInitialization" is True, this
        # is where the points are broken up between processors as well.
        NodeGeneratorBase.__init__(self, serialInitialization,
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
        assert i >= 0 and i < len(self.rho)
        return self.rho[i]

    #-------------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #-------------------------------------------------------------------------------
    def localHtensor(self, i):
        assert i >= 0 and i < len(self.H)
        return self.H[i]

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
                            origin,
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

        imin, imax = self.globalIDRange(nx*ny*nz)
        for iglobal in xrange(imin, imax):
            i = iglobal % nx
            j = (iglobal // nx) % ny
            k = iglobal // (nx*ny)
            xx = xmin[0] + (i + 0.5)*dx
            yy = xmin[1] + (j + 0.5)*dy
            zz = xmin[2] + (k + 0.5)*dz
            r2 = (xx - origin[0])**2 + (yy - origin[1])**2 + (zz - origin[2])**2
            if ((rmin is None or r2 >= rmin**2) and
                (rmax is None or r2 <= rmax**2)):
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
                mring = 0.5*(rOuter**2 - rInner**2)*Dtheta*dz*rho
                ## if ir > 0:
##                     mring = 0.5*(rOuter**2 - rInner**2)*Dtheta*dz*rho
##                 else:
##                     mring = 0.5*(rOuter**2)*Dtheta*dz*rho
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

    #-------------------------------------------------------------------------------
    # Seed positions/masses on an hcp lattice
    #-------------------------------------------------------------------------------
    def hcpDistribution(self, nx, ny, nz, rho,
                        xmin,
                        xmax,
                        rmin,
                        rmax,
                        origin,
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

        imin, imax = self.globalIDRange(nx*ny*nz)
        for iglobal in xrange(imin, imax):
            i = iglobal % nx
            j = (iglobal // nx) % ny
            k = iglobal // (nx*ny)
            xx = xmin[0] + (i + 0.5*((j % 2) + (k % 2)))*dx
            yy = xmin[1] + (j + 0.5*(k % 2))*dy
            zz = xmin[2] + (k + 0.5)*dz

            # xx = xmin[0] + 0.5*dx*(2*i + (j % 2) + (k % 2))
            # yy = xmin[1] + 0.5*dy*sqrt(3.0)/3.0*(j + (k % 2))
            # zz = xmin[2] + dz*sqrt(6.0)/3.0*k
            r2 = (xx - origin[0])**2 + (yy - origin[1])**2 + (zz - origin[2])**2
            if ((rmin is None or r2 >= rmin**2) and
                (rmax is None or r2 <= rmax**2)):
                x.append(xx)
                y.append(yy)
                z.append(zz)
                m.append(m0)
                H.append(H0)

        return x, y, z, m, H

#--------------------------------------------------------------------------------
# Create cylindrical distributions similar to the RZ generator.
#--------------------------------------------------------------------------------
from GenerateNodeDistribution2d import GenerateNodeDistributionRZ
from Spheral import generateCylDistributionFromRZ
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
        self.H = mpi.allreduce(self.H[:], mpi.SUM)
        n = len(self.x)
        self.z = [0.0]*n
        self.globalIDs = [0]*n

        # Convert the 2-D H tensors to 3-D, and correct the masses.
        for i in xrange(n):
            xi = self.x[i]
            yi = self.y[i]
            H2d = self.H[i]
            H2dinv = H2d.Inverse()

            hxy0 = 0.5*(H2dinv.Trace())
            dphi = CylindricalBoundary.angularSpacing(yi, hxy0, nNodePerh, 2.0)
            assert dphi > 0.0
            nsegment = max(1, int(phi/dphi + 0.5))
            dphi = phi/nsegment

            hz = dphi*yi*nNodePerh
            self.H[i] = SymTensor3d(H2d.xx, H2d.xy, 0.0,
                                    H2d.yx, H2d.yy, 0.0,
                                    0.0,    0.0,    1.0/hz)
            if SPH:
                h0 = self.H[-1].Determinant()**(1.0/3.0)
                self.H[-1] = SymTensor3d.one * h0

            # Convert the mass to the full hoop mass, which will then be used in
            # generateCylDistributionFromRZ to compute the actual nodal masses.
            mi = self.m[i]
            circ = 2.0*pi*yi
            mhoop = mi*circ
            self.m[i] = mhoop

        assert len(self.m) == n
        assert len(self.H) == n

        # Duplicate the nodes from the xy-plane, creating rings of nodes about
        # the x-axis.  We use a C++ helper method for the sake of speed.
        kernelExtent = 2.0
        extras = []
        xvec = self.vectorFromList(self.x, vector_of_double)
        yvec = self.vectorFromList(self.y, vector_of_double)
        zvec = self.vectorFromList(self.z, vector_of_double)
        mvec = self.vectorFromList(self.m, vector_of_double)
        Hvec = self.vectorFromList(self.H, vector_of_SymTensor3d)
        globalIDsvec = self.vectorFromList(self.globalIDs, vector_of_int)
        extrasVec = vector_of_vector_of_double()
        for extra in extras:
            extrasVec.append(self.vectorFromList(extra, vector_of_double))
        generateCylDistributionFromRZ(xvec, yvec, zvec, mvec, Hvec, globalIDsvec,
                                      extrasVec,
                                      nNodePerh, kernelExtent, phi,
                                      procID, nProcs)
        self.x = [x for x in xvec]
        self.y = [x for x in yvec]
        self.z = [x for x in zvec]
        self.m = [x for x in mvec]
        self.H = [SymTensor3d(x) for x in Hvec]
        self.globalIDs = [x for x in globalIDsvec]
        for i in xrange(len(extras)):
            extras[i] = [x for x in extrasVec[i]]

        return

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localPosition(self, i):
        from Spheral import Vector3d
        return Vector3d(self.x[i], self.y[i], self.z[i])

    #---------------------------------------------------------------------------
    # Get the mass for the given node index.
    #---------------------------------------------------------------------------
    def localMass(self, i):
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
        return self.H[i]

#-------------------------------------------------------------------------------
# Specialized version that generates a variable radial stepping to try
# and match a given density profile with (nearly) constant mass nodes.  This
# only supports random node placement in each shell.
#-------------------------------------------------------------------------------
class GenerateRandomNodesMatchingProfile3d(NodeGeneratorBase):

    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def __init__(self, n, densityProfileMethod,
                 rmin = 0.0,
                 rmax = 1.0,
                 thetaMin = 0.0,
                 thetaMax = pi,
                 phiMin = 0.0,
                 phiMax = 2.0*pi,
                 nNodePerh = 2.01):
        
        assert n > 0
        assert rmin < rmax
        assert thetaMin < thetaMax
        assert thetaMin >= 0.0 and thetaMin <= 2.0*pi
        assert thetaMax >= 0.0 and thetaMax <= 2.0*pi
        assert phiMin < phiMax
        assert phiMin >= 0.0 and phiMin <= 2.0*pi
        assert phiMax >= 0.0 and phiMax <= 2.0*pi
        assert nNodePerh > 0.0
        
        self.n = n
        self.rmin = rmin
        self.rmax = rmax
        self.thetaMin = thetaMin
        self.thetaMax = thetaMax
        self.phiMin = phiMin
        self.phiMax = phiMax
        self.nNodePerh = nNodePerh
        
        # If the user provided a constant for rho, then use the constantRho
        # class to provide this value.
        if type(densityProfileMethod) == type(1.0):
            self.densityProfileMethod = ConstantRho(densityProfileMethod)
        else:
            self.densityProfileMethod = densityProfileMethod
        
        # Determine how much total mass there is in the system.
        self.totalMass = self.integrateTotalMass(self.densityProfileMethod,
                                                 rmin, rmax,
                                                 thetaMin, thetaMax,
                                                 phiMin, phiMax)
        print "Total mass of %g in the range r = (%g, %g), theta = (%g, %g), phi = (%g, %g)" % \
            (self.totalMass, rmin, rmax, thetaMin, thetaMax, phiMin, phiMax)

        # Now set the nominal mass per node.
        self.m0 = self.totalMass/(4.0/3.0*pi*pow(self.n,3))
        assert self.m0 > 0.0
        print "Nominal mass per node of %g." % self.m0

        # OK, we now know enough to generate the node positions.
        from Spheral import SymTensor3d
        import random
        self.x = []
        self.y = []
        self.z = []
        self.m = []
        self.H = []
        ri = rmax
        
        while ri > rmin:
            rhoi = densityProfileMethod(ri)
            dr = pow(self.m0/(4.0/3.0*pi*rhoi),1.0/3.0)

            hi = nNodePerh*dr
            Hi = SymTensor3d(1.0/hi, 0.0, 0.0,
                             0.0, 1.0/hi, 0.0,
                             0.0, 0.0, 1.0/hi)
            mshell  = rhoi * 4.0*pi*ri*ri*dr
            nshell  = int(mshell / self.m0)
            mi      = self.m0 * (mshell/(nshell*self.m0))
            
            for n in xrange(nshell):
                random.seed()
                u       = random.random()
                v       = random.random()
                w       = random.random()
                theta   = 2.0*pi*u
                phi     = acos(2.0*v-1)
                deltar  = w*(0.45*dr)
                self.x.append((ri-deltar)*cos(theta)*sin(phi))
                self.y.append((ri-deltar)*sin(theta)*sin(phi))
                self.z.append((ri-deltar)*cos(phi))
                self.m.append(mi)
                self.H.append(Hi)
                    
            # Decrement to the next radial bin inward.
            ri = max(0.0, ri - dr)


        print "Generated a total of %i nodes." % len(self.x)

        # Make sure the total mass is what we intend it to be, by applying
        # a multiplier to the particle masses.
        sumMass = 0.0
        for m in self.m:
            sumMass += m
        assert sumMass > 0.0
        massCorrection = self.totalMass/sumMass
        for i in xrange(len(self.m)):
            self.m[i] *= massCorrection
        print "Applied a mass correction of %f to ensure total mass is %f." % (massCorrection, self.totalMass)

        # Have the base class break up the serial node distribution
        # for parallel cases.
        NodeGeneratorBase.__init__(self, True,
                                self.x, self.y, self.z, self.m, self.H)
        return

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localPosition(self, i):
        assert i >= 0 and i < len(self.x)
        assert len(self.x) == len(self.y) == len(self.z)
        return Vector3d(self.x[i], self.y[i], self.z[i])

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
        return self.densityProfileMethod(self.localPosition(i).magnitude())

    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        assert i >= 0 and i < len(self.H)
        return self.H[i]


    #---------------------------------------------------------------------------
    # Numerically integrate the given density profile to determine the total
    # enclosed mass.
    #---------------------------------------------------------------------------
    def integrateTotalMass(self, densityProfileMethod,
                           rmin, rmax,
                           thetaMin, thetaMax,
                           phiMin, phiMax,
                           nbins = 10000):
        assert nbins > 0
        assert nbins % 2 == 0
        
        result = 0
        dr = (rmax-rmin)/nbins
        for i in xrange(1,nbins):
            r1 = rmin + (i-1)*dr
            r2 = rmin + i*dr
            result += 0.5*dr*(r2*r2*densityProfileMethod(r2)+r1*r1*densityProfileMethod(r1))
        result = result * (phiMax-phiMin) * (cos(thetaMin)-cos(thetaMax))
        return result


#-------------------------------------------------------------------------------
# Specialized version that generates a variable radial stepping to try
# and match a given density profile with (nearly) constant mass nodes.  This
# only supports the equivalent of the constant DthetaDphi method.
#-------------------------------------------------------------------------------
class GenerateNodesMatchingProfile3d(NodeGeneratorBase):

    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def __init__(self, n, densityProfileMethod,
                 rmin = 0.0,
                 rmax = 1.0,
                 thetaMin = 0.0,
                 thetaMax = pi,
                 phiMin = 0.0,
                 phiMax = 2.0*pi,
                 nNodePerh = 2.01):
        
        assert n > 0
        assert rmin < rmax
        assert thetaMin < thetaMax
        assert thetaMin >= 0.0 and thetaMin <= 2.0*pi
        assert thetaMax >= 0.0 and thetaMax <= 2.0*pi
        assert phiMin < phiMax
        assert phiMin >= 0.0 and phiMin <= 2.0*pi
        assert phiMax >= 0.0 and phiMax <= 2.0*pi
        assert nNodePerh > 0.0

        self.n = n
        self.rmin = rmin
        self.rmax = rmax
        self.thetaMin = thetaMin
        self.thetaMax = thetaMax
        self.phiMin = phiMin
        self.phiMax = phiMax
        self.nNodePerh = nNodePerh

        # If the user provided a constant for rho, then use the constantRho
        # class to provide this value.
        if type(densityProfileMethod) == type(1.0):
            self.densityProfileMethod = ConstantRho(densityProfileMethod)
        else:
            self.densityProfileMethod = densityProfileMethod

        # Determine how much total mass there is in the system.
        self.totalMass = self.integrateTotalMass(self.densityProfileMethod,
                                                 rmin, rmax,
                                                 thetaMin, thetaMax,
                                                 phiMin, phiMax)
        print "Total mass of %g in the range r = (%g, %g), theta = (%g, %g), phi = (%g, %g)" % \
            (self.totalMass, rmin, rmax, thetaMin, thetaMax, phiMin, phiMax)
                
        # Now set the nominal mass per node.
        self.m0 = self.totalMass/(4.0/3.0*pi*pow(self.n,3))
        assert self.m0 > 0.0
        print "Nominal mass per node of %g." % self.m0
            
        # OK, we now know enough to generate the node positions.
        self.x, self.y, self.z, self.m, self.H = \
            self.constantDThetaDPhiDistribution(self.densityProfileMethod,
                                                self.m0, self.n,
                                                rmin, rmax,
                                                thetaMin, thetaMax,
                                                phiMin, phiMax,
                                                nNodePerh)
        print "Generated a total of %i nodes." % len(self.x)
            
        # Make sure the total mass is what we intend it to be, by applying
        # a multiplier to the particle masses.
        sumMass = 0.0
        for m in self.m:
            sumMass += m
        assert sumMass > 0.0
        massCorrection = self.totalMass/sumMass
        for i in xrange(len(self.m)):
            self.m[i] *= massCorrection
        print "Applied a mass correction of %f to ensure total mass is %f." % (massCorrection, self.totalMass)

        # Have the base class break up the serial node distribution
        # for parallel cases.
        NodeGeneratorBase.__init__(self, True,
                                   self.x, self.y, self.z, self.m, self.H)
        return
            
    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localPosition(self, i):
        assert i >= 0 and i < len(self.x)
        assert len(self.x) == len(self.y) == len(self.z)
        return Vector3d(self.x[i], self.y[i], self.z[i])
    
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
        return self.densityProfileMethod(self.localPosition(i).magnitude())
    
    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        assert i >= 0 and i < len(self.H)
        return self.H[i]


    #---------------------------------------------------------------------------
    # Numerically integrate the given density profile to determine the total
    # enclosed mass.
    #---------------------------------------------------------------------------
    def integrateTotalMass(self, densityProfileMethod,
                           rmin, rmax,
                           thetaMin, thetaMax,
                           phiMin, phiMax,
                           nbins = 10000):
        assert nbins > 0
        assert nbins % 2 == 0
        
        result = 0
        dr = (rmax-rmin)/nbins
        for i in xrange(1,nbins):
            r1 = rmin + (i-1)*dr
            r2 = rmin + i*dr
            result += 0.5*dr*(r2*r2*densityProfileMethod(r2)+r1*r1*densityProfileMethod(r1))
        result = result * (phiMax-phiMin) * (cos(thetaMin)-cos(thetaMax))
        return result

    #---------------------------------------------------------------------------
    # Seed positions/masses for circular symmetry
    # This version tries to seed nodes in circular rings with constant spacing.
    #---------------------------------------------------------------------------
    def constantDThetaDPhiDistribution(self, densityProfileMethod,
                                       m0, nr,
                                       rmin, rmax,
                                       thetaMin, thetaMax,
                                       phiMin, phiMax,
                                       nNodePerh = 2.01):
        
        from Spheral import SymTensor3d
        
        # Return lists for positions, masses, and H's.
        x = []
        y = []
        z = []
        m = []
        H = []
        
        # Start at the outermost radius, and work our way inward.
        arcTheta = thetaMax - thetaMin
        arcPhi = phiMax - phiMin
        ri = rmax
        
        while ri > rmin:
            
            # Get the nominal delta r, delta theta, delta phi, number of nodes,
            # and mass per node at this radius.
            rhoi = densityProfileMethod(ri)
            dr = pow(m0/(4.0/3.0*pi*rhoi),1.0/3.0)
            dTheta = 2.0*dr
            dPhi = 2.0*dr
            arclength = arcTheta*ri
            nTheta = max(1,int(arclength/dTheta))
            dTheta = arcTheta/nTheta
            hi = nNodePerh*0.5*(dr + ri*dTheta)
            Hi = SymTensor3d(1.0/hi, 0.0, 0.0,
                             0.0, 1.0/hi, 0.0,
                             0.0, 0.0, 1.0/hi)
                             
            
                             
            # Now assign the nodes for this radius.
            for i in xrange(nTheta):
                thetai = thetaMin + (i + 0.5)*dTheta
                rp = ri*sin(thetai)
                arclength = arcPhi * rp
                nPhi = max(1,int(arclength/dPhi))
                dPhij = arcPhi/nPhi
                #print "at r=%g theta=%g rp=%g nTheta=%d nPhi=%d" % (ri,thetai,rp,nTheta,nPhi)
                for j in xrange(nPhi):
                    phij = phiMin + (j+0.5)*dPhij
                    x.append(ri*cos(thetai)*sin(phij))
                    y.append(ri*sin(thetai)*sin(phij))
                    z.append(ri*cos(phij))
                    m.append(m0)
                    H.append(Hi)

            # Decrement to the next radial bin inward.
            ri = max(0.0, ri - dr)
        
        return x, y, z, m, H

#-------------------------------------------------------------------------------
# Specialized version that generates a variable radial stepping to try
# and match a given density profile with (nearly) constant mass nodes
# on a variable resolution Icosahedron.
# It is recommended for now that you use 0-pi and 0-2pi for theta,phi.
#-------------------------------------------------------------------------------
class GenerateIcosahedronMatchingProfile3d(NodeGeneratorBase):
    
    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def __init__(self, n, densityProfileMethod,
                 rmin = 0.0,
                 rmax = 1.0,
                 thetaMin = 0.0,
                 thetaMax = pi,
                 phiMin = 0.0,
                 phiMax = 2.0*pi,
                 nNodePerh = 2.01):
        
        assert n > 0
        assert rmin < rmax
        assert thetaMin < thetaMax
        assert thetaMin >= 0.0 and thetaMin <= 2.0*pi
        assert thetaMax >= 0.0 and thetaMax <= 2.0*pi
        assert phiMin < phiMax
        assert phiMin >= 0.0 and phiMin <= 2.0*pi
        assert phiMax >= 0.0 and phiMax <= 2.0*pi
        assert nNodePerh > 0.0
        
        self.n = n
        self.rmin = rmin
        self.rmax = rmax
        self.thetaMin = thetaMin
        self.thetaMax = thetaMax
        self.phiMin = phiMin
        self.phiMax = phiMax
        self.nNodePerh = nNodePerh
        
        # If the user provided a constant for rho, then use the constantRho
        # class to provide this value.
        if type(densityProfileMethod) == type(1.0):
            self.densityProfileMethod = ConstantRho(densityProfileMethod)
        else:
            self.densityProfileMethod = densityProfileMethod
        
        # Determine how much total mass there is in the system.
        self.totalMass = self.integrateTotalMass(self.densityProfileMethod,
                                                 rmin, rmax,
                                                 thetaMin, thetaMax,
                                                 phiMin, phiMax)
        print "Total mass of %g in the range r = (%g, %g), theta = (%g, %g), phi = (%g, %g)" % \
            (self.totalMass, rmin, rmax, thetaMin, thetaMax, phiMin, phiMax)

        # Now set the nominal mass per node.
        self.m0 = self.totalMass/(4.0/3.0*pi*pow(self.n,3))
        assert self.m0 > 0.0
        print "Nominal mass per node of %g." % self.m0
        
        # Create the rotation matrices
        self.R = []
        for i in xrange(20):
            self.R.append([[0 for x in xrange(3)] for x in xrange(3)])
        self.compute_matrices()
        self.v = []
        for i in xrange(12):
            self.v.append([0,0,0])
        self.compute_corners()
        
        # pixel count goes as 40*res*(res-1)+12
        # res = 1 -> 0  pixels per face,  12 pixels total
        # res = 2 -> 4  pixels per face,  92 pixels total
        # res = 3 -> 12 pixels per face, 252 pixels total
        # res = 4 -> 24 pixels per face, 492 pixels total
        
        from Spheral import SymTensor3d
        self.x = []
        self.y = []
        self.z = []
        self.m = []
        self.H = []
        ri = rmax
        
        while ri > rmin:
            
            # Get the nominal delta r, number of nodes,
            # and mass per node at this radius.
            rhoi = densityProfileMethod(ri)
            dr = pow(self.m0/(4.0/3.0*pi*rhoi),1.0/3.0)
            mshell = rhoi * 4.0*pi*ri*ri*dr
            nshell = int(mshell / self.m0)
            if (nshell<12):
                nshell = 12
            res = int((1+sqrt(1+4*(nshell-12)/40))/2)
            nc = (40*res*(res-1)+12)
            print "@r = %g resolution = %d" %(ri,res)
            print "nominal shell node count = %d" % nc
            print "resulting in m0 scaling by %2.2f" % (mshell/(nc*self.m0))
            mi = self.m0 * (mshell/(nc*self.m0))
            hi = nNodePerh*(dr)
            Hi = SymTensor3d(1.0/hi, 0.0, 0.0,
                             0.0, 1.0/hi, 0.0,
                             0.0, 0.0, 1.0/hi)
                
            for n in xrange(nc):
                v1 = self.pixel2vector(n,res)
                pix = self.vector2pixel(v1,res)
                #print "[%g %g %g]" % (v1[0],v1[1],v1[2])
                self.x.append(ri*v1[0])
                self.y.append(ri*v1[1])
                self.z.append(ri*v1[2])
                self.m.append(mi)
                self.H.append(Hi)
            
            ri = max(0.0, ri - dr)
                
        print "Generated a total of %i nodes." % len(self.x)
        NodeGeneratorBase.__init__(self, True,
                                   self.x, self.y, self.z, self.m, self.H)
        return

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localPosition(self, i):
        assert i >= 0 and i < len(self.x)
        assert len(self.x) == len(self.y) == len(self.z)
        return Vector3d(self.x[i], self.y[i], self.z[i])

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
        return self.densityProfileMethod(self.localPosition(i).magnitude())

    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        assert i >= 0 and i < len(self.H)
        return self.H[i]

    #---------------------------------------------------------------------------
    # Numerically integrate the given density profile to determine the total
    # enclosed mass.
    #---------------------------------------------------------------------------
    def integrateTotalMass(self, densityProfileMethod,
                           rmin, rmax,
                           thetaMin, thetaMax,
                           phiMin, phiMax,
                           nbins = 10000):
        assert nbins > 0
        assert nbins % 2 == 0
        
        result = 0
        dr = (rmax-rmin)/nbins
        for i in xrange(1,nbins):
            r1 = rmin + (i-1)*dr
            r2 = rmin + i*dr
            result += 0.5*dr*(r2*r2*densityProfileMethod(r2)+r1*r1*densityProfileMethod(r1))
        result = result * (phiMax-phiMin) * (cos(thetaMin)-cos(thetaMax))
        return result

    

    def compute_matrices(self):
        self.A = [[0 for x in xrange(3)] for x in xrange(3)]
        self.B = [[0 for x in xrange(3)] for x in xrange(3)]
        self.C = [[0 for x in xrange(3)] for x in xrange(3)]
        self.D = [[0 for x in xrange(3)] for x in xrange(3)]
        self.E = [[0 for x in xrange(3)] for x in xrange(3)]

        x = 2.0/5.0*pi
        cs = cos(x)
        sn = sin(x)
        self.A[0][0] = cs
        self.A[0][1] = -sn
        self.A[1][0] = sn
        self.A[1][1] = cs
        self.A[2][2] = 1.0
        # A rotates 72 degrees around the z-axis
        x = pi/5.0
        ct = cos(x)/sin(x)
        cs = ct/sqrt(3.0)
        sn = sqrt(1.0-ct*ct/3.0)
        self.C[0][0] = 1.0
        self.C[1][1] = cs
        self.C[1][2] = -sn
        self.C[2][1] = sn
        self.C[2][2] = cs
        # C rotates around the x-axis so the north pole
        # ends up at the center of face 1
        cs = -0.5
        sn = sqrt(3.0)/2.0
        self.D[0][0] = cs
        self.D[0][1] = -sn
        self.D[1][0] = sn
        self.D[1][1] = cs
        self.D[2][2] = 1.0
        # D rotates 120 degrees around z-axis
        self.E = self.matmul1(self.C,self.D)
        self.B = self.matmul2(self.E,self.C)
        # B rotates faces 1 by 120 degrees
        for i in xrange(3):
            for j in xrange(3):
                self.E[i][j] = 0.0
            self.E[i][i] = 1.0
        # Now E is the identity matrix
        self.putMatrix(0,self.E)
        self.E = self.matmul1(self.B,self.A)
        self.E = self.matmul1(self.B,self.E)
        self.putMatrix(5,self.E)
        self.E = self.matmul1(self.E,self.A)
        self.putMatrix(10,self.E)
        self.E = self.matmul1(self.E,self.B)
        self.E = self.matmul1(self.E,self.B)
        self.E = self.matmul1(self.E,self.A)
        self.putMatrix(15,self.E)
        for n in xrange(4):
            nn = n*5
            self.E = self.getMatrix(nn)
            for i in xrange(1,5):
                self.E = self.matmul1(self.A,self.E)
                self.putMatrix(nn+i,self.E)
            
        for n in xrange(20):
            self.E = self.getMatrix(n)
            self.E = self.matmul1(self.E,self.C)
            self.putMatrix(n,self.E)
        return
            
    def compute_corners(self):
        dph = 2.0/5.0*pi
        # First corner is the north pole
        self.v[0][0] = 0
        self.v[0][1] = 0
        self.v[0][2] = 1
        # Next 5 lie on a circle with one on the y-axis
        z = 0.447213595 # This is 1/(2 sin^2(pi/5)) - 1
        rho = sqrt(1.0-z*z)
        for i in xrange(5):
            self.v[1+i][0] = -rho*sin(i*dph)
            self.v[1+i][1] = rho*cos(i*dph)
            self.v[1+i][2] = z
        # Second half are opposite of first half
        for i in xrange(6):
            for j in xrange(3):
                self.v[6+i][j] = -self.v[i][j]
        return

    def vector2pixel(self,v1,res):
        #assert res > 0
        pixperface = 2*res*(res-1)
        face    = self.find_face(v1)
        A       = self.getMatrix(face)
        v2      = self.vecmatmul2(A,v1)
        x       = v2[0]/v2[2]
        y       = v2[1]/v2[2]
        x,y     = self.adjust(x,y)
        ifail,pix   = self.tangentPlanePixel(res,x,y)
        if (ifail > 0):
            # try the runner-up face
            face    = self.find_another_face(v1,face)
            A       = self.getMatrix(face)
            v2      = self.vecmatmul2(A,v1)
            x       = v2[0]/v2[2]
            y       = v2[1]/v2[2]
            x,y     = self.adjust(x,y)
            ifail,pix  = self.tangentPlanePixel(res,x,y)
        pixel = face*pixperface + pix
        if (ifail > 0):
            # the pixel wasn't either of those, so must be a corner
            pix = self.find_corner(v1)
            pixel = 20*pixperface + pix
        return pixel
            
    def pixel2vector(self,pixel,res):
        # returns a unit vector pointing to pixel
        v1 = []
        pixperface = 2*res*(res-1)
        if (pixperface>0):
            face = pixel/pixperface
            if (face>20):
                face = 20
        else:
            print "it's all corners!"
            face = 20
        pix = pixel - face*pixperface
        print "passed in %d -> %d" % (pixel,pix)
        if (face<20):
            # This pixel is on a face
            x,y     = self.tangentPlaneVector(res,pix)
            x,y     = self.unadjust(x,y)
            norm    = sqrt(x*x + y*y + 1.0)
            v1.append(x/norm)
            v1.append(y/norm)
            v1.append(1.0/norm)
            A       = self.getMatrix(face)
            v1      = self.vecmatmul1(A,v1)
        else:
            # This is a corner pixel
            v1.append(self.v[pix][0])
            v1.append(self.v[pix][1])
            v1.append(self.v[pix][2])
        return v1
    #---------------------------------------------------------------------------
    # Lots of Linear Algebra follows. You have been warned!!!
    #---------------------------------------------------------------------------
    def copymatrix(self,M1,M2):
        for i in xrange(3):
            for j in xrange(3):
                M2[i][j] = M1[i][j]
        return
                
    def copyvector(self,v1,v2):
        for i in xrange(3):
            v2[i] = v1[i]
        return

    def getMatrix(self,n):
        M1 = [[0 for x in xrange(3)] for x in xrange(3)]
        for i in xrange(3):
            for j in xrange(3):
                M1[i][j] = self.R[n][i][j]
        return M1

    def putMatrix(self,n,M1):
        for i in xrange(3):
            for j in xrange(3):
                self.R[n][i][j] = M1[i][j]
        return

    def matmul1(self,M1,M2):
        # M3 = M1*M2
        M4 = [[0 for x in xrange(3)] for x in xrange(3)]
        sum = 0
        for i in xrange(3):
            for j in xrange(3):
                sum = 0
                for k in xrange(3):
                    sum = sum + M1[i][k] * M2[k][j]
                M4[i][j] = sum
        return M4

    def matmul2(self,M1,M2):
        # M3 = M1*M2^t
        M4 = [[0 for x in xrange(3)] for x in xrange(3)]
        sum = 0
        for i in xrange(3):
            for j in xrange(3):
                sum = 0
                for k in xrange(3):
                    sum = sum + M1[i][k] * M2[j][k]
                M4[i][j] = sum
        return M4

    def matmul3(self,M1,M2):
        # M3 = M1^t *M2
        M4 = [[0 for x in xrange(3)] for x in xrange(3)]
        sum = 0
        for i in xrange(3):
            for j in xrange(3):
                sum = 0
                for k in xrange(3):
                    sum = sum + M1[k][i] * M2[k][j]
                M4[i][j] = sum
        return M4

    def vecmatmul1(self,M1,v1):
        # v2 = M1*v1
        sum = 0
        v2 = []
        for i in xrange(3):
            sum = 0
            for j in xrange(3):
                sum = sum + M1[i][j]*v1[j]
            v2.append(sum)
        return v2

    def vecmatmul2(self,M1,v1):
        # v2 = M1^t * v1
        sum = 0
        v2 = []
        for i in xrange(3):
            sum = 0
            for j in xrange(3):
                sum = sum + M1[j][i]*v1[j]
            v2.append(sum)
        return v2

    def find_face(self,v1):
        # Locates the face to which the vector v1 points
        max,face = 0,0
        for n in xrange(20):
            dot = 0
            for i in xrange(3):
                dot = dot + self.R[n][i][2] * v1[i]
            if (dot > max):
                max = dot
                face = n
        return face

    def find_another_face(self,v1,face):
        # Computes the dot product with the vectors
        # pointing to the center of each face and picks the
        # largest one other than face.
        max = -17
        facetoavoid = face
        for n in xrange(20):
            if (n!=facetoavoid):
                dot = 0
                for i in xrange(3):
                    dot = dot + self.R[n][i][2] * v1[i]
                if (dot>max):
                    face = n
                    max = dot
        return face

    def find_corner(self,v1):
    	# Locates the corner to which v1 points.
        # Computes the dot product with the vectors
        # pointing to each corner and picks the
        # largest one.
        max = -17
        corner = 0
        for n in xrange(12):
            dot = 0
            for i in xrange(3):
                dot = dot + self.v[n][i] * v1[i]
            if (dot > max):
                corner = n
                max = dot
        return corner

    def find_mn(self,pixel,res):
        # Computes the integer coordinates (m,n) of the pixel
        # numbered pix on the basic triangle.
        found = 0
        pix, m, n = 0,0,0
        interiorpix, pixperedge = 0,0
        pix 	    = pixel
        interiorpix = (2*res-3)*(res-1)
        pixperedge  = (res)-1
        if (pix < interiorpix) and (found == 0):
            #pixel %d lies in the interior of the triangle." % pixel
            m = (sqrt(1.+8.*pix)-1.)/2. + 0.5/res
            # 0.5/resolution was added to avoid problems with
            # rounding errors for the case when n=0.
            # As long as you don't add more than 2/m, you're OK.
            n = pix - m*(m+1)/2
            m = m + 2
            n = n + 1
            found = 1
        pix = pix - interiorpix
        if (pix < pixperedge) and (found == 0):
            #pixel %d lies on the bottom edge." % pixel
            m = 2*res-1
            n = pix+1
            found = 1
        pix = pix - pixperedge
        if (pix < pixperedge) and (found == 0):
            #pixel %d lies on the right edge." % pixel
            m = 2*res-(pix+2)
            n = m
            found = 1
        pix = pix - pixperedge
        #pixel %d lies on the left edge." % pixel
        if (found == 0):
            m = pix+1
            n = 0
        return m,n

    def tangentPlanePixel(self,res,x,y):
        # Finds the hexagon in which the point (x,y) lies
        # and computes the corresponding pixel number pix.
        # if (x,y) lies on face, ifail = 0
        c=0.866025404
        edgelength=1.3231690765
        # The edge length of the icosahedron is
        # sqrt(9 tan^2(pi/5) - 3) when scaled so that
        # it circumscribes the unit sphere.
        r2	= 2*res
        a 	= 0.5*x
        b 	= c*y
        d 	= 0.5*edgelength/(r2-1)
        i 	= x/d 	  + r2
        j 	= (a+b)/d + r2
        k 	= (a-b)/d + r2
        m 	= (r2+r2-j+k-1)/3
        n 	= (i+k+1-r2)/3
        pix = (m-2)*(m-1)/2 + (n-1)
        ifail = 0
        correct = 0
        if (m  == (r2-1)) and (correct == 0):
            # On bottom row
            if ((n < 0) or (n > res)):
                ifail=1
            correct = 1
        if (n == m) and (correct == 0):
            # On right edge
            k = (r2-1) - m
            if ((k < 0) or (k> res)):
                ifail = 1
            else:
                pix = (r2-2)*(resolution-1) + k - 1
            correct = 1
        if (n == 0) and (correct == 0):
            # On left edge
            if ((m < 0) or (m > res)):
                ifail = 1
            else:
                pix = (r2-1)*(res-1) + m - 1
        return ifail, pix

    def tangentPlaneVector(self,res,pix):
        # Computes the coords(x,y) of pixel(pix)
        c1 = 0.577350269
        c2 = 0.866025404
        edgelength=1.3231690765
        m,n = self.find_mn(pix,res)
        x	= edgelength*(n-0.5*m)/(2*res-1)
        y 	= edgelength*(c1-(c2/(2*res-1))*m)
        return x,y

    def find_sixth(self,x,y):
        # Find which sixth of the basic triangle the point (x,y)
        # lies in.
        c=1.73205081
        d = c*y
        rot, flip = 0, 0
        if (x>0):
            if (x<-d):
                rot  = 0
                flip = 0
            else:
                if (x>d):
                    rot  = 2
                    flip = 1
                else:
                    rot  = 2
                    flip = 0
        else:
            if (x>-d):
                rot  = 1
                flip = 1
            else:
                if (x<d):
                    rot  = 1
                    flip = 0
                else:
                    rot  = 0
                    flip = 1
        return rot,flip

    def rotate_and_flip(self,rot,flip,x,y):
        cs = -0.5
        c  = 0.866025404
        sn = 0
        if (rot>0):
            if (rot == 1):
                # rotate by 120 degrees counter-clockwise
                sn = c
            else:
                sn = -c
            x1 = x
            x  = cs*x1 - sn*y
            y  = sn*x1 + cs*y
        if (flip>0):
            x = -x
        return x,y

    def adjust(self,x,y):
        # Maps the triangle onto the sphere (sort of)
        rot,flip    = self.find_sixth(x,y)
        x,y         = self.rotate_and_flip(rot,flip,x,y)
        x,y         = self.adjust_sixth(x,y)
        # Now rotate and flip the sixth back to original position
        if ((flip ==0) and (rot > 0)):
            x,y = self.rotate_and_flip(3-rot,flip,x,y)
        else:
            x,y = self.rotate_and_flip(rot,flip,x,y)
        return x,y

    def unadjust(self,x,y):
        # Does the inverse of adjust
        rot,flip    = self.find_sixth(x,y)
        x,y         = self.rotate_and_flip(rot,flip,x,y)
        x,y         = self.unadjust_sixth(x,y)
        # Now rotate and flip the sixth back to original position
        if ((flip ==0) and (rot > 0)):
            x,y = self.rotate_and_flip(3-rot,flip,x,y)
        else:
            x,y = self.rotate_and_flip(rot,flip,x,y)
        return x,y

    def adjust_sixth(self,x,y):
        # Maps the basic right triangle (the sixth of the face that
        # is in the lower right corner) onto itself in such a way
        # that pixels will have equal area when mapped onto the sphere.
        eps     = 1.e-14
        scale   = 1.09844
        g       = 1.7320508075689
        u 		= x  + eps
        v		= -y + eps
        v2		= v*v
        root	= sqrt(1.+4.*v2)
        trig	= atan((g*root-g)/(root+3.))
        y		= sqrt(trig*2./g)
        x		= sqrt((1.+4.*v2)/(1.+u*u+v2))*u*y/v
        x		= scale*x
        y 		= -scale*y	
        return x,y

    def unadjust_sixth(self,x,y):
        # Performs the inverse of what adjust_sixth does.
        eps     = 1.e-14
        scale   = 1.09844
        g       = 1.7320508075689
        u 		=  x/scale + eps
        v		= -y/scale + eps
        v2		= v*v
        trig	= tan(g*v2/2.)
        tmp		= (g+3.*trig)/(g-trig)
        y2		= (tmp*tmp-1.)/4.
        y		= sqrt(y2)
        tmp		= v2*(1.+4.*y2) - u*u*y2
        x	 	= u*y*sqrt((1.+y2)/tmp)
        y 		= -y
        return x,y

