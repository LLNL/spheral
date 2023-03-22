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
                 offset = None,
                 nNodePerh = 2.01,
                 SPH = False,
                 rejecter = None):

        # Check that the input parameters are good.
        assert n1 > 0
        assert n2 > 0
        assert n3 > 0 or distributionType == "cylindrical"
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
        assert offset is None or len(offset) == 3

        # If the user provided a constant for rho, then use the constantRho
        # class to provide this value.
        if type(rho) == type(1.0):
            self.rhofunc = ConstantRho(rho)
        else:
            self.rhofunc = rho

        self.n1 = n1
        self.n2 = n2
        self.n3 = n3
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
                                                      self.rhofunc,
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
                                                         self.rhofunc,
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
                                                         self.rhofunc,
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
                                         self.rhofunc,
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
                                 self.rhofunc,
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
                                    self.rhofunc,
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
                                      self.rhofunc,
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
                                     self.rhofunc,
                                     self.xmin,    # (xmin, ymin, zmin)
                                     self.xmax,    # (xmax, ymax, zmax)
                                     self.rmin,
                                     self.rmax,
                                     self.origin,
                                     self.nNodePerh)

        # If requested, shift the nodes.
        if offset:
            for i in range(len(self.x)):
                self.x[i] += offset[0]
                self.y[i] += offset[1]
                self.z[i] += offset[2]

        # If the user provided a "rejecter", give it a pass
        # at the nodes.
        if rejecter:
            self.x, self.y, self.z, self.m, self.H = rejecter(self.x,
                                                              self.y,
                                                              self.z,
                                                              self.m,
                                                              self.H)

        # Make rho a list
        self.rho = [self.rhofunc(Vector3d(self.x[i], self.y[i], self.z[i])) for i in range(len(self.m))]

        # Initialize the base class.  If "serialInitialization" is True, this
        # is where the points are broken up between processors as well.
        NodeGeneratorBase.__init__(self, serialInitialization,
                                   self.x, self.y, self.z, self.m, self.rho, self.H)

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
    def optimalSphericalDistribution(self, nRadial, nTheta, rhofunc,
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

        print('optimalSphericalDistribution: Cutoff radius is ', r0)
        nRadial0 = max(0, min(nRadial, int((r0 - rmin)/dr)))
        nRadial1 = nRadial - nRadial0
        r0 = rmin + nRadial0*dr
        print('Shifted to ', r0)
        print(nRadial0, nRadial1)

        if nRadial0 and not nRadial1:
            # Only use constant spacing for the nodes.
            x, y, m, H = self.constantDThetaSphericalDistribution(nRadial0, rhofunc,
                                                                  rmin,
                                                                  rmax,
                                                                  nNodePerh,
                                                                  theta)

        elif not nRadial0 and nRadial1:
            # Only use constant nodes per radial bin
            x, y, m, H = self.constantNThetaSphericalDistribution(nRadial1,
                                                                  nTheta,
                                                                  rhofunc,
                                                                  rmin,
                                                                  rmax,
                                                                  nNodePerh,
                                                                  theta)

        else:
            # Combine the two schemes.
            x0, y0, m0, H0 = self.constantDThetaSphericalDistribution(nRadial0, rhofunc,
                                                                      rmin,
                                                                      r0,
                                                                      nNodePerh,
                                                                      theta)
            x1, y1, m1, H1 = self.constantNThetaSphericalDistribution(nRadial1, nTheta, rhofunc,
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
    def constantMassShells(self, nNodes, nRadialShells, rhofunc,
                           rmin, rmax,
                           thetamin, thetamax,
                           phimin, phimax,
                           nNodePerh):

        assert nNodes > 0
        assert nRadialShells > 0 and nRadialShells < nNodes
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
        assert isinstance(rhofunc, ConstantRho)
        rho = rhofunc(Vector3d.zero)
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
        for ishell in range(nRadial):

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
        for i in range(n):
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
            print('Iteration %i for shell @ r=%f' % (iter, r))
            maxDisplacement = 0.0
            iter += 1

            # Loop over each point in the shell.
            for i in range(n):
                xi = r*theta[i]
                yi = r*phi[i]

                # Loop over every other point, and push it away if it's
                # within the mean seperation from this point.
                for j in range(n):
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
        for i in range(n):
            x.append(r*sin(phi[i])*cos(theta[i]))
            y.append(r*sin(phi[i])*sin(theta[i]))
            z.append(r*cos(phi[i]))

        return x, y, z

    #-------------------------------------------------------------------------------
    # Seed positions/masses for circular symmetry.
    # This version seeds a constant number of nodes per radial bin.
    #-------------------------------------------------------------------------------
    def constantNThetaSphericalDistribution(self, nRadial, nTheta, rhofunc,
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

        assert isinstance(rhofunc, ConstantRho)
        rho = rhofunc(Vector3d.zero)

        for i in range(0, nRadial):
            rInner = rmin + i*dr
            rOuter = rmin + (i + 1)*dr
            ri = rmin + (i + 0.5)*dr

            mRing = (rOuter**2 - rInner**2) * theta/2.0 * rho
            mi = mRing/nTheta

            hTheta = 1.0/(nNodePerh*ri*dTheta)

            for j in range(nTheta):
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
    def latticeDistribution(self, nx, ny, nz, rhofunc,
                            xmin,
                            xmax,
                            rmin,
                            rmax,
                            origin,
                            nNodePerh = 2.01):

        assert nx > 0
        assert ny > 0
        assert nz > 0

        dx = (xmax[0] - xmin[0])/nx
        dy = (xmax[1] - xmin[1])/ny
        dz = (xmax[2] - xmin[2])/nz

        n = nx*ny*nz
        volume = ((xmax[0] - xmin[0])*
                  (xmax[1] - xmin[1])*
                  (xmax[2] - xmin[2]))

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
        for iglobal in range(imin, imax):
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
                m.append(dx*dy*dz * rhofunc(Vector3d(xx, yy, zz)))
                H.append(H0)

        return x, y, z, m, H

    #-------------------------------------------------------------------------------
    # Seed positions/masses on a line
    #-------------------------------------------------------------------------------
    def lineDistribution(self, nx, ny, nz, rhofunc,
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

        dx = (xmax[0] - xmin[0])/nx
        dy = (xmax[1] - xmin[1])/ny
        dz = (xmax[2] - xmin[2])/nz
        ns = max([nx, ny, nz])
        ds = min([dx, dy, dz])

        assert isinstance(rhofunc, ConstantRho)
        rho = rhofunc(Vector3d.zero)

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

        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    xx = xmin[0] + (i + 0.5)*dx
                    yy = xmin[1] + (j + 0.5)*dy
                    zz = xmin[2] + (k + 0.5)*dz
                    r = sqrt(xx*xx + yy*yy + zz*zz)
                    if ((rmin is None or r >= rmin) and
                        (rmax is None or r <= rmax)):
                        x.append(xx)
                        y.append(yy)
                        z.append(zz)
                        m.append(m0)
                        H.append(H0)

        return x, y, z, m, H

    #-------------------------------------------------------------------------------
    # Seed positions/masses on a cylindrical section about the z-axis.
    #-------------------------------------------------------------------------------
    def cylindrical(self, nr, nz, rhofunc,
                    rmin,
                    rmax,
                    thetamin,
                    thetamax,
                    zmin,
                    zmax,
                    nNodePerh = 2.01):

        assert nr > 0
        assert nz > 0
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

        assert isinstance(rhofunc, ConstantRho)
        rho = rhofunc(Vector3d.zero)

        for iz in range(nz):
            zi = zmin + (iz + 0.5)*dz
            for ir in range(nr):
                rInner = rmin + ir*dr
                rOuter = rInner + dr
                ri = 0.5*(rInner + rOuter)
                circumference = ri*Dtheta
                nominalNtheta = int(circumference/max(dr, dz))
                ntheta = max(minNtheta, nominalNtheta)
                dtheta = Dtheta/ntheta
                mring = 0.5*(rOuter**2 - rInner**2)*Dtheta*dz*rho
                mi = mring/ntheta
                htheta = 1.0/(nNodePerh*dtheta*ri)
                for itheta in range(ntheta):
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
    def constantNTheta(self, nr, ntheta, nz, rhofunc,
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
        assert rmin < rmax
        assert zmin < zmax

        dr = (rmax - rmin)/nr
        Dtheta = thetamax - thetamin
        dtheta = Dtheta/ntheta
        dz = (zmax - zmin)/nz

        hr = 1.0/(nNodePerh*dr)
        hz = 1.0/(nNodePerh*dz)

        assert isinstance(rhofunc, ConstantRho)
        rho = rhofunc(Vector3d.zero)

        x = []
        y = []
        z = []
        m = []
        H = []

        for iz in range(nz):
            zi = zmin + (iz + 0.5)*dz
            for ir in range(nr):
                rInner = rmin + ir*dr
                rOuter = rInner + dr
                ri = 0.5*(rInner + rOuter)
                circumference = ri*Dtheta
                mring = 0.5*(rOuter**2 - rInner**2)*Dtheta*dz*rho
                mi = mring/ntheta
                htheta = 1.0/(nNodePerh*dtheta*ri)
                for itheta in range(ntheta):
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
    def hcpDistribution(self, nx, ny, nz, rhofunc,
                        xmin,
                        xmax,
                        rmin,
                        rmax,
                        origin,
                        nNodePerh = 2.01):

        assert nx > 0
        assert ny > 0
        assert nz > 0

        dx = (xmax[0] - xmin[0])/nx
        dy = (xmax[1] - xmin[1])/ny
        dz = (xmax[2] - xmin[2])/nz

        n = nx*ny*nz
        volume = ((xmax[0] - xmin[0])*
                  (xmax[1] - xmin[1])*
                  (xmax[2] - xmin[2]))

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
        for iglobal in range(imin, imax):
            i = iglobal % nx
            j = (iglobal // nx) % ny
            k = iglobal // (nx*ny)
            xx = xmin[0] + (i + 0.5*((j % 2) + (k % 2)))*dx
            yy = xmin[1] + (j + 0.5*(k % 2))*dy
            zz = xmin[2] + (k + 0.5)*dz

            r2 = (xx - origin[0])**2 + (yy - origin[1])**2 + (zz - origin[2])**2
            if ((rmin is None or r2 >= rmin**2) and
                (rmax is None or r2 <= rmax**2)):
                x.append(xx)
                y.append(yy)
                z.append(zz)
                m.append(dx*dy*dz*rhofunc(Vector3d(xx, yy, zz)))
                H.append(H0)

        return x, y, z, m, H

#--------------------------------------------------------------------------------
# Create cylindrical distributions similar to the RZ generator.
#--------------------------------------------------------------------------------
from GenerateNodeDistribution2d import GenerateNodeDistribution2d, RZGenerator
from Spheral import generateCylDistributionFromRZ

class CylindricalSpunGenerator3d(NodeGeneratorBase):

    def __init__(self,
                 gen2d,
                 rho,
                 nNodePerh = 2.01,
                 SPH = False,
                 rejecter = None,
                 phi = 2.0*pi):

        from Spheral import Vector3d, CylindricalBoundary

        genRZ = RZGenerator(gen2d)

        # If the user provided a constant for rho, then use the constantRho
        # class to provide this value.
        if type(rho) == type(1.0):
            self.rho = ConstantRho(rho)
        else:
            self.rho = rho

        # The 2D generator already split the nodes up between processors, but
        # we want to handle that ourselves.  Distribute the full set of RZ
        # nodes to every process, then redecompose them below.
        self.x = mpi.allreduce(genRZ.x[:], mpi.SUM)
        self.y = mpi.allreduce(genRZ.y[:], mpi.SUM)
        self.m = mpi.allreduce(genRZ.m[:], mpi.SUM)
        self.H = mpi.allreduce(genRZ.H[:], mpi.SUM)
        n = len(self.x)
        self.z = [0.0]*n
        self.globalIDs = [0]*n

        # Convert the 2-D H tensors to 3-D, and correct the masses.
        for i in range(n):
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

        assert len(self.m) == n
        assert len(self.H) == n

        # Duplicate the nodes from the xy-plane, creating rings of nodes about
        # the x-axis.  We use a C++ helper method for the sake of speed.
        kernelExtent = 2.0
        xvec = self.vectorFromList(self.x, vector_of_double)
        yvec = self.vectorFromList(self.y, vector_of_double)
        zvec = self.vectorFromList(self.z, vector_of_double)
        mvec = self.vectorFromList(self.m, vector_of_double)
        Hvec = self.vectorFromList(self.H, vector_of_SymTensor3d)
        globalIDsvec = self.vectorFromList(self.globalIDs, vector_of_int)
        extrasVec = vector_of_vector_of_double()
        generateCylDistributionFromRZ(xvec, yvec, zvec, mvec, Hvec, globalIDsvec,
                                      extrasVec,
                                      nNodePerh, kernelExtent, phi,
                                      procID, nProcs)

        # Allow some 3D rejecter logic.
        if rejecter:
            self.x, self.y, self.z, self.m, self.H, self.globalIDs = [], [], [], [], [], []
            for i in range(len(xvec)):
                if rejecter.accept(xvec[i],yvec[i],zvec[i]):
                    self.x.append(xvec[i])
                    self.y.append(yvec[i])
                    self.z.append(zvec[i])
                    self.m.append(mvec[i])
                    self.H.append(SymTensor3d(Hvec[i]))
                    self.globalIDs.append(globalIDsvec[i])

        else:
            self.x = [x for x in xvec]
            self.y = [x for x in yvec]
            self.z = [x for x in zvec]
            self.m = [x for x in mvec]
            self.H = [SymTensor3d(x) for x in Hvec]
            self.globalIDs = [x for x in globalIDsvec]

        # Initialize the base.
        NodeGeneratorBase.__init__(self, False,
                                   self.x, self.y, self.z, self.m, self.H)
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
# Backwards compatible version of our 3D spun generator.
#-------------------------------------------------------------------------------
class GenerateCylindricalNodeDistribution3d(CylindricalSpunGenerator3d):

    def __init__(self, nRadial, nTheta, rho,
                 distributionType = 'optimal',
                 xmin = None,
                 xmax = None,
                 rmin = None,
                 rmax = None,
                 nNodePerh = 2.01,
                 theta = pi/2.0,
                 azimuthalOffsetFraction = 0.0,
                 SPH = False,
                 rotation = 0.0,
                 offset = None,
                 xminreject = None,
                 xmaxreject = None,
                 rreject = None,
                 originreject = None,
                 reversereject = False,
                 relaxation = None,
                 rejecter = None,
                 rejecter3d = None,
                 phi = 2.0*pi):
        gen2d = GenerateNodeDistribution2d(nRadial = nRadial,
                                           nTheta = nTheta,
                                           rho = rho,
                                           distributionType = distributionType,
                                           xmin = xmin,
                                           xmax = xmax,
                                           rmin = rmin,
                                           rmax = rmax,
                                           nNodePerh = nNodePerh,
                                           theta = theta,
                                           azimuthalOffsetFraction = azimuthalOffsetFraction,
                                           SPH = SPH,
                                           rotation = rotation,
                                           offset = offset,
                                           xminreject = xminreject,
                                           xmaxreject = xmaxreject,
                                           rreject = rreject,
                                           originreject = originreject,
                                           reversereject = reversereject,
                                           relaxation = relaxation,
                                           rejecter = rejecter)
        CylindricalSpunGenerator3d.__init__(self,
                                            gen2d = gen2d,
                                            rho = rho,
                                            nNodePerh = nNodePerh,
                                            SPH = SPH,
                                            rejecter = rejecter3d,
                                            phi = phi)
        return

#-------------------------------------------------------------------------------
# Specialized version that generates a sphere matching a radial profile using
# the monte carlo method
#-------------------------------------------------------------------------------
class GenerateMonteCarloProfile3d(NodeGeneratorBase):
    
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
                 nNodePerh = 2.01,
                 m0 = None):
        
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
        
        self.totalMass = self.integrateTotalMass(self.densityProfileMethod,
                                                 rmin, rmax,
                                                 0.0, pi,
                                                 0.0, 2.0*pi)
        
        if m0 is None:
            # Now set the nominal mass per node.
            self.m0 = self.totalMass/n # In this case, we're setting the total number of particles
        else:
            self.m0 = m0
            n = self.totalMass/m0
        assert self.m0 > 0.0
        print("Nominal mass per node of %g for %d nodes" % (self.m0,n))
        
        # Get the normalization constant for radial probability distribution
        norm = self.integrateDensityFunction(self.densityProfileMethod,rmin,rmax)
        pMax = 0
        
        dr = (rmax-rmin)/10000.0
        for i in range(10000):
            ri = rmin + dr*i
            pMax = max(self.densityProfileMethod(ri)/norm,pMax)

        # OK, we now know enough to generate the node positions.
        from Spheral import SymTensor3d
        import random
        self.x = []
        self.y = []
        self.z = []
        self.m = []
        self.H = []

        # Below was an attempt at a sophisticated probability sampling, using accept/reject instead
        '''
        while (len(self.x)<n):
            random.seed()
            u       = random.random()
            v       = random.random()
            w       = random.random()
            theta   = 2.0*pi*u
            phi     = acos(2.0*v-1)
            r       = self.integrateDensityFunction(self.densityProfileMethod,rmin,rmax,w,norm)
            print "r = %f" % r
            rhoi    = self.densityProfileMethod(r)
            voli    = self.m0/rhoi
            hi      = 2.0*nNodePerh*pow(3.0/(4.0*pi)*voli,1.0/3.0)
            Hi      = SymTensor3d(1.0/hi, 0.0, 0.0,
                                  0.0, 1.0/hi, 0.0,
                                  0.0, 0.0, 1.0/hi)

            if ((theta >= thetaMin) and (theta <= thetaMax) and (phi <= phiMax) and (phi >= phiMin) and (r <= rmax) and (r >= rmin)):
                self.x.append((r)*cos(theta)*sin(phi))
                self.y.append((r)*sin(theta)*sin(phi))
                self.z.append((r)*cos(phi))
                self.m.append(self.m0)
                self.H.append(Hi)
        '''
        
        while (len(self.x)<n):
            random.seed()
            x       = random.uniform(-rmax,rmax)
            y       = random.uniform(-rmax,rmax)
            z       = random.uniform(-rmax,rmax)
            r       = sqrt(x*x+y*y+z*z)
            u       = random.random()*pMax
            
            if ((r<=rmax) and (r>=rmin)):
                if (self.densityProfileMethod(ri)/norm >= u):
                    # accept
                    rhoi    = self.densityProfileMethod(r)
                    voli    = self.m0/rhoi
                    hi      = 2.0*nNodePerh*pow(3.0/(4.0*pi)*voli,1.0/3.0)
                    Hi      = SymTensor3d(1.0/hi, 0.0, 0.0,
                                          0.0, 1.0/hi, 0.0,
                                          0.0, 0.0, 1.0/hi)
                    phi     = atan2(y,x)+pi
                    theta   = acos(z/r)
                    if ((theta >= thetaMin) and (theta <= thetaMax) and (phi <= phiMax) and (phi >= phiMin)):
                        self.x.append(x)
                        self.y.append(y)
                        self.z.append(z)
                        self.m.append(self.m0)
                        self.H.append(Hi)


        print("Generated a total of %i nodes." % len(self.x))

        # Make sure the total mass is what we intend it to be, by applying
        # a multiplier to the particle masses.
        sumMass = 0.0
        for m in self.m:
            sumMass += m
        assert sumMass > 0.0

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
                           nbins = 1000):
        assert nbins > 0
        assert nbins % 2 == 0
        
        result = 0
        dr = (rmax-rmin)/nbins
        for i in range(1,nbins):
            r1 = rmin + (i-1)*dr
            r2 = rmin + i*dr
            result += 0.5*dr*(r2*r2*densityProfileMethod(r2)+r1*r1*densityProfileMethod(r1))
        result = result * (phiMax-phiMin) * (cos(thetaMin)-cos(thetaMax))
        return result

    def integrateDensityFunction(self,densityProfileMethod,rmin,rmax,
                                 lam = -1.0,
                                 norm = 1.0,
                                 nbins = 1000):
        assert nbins > 0
        assert nbins % 2 == 0

        result = 0
        dr = (rmax-rmin)/nbins
        for i in range(1,nbins):
            r1 = rmin + (i-1)*dr
            r2 = rmin + i*dr
            result += 0.5*dr*norm*(densityProfileMethod(r2)+densityProfileMethod(r1))
            if (lam >=0 and result >= lam):
                result = r1
                break
        return result



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
                 nNodePerh = 2.01,
                 m0 = -1.0):
        
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
        print("Total mass of %g in the range r = (%g, %g), theta = (%g, %g), phi = (%g, %g)" % \
            (self.totalMass, rmin, rmax, thetaMin, thetaMax, phiMin, phiMax))

        # Now set the nominal mass per node.
        self.m0 = self.totalMass/(4.0/3.0*pi*pow(self.n,3))
        if m0 > 0:
            self.m0 = m0
        assert self.m0 > 0.0
        print("Nominal mass per node of %g." % self.m0)

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
            #print "nshell=%d mshell=%f m0=%f" % (nshell,mshell,self.m0)
            #mi      = self.m0 * (mshell/(nshell*self.m0))
            
            for n in range(nshell):
                random.seed()
                u       = random.random()
                v       = random.random()
                w       = random.random()
                theta   = 2.0*pi*u
                phi     = acos(2.0*v-1)
                #deltar  = w*(0.45*dr)
                deltar  = 0.0
                self.x.append((ri-deltar)*cos(theta)*sin(phi))
                self.y.append((ri-deltar)*sin(theta)*sin(phi))
                self.z.append((ri-deltar)*cos(phi))
                self.m.append(self.m0)
                self.H.append(Hi)
                    
            # Decrement to the next radial bin inward.
            ri = max(0.0, ri - dr)


        print("Generated a total of %i nodes." % len(self.x))

        # Make sure the total mass is what we intend it to be, by applying
        # a multiplier to the particle masses.
        sumMass = 0.0
        for m in self.m:
            sumMass += m
        assert sumMass > 0.0
        massCorrection = self.totalMass/sumMass
        for i in range(len(self.m)):
            self.m[i] *= massCorrection
        print("Applied a mass correction of %f to ensure total mass is %f." % (massCorrection, self.totalMass))

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
        for i in range(1,nbins):
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
class GenerateLongitudinalNodesMatchingProfile3d(NodeGeneratorBase):

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
        print("Total mass of %g in the range r = (%g, %g), theta = (%g, %g), phi = (%g, %g)" % \
            (self.totalMass, rmin, rmax, thetaMin, thetaMax, phiMin, phiMax))
                
        # Now set the nominal mass per node.
        self.m0 = self.totalMass/(4.0/3.0*pi*pow(self.n,3))
        assert self.m0 > 0.0
        print("Nominal mass per node of %g." % self.m0)
            
        # OK, we now know enough to generate the node positions.
        self.x, self.y, self.z, self.m, self.H = \
            self.constantDThetaDPhiDistribution(self.densityProfileMethod,
                                                self.m0, self.n,
                                                rmin, rmax,
                                                thetaMin, thetaMax,
                                                phiMin, phiMax,
                                                nNodePerh)
        print("Generated a total of %i nodes." % len(self.x))
            
        # Make sure the total mass is what we intend it to be, by applying
        # a multiplier to the particle masses.
        sumMass = 0.0
        for m in self.m:
            sumMass += m
        assert sumMass > 0.0
        massCorrection = self.totalMass/sumMass
        for i in range(len(self.m)):
            self.m[i] *= massCorrection
        print("Applied a mass correction of %f to ensure total mass is %f." % (massCorrection, self.totalMass))

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
        for i in range(1,nbins):
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
            for i in range(nTheta):
                thetai = thetaMin + (i + 0.5)*dTheta
                rp = ri*sin(thetai)
                arclength = arcPhi * rp
                nPhi = max(1,int(arclength/dPhi))
                dPhij = arcPhi/nPhi
                #print "at r=%g theta=%g rp=%g nTheta=%d nPhi=%d" % (ri,thetai,rp,nTheta,nPhi)
                for j in range(nPhi):
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
# in a disk described by a density profile in r and z.
#-------------------------------------------------------------------------------
class GenerateIdealDiskMatchingProfile3d(NodeGeneratorBase):
    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def __init__(self, n, densityProfileMethod,
                 rmin = 0.0,
                 rmax = 1.0,
                 zmax = 0.1,
                 nNodePerh = 2.01,
                 offset=None):
        
        assert n > 0
        assert rmin < rmax
        assert zmax >= 0
        assert nNodePerh > 0.0
        assert offset is None or len(offset)==3

        self.n = n
        self.rmin = rmin
        self.rmax = rmax
        self.zmax = zmax
        self.nNodePerh = nNodePerh
        self.densityProfileMethod = densityProfileMethod
    
        # Determine how much total mass there is in the system.
        self.totalMass = self.integrateTotalMass(self.densityProfileMethod,
                                                 rmin, rmax,
                                                 zmax)
        print("Total mass of %g in the range r = (%g, %g), z = (0, %g)" % \
            (self.totalMass, rmin, rmax, zmax))
        
        self.laminarMass = self.integrateLaminarMass(self.densityProfileMethod,rmin,rmax)

        self.m0 = self.laminarMass/(self.n*self.n*pi)
            
        # Return lists for positions, masses, and H's.
        self.x = []
        self.y = []
        self.z = []
        self.m = []
        self.H = []

        zi = 0
        while zi <= zmax:
            ri = rmax
            while ri > rmin:
                rhoi        = densityProfileMethod(ri,zi)
                dr          = sqrt(self.m0/rhoi)
                arclength   = 2.0*pi*ri
                arcmass     = arclength*dr*rhoi
                nTheta      = max(1,int(arcmass/self.m0))
                dTheta      = 2.0*pi/nTheta
                mi          = arcmass/nTheta
                hi          = nNodePerh*0.5*(dr + ri*dTheta)
                Hi          = SymTensor3d(1.0/hi,0.0,0.0,
                                          0.0,1.0/hi,0.0,
                                          0.0,0.0,1.0/hi)
                                          
                # Now assign the nodes for this radius.
                for i in range(nTheta):
                    thetai = (i + 0.5)*dTheta
                    self.x.append(ri*cos(thetai))
                    self.y.append(ri*sin(thetai))
                    self.z.append(zi)
                    self.m.append(mi)
                    self.H.append(Hi)
                    #now same for negative z
                    if (zi>0):
                        thetai = (i + 0.5)*dTheta
                        self.x.append(ri*cos(thetai))
                        self.y.append(ri*sin(thetai))
                        self.z.append(-zi)
                        self.m.append(mi)
                        self.H.append(Hi)
                #move inward
                ri          = max(0.0, ri - dr)
            rho0 = densityProfileMethod(rmin, zi)
            dz = sqrt(self.m0/rho0)
            zi += dz
            
        NodeGeneratorBase.__init__(self, True,
                                   self.x, self.y, self.z, self.m, self.H)

        return



    #---------------------------------------------------------------------------
    # Numerically integrate the given density profile to determine the total
    # enclosed mass.
    #---------------------------------------------------------------------------
    def integrateTotalMass(self, densityProfileMethod,
                           rmin, rmax,
                           zmax,
                           nrbins = 10000,
                           nzbins = 1000):
        assert nrbins > 0
        assert nrbins % 2 == 0
        assert nzbins > 0
        assert nzbins % 2 == 0
        
        result = 0
        dr = (rmax-rmin)/nrbins
        dz = zmax/nzbins
        I = []
        
        for i in range(1,nrbins):
            r1 = rmin + (i-1)*dr
            r2 = rmin + i*dr
            ij = 0
            for j in range(1,nzbins):
                z1 = (j-1)*dz
                z2 = j*dz
                ij += 0.5*r1*dz*(densityProfileMethod(r1,z1)+densityProfileMethod(r1,z2))
            I.append(ij)
            if (i>1):
                result += 0.5*dr*(I[i-1]+I[i-2])
    
        return result

    def integrateLaminarMass(self, densityProfileMethod,rmin,rmax,
                             nrbins = 10000):
        
        assert nrbins > 0
        assert nrbins % 2 == 0
        
        h = (rmax - rmin)/nrbins
        result = (rmin*densityProfileMethod(rmin,0) +
                  rmax*densityProfileMethod(rmax,0))
        for i in range(1, nrbins):
            ri = rmin + i*h
            if i % 2 == 0:
                result += 4.0*ri*densityProfileMethod(ri,0)
            else:
                result += 2.0*ri*densityProfileMethod(ri,0)

        result *= 2*pi*h/3.0
        return result

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
        return self.densityProfileMethod(sqrt(self.localPosition(i)[0]**2+self.localPosition(i)[1]**2),self.localPosition(i)[2])

    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        assert i >= 0 and i < len(self.H)
        return self.H[i]


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
                 nNodePerh = 2.01,
                 offset=None,
                 rejecter=None,
                 m0 = 0.0):
        
        print("This generator has been deprecated and will always fail.")
        print("Please import * from RecursivePrimitiveRefinementGenerator.py and use RPRPSGenerator()")
        return

class GenerateSpiralMatchingProfile3d(NodeGeneratorBase):
    
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
                 nNodePerh = 2.01,
                 offset=None,
                 rMaxForMassMatching=None):
        
        assert n > 0
        assert rmin < rmax
        assert thetaMin < thetaMax
        assert thetaMin >= 0.0 and thetaMin <= 2.0*pi
        assert thetaMax >= 0.0 and thetaMax <= 2.0*pi
        assert phiMin < phiMax
        assert phiMin >= 0.0 and phiMin <= 2.0*pi
        assert phiMax >= 0.0 and phiMax <= 2.0*pi
        assert nNodePerh > 0.0
        assert offset is None or len(offset)==3
        
        import random
        
        if offset is None:
            self.offset = Vector3d(0,0,0)
        else:
            self.offset = Vector3d(offset[0],offset[1],offset[2])
        
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
        if rMaxForMassMatching is None:
            self.totalMass = self.integrateTotalMass(self.densityProfileMethod,
                                                     rmin, rmax,
                                                     thetaMin, thetaMax,
                                                     phiMin, phiMax)
        else:
            self.totalMass = self.integrateTotalMass(self.densityProfileMethod,
                                                     0.0, rMaxForMassMatching,
                                                     thetaMin, thetaMax,
                                                     phiMin, phiMax)
        print("Total mass of %g in the range r = (%g, %g), theta = (%g, %g), phi = (%g, %g)" % \
            (self.totalMass, rmin, rmax, thetaMin, thetaMax, phiMin, phiMax))

        # Now set the nominal mass per node.
        self.m0 = self.totalMass/(4.0/3.0*pi*pow(self.n,3))
        assert self.m0 > 0.0
        print("Nominal mass per node of %g." % self.m0)

        from Spheral import SymTensor3d
        self.x = []
        self.y = []
        self.z = []
        self.m = []
        self.H = []
        ri = rmax
        
        while ri > rmin:
            # create the database of faces and positions
            self.positions      = []     # [index,[point]]
            self.middlePoints   = []  # [i,[key,index]]
            self.faces          = []
            self.index          = 0
            
            # Get the nominal delta r, number of nodes,
            # and mass per node at this radius.
            rhoi    = self.densityProfileMethod(ri)
            dr      = pow(self.m0/(rhoi),1.0/3.0)
            dr      = min(dr,ri-rmin)
            #mshell  = rhoi * 4.0*pi*ri*ri*dr
            mshell  = self.integrateTotalMass(self.densityProfileMethod,
                                              ri-dr, ri,
                                              0, pi,
                                              0, 2*pi)
            nshell  = int(mshell / self.m0+0.5)
            nshell  = max(nshell,1)
            nr      = 0
            ver     = 0
            counts  = []

            hi = nNodePerh*(dr)
            Hi = SymTensor3d(1.0/hi, 0.0, 0.0,
                           0.0, 1.0/hi, 0.0,
                           0.0, 0.0, 1.0/hi)
                           
            mi  = mshell / float(nshell)
            
            if (nshell > 1):
                p = 0
                for i in range(1,nshell+1):
                    h = -1.0+(2.0*(i-1.0)/(nshell-1.0))
                    t = acos(h)
                    
                    if (i>1 and i<nshell):
                        p = (p + 3.8/sqrt(nshell)*1.0/sqrt(1.0-h*h)) % (2.0*pi)
                    elif (i==nshell):
                        p = 0
                    
                    x = sin(t)*cos(p)
                    y = sin(t)*sin(p)
                    z = cos(t)

                    self.positions.append([x,y,z])
            elif (nshell==1 and mi > 0.5 * self.m0):
                self.positions.append([0,0,0])

            #mi = self.m0 * (float(nshell)/float(len(self.positions)))
            
            rii = ri - 0.5*dr
            print("at r=%g, wanted %d; computed %d total nodes with mass=%g" %(rii,nshell,len(self.positions),mi))
            for n in range(len(self.positions)):
                x       = rii*self.positions[n][0]
                y       = rii*self.positions[n][1]
                z       = rii*self.positions[n][2]
                
                
                random.seed(nshell)
                dt = random.random()*pi
                dt2 = random.random()*pi
                
                rot = [[1.0,0.0,0.0],[0.0,cos(dt),-sin(dt)],[0.0,sin(dt),cos(dt)]]
                rot2 = [[cos(dt2),0.0,sin(dt2)],[0.0,1.0,0.0],[-sin(dt2),0.0,cos(dt2)]]
                pos = [x,y,z]
                posp= [0,0,0]
                for k in range(3):
                    for j in range(3):
                        posp[k] += pos[j]*rot[k][j]
                
                x = posp[0]
                y = posp[1]
                z = posp[2]
                
                pos = [x,y,z]
                posp= [0,0,0]
                for k in range(3):
                    for j in range(3):
                        posp[k] += pos[j]*rot2[k][j]
                x = posp[0]
                y = posp[1]
                z = posp[2]
                
                if(nshell>1):
                    theta   = acos(z/sqrt(x*x+y*y+z*z))
                    phi     = atan2(y,x)
                    if (phi<0.0):
                        phi = phi + 2.0*pi
                else:
                    theta = (thetaMax - thetaMin)/2.0
                    phi = (phiMax - phiMin)/2.0
                if (theta<=thetaMax and theta>=thetaMin) and (phi<=phiMax and phi>=phiMin):
                    self.x.append(x)
                    self.y.append(y)
                    self.z.append(z)
                    self.m.append(mi)
                    self.H.append(SymTensor3d.one*(1.0/hi))
                
            ri = max(rmin, ri - dr)
                
        # If requested, shift the nodes.
        if offset:
            for i in range(len(self.x)):
                self.x[i] += offset[0]
                self.y[i] += offset[1]
                self.z[i] += offset[2]
        
        print("Generated a total of %i nodes." % len(self.x))
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
        loc = Vector3d(0,0,0)
        loc = self.localPosition(i) - self.offset
        return self.densityProfileMethod(loc.magnitude())

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
        for i in range(1,nbins):
            r1 = rmin + (i-1)*dr
            r2 = rmin + i*dr
            result += 0.5*dr*(r2*r2*densityProfileMethod(r2)+r1*r1*densityProfileMethod(r1))
        result = result * (phiMax-phiMin) * (cos(thetaMin)-cos(thetaMax))
        return result
