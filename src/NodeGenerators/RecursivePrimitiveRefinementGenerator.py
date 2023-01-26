from math import *

from NodeGeneratorBase import *

from Spheral import Vector3d
from Spheral import Tensor3d
from Spheral import SymTensor3d

from Spheral import vector_of_int, vector_of_double, vector_of_SymTensor3d, vector_of_vector_of_double

import mpi
rank = mpi.rank
procs = mpi.procs

class RPRPSGenerator3d(NodeGeneratorBase):
    def __init__(self,n,densityProfileMethod,
                 rmin = 0.0,
                 rmax = 0.0,
                 thetaMin = 0.0,
                 thetaMax = pi,
                 phiMin = 0.0,
                 phiMax = 2.0*pi,
                 nNodePerh = 2.01,
                 offset = None,
                 rejecter = None,
                 m0 = 0.0):

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
        
        self.rejecter = None
        if rejecter:
            self.rejecter = rejecter
        
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
        self.totalMass = self.integrateTotalMass(self.densityProfileMethod,
                                                     rmin, rmax,
                                                     thetaMin, thetaMax,
                                                     phiMin, phiMax)

        print("Total mass of %g in the range r = (%g, %g), theta = (%g, %g), phi = (%g, %g)" % \
            (self.totalMass, rmin, rmax, thetaMin, thetaMax, phiMin, phiMax))

        # Now set the nominal mass per node.
        if (m0 == 0.0):
            self.m0 = self.totalMass/n
        else:
            self.m0 = m0
            n = int(self.totalMass/self.m0)
        assert self.m0 > 0.0
        print("Nominal mass per node of %g for %d nodes." % (self.m0,n))
        
        from Spheral import SymTensor3d
        self.x = []
        self.y = []
        self.z = []
        self.m = []
        self.H = []
        ri = rmax

        # new formula for calculating number of points for a given subdivision level
        # (Nf * Np(n) - Ne * Npe(n) + Nc)
        # Nf = Number of faces of primitive shape
        # Np(n) = Number of points in a triangle subdivided n times
        #       2^(2n-1) + 3*2^(n-1) + 1
        # Ne = Number of edges of primitive shape
        # Npe(n) = Number of points along an edge of primitive shape subdivided n times
        #       2^n + 1
        # Nc = Number of corners
        
        # shapeData = [Nf,Ne,Nc]
        
        shapeData = [[ 6, 9, 5],
                     [ 8,12, 6],
                     [12,18, 8],
                     [20,30,12]]
        
        # first column is total number of shell points
        # second column is number of refinements to reach that shell count
        # third column is shape choice that reaches that shell count
        resolution = [[5,0,0],
                      [6,0,1],
                      [8,0,2],
                      [12,0,3],
                      [14,1,0],
                      [18,1,1],
                      [26,1,2],
                      [42,1,3],
                      [50,2,0],
                      [66,2,1],
                      [98,2,2],
                      [162,2,3],
                      [194,3,0],
                      [258,3,1],
                      [386,3,2],
                      [642,3,3],
                      [770,4,0],
                      [1026,4,1],
                      [1538,4,2],
                      [2562,4,3],
                      [3074,5,0],
                      [4098,5,1],
                      [6146,5,2],
                      [10242,5,3],
                      [12290,6,0],
                      [16386,6,1],
                      [24578,6,2],
                      [40962,6,3],
                      [49154,7,0],
                      [65538,7,1],
                      [98306,7,2],
                      [163842,7,3],
                      [196610,8,0],
                      [262146,8,1],
                      [393218,8,2]]
        
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
            #dr      = min(dr,ri-rmin)
            rii = abs(ri - 0.5*dr)
            # now compute a new dr based on rii
            # this should in theory protect against the half bin radius being
            # below rmin while not sacrificing the mass of the entire shell
            # with a simple if condition
            rhoi    = self.densityProfileMethod(rii)
            dr      = pow(self.m0/rhoi,1.0/3.0)
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
            
            random.seed(nshell)
            dt = random.random()*pi
            dt2 = random.random()*pi
            
            rot = [[1.0,0.0,0.0],[0.0,cos(dt),-sin(dt)],[0.0,sin(dt),cos(dt)]]
            rot2 = [[cos(dt2),0.0,sin(dt2)],[0.0,1.0,0.0],[-sin(dt2),0.0,cos(dt2)]]
            
            if (nshell > 4 and nshell<163):
                if (mpi.rank == 0):
                    for i in range(len(shapeData)):
                        nc  = 0
                        nco = 0
                        nrf = 0
                        while (nc < nshell):
                            nrf += 1
                            nco = nc
                            nc = self.shapeCount(nrf,shapeData[i])
                        counts.append([i,nrf-1,nco])
                        counts.append([i,nrf,nc])

                    diff = 1e13
                    for i in range(len(counts)):
                        dd = abs(counts[i][2] - nshell)
                        if (dd < diff):
                            diff = dd
                            ver = counts[i][0]
                            nr = counts[i][1]

                    if (nr<0):
                        nr = 0

                    if (ver==0):
                        self.createHexaSphere(nr)
                    elif (ver==1):
                        self.createOctaSphere(nr)
                    elif (ver==2):
                        self.createCubicSphere(nr)
                    else:
                        self.createIcoSphere(nr)

                    for n in range(len(self.positions)):
                        self.positions[n] = self.rotater(self.positions[n],rot,rot2)
            elif(nshell==1 and mi> 0.5 * self.m0):
                if (mpi.rank == 0):
                    if rejecter:
                        if rejecter.accept(0,0,0):
                            self.positions.append([0,0,0])
                    else:
                        self.positions.append([0,0,0])
            elif(nshell==2):
                if (mpi.rank == 0):
                    position1 = self.rotater([0,0,1],rot,rot2)
                    position2 = self.rotater([0,0,-1],rot,rot2)
                    if rejecter:
                        if rejecter.accept(rii*position1[0],rii*position1[1],rii*position1[2]):
                            self.positions.append(position1)
                        if rejecter.accept(rii*position2[0],rii*position2[1],rii*position2[2]):
                            self.positions.append(position2)
                    else:
                        self.positions.append(position1)
                        self.positions.append(position2)
            elif(nshell==3):
                if (mpi.rank == 0):
                    t = sqrt(3)/2.0
                    position1 = self.rotater([0,1,0],rot,rot2)
                    position2 = self.rotater([t,-0.5,0],rot,rot2)
                    position3 = self.rotater([-t,-0.5,0],rot,rot2)
                    if rejecter:
                        if rejecter.accept(rii*position1[0],rii*position1[1],rii*position1[2]):
                            self.positions.append(position1)
                        if rejecter.accept(rii*position2[0],rii*position2[1],rii*position2[2]):
                            self.positions.append(position2)
                        if rejecter.accept(rii*position3[0],rii*position3[1],rii*position3[2]):
                            self.positions.append(position3)
                    else:
                        self.positions.append(position1)
                        self.positions.append(position2)
                        self.positions.append(position3)
            elif(nshell==4):
                if (mpi.rank == 0):
                    t = sqrt(3.0)/3.0
                    position1 = self.rotater([t,t,t],rot,rot2)
                    position2 = self.rotater([t,-t,-t],rot,rot2)
                    position3 = self.rotater([-t,-t,t],rot,rot2)
                    position4 = self.rotater([-t,t,-t],rot,rot2)
                    if rejecter:
                        if rejecter.accept(rii*position1[0],rii*position1[1],rii*position1[2]):
                            self.positions.append(position1)
                        if rejecter.accept(rii*position2[0],rii*position2[1],rii*position2[2]):
                            self.positions.append(position2)
                        if rejecter.accept(rii*position3[0],rii*position3[1],rii*position3[2]):
                            self.positions.append(position3)
                        if rejecter.accept(rii*position4[0],rii*position4[1],rii*position4[2]):
                            self.positions.append(position4)
                    else:
                        self.positions.append(position1)
                        self.positions.append(position2)
                        self.positions.append(position3)
                        self.positions.append(position4)
            elif(nshell>=163):
                if (nshell > mpi.procs and mpi.procs > 1):
                    p = 0
                    npp = 0
                    if mpi.procs > 2:
                        npp = nshell/(mpi.procs -1)
                    else:
                        npp = (nshell/2) if (rank == 0) else (nshell - nshell/2)
                    print("npp = %d"%npp)
                    
                    if(rank>0 and rank*npp!=nshell):
                        imax = rank*npp + 1
                        for i in range(1,imax):
                            h = -1.0+(2.0*(i-1.0)/(nshell-1.0))
                            if(i>1 and i<nshell):
                                p = (p+3.8/sqrt(nshell)*1.0/sqrt(1.0-h*h))%(2.0*pi)
                    rankmin = rank*npp + 1
                    rankmax = ((rank+1)*npp + 1) if (rank != procs -1) else (nshell + 1)
                    for i in range(rankmin, rankmax):
                        h = -1.0+(2.0*(i-1.0)/(nshell-1.0))
                        t = acos(h)
                        if(i>1 and i<nshell):
                            p = (p+3.8/sqrt(nshell)*1.0/sqrt(1.0-h*h))%(2.0*pi)
                        else:
                            p = 0
                        x = sin(t)*cos(p)
                        y = sin(t)*sin(p)
                        z = cos(t)
                        
                        position = self.rotater([x,y,z],rot,rot2)
                        x = position[0]
                        y = position[1]
                        z = position[2]

                        if rejecter:
                            if rejecter.accept(rii*x,rii*y,rii*z):
                                self.positions.append([x,y,z])
                        else:
                            self.positions.append([x,y,z])
                else:
                    # let rank 0 do all the work
                    p = 0
                    if (mpi.rank == 0):
                        for i in range(1, nshell+1):
                            h = -1.0+(2.0*(i-1.0)/(nshell-1.0))
                            t = acos(h)
                            if(i>1 and i<nshell):
                                p = (p+3.8/sqrt(nshell)*1.0/sqrt(1.0-h*h))%(2.0*pi)
                            else:
                                p = 0
                            x = sin(t)*cos(p)
                            y = sin(t)*sin(p)
                            z = cos(t)
                            position = self.rotater([x,y,z],rot,rot2)
                            x = position[0]
                            y = position[1]
                            z = position[2]

                            if rejecter:
                                if rejecter.accept(rii*x,rii*y,rii*z):
                                    self.positions.append([x,y,z])
                            else:
                                self.positions.append([x,y,z])

            # now reduce some lengths for output
            numNodes = mpi.allreduce(len(self.positions),mpi.SUM)

            print("at r=%3.4g\t wanted %d;\t computed %d total nodes with\t mass=%3.4g" %(rii,nshell,numNodes,mi))
            for n in range(len(self.positions)):
                x       = rii*self.positions[n][0]
                y       = rii*self.positions[n][1]
                z       = rii*self.positions[n][2]
                

                
                if(nshell>1):
                    theta   = acos(z/sqrt(x*x+y*y+z*z))
                    phi     = atan2(y,x)
                    if (phi<0.0):
                        phi = phi + 2.0*pi
                else:
                    theta = (thetaMax - thetaMin)/2.0
                    phi = (phiMax - phiMin)/2.0
                if (theta<=thetaMax and theta>=thetaMin) and (phi<=phiMax and phi>=phiMin):
                    # run a final pass on the rejecter
                    if rejecter:
                        if rejecter.accept(x,y,z):
                            self.x.append(x)
                            self.y.append(y)
                            self.z.append(z)
                            self.m.append(mi)
                            self.H.append(SymTensor3d.one*(1.0/hi))
                    else:
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
            
        print("Generated a total of %i nodes." % mpi.allreduce(len(self.x),mpi.SUM))
        NodeGeneratorBase.__init__(self, False,
                                   self.x, self.y, self.z, self.m, self.H)
        return

    def rotater(self,pos,rot1,rot2):
        posp = [0,0,0]
        for k in range(3):
            for j in range(3):
                posp[k] += pos[j]*rot1[k][j]
                
        x = posp[0]
        y = posp[1]
        z = posp[2]
        
        pos = [x,y,z]
        posp= [0,0,0]
        for k in range(3):
            for j in range(3):
                posp[k] += pos[j]*rot2[k][j]
        return posp

    #---------------------------------------------------------------------------
    # Compute the number of vertices for a given shape at a specific refinement
    # level.
    #  new formula for calculating number of points for a given subdivision level
    #  (Nf * Np(n) - Ne * Npe(n) + Nc)
    #  Nf = Number of faces of primitive shape
    #  Np(n) = Number of points in a triangle subdivided n times
    #       2^(2n-1) + 3*2^(n-1) + 1
    #  Ne = Number of edges of primitive shape
    #  Npe(n) = Number of points along an edge of primitive shape subdivided n times
    #       2^n + 1
    #  Nc = Number of corners
    #---------------------------------------------------------------------------
    def shapeCount(self, refinement, shape):
        Nf  = shape[0]
        Ne  = shape[1]
        Nc  = shape[2]
        n   = refinement
    
        Npe = 2**n + 1
        Np  = 2**(2*n-1) + 3*(2**(n-1)) + 1
        return (Nf * Np - Ne * Npe + Nc)
    
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
        nbp = nbins/procs
        binmin = nbp*rank if (rank!=0) else (1)
        binmax = nbp*(rank+1) if (rank!=procs-1) else (nbins)
        for i in range(binmin,binmax):
            r1 = rmin + (i-1)*dr
            r2 = rmin + i*dr
            result += 0.5*dr*(r2*r2*densityProfileMethod(r2)+r1*r1*densityProfileMethod(r1))
        result = result * (phiMax-phiMin) * (cos(thetaMin)-cos(thetaMax))
        result = mpi.allreduce(result,mpi.SUM)
        return result

    #---------------------------------------------------------------------------
    # Mechanics for creating and refining the icosahedron
    #---------------------------------------------------------------------------
    def addVertex(self,point):
        length = sqrt(point[0]*point[0] + point[1]*point[1] + point[2]*point[2])
        self.positions.append([point[0]/length,point[1]/length,point[2]/length])
        self.index = self.index + 1
        return self.index
        
    def checkMiddlePoint(self,key):
        exists  = 0
        myidx   = 0
        for i in range(len(self.middlePoints)):
            if (self.middlePoints[i][0] == key):
                exists = 1
                myidx = self.middlePoints[i][1]
        return exists, myidx
    
    def getMiddlePoint(self,p1,p2):
        firstIsSmaller = (p1<p2)
        smallerIndex = 0
        greaterIndex = 0
        if firstIsSmaller:
            smallerIndex = p1
            greaterIndex = p2
        else:
            smallerIndex = p2
            greaterIndex = p1
        key = smallerIndex * (1e10) + greaterIndex # some giant number
        
        # check if this key already exists in middlepoints
        exists, idx = self.checkMiddlePoint(key)
        if (exists):
            return idx
        
        # otherwise, not already cached, time to add one
        point1 = self.positions[p1]
        point2 = self.positions[p2]
        middle = [(point1[0]+point2[0])/2.0,(point1[1]+point2[1])/2.0,(point1[2]+point2[2])/2.0]
        
        idx = self.addVertex(middle)
        self.middlePoints.append([key,idx-1])
        
        return idx-1
    
    def createIcoSphere(self,np):
        n = 0
        t = (1.0+sqrt(5.0))/2.0
        # create 12 vertices of an icosahedron
        self.addVertex([-1, t, 0])
        self.addVertex([ 1, t, 0])
        self.addVertex([-1,-t, 0])
        self.addVertex([ 1,-t, 0])
        
        self.addVertex([ 0,-1, t])
        self.addVertex([ 0, 1, t])
        self.addVertex([ 0,-1,-t])
        self.addVertex([ 0, 1,-t])
        
        self.addVertex([ t, 0,-1])
        self.addVertex([ t, 0, 1])
        self.addVertex([-t, 0,-1])
        self.addVertex([-t, 0, 1])
        
        # create the 20 initial faces
        # 5 faces around point 0
        self.faces.append([ 0,11, 5])
        self.faces.append([ 0, 5, 1])
        self.faces.append([ 0, 1, 7])
        self.faces.append([ 0, 7,10])
        self.faces.append([ 0,10,11])
        # 5 adjacent faces
        self.faces.append([ 1, 5, 9])
        self.faces.append([ 5,11, 4])
        self.faces.append([11,10, 2])
        self.faces.append([10, 7, 6])
        self.faces.append([ 7, 1, 8])
        # 5 faces around point 3
        self.faces.append([ 3, 9, 4])
        self.faces.append([ 3, 4, 2])
        self.faces.append([ 3, 2, 6])
        self.faces.append([ 3, 6, 8])
        self.faces.append([ 3, 8, 9])
        # 5 adjacent faces
        self.faces.append([ 4, 9, 5])
        self.faces.append([ 2, 4,11])
        self.faces.append([ 6, 2,10])
        self.faces.append([ 8, 6, 7])
        self.faces.append([ 9, 8, 1])
        
        # now refine triangles until you're done
        for i in range(np):
            faces2 = []
            for j in range(len(self.faces)):
                x,y,z = self.faces[j][0], self.faces[j][1], self.faces[j][2]
                a = self.getMiddlePoint(x,y)
                b = self.getMiddlePoint(y,z)
                c = self.getMiddlePoint(z,x)
                
                faces2.append([x,a,c])
                faces2.append([y,b,a])
                faces2.append([z,c,b])
                faces2.append([a,b,c])
            self.faces = faces2
            n = len(self.positions)

    def createOctaSphere(self,np):
        n = 0
        t = sqrt(2.0)/2.0
        # create the 6 vertices of the octahedron
        self.addVertex([ 0, 0, 1])
        self.addVertex([ t, t, 0])
        self.addVertex([ t,-t, 0])
        self.addVertex([-t,-t, 0])
        self.addVertex([-t, t, 0])
        self.addVertex([ 0, 0,-1])
        
        # create the 8 initial faces
        # 4 faces around point 0
        self.faces.append([ 0, 1, 2])
        self.faces.append([ 0, 2, 3])
        self.faces.append([ 0, 3, 4])
        self.faces.append([ 0, 4, 1])
        # 4 faces around point 5
        self.faces.append([ 5, 2, 1])
        self.faces.append([ 5, 3, 2])
        self.faces.append([ 5, 4, 3])
        self.faces.append([ 5, 1, 4])
        
        # now refine triangles until you're done
        for i in range(np):
            faces2 = []
            for j in range(len(self.faces)):
                x,y,z = self.faces[j][0], self.faces[j][1], self.faces[j][2]
                a = self.getMiddlePoint(x,y)
                b = self.getMiddlePoint(y,z)
                c = self.getMiddlePoint(z,x)
                
                faces2.append([x,a,c])
                faces2.append([y,b,a])
                faces2.append([z,c,b])
                faces2.append([a,b,c])
            self.faces = faces2
            n = len(self.positions)

    def createHexaSphere(self,np):
        n = 0
        t = sqrt(3.0)/2.0
        # create the 5 vertices of the hexahedron
        self.addVertex([ 0, 0, 1])
        self.addVertex([ 0, 1, 0])
        self.addVertex([ t,-0.5,0])
        self.addVertex([-t,-0.5,0])
        self.addVertex([ 0, 0,-1])
        
        # create the 6 initial faces
        # 3 faces around point 0
        self.faces.append([ 0, 1, 2])
        self.faces.append([ 0, 2, 3])
        self.faces.append([ 0, 3, 1])
        # 3 faces around point 4
        self.faces.append([ 4, 2, 1])
        self.faces.append([ 4, 3, 2])
        self.faces.append([ 4, 1, 3])
        
        # now refine triangles until you're done
        for i in range(np):
            faces2 = []
            for j in range(len(self.faces)):
                x,y,z = self.faces[j][0], self.faces[j][1], self.faces[j][2]
                a = self.getMiddlePoint(x,y)
                b = self.getMiddlePoint(y,z)
                c = self.getMiddlePoint(z,x)
                
                faces2.append([x,a,c])
                faces2.append([y,b,a])
                faces2.append([z,c,b])
                faces2.append([a,b,c])
            self.faces = faces2
            n = len(self.positions)

    def createCubicSphere(self,np):
        n = 0
        t = sqrt(3.0)/3.0
        # create the 8 vertices of the cube
        self.addVertex([-t, t, t])
        self.addVertex([-t,-t, t])
        self.addVertex([ t,-t, t])
        self.addVertex([ t, t, t])
        self.addVertex([ t, t,-t])
        self.addVertex([ t,-t,-t])
        self.addVertex([-t,-t,-t])
        self.addVertex([-t, t,-t])
        
        # create the 6 initial faces
        # 5 faces around point 0
        self.faces.append([ 0, 4, 7])
        self.faces.append([ 0, 1, 7])
        self.faces.append([ 0, 1, 2])
        self.faces.append([ 0, 2, 3])
        self.faces.append([ 0, 3, 4])
        # 5 faces around point 5
        self.faces.append([ 5, 2, 3])
        self.faces.append([ 5, 3, 4])
        self.faces.append([ 5, 4, 7])
        self.faces.append([ 5, 6, 7])
        self.faces.append([ 5, 2, 6])
        # 2 faces around point 1
        self.faces.append([ 1, 6, 7])
        self.faces.append([ 1, 2, 6])
        
        # now refine triangles until you're done
        for i in range(np):
            faces2 = []
            for j in range(len(self.faces)):
                x,y,z = self.faces[j][0], self.faces[j][1], self.faces[j][2]
                a = self.getMiddlePoint(x,y)
                b = self.getMiddlePoint(y,z)
                c = self.getMiddlePoint(z,x)
                
                faces2.append([x,a,c])
                faces2.append([y,b,a])
                faces2.append([z,c,b])
                faces2.append([a,b,c])
            self.faces = faces2
            n = len(self.positions)
