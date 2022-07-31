from math import *

from NodeGeneratorBase import *

from Spheral import Vector1d, Vector2d, Vector3d, SymTensor1d, SymTensor2d, SymTensor3d
from SpheralTestUtilities import fuzzyEqual

#-------------------------------------------------------------------------------
# Wrapper Generator for DEM based on SPH generators
#-------------------------------------------------------------------------------
class GenerateDEMfromSPHGenerator1d(NodeGeneratorBase):
    
    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def __init__(self,
                 W,
                 SPHGenerator,
                 DEMParticleGenerator=None,
                 particleRadius=None,
                 nPerh=None,
                 kernelExtentOverParticleRadius=2.05): 

        assert (particleRadius or nPerh), "a constant particle radius or the nPerh must be specified"
        assert (kernelExtentOverParticleRadius >= 2.0)

        hOverR = kernelExtentOverParticleRadius/W.kernelExtent

        # set up our initial radius
        #--------------------------------------------------
        initialParticleRadius=[]
        if not particleRadius:
            for i in range(SPHGenerator.localNumNodes()):
                hi = 1.0/SPHGenerator.H[i](0,0)
                di = hi/nPerh
                initialParticleRadius.append(di/2.0)
        elif isfloat(particleRadius):
            for i in range(SPHGenerator.localNumNodes()):
                initialParticleRadius.append(particleRadius)
        else:
            raise ValueError, "particleRadius must be none or float"

        # create composite particle distribution from base
        #--------------------------------------------------
        if DEMParticleGenerator:
            self.x = []
            self.m = []
            self.particleRadius = []
            self.compositeParticleIndex = []
            self.H = []
            for i in range(SPHGenerator.localNumNodes()):

                xi,mi,Ri = DEMParticleGenerator(SPHGenerator.x[i],
                                                SPHGenerator.H[i],
                                                SPHGenerator.m[i],
                                                initialParticleRadius[i])

                pIDi = [SPHGenerator.globalIDs[i] for entry in Ri]
                self.x.extend(xi)
                self.m.extend(mi)
                self.particleRadius.extend(Ri)
                self.compositeParticleIndex.extend(pIDi)
                
        else:
            self.x = SPHGenerator.x
            self.m = SPHGenerator.m
            self.particleRadius = initialParticleRadius
            self.compositeParticleIndex = SPHGenerator.globalIDs
                
        self.H = [SymTensor1d(1.0/(Rj*hOverR)) for Rj in self.particleRadius]

        NodeGeneratorBase.__init__(self, False,
                                   self.x, self.m, self.H, self.particleRadius, self.compositeParticleIndex)
        return

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localParticleRadius(self, i):
        return self.particleRadius[i]

    def globalParticleRadius(self, i):
        return self._globalValue(i, self.localParticleRadius)

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localCompositeParticleIndex(self, i):
        return self.compositeParticleIndex[i]

    def globalCompositeParticleIndex(self, i):
        return self._globalValue(i, self.localCompositeParticleIndex)

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localPosition(self, i):
        return Vector1d(self.x[i])
    
    #---------------------------------------------------------------------------
    # Get the mass for the given node index.
    #---------------------------------------------------------------------------
    def localMass(self, i):
        return self.m[i]
    
    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        return self.H[i]
    

#-------------------------------------------------------------------------------
# Wrapper Generator for DEM based on SPH generators
#-------------------------------------------------------------------------------
class GenerateDEMfromSPHGenerator2d(NodeGeneratorBase):
    
    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def __init__(self,
                 W,
                 SPHGenerator,
                 DEMParticleGenerator=None,
                 particleRadius=None,
                 nPerh=None,
                 kernelExtentOverParticleRadius=2.05): 

        assert (particleRadius or nPerh), "a constant particle radius or the nPerh must be specified"
        assert (kernelExtentOverParticleRadius >= 2.0)

        hOverR = kernelExtentOverParticleRadius/W.kernelExtent

        # set up our initial radius
        #--------------------------------------------------
        initialParticleRadius=[]
        if not particleRadius:
            SPHGenerator.makeHround()
            for i in range(SPHGenerator.localNumNodes()):
                hi = 1.0/SPHGenerator.H[i](0,0)
                di = hi/nPerh
                initialParticleRadius.append(di/2.0)
        elif isfloat(particleRadius):
            for i in range(SPHGenerator.localNumNodes()):
                initialParticleRadius.append(particleRadius)
        else:
            raise ValueError, "particleRadius must be none or float"

        # create composite particle distribution from base
        #--------------------------------------------------
        if DEMParticleGenerator:
            self.x = []
            self.y = []
            self.m = []
            self.particleRadius = []
            self.compositeParticleIndex = []
            for i in range(SPHGenerator.localNumNodes()):

                xi,yi,mi,Ri = DEMParticleGenerator(SPHGenerator.x[i],
                                                   SPHGenerator.y[i],
                                                   SPHGenerator.H[i],
                                                   SPHGenerator.m[i],
                                                   initialParticleRadius[i])

                pIDi = [SPHGenerator.globalIDs[i] for entry in Ri]
                self.x.extend(xi)
                self.y.extend(yi)
                self.m.extend(mi)
                self.particleRadius.extend(Ri)
                self.compositeParticleIndex.extend(pIDi)
                
        else:
            self.x = SPHGenerator.x
            self.y = SPHGenerator.y
            self.m = SPHGenerator.m
            self.particleRadius = initialParticleRadius
            self.compositeParticleIndex = SPHGenerator.globalIDs
    
        self.H = [SymTensor2d(1.0/(Rj*hOverR), 0.0, 0.0, 1.0/(Rj*hOverR)) for Rj in self.particleRadius]

        NodeGeneratorBase.__init__(self, False,
                                   self.x, self.y, self.m, self.H, self.particleRadius, self.compositeParticleIndex)
        return

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localParticleRadius(self, i):
        return self.particleRadius[i]

    def globalParticleRadius(self, i):
        return self._globalValue(i, self.localParticleRadius)

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localCompositeParticleIndex(self, i):
        return self.compositeParticleIndex[i]

    def globalCompositeParticleIndex(self, i):
        return self._globalValue(i, self.localCompositeParticleIndex)

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localPosition(self, i):
        return Vector2d(self.x[i],self.y[i])
    
    #---------------------------------------------------------------------------
    # Get the mass for the given node index.
    #---------------------------------------------------------------------------
    def localMass(self, i):
        return self.m[i]
    
    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        return self.H[i]
    


#-------------------------------------------------------------------------------
# Wrapper Generator for DEM based on SPH generators
#-------------------------------------------------------------------------------
class GenerateDEMfromSPHGenerator3d(NodeGeneratorBase):
    
    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def __init__(self,
                 W,
                 SPHGenerator,
                 DEMParticleGenerator=None,
                 particleRadius=None,
                 nPerh=None,
                 kernelExtentOverParticleRadius=2.05): 

        assert (particleRadius or nPerh), "a constant particle radius or the nPerh must be specified"
        assert (kernelExtentOverParticleRadius >= 2.0)

        hOverR = kernelExtentOverParticleRadius/W.kernelExtent

        # set up our initial radius
        #--------------------------------------------------
        initialParticleRadius=[]
        if not particleRadius:
            SPHGenerator.makeHround()
            for i in range(SPHGenerator.localNumNodes()):
                hi = 1.0/SPHGenerator.H[i](0,0)
                di = hi/nPerh
                initialParticleRadius.append(di/2.0)
        elif isfloat(particleRadius):
            for i in range(SPHGenerator.localNumNodes()):
                initialParticleRadius.append(particleRadius)
        else:
            raise ValueError, "particleRadius must be none or float"

        # create composite particle distribution from base
        #--------------------------------------------------
        if DEMParticleGenerator:
            self.x = []
            self.y = []
            self.z = []
            self.m = []
            self.particleRadius = []
            self.compositeParticleIndex = []
            self.H = []
            for i in range(SPHGenerator.localNumNodes()):

                xi,yi,zi,mi,Ri = DEMParticleGenerator(SPHGenerator.x[i],
                                                   SPHGenerator.y[i],
                                                   SPHGenerator.z[i],
                                                   SPHGenerator.H[i],
                                                   SPHGenerator.m[i],
                                                   initialParticleRadius[i])

                pIDi = [SPHGenerator.globalIDs[i] for entry in Ri]
                self.x.extend(xi)
                self.y.extend(yi)
                self.z.extend(zi)
                self.m.extend(mi)
                self.particleRadius.extend(Ri)
                self.compositeParticleIndex.extend(pIDi)
                
        else:
            self.x = SPHGenerator.x
            self.y = SPHGenerator.y
            self.z = SPHGenerator.z
            self.m = SPHGenerator.m
            self.particleRadius = initialParticleRadius
            self.compositeParticleIndex = SPHGenerator.globalIDs
                
        self.H = [SymTensor3d(1.0/(Rj*hOverR), 0.0, 0.0, 0.0, 1.0/(Rj*hOverR), 0.0, 0.0, 0.0, 1.0/(Rj*hOverR)) for Rj in self.particleRadius]

        NodeGeneratorBase.__init__(self, False,
                                   self.x, self.y, self.z, self.m, self.H, self.particleRadius, self.compositeParticleIndex)
        return

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localParticleRadius(self, i):
        return self.particleRadius[i]

    def globalParticleRadius(self, i):
        return self._globalValue(i, self.localParticleRadius)

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localCompositeParticleIndex(self, i):
        return self.compositeParticleIndex[i]

    def globalCompositeParticleIndex(self, i):
        return self._globalValue(i, self.localCompositeParticleIndex)

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localPosition(self, i):
        return Vector3d(self.x[i],self.y[i],self.z[i])
    
    #---------------------------------------------------------------------------
    # Get the mass for the given node index.
    #---------------------------------------------------------------------------
    def localMass(self, i):
        return self.m[i]
    
    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        return self.H[i]
    
