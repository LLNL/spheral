from math import *

import loadmpi
mpi, rank, procs = loadmpi.loadmpi()

#-------------------------------------------------------------------------------
# NodeGeneratorBase.
#
# Base class for Spheral++ NodeGenerators.  Implements the default methods
# for getting node info by global IDs and such.
#-------------------------------------------------------------------------------
class NodeGeneratorBase:

    def __init__(self,
                 serialInitialization,
                 *vars):

        if serialInitialization:
            ntot = len(vars[0])
            minGlobalID, maxGlobalID = self.globalIDRange(ntot)
            self.globalIDs = range(minGlobalID, maxGlobalID)
            self._cullVars(minGlobalID, maxGlobalID, *vars)

        return

    #---------------------------------------------------------------------------
    # Simple minded method to compute and assign unique global ID ranges to
    # each processor.
    #---------------------------------------------------------------------------
    def globalIDRange(self, ntot):
        ndomain0 = ntot/procs
        remainder = ntot % procs
        assert remainder < procs
        ndomain = ndomain0
        if rank < remainder:
            ndomain += 1
        minGlobalID = rank*ndomain0 + min(rank, remainder)
        maxGlobalID = minGlobalID + ndomain

        assert mpi.allreduce(minGlobalID, mpi.MIN) == 0
        assert mpi.allreduce(maxGlobalID, mpi.MAX) == ntot

        return minGlobalID, maxGlobalID

    #---------------------------------------------------------------------------
    # Cull the given set of variables to the given indicies.
    #---------------------------------------------------------------------------
    def _cullVars(self, minID, maxID, *args):
        for arg in args:
            del arg[maxID:]
            del arg[:minID]
        return

    #---------------------------------------------------------------------------
    # Find the domain that possesses the given globalID.
    #---------------------------------------------------------------------------
    def domainForID(self, i):
        result = -1
        if i in self.globalIDs:
            result = rank
        mpi.allreduce(result, mpi.MAX)
        return result

    #---------------------------------------------------------------------------
    # Get the position for the given global node index.
    #---------------------------------------------------------------------------
    def globalPosition(self, i):
        return self._globalValue(i, self.localPosition)

    #-------------------------------------------------------------------------------
    # Get the mass for the given global node index.
    #-------------------------------------------------------------------------------
    def globalMass(self, i):
        return self._globalValue(i, self.localMass)

    #-------------------------------------------------------------------------------
    # Get the mass density for the given global node index.
    #-------------------------------------------------------------------------------
    def globalMassDensity(self, i):
        return self._globalValue(i, self.localMassDensity)

    #-------------------------------------------------------------------------------
    # Get the H tensor for the given global node index.
    #-------------------------------------------------------------------------------
    def globalHtensor(self, i):
        return self._globalValue(i, self.localHtensor)

    #-------------------------------------------------------------------------------
    # Get the local number of nodes.
    #-------------------------------------------------------------------------------
    def localNumNodes(self):
        return len(self.m)

    #-------------------------------------------------------------------------------
    # Get the total number of nodes.
    #-------------------------------------------------------------------------------
    def globalNumNodes(self):
        return mpi.allreduce(len(self.m), mpi.SUM)

    #---------------------------------------------------------------------------
    # Generic method to get the global value.
    #---------------------------------------------------------------------------
    def _globalValue(self, i, method):
        dom = self.domainForID(i)
        result = None
        if rank == dom:
            j = self.globalIDs[i]
            result = method(j)
        mpi.bcast(result, dom)
        return result

#-------------------------------------------------------------------------------
# Helper class for providing default constant density behaviour.
#-------------------------------------------------------------------------------
class ConstantRho:
    def __init__(self, rho0):
        self.rho0 = rho0
        return
    def rho(self, r):
        return self.rho0
    def __call__(self, r):
        return self.rho(r)

