#-------------------------------------------------------------------------------
# Dumper class for outputting Spheral++ format data for visualization with
# Visit.
# ftp://ftp.llnl.gov/pub/pub/visit
#
# Format version:
# 1.0 -- original release
# 1.1 -- replacing all spaces in strings with underscores, so Visit can
#        manipulate them.
#-------------------------------------------------------------------------------
from SpheralModules import *
from SpheralModules.Spheral import *
from SpheralModules.Spheral.NodeSpace import *
from SpheralModules.Spheral.FieldSpace import *
from SpheralModules.Spheral.DataBaseSpace import *
from SpheralModules.Spheral.FileIOSpace import *
from SpheralModules.Spheral.ArtificialViscositySpace import *
from SpheralModules.Spheral.DataOutput import *
from SpheralModules.Spheral.KernelSpace import *
from SpheralModules.Spheral.NeighborSpace import *
from SpheralModules.Spheral.Material import *
from SpheralModules.Spheral.BoundarySpace import *
from SpheralModules.Spheral.PhysicsSpace import *
from SpheralModules.Spheral.GravitySpace import *
from SpheralModules.Spheral.IntegratorSpace import *

import os, time, mpi
from siloPointmeshDump import siloPointmeshDump

#-------------------------------------------------------------------------------
# Dump out all the Fields in a State object.
# You can pass any of the following for stateThingy:
#     Integrator
#     State
#-------------------------------------------------------------------------------
def dumpPhysicsState(stateThingy,
                     baseFileName,
                     baseDirectory = ".",
                     fields = None,
                     fieldLists = None,
                     currentTime = None,
                     currentCycle = None,
                     dumpGhosts = False,
                     dumpDerivatives = False,
                     boundaries = None):

    # What did we get passed?
    if max([isinstance(stateThingy, x) for x in [Integrator1d, Integrator2d, Integrator3d]]):
        integrator = stateThingy
        dataBase = integrator.dataBase()
        state = eval("State%id(integrator.dataBase(), integrator.physicsPackages())" % integrator.dataBase().nDim)
        derivs = None
        if dumpDerivatives:
            derivs = eval("StateDerivatives%id(integrator.dataBase(), integrator.physicsPackages())" % integrator.dataBase().nDim)
        currentTime = integrator.currentTime
        currentCycle = integrator.currentCycle
    elif max([isinstance(stateThingy, x) for x in [State1d, State2d, State3d]]):
        integrator = None
        dataBase = None
        state = stateThingy
        derivs = None
        assert currentTime is not None
        assert currentCycle is not None
    assert state is not None

    # Did the user specify any data to be dumped?
    if not fields:
        fields = []
    if not fieldLists:
        fieldLists = []

    # Build up the list of fields in the state object.
    fields += [x for x in state.allIntFields()]
    fields += [x for x in state.allScalarFields()]
    fields += [x for x in state.allVectorFields()]
    fields += [x for x in state.allTensorFields()]
    fields += [x for x in state.allSymTensorFields()]

    # Are we also dumping the derivative fields?
    if not derivs is None:
        fields += [x for x in derivs.allIntFields()]
        fields += [x for x in derivs.allScalarFields()]
        fields += [x for x in derivs.allVectorFields()]
        fields += [x for x in derivs.allTensorFields()]
        fields += [x for x in derivs.allSymTensorFields()]

    # If available, add the work, H inverse and hmin, hmax, & hmin_hmax_ratio by default.
    if dataBase:
        work = dataBase.globalWork
        fieldLists.append(work)
        Hfl = dataBase.fluidHfield
        Hi = dataBase.newGlobalSymTensorFieldList()
        dataBase.fluidHinverse(Hi)
        fieldLists.append(Hi)
        hmin = dataBase.newGlobalScalarFieldList()
        hmax = dataBase.newGlobalScalarFieldList()
        hminhmax = dataBase.newGlobalScalarFieldList()
        for H, fmin, fmax, fratio in zip(Hfl,
                                         hmin,
                                         hmax,
                                         hminhmax):
            fmin.name = "hmin"
            fmax.name = "hmax"
            fratio.name = "hmin_hmax_ratio"
            n = H.nodeList().numInternalNodes
            for i in xrange(n):
                ev = H[i].eigenValues()
                fmin[i] = 1.0/ev.maxElement()
                fmax[i] = 1.0/ev.minElement()
                fratio[i] = ev.minElement()/ev.maxElement()
        fieldLists.append(hmin)
        fieldLists.append(hmax)
        fieldLists.append(hminhmax)

    # Add a domain decomposition tag (if we're parallel).
    try:
        import mpi
        domains = dataBase.newGlobalScalarFieldList()
        for f in domains:
            f.name = "Domains"
            for i in xrange(f.nodeList().numInternalNodes):
                f[i] = mpi.rank
        fieldLists.append(domains)
    except:
        pass

    # Dump the sucker.
    fullBaseName = baseFileName + "-time=%g-cycle=%i" % (currentTime, currentCycle)
    siloPointmeshDump(baseName = fullBaseName,
                      baseDirectory = baseDirectory,
                      fields = fields,
                      fieldLists = fieldLists,
                      time = currentTime,
                      cycle = currentCycle)

    # Add to the master file.
    masterFileName = os.path.join(baseDirectory, baseFileName + ".visit")
    mf = open(masterFileName, "a")
    mf.write("%s\n" % (fullBaseName + ".silo"))
    mf.close()

    return

