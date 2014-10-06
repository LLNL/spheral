#-------------------------------------------------------------------------------
# Package for outputting Spheral++ data in silo files on automatically generated
# Voronoi meshes.
# This data can be easily visualized using Visit:
# ftp://ftp.llnl.gov/pub/pub/visit
#
# Format version:
# 1.0 -- original release
#-------------------------------------------------------------------------------
import os, time, gc, mpi

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

from PolytopeModules import polytope
from siloMeshDump import *

class SpheralVoronoiSiloDump:

    #---------------------------------------------------------------------------
    # Constructor.
    #---------------------------------------------------------------------------
    def __init__(self,
                 baseFileName,
                 baseDirectory = ".",
                 listOfFields = [],
                 listOfFieldLists = [],
                 boundaries = []):

        # Store the file name template.
        self.baseFileName = baseFileName
        self.masterFileName = baseFileName + ".visit"

        # Store the base directory name which will be used to create files
        # under.
        self.baseDirectory = baseDirectory
        if self.baseDirectory[-1] == "":
            self.baseDirectory = "."
        if self.baseDirectory[-1] == "/":
            self.baseDirectory = self.baseDirectory[:-1]

        # The internal list of Fields we're dumping.
        self._fields = []

        # If the user volunteered any Fields to be stored, stash them.
        for field in listOfFields:
            self.addField(field)

        # If the user requested any FieldLists be stored, stash their
        # component fields.
        for fieldList in listOfFieldLists:
            self.addFieldList(fieldList)

        # Build the set of NodeLists for all Fields we're dumping.
        self._nodeLists = self.uniqueNodeLists(self._fields)
        assert len(self._nodeLists) > 0 and len(self._nodeLists) <= len(self._fields)

        # Store the boundaries.
        self._boundaries = boundaries

        # Determine our dimensionality.
        assert len(self._nodeLists) > 0
        if "List2d" in str(self._nodeLists[0]):
            self.dimension = "2d"
        elif "List3d" in str(self._nodeLists[0]):
            self.dimension = "3d"
        else:
            assert False

        # Version of the format written by this class.
        self.version = "1.0"

        return
    
    #---------------------------------------------------------------------------
    # Dump the set of Fields to a dump file.
    #---------------------------------------------------------------------------
    def dump(self, simulationTime, cycle,
             format = "ascii"):

        # Build the name of the directory we will be stuffing the viz file
        # into.
        outputdir = os.path.join(self.baseDirectory, self.baseFileName)

        # Make sure the output directory exists.
        if mpi.rank == 0:
            if not os.path.exists(outputdir):
                try:
                    os.makedirs(outputdir)
                except:
                    raise ValueError, "Cannot create output directory %s" % outputdir
        mpi.barrier()

        # Now build the output file name, including directory.  Make sure
        # the file does not already exist -- if it does we default to overwriting.
        filename = os.path.join(outputdir,
                                self.baseFileName + "-time=%g-cycle=%i" % (simulationTime, cycle))
##         if os.path.exists(filename):
##             raise ValueError, "File %s already exists!  Aborting." % filename

        # Build the set of generators from our points.
        gens = vector_of_double()
        nDim = eval("Vector%s.nDimensions" % self.dimension)
        xmin = vector_of_double(nDim,  1e100)
        xmax = vector_of_double(nDim, -1e100)
        for nodes in self._nodeLists:
            pos = nodes.positions()
            for i in xrange(nodes.numInternalNodes):
                for j in xrange(nDim):
                    gens.append(pos[i][j])
                    xmin[j] = min(xmin[j], pos[i][j])
                    xmax[j] = max(xmax[j], pos[i][j])

        # Check the boundaries for any additional points we want to use for the bounding box.
        for bound in self._boundaries:
            try:
                pb = dynamicCastBoundaryToPlanarBoundary2d(bound)
                for p in (pb.enterPlane.point, pb.exitPlane.point):
                    for j in xrange(nDim):
                        xmin[j] = min(xmin[j], p[j])
                        xmax[j] = max(xmax[j], p[j])
            except:
                pass

        xmin[0] = mpi.allreduce(xmin[0], mpi.MIN)
        xmin[1] = mpi.allreduce(xmin[1], mpi.MIN)
        xmax[0] = mpi.allreduce(xmax[0], mpi.MAX)
        xmax[1] = mpi.allreduce(xmax[1], mpi.MAX)

        # Build the PLC.
        plc = polytope.PLC2d()
        plc.facets.resize(4)
        for i in xrange(4):
            plc.facets[i].resize(2)
            plc.facets[i][0] = i
            plc.facets[i][1] = i % 4
        plccoords = vector_of_double(8)
        plccoords[0] = xmin[0]
        plccoords[1] = xmin[1]
        plccoords[2] = xmax[0]
        plccoords[3] = xmin[1]
        plccoords[4] = xmax[0]
        plccoords[5] = xmax[1]
        plccoords[6] = xmin[0]
        plccoords[7] = xmax[1]

        # Build the tessellation.
        if self.dimension == "2d":
            mesh = polytope.Tessellation2d()
            if "TriangleTessellator2d" in dir(polytope):
                serial_tessellator = polytope.TriangleTessellator2d()
            else:
                assert "BoostTessellator2d" in dir(polytope)
                serial_tessellator = polytope.BoostTessellator2d()
        else:
            assert self.dimension == "3d"
            raise RuntimeError, "Sorry: 3D tessellation silo dumps are not supported yet."
        if mpi.procs > 1:
            tessellator = eval("polytope.DistributedTessellator%s(serial_tessellator, False, True)" % self.dimension)
        else:
            tessellator = serial_tessellator
        tessellator.tessellate(gens, plccoords, plc, mesh)

        # Figure out how many of each type of field we're dumping.
        scalarFields = [x for x in self._fields if isinstance(x, eval("ScalarField%s" % self.dimension))]
        vectorFields = [x for x in self._fields if isinstance(x, eval("VectorField%s" % self.dimension))]
        tensorFields = [x for x in self._fields if isinstance(x, eval("TensorField%s" % self.dimension))]
        symTensorFields = [x for x in self._fields if isinstance(x, eval("SymTensorField%s" % self.dimension))]

        # For tensor fields we like to dump out some extra info.
        for f in (tensorFields + symTensorFields):
            n = f.nodeList()
            tr = eval("ScalarField%s('%s_trace', n)" % (self.dimension, f.name))
            det = eval("ScalarField%s('%s_determinant', n)" % (self.dimension, f.name))
            mineigen = eval("ScalarField%s('%s_eigen_min', n)" % (self.dimension, f.name))
            maxeigen = eval("ScalarField%s('%s_eigen_max', n)" % (self.dimension, f.name))
            for i in xrange(n.numInternalNodes):
                tr[i] = f[i].Trace()
                det[i] = f[i].Determinant()
                eigen = f[i].eigenValues()
                mineigen[i] = eigen.minElement()
                maxeigen[i] = eigen.maxElement()
            scalarFields += [tr, det, mineigen, maxeigen]

        # Write the output.
        timeslice = siloMeshDump(filename, mesh,
                                 nodeLists = self._nodeLists,
                                 time = simulationTime,
                                 cycle = cycle,
                                 scalarFields = scalarFields,
                                 vectorFields = vectorFields,
                                 tensorFields = tensorFields,
                                 symTensorFields = symTensorFields)

        # Write the master file listing all the time slices.
        if mpi.rank == 0:
            mastername = os.path.join(self.baseDirectory, self.baseFileName, self.masterFileName)
            mf = open(mastername, "a")
            mf.write("%s\n" % timeslice)
            mf.close()
        mpi.barrier()

        # That's it.
        del mesh
        while gc.collect():
            pass
        return

    #---------------------------------------------------------------------------
    # Add the given Field to the list of Fields that will be dumped to file.
    #---------------------------------------------------------------------------
    def addField(self, field):
        if ((field.name, field.nodeList().name) not in
            [(f.name, f.nodeList().name) for f in self._fields]):
            self._fields.append(field)
        return

    #---------------------------------------------------------------------------
    # Add the Fields from a FieldList to the set of dumped Fields.
    #---------------------------------------------------------------------------
    def addFieldList(self, fieldList):
        for field in fieldList:
            self.addField(field)
        return

    #---------------------------------------------------------------------------
    # Return the unique set of NodeLists based on the Fields we're dumping.
    #---------------------------------------------------------------------------
    def uniqueNodeLists(self, listOfFields):
        redundantSet = [(f.nodeList().name, f.nodeList()) for f in listOfFields]
        redundantSet.sort()
        assert len(redundantSet) > 0
        result = [redundantSet[0][1]]
        for i in xrange(1, len(redundantSet)):
            if redundantSet[i][0] != redundantSet[i - 1][0]:
                result.append(redundantSet[i][1])
        names = [x.name for x in result]
        assert max([names.count(x) for x in names]) == 1
        return result

    #---------------------------------------------------------------------------
    # Check if the given NodeList is already listed as being dumped.
    #---------------------------------------------------------------------------
    def haveNodeList(self, nodeList):
        return nodeList.name in [n.name for n in self._nodeLists]

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
    assert not dumpGhosts

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
        fields += [x for x in state.allIntFields()]
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

        # We also like to dump the moments of the local point distribution.
        zerothMoment = eval("ScalarFieldList%id(Copy)" % dataBase.nDim)
        firstMoment = eval("VectorFieldList%id(Copy)" % dataBase.nDim)
        W = eval("TableKernel%id(BSplineKernel%id(), 1000)" % (dataBase.nDim, dataBase.nDim))
        zerothAndFirstNodalMoments(dataBase.nodeLists(), W, True, zerothMoment, firstMoment)
        fieldLists += [zerothMoment, firstMoment]

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

    # Now build the visit dumper.
    dumper = SpheralVoronoiSiloDump(baseFileName,
                                    baseDirectory,
                                    fields,
                                    fieldLists,
                                    boundaries)

    # Dump the sucker.
    dumper.dump(currentTime, currentCycle)

    # Try to clean up to free memory.
    del dumper, fields, fieldLists
    while gc.collect():
        pass

    return
