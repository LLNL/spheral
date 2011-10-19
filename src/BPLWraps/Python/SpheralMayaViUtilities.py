#-------------------------------------------------------------------------------
# This is a collection of utilities useful for visualizing Spheral++ results
# using the free MayaVi scientific viz tool ("http://mayavi.sourceforge.net").
# It relies upon PyVTK from Pearu Peterson, available at
# "http://cens.ioc.ee/projects/pyvtk/".
#
# Created by J. Michael Owen, Wed Sep 11 17:55:10 PDT 2002
#-------------------------------------------------------------------------------

import pyvtk
import Spheral
import FieldFunctions
from SpheralGnuPlotUtilities import loadmpi

#-------------------------------------------------------------------------------
# Sample a given FieldListSet to a grid of points in 2-D, returning a new
# FieldListSet of the results and a NodeList defining the sampling grid
# positions.
#-------------------------------------------------------------------------------
def sampleFieldListSetToGrid2d(fieldListSet,
                               dataBase,
                               topGridCellSize,
                               nx = 100,
                               ny = 100,
                               rmin = Spheral.Vector2d(0,0),
                               rmax = Spheral.Vector2d(1,1),
                               W = Spheral.BSplineKernel2d(),
                               boundaryList = []):

    assert nx > 0
    assert ny > 0
    assert rmax.x > rmin.x
    assert rmax.y > rmin.y

    # If this is a parallel run, we have to create a DistributedBoundary
    # and divy up the sampling across processors.
    from SpheralController import insertDistributedBoundary
    domainbc, bcList = insertDistributedBoundary(boundaryList, '2d')

    # Begin by creating a new NodeList defining the positions we are going to
    # sample to.
    # At the same time initialize the weight and H fields for the sampling
    # nodes.
    print 'Generating node distribution.'
    import GenerateNodeDistribution2d
    import DistributeNodes
    nodes = Spheral.AsphNodeList2d(0)
    generator = GenerateNodeDistribution2d. \
                GenerateNodeDistribution2d(nx, ny, 1.0,
                                           'lattice',
                                           rmin = (rmin.x, rmin.y),
                                           rmax = (rmax.x, rmax.y),
                                           nNodePerh = 2.01)
    n = generator.globalNumNodes()
    nodeInfo = DistributeNodes.distributeNodes2d([(nodes, n, generator)])
    assert len(nodeInfo[nodes]['globalNodeListID']) == nodes.numInternalNodes
    for i in xrange(nodes.numInternalNodes):
        assert (nodes.positions[i].x > rmin.x and
                nodes.positions[i].x < rmax.x)
        assert (nodes.positions[i].y > rmin.y and
                nodes.positions[i].y < rmax.y)
        globalID = nodeInfo[nodes]['globalNodeListID'][i]
        nodes.mass[i] = 1.0
        nodes.weight[i] = 1.0
        nodes.massDensity[i] = 1.0
        nodes.Hfield[i] = generator.Htensor(globalID)
    print 'Done.'

    # Construct a neighbor object for this NodeList.
    neighbor0 = dataBase.fluidNodeLists[0].neighbor
    neighborSearchType = 2 # Gather neighbor0.searchType
    numGridLevels = 20
    kernelExtent = W.kernelExtent
    neighbor = Spheral.NestedGridNeighbor2d(nodes,
                                            neighborSearchType,
                                            numGridLevels,
                                            topGridCellSize,
                                            rmin,
                                            kernelExtent)
    nodes.neighbor = neighbor

    # Get FieldLists of the input nodes positions, weights, and Hfield.
    nodePosition = dataBase.fluidPosition
    nodeWeight = dataBase.fluidWeight
    nodeHfield = dataBase.fluidHfield

    # Get FieldList versions of the sampling positions, weights, and Hfield.
    samplePosition = Spheral.VectorFieldList2d(1)
    sampleWeight = Spheral.ScalarFieldList2d(1)
    sampleHfield = Spheral.SymTensorFieldList2d(1)
    samplePosition.appendField(nodes.positions)
    sampleWeight.appendField(nodes.weight)
    sampleHfield.appendField(nodes.Hfield)

    # Make a TableKernel representation of the sampling kernel.
    WT = Spheral.TableKernel2d()
    WT.setTableData(W, 100)

    # Add our sampling NodeList to the DataBase, and clear out any prior boundary info.
    dataBase.appendNodeList(nodes)
    for nodeList in dataBase.fluidNodeLists:
        nodeList.numGhostNodes = 0
        nodeList.neighbor.updateNodes()

    # Set the ghost nodes for any boundary conditions we're enforcing (including
    # the parallel distributed boundary).
    print 'Applying boundary conditions.'
    for boundary in bcList:
        boundary.setGhostNodes(dataBase)
        for nodeList in dataBase.fluidNodeLists:
            boundary.updateGhostNodes(nodeList)
            boundary.applyScalarGhostBoundary(nodeList.weight)
        for fieldList in fieldListSet.ScalarFieldLists:
            boundary.applyScalarFieldListGhostBoundary(fieldList)
        for fieldList in fieldListSet.VectorFieldLists:
            boundary.applyVectorFieldListGhostBoundary(fieldList)
        for fieldList in fieldListSet.TensorFieldLists:
            boundary.applyTensorFieldListGhostBoundary(fieldList)
        for fieldList in fieldListSet.SymTensorFieldLists:
            boundary.applySymTensorFieldListGhostBoundary(fieldList)
    print 'Done.'

    # Sample the FieldListSet to the new positions.
    print 'Sampling data to regular lattice.'
    funcs = FieldFunctions.FieldFunctions()
    sampleFieldListSet = funcs.sampleMultipleFieldsMash2d(fieldListSet,
                                                          nodePosition,
                                                          nodeWeight,
                                                          nodeHfield,
                                                          WT,
                                                          samplePosition,
                                                          sampleWeight,
                                                          sampleHfield)
    print 'Done.'

    # Remove the sampling nodes from the DataBase.
    dataBase.deleteNodeList(nodes)

    # That's it, return the FieldListSet and sampling NodeList.
    return nodes, neighbor, nodeInfo, sampleFieldListSet

#-------------------------------------------------------------------------------
# A helpers class for collecting data from a set of FieldLists, and writing
# them out to a VTK formatted file for post-processing.
#-------------------------------------------------------------------------------
class VtkFile:

    #---------------------------------------------------------------------------
    # Constructor.
    #---------------------------------------------------------------------------
    def __init__(self, fieldListTuples, dataBase, topGridCellSize,
                 fieldListLabels = {},
                 header = 'Spheral VTK file',
                 masterProc = 0,
                 nx = 100,
                 ny = 100,
                 rmin = Spheral.Vector2d(0,0),
                 rmax = Spheral.Vector2d(1,1),
                 W = Spheral.BSplineKernel2d(),
                 boundaryList = []):
        self.masterProc = masterProc
        self.fieldListTuples = fieldListTuples
        self.dataBase = dataBase
        self.topGridCellSize = topGridCellSize
        self.fieldListLabels = fieldListLabels
        self.header = header
        self.nx = nx
        self.ny = ny
        self.rmin = rmin
        self.rmax = rmax
        self.W = W
        self.boundaryList = boundaryList

        # Put the fieldLists into a FieldListSet.  At the same time, build up
        # a dictionary relating the provided names to the FieldLists in the
        # FieldListSet.
        self.fieldListSet = Spheral.FieldListSet2d()
        self.fieldListLabels = []
        for pair in fieldListTuples:
            fieldList = pair[0]
            name = pair[1]
            if type(fieldList) == type(Spheral.ScalarFieldList2d()):
                self.fieldListSet.ScalarFieldLists.append(fieldList)
            elif type(fieldList) == type(Spheral.VectorFieldList2d()):
                self.fieldListSet.VectorFieldLists.append(fieldList)
            elif type(fieldList) == type(Spheral.TensorFieldList2d()):
                self.fieldListSet.TensorFieldLists.append(fieldList)
            elif type(fieldList) == type(Spheral.SymTensorFieldList2d()):
                self.fieldListSet.SymTensorFieldLists.append(fieldList)
            else:
                raise "Unknown FieldListType ", fieldList
            self.fieldListLabels.append(name)

        # Sample the given FieldListSet to a regular grid of points.
        self.sampleNodes, \
                          self.sampleNeighbor,\
                          self.nodeInfo, \
                          self.sampleFieldListSet = \
                          sampleFieldListSetToGrid2d(self.fieldListSet,
                                                     dataBase,
                                                     topGridCellSize,
                                                     nx, ny,
                                                     rmin, rmax,
                                                     W, boundaryList)

        # Determine the distribution of global nodes across processors.
        self.localToGlobalMap = self.determineGlobalNodeDistribution(self.nodeInfo)

        # Collect the FieldList data to the master process
        self.sampleScalars = []
        self.sampleVectors = []
        self.sampleTensors = []
        self.sampleSymTensors = []
        for fieldList in self.sampleFieldListSet.ScalarFieldLists:
            self.sampleScalars.append(self.collectFieldList(fieldList))
        for fieldList in self.sampleFieldListSet.VectorFieldLists:
            self.sampleVectors.append(self.collectFieldList(fieldList))
        for fieldList in self.sampleFieldListSet.TensorFieldLists:
            self.sampleTensors.append(self.collectFieldList(fieldList))
        for fieldList in self.sampleFieldListSet.SymTensorFieldLists:
            self.sampleSymTensors.append(self.collectFieldList(fieldList))

        # That's it, we should now be ready to write data.
        return

    #---------------------------------------------------------------------------
    # Create a map of all global node IDs for local node IDs per processor.
    #---------------------------------------------------------------------------
    def determineGlobalNodeDistribution(self, nodeInfo):

        # Create nested dictionaries, which will relate
        # procID -> localNodeID -> global NodeIDs.
        localToGlobalIDMap = {}

        try:
            import mpi
            procID = mpi.rank
        except:
            procID = 0

        # Set up this processors portion of the map.
        localIDMap = {procID: nodeInfo[self.sampleNodes]['globalNodeID'][:]}
##        localID = 0
##        for globalNodeID in nodeInfo[self.sampleNodes]['globalNodID']:
##            localIDMap[procID][localID] = globalNodeID
##            localID += 1

        # Now collect the complete map of procID ->localID -> globalID to the
        # master process.
        globalIDMap = localIDMap
        if mpi:
            if mpi.rank == self.masterProc:
                recvRange = range(0, self.masterProc) + \
                            range(self.masterProc + 1, mpi.procs)
                for recvID in recvRange:
                    globalIDMap[recvID] = mpi.recv(recvID)[0]
##                    recvDict = mpi.recv(recvID)[0]
##                    for key in recvDict:
##                        globalIDMap[recvID][key] = recvDict[key]
            else:
                mpi.send(localIDMap[procID], self.masterProc)

        return globalIDMap

    #-------------------------------------------------------------------------------
    # Convert a list of Vector2d to a list of three tuples.
    #-------------------------------------------------------------------------------
    def convertVector2dElements(self, listOfVectors):
        for i in xrange(len(listOfVectors)):
            listOfVectors[i] = tuple(listOfVectors[i][:] + [0.0])

    #-------------------------------------------------------------------------------
    # Convert a list of Tensor2d to a list of 3x3 tuples.
    #-------------------------------------------------------------------------------
    def convertTensor2dElements(self, listOfTensors):
        for i in xrange(len(listOfTensors)):
            listOfTensors[i] = ((listOfTensors[i].xx, listOfTensors[i].xy, 0.0),
                                (listOfTensors[i].yx, listOfTensors[i].yy, 0.0),
                                (0.0, 0.0, 1.0))

    #-------------------------------------------------------------------------------
    # Collect the values from the the given FieldList to processor zero as a single
    # contiguous list.
    #-------------------------------------------------------------------------------
    def collectFieldList(self, fieldList,
                         targetProc = 0):

        try:
            import mpi
            procID = mpi.rank
        except:
            procID = 0

        # Make a list out of the FieldList on each processor.
        listLocal = []
        for field in fieldList:
            listLocal += field[:field.nodeList.numInternalNodes]
        assert len(listLocal) == len(self.localToGlobalMap[procID])

        # If this is a Vector or Tensor FieldList, convert the thing to
        # a list of tuples of the values.
        if listLocal:
            if type(listLocal[0]) == type(Spheral.Vector2d()):
                self.convertVector2dElements(listLocal)
            elif (type(listLocal[0]) == type(Spheral.Tensor2d()) or
                  type(listLocal[0]) == type(Spheral.SymTensor2d())):
                self.convertTensor2dElements(listLocal)

        # The master process reserves a list of the appropriate size
        # to hold all the field values.
        if procID == targetProc:
            nGlobal = 0
            for proc in self.localToGlobalMap:
                nGlobal += len(self.localToGlobalMap[proc])
            listGlobal = [None]*nGlobal

            # Go ahead and fill the local values into the appropriate slots
            # in the global list.
            for localID in xrange(len(listLocal)):
                globalID = self.localToGlobalMap[procID][localID]
                assert globalID >= 0 and globalID < len(listGlobal)
                listGlobal[globalID] = listLocal[localID]

            # Receive the other processors values for the field values.
            if mpi:
                recvRange = range(0, self.masterProc) + \
                            range(self.masterProc + 1, mpi.procs)
                for recvID in recvRange:
                    recvList = mpi.recv(recvID)[0]
                    for localID in xrange(len(recvList)):
                        globalID = self.localToGlobalMap[recvID][localID]
                        assert globalID >= 0 and globalID < len(listGlobal)
                        listGlobal[globalID] = recvList[localID]

        # Have all other processors send their info to the master.
        else:
            if mpi:
                mpi.send(listLocal, targetProc)

        # That's it, return the resulting list.
        if procID == targetProc:
            assert None not in listGlobal
            return listGlobal
        else:
            return None

    #-------------------------------------------------------------------------------
    # Write the given set of FieldLists to a VTK file.
    #-------------------------------------------------------------------------------
    def writeVTKFieldLists(self, fileName,
                           format = 'ascii'):

        # In the parallel case only the "master" process should write the file.
        try:
            import mpi
            procID = mpi.rank
        except:
            procID = 0
        if procID == self.masterProc:

            # Create a VTK structure (in this case structured points).
            assert self.nx > 1
            assert self.ny > 1
            assert self.rmin.x < self.rmax.x
            assert self.rmin.y < self.rmax.y
            dim = (self.nx, self.ny, 1)
            dr = ((self.rmax.x - self.rmin.x)/dim[0],
                  (self.rmax.y - self.rmin.y)/dim[1],
                  1.0)
            structure = pyvtk.StructuredPoints(dim,
                                               origin = (self.rmin.x, self.rmin.y, 0.0),
                                               spacing = dr)

            # Create PointData data fields for each of the fields we're storing.
            vtkPointData = []
            allData = self.sampleScalars + \
                      self.sampleVectors + \
                      self.sampleTensors + \
                      self.sampleSymTensors
            i = 0
            for fieldData in self.sampleScalars:
                assert i < len(self.fieldListLabels)
                vtkPointData.append(pyvtk.Scalars(fieldData,
                                                  name = self.fieldListLabels[i]))
                i += 1
            for fieldData in self.sampleVectors:
                assert i < len(self.fieldListLabels)
                vtkPointData.append(pyvtk.Vectors(fieldData,
                                                  name = self.fieldListLabels[i]))
                i += 1
            for fieldData in self.sampleTensors:
                assert i < len(self.fieldListLabels)
                vtkPointData.append(pyvtk.Tensors(fieldData,
                                                  name = self.fieldListLabels[i]))
                i += 1
            for fieldData in self.sampleSymTensors:
                assert i < len(self.fieldListLabels)
                vtkPointData.append(pyvtk.Tensors(fieldData,
                                                  name = self.fieldListLabels[i]))
                i += 1

            # Create a VTK file object, and write that file out.
            file = pyvtk.VtkData(structure,
                                 header = self.header)
            for data in vtkPointData:
                file.point_data.append(data)
            file.tofile(fileName, format)

        # That's it.
        return
