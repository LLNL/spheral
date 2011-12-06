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
import Spheral
import time

# Load the mpi module if we're parallel.
import loadmpi
mpi, rank, procs = loadmpi.loadmpi()

class SpheralVisitDump:

    #---------------------------------------------------------------------------
    # Constructor.
    #---------------------------------------------------------------------------
    def __init__(self,
                 baseFileName,
                 baseDirectory = ".",
                 listOfFields = [],
                 listOfFieldLists = [],
                 dumpGhosts = False):

        # Should we dump ghost values.
        self.dumpGhosts = dumpGhosts

        # Store the file name template.
        self.extension = ".sv"
        self.baseFileName = baseFileName + "-time=%g-cycle=%i"
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

        # A dictionary of the valid keywords in .sv files.
        self.keyWords = {"nodelist": "!NodeList",
                         "field": "!Field",
                         "header": "!Header",
                         "endheader": "!EndHeader",
                         "time": "!Time",
                         "cycle": "!Cycle",
                         "filelist": "!FileList",
                         "ascii": "!ASCIIData",
                         "binary": "!BinaryData",
                         "vector": "Vector",
                         "tensor": "Tensor",
                         "symtensor": "SymTensor",
                         }

        # Build the internal map of Spheral data types to descriptors.
        self.dataType = {type(Spheral.ScalarField1d("thpt")): "Scalar 1",
                         type(Spheral.VectorField1d("thpt")): "Vector 1",
                         type(Spheral.TensorField1d("thpt")): "Tensor 1 1",
                         type(Spheral.SymTensorField1d("thpt")): "SymTensor 1 1",
                         type(Spheral.IntField1d("thpt")): "Scalar 1",
                         type(Spheral.ScalarField2d("thpt")): "Scalar 1",
                         type(Spheral.VectorField2d("thpt")): "Vector 2",
                         type(Spheral.TensorField2d("thpt")): "Tensor 2 2",
                         type(Spheral.SymTensorField2d("thpt")): "SymTensor 2 2",
                         type(Spheral.IntField2d("thpt")): "Scalar 1",
                         type(Spheral.ScalarField3d("thpt")): "Scalar 1",
                         type(Spheral.VectorField3d("thpt")): "Vector 3",
                         type(Spheral.TensorField3d("thpt")): "Tensor 3 3",
                         type(Spheral.SymTensorField3d("thpt")): "SymTensor 3 3",
                         type(Spheral.IntField3d("thpt")): "Scalar 1",
                         }

        # Build the map of data types to element writers.
        self.elementDumper = {type(Spheral.ScalarField1d("thpt")): self._dumpFloat,
                              type(Spheral.VectorField1d("thpt")): self._dumpVector,
                              type(Spheral.TensorField1d("thpt")): self._dumpTensor,
                              type(Spheral.SymTensorField1d("thpt")): self._dumpTensor,
                              type(Spheral.IntField1d("thpt")): self._dumpFloat,
                              type(Spheral.ScalarField2d("thpt")): self._dumpFloat,
                              type(Spheral.VectorField2d("thpt")): self._dumpVector,
                              type(Spheral.TensorField2d("thpt")): self._dumpTensor,
                              type(Spheral.SymTensorField2d("thpt")): self._dumpTensor,
                              type(Spheral.IntField2d("thpt")): self._dumpFloat,
                              type(Spheral.ScalarField3d("thpt")): self._dumpFloat,
                              type(Spheral.VectorField3d("thpt")): self._dumpVector,
                              type(Spheral.TensorField3d("thpt")): self._dumpTensor,
                              type(Spheral.SymTensorField3d("thpt")): self._dumpTensor,
                              type(Spheral.IntField3d("thpt")): self._dumpFloat,
                         }

        # Version of the format written by this class.
        self.version = "1.1"

        return
    
    #---------------------------------------------------------------------------
    # Dump the set of Fields to a dump file.
    #---------------------------------------------------------------------------
    def dump(self, simulationTime, cycle,
             format = "ascii"):

        # Build the name of the directory we will be stuffing the viz file
        # into.
        outputdir = self.baseDirectory
        if mpi.procs > 1:
            outputdir += "/domain%04i" % mpi.rank

        # Make sure the output directory exists.
        import os
        if not os.path.exists(outputdir):
            try:
                os.makedirs(outputdir)
            except:
                raise "Cannot create output directory %s" % outputdir
        mpi.barrier()

        # Now build the output file name, including directory.  Make sure
        # the file does not already exist.
        filename = os.path.join(outputdir,
                                self.baseFileName % (simulationTime, cycle) + 
                                self.extension)
        if os.path.exists(filename):
            raise "File %s already exists!  Aborting." % filename

        # Open the output file for writing.
        try:
            f = open(filename, "w")
        except:
            raise "Unable to open file %s for output" % filename

        # Build the header info for this file, and write it.
        header = self.constructHeader(simulationTime, cycle)
        f.writelines(header)

        # Write whether we're dumping text or binary data.
        if format == "binary":
            f.write(self.keyWords["binary"] + "\n")

            # Close the file, and reopen it in binary mode.
            f.close()
            f = open(filename, "ab")

        else:
            assert format == "ascii"
            f.write(self.keyWords["ascii"] + "\n")

        # Loop over each NodeList.
        for nodeList in self._nodeLists:

            # Write the NodeList declarator.
            if self.dumpGhosts:
                n = nodeList.numNodes
            else:
                n = nodeList.numInternalNodes
            f.write(self.keyWords["nodelist"] + ' "' +
                    nodeList.name() + '" %i\n' % n)

            # Loop over the Fields, and call the private method to write
            # the Field data.
            for field in self.fields(nodeList):
                self._dumpField(field, f)

        # Close the data file.
        f.close()

        # Now, if this is a multidomain run we need to write a root
        # file linking together the individiual domain files.
        if mpi.procs > 1 and mpi.rank == 0:

            # Build the root file name, and open the file.
            rootname = self.baseFileName % (simulationTime, cycle) + "-root" + self.extension
            rootfile = os.path.join(self.baseDirectory, rootname)
            if os.path.exists(rootfile):
                raise "File %s already exists!  Aborting." % rootfile
            try:
                f = open(rootfile, "w")
            except:
                raise "Unable to open root file %s for output" % rootfile

            # Write the standard header info.
            f.writelines(header)

            # Now write the individual domain file names.
            f.write(self.keyWords["filelist"] + " %i\n" % mpi.procs)
            for domain in xrange(mpi.procs):
                filename = "domain%04i" % domain + "/" + \
                           self.baseFileName % (simulationTime, cycle) + \
                           self.extension
                f.write(filename + "\n")

            # That's it, close the root file.
            f.close()

            # Add this root file to the master file.
            mastername = os.path.join(self.baseDirectory, self.masterFileName)
            mf = open(mastername, "a")
            mf.write("%s\n" % rootname)
            mf.close()

        elif mpi.procs == 1:

            # In serial we need to provide the master file.
            mastername = os.path.join(self.baseDirectory, self.masterFileName)
            mf = open(mastername, "a")
            mf.write("%s\n" % os.path.basename(filename))
            mf.close()

        mpi.barrier()
        return

    #---------------------------------------------------------------------------
    # Add the given Field to the list of Fields that will be dumped to file.
    #---------------------------------------------------------------------------
    def addField(self, field):
        if ((field.name(), field.nodeList().name()) not in
            [(f.name(), f.nodeList().name()) for f in self._fields]):
            self._fields.append(field)
        return

    #---------------------------------------------------------------------------
    # Add the Fields from a FieldList to the set of dumped Fields.
    #---------------------------------------------------------------------------
    def addFieldList(self, fieldList):
        for field in fieldList.fields():
            self.addField(field)
        return

    #---------------------------------------------------------------------------
    # Return the unique set of NodeLists based on the Fields we're dumping.
    #---------------------------------------------------------------------------
    def uniqueNodeLists(self, listOfFields):
        redundantSet = [(f.nodeList().name(), f.nodeList()) for f in listOfFields]
        redundantSet.sort()
        assert len(redundantSet) > 0
        result = [redundantSet[0][1]]
        for i in xrange(1, len(redundantSet)):
            if redundantSet[i][0] != redundantSet[i - 1][0]:
                result.append(redundantSet[i][1])
        names = [x.name() for x in result]
        assert max([names.count(x) for x in names]) == 1
        return result

    #---------------------------------------------------------------------------
    # Check if the given NodeList is already listed as being dumped.
    #---------------------------------------------------------------------------
    def haveNodeList(self, nodeList):
        return nodeList.name() in [n.name() for n in self._nodeLists]

    #---------------------------------------------------------------------------
    # Build a header, as a list of strings.
    #---------------------------------------------------------------------------
    def constructHeader(self, simulationTime, cycle):

        # List to build up the result in.
        header = []

        # First put some comments in, describing the file.
        header.append("# Spheral++ Visit dump file %s\n" % self.version)
        header.append("# time  = %g\n" % simulationTime)
        header.append("# cycle = %i\n" % cycle)

        # Put in the begin header descriptor.
        header.append(self.keyWords["header"] + "\n")

        # Insert the time and cycle descriptors.
        header.append(self.keyWords["time"] + " %g\n" % simulationTime)
        header.append(self.keyWords["cycle"] + " %i\n" % cycle)

        # Now loop over each NodeList.
        for nodeList in self._nodeLists:

            # If so, then write out the NodeList descriptor, and any
            # associated Field descriptors.
            header.append(self.keyWords["nodelist"] + ' "' +
                          nodeList.name() + '"\n')
            for field in self.fields(nodeList):
                header.append(self.keyWords["field"] + ' "' +
                              self._processName(field.name()) + '" ' +
                              self.dataType[type(field)] + "\n")

        # Terminate the header lines.
        header.append(self.keyWords["endheader"] + "\n")

        # That's, return the list of header lines.
        return header

    #---------------------------------------------------------------------------
    # Return the list of Fields for the given NodeList.
    #---------------------------------------------------------------------------
    def fields(self, nodeList):
        result = [f for f in self._fields if nodeList.haveField(f)]
        if result != []:
            result = [nodeList.positions()] + result
        return result

    #---------------------------------------------------------------------------
    # Private method dump a field to the given file.
    #---------------------------------------------------------------------------
    def _dumpField(self, field, file):
        import string

        # Write the field description.
        assert type(field) in self.dataType.keys()
        file.write(self.keyWords["field"] +
                   ' "' +
                   self._processName(field.name()) +
                   '" ' +
                   self.dataType[type(field)] + "\n")

        # Get the appropriate write method for the field DataType.
        elementDumpMethod = self.elementDumper[type(field)]

        # Loop over the field elements, and dump them.
        if self.dumpGhosts:
            for element in field.allValues():
                elementDumpMethod(element, file)
        else:
            for element in field.internalValues():
                elementDumpMethod(element, file)

        return

    #---------------------------------------------------------------------------
    # Private method to process a Field name into a Visit palatable label
    #---------------------------------------------------------------------------
    def _processName(self, name):
        import string
        return string.replace(str(name), " ", "_")

    #---------------------------------------------------------------------------
    # Private method to write an atomic DataType element to a file.
    #---------------------------------------------------------------------------
    def _dumpFloat(self, element, file):
        file.write("%g\n" % element)
        return

    #---------------------------------------------------------------------------
    # Private method to write an array DataType element to a file.
    #---------------------------------------------------------------------------
    def _dumpArray(self, element, file):
        buf = ""
        for x in element:
            buf += "%g " % x
        buf += "\n"
        file.write(buf)
        return

    #---------------------------------------------------------------------------
    # Private method to write a Vector element to a file.
    #---------------------------------------------------------------------------
    def _dumpVector(self, element, file):
        buf = ""
        for i in xrange(element.nDimensions):
            buf += "%g " % element(i)
        buf += "\n"
        file.write(buf)
        return

    #---------------------------------------------------------------------------
    # Private method to write a Tensor element to a file.
    #---------------------------------------------------------------------------
    def _dumpTensor(self, element, file):
        buf = ""
        for i in xrange(element.nDimensions):
            for j in xrange(element.nDimensions):
                buf += "%g " % element(i,j)
        buf += "\n"
        file.write(buf)
        return

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
                     dumpDerivatives = False):

    # What did we get passed?
    if max([isinstance(stateThingy, x) for x in [Spheral.Integrator1d, Spheral.Integrator2d, Spheral.Integrator3d]]):
        integrator = stateThingy
        dataBase = integrator.dataBase()
        state = eval("Spheral.State%id(integrator.dataBase(), integrator.physicsPackages())" % integrator.dataBase().nDim)
        derivs = None
        if dumpDerivatives:
            derivs = eval("Spheral.StateDerivatives%id(integrator.dataBase(), integrator.physicsPackages())" % integrator.dataBase().nDim)
        currentTime = integrator.currentTime
        currentCycle = integrator.currentCycle
    elif max([isinstance(stateThingy, x) for x in [Spheral.State1d, Spheral.State2d, Spheral.State3d]]):
        integrator = None
        dataBase = None
        state = stateThingy
        derivs = None
        assert currentTime is not None
        assert currentCycle is not None
    assert state is not None

    # Did the user specify any data to be dumped?
    if fields is None:
        fields = []
    if fieldLists is None:
        fieldLists = []

    # Build up the list of field lists in the state object.
    fieldLists += [state.scalarFields(name) for name in state.scalarFieldNames()]
    fieldLists += [state.vectorFields(name) for name in state.vectorFieldNames()]
    fieldLists += [state.tensorFields(name) for name in state.tensorFieldNames()]
    fieldLists += [state.symTensorFields(name) for name in state.symTensorFieldNames()]

    # Are we also dumping the derivative fields?
    if not derivs is None:
        fieldLists += [derivs.scalarFields(name) for name in derivs.scalarFieldNames()]
        fieldLists += [derivs.vectorFields(name) for name in derivs.vectorFieldNames()]
        fieldLists += [derivs.tensorFields(name) for name in derivs.tensorFieldNames()]
        fieldLists += [derivs.symTensorFields(name) for name in derivs.symTensorFieldNames()]

    # If available, add the work, H inverse and hmin, hmax, & hmin_hmax_ratio by default.
    if dataBase:
        work = dataBase.globalWork
        fieldLists.append(work)
        Hi = dataBase.fluidHinverse
        fieldLists.append(Hi)
        hmin = dataBase.newGlobalScalarFieldList()
        hmax = dataBase.newGlobalScalarFieldList()
        hminhmax = dataBase.newGlobalScalarFieldList()
        for fmin, fmax, fratio in zip(hmin.fields(),
                                      hmax.fields(),
                                      hminhmax.fields()):
            fmin.setName("hmin")
            fmax.setName("hmax")
            fratio.setName("hmin_hmax_ratio")
            n = fmin.nodeListPtr()
            H = n.Hfield()
            for i in xrange(n.numInternalNodes):
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
        for f in domains.fields():
            f.setName("Domains")
            for i in xrange(f.nodeList().numInternalNodes):
                f[i] = mpi.rank
        fieldLists.append(domains)
    except:
        pass

    # Now build the visit dumper.
    dumper = SpheralVisitDump(baseFileName,
                              baseDirectory,
                              fields,
                              fieldLists,
                              dumpGhosts = dumpGhosts)

    # Dump the sucker.
    dumper.dump(currentTime, currentCycle)
    return

#-------------------------------------------------------------------------------
# Dump out all the Fields in a State object by sampling to a regular lattice.
# You can pass any of the following for stateThingy:
#     Integrator
#     State
#-------------------------------------------------------------------------------
def dumpPhysicsState2Lattice(stateThingy,
                             nsample,
                             xmin,
                             xmax,
                             W,
                             baseFileName,
                             baseDirectory = ".",
                             mask = None,
                             currentTime = None,
                             currentCycle = None,
                             binary = True,
                             dumpGhosts = True,
                             scalarNames = "all",
                             vectorNames = "all",
                             tensorNames = [],
                             symTensorNames = []):

    assert (isinstance(stateThingy, Spheral.Integrator2d) or
            isinstance(stateThingy, Spheral.Integrator3d) or
            isinstance(stateThingy, Spheral.State2d) or
            isinstance(stateThingy, Spheral.State3d))

    import mpi

    # Prepare to time how long this takes.
    t0 = time.clock()

    # What did we get passed?
    if max([isinstance(stateThingy, x) for x in [Spheral.Integrator2d, Spheral.Integrator3d]]):
        integrator = stateThingy
        dataBase = integrator.dataBase()
        state = eval("Spheral.State%id(integrator.dataBase(), integrator.physicsPackages())" % integrator.dataBase().nDim)
        currentTime = integrator.currentTime
        currentCycle = integrator.currentCycle
    elif max([isinstance(stateThingy, x) for x in [Spheral.State2d, Spheral.State3d]]):
        integrator = None
        dataBase = None
        state = stateThingy
        assert currentTime is not None
        assert currentCycle is not None
    assert state is not None

    # Determine the dimensionality.
    if isinstance(state, Spheral.State2d):
        nDim = 2
        fieldListSet = Spheral.FieldListSet2d()
        nsample = (nsample[0], nsample[1], 1)
        xmin = Spheral.Vector3d(xmin.x, xmin.y, 0.0)
        xmax = Spheral.Vector3d(xmax.x, xmax.y, 0.0)
    else:
        nDim = 3
        fieldListSet = Spheral.FieldListSet3d()

    # If requested, set ghost node info.
    if dumpGhosts and not integrator is None:
        derivs = eval("Spheral.StateDerivatives%id(integrator.dataBase(), integrator.physicsPackages())" % nDim)
        integrator.setGhostNodes()
        integrator.applyGhostBoundaries(state, derivs)

    # Determine the set of field names we're going to extract.
    if scalarNames == "all":
        scalarNames = state.scalarFieldNames()
    if vectorNames == "all":
        vectorNames = state.vectorFieldNames()
    if tensorNames == "all":
        tensorNames = state.tensorFieldNames()
    if symTensorNames == "all":
        symTensorNames = state.symTensorFieldNames()

    # Check that the requested names actually exist.
    for name in scalarNames:
        assert name in state.scalarFieldNames()
    for name in vectorNames:
        assert name in state.vectorFieldNames()
    for name in tensorNames:
        assert name in state.tensorFieldNames()
    for name in symTensorNames:
        assert name in state.symTensorFieldNames()

    # Build up the list of field lists in the state object.
    [fieldListSet.ScalarFieldLists.append(state.scalarFields(name)) for name in scalarNames]
    [fieldListSet.VectorFieldLists.append(state.vectorFields(name)) for name in vectorNames]
    [fieldListSet.TensorFieldLists.append(state.tensorFields(name)) for name in tensorNames]
    [fieldListSet.SymTensorFieldLists.append(state.symTensorFields(name)) for name in symTensorNames]

    # Build the sampling coordinates.
    assert len(nsample) == 3
    assert min(nsample) > 0
    assert xmax > xmin
    dx = (xmax.x - xmin.x)/nsample[0]
    dy = (xmax.y - xmin.y)/nsample[1]
    dz = (xmax.z - xmin.z)/nsample[2]
    coords = ([xmin.x + i*dx for i in xrange(nsample[0] + 1)],
              [xmin.y + i*dy for i in xrange(nsample[1] + 1)],
              [xmin.z + i*dz for i in xrange(nsample[2] + 1)])
    
    # Get the state fields.
    r = state.vectorFields("position")
    w = state.scalarFields("weight")
    H = state.symTensorFields("H")

    # Did the user provide a mask?
    if mask is None:
        mask = eval("Spheral.IntFieldList%id()" % nDim)
        mask.copyFields()
        for field in r.fields():
            mask.appendField(eval("Spheral.Intfield%id('mask', field.nodeList(), 1)" % nDim))

    # Make sure the mask covers all our fields!
    for field in r.fields():
        assert(mask.haveNodeList(field.nodeList()))

    # Sample the fields.
    if nDim == 2:
        sampledFields = Spheral.sampleMultipleFields2LatticeMash(fieldListSet,
                                                                 r,
                                                                 w,
                                                                 H,
                                                                 mask,
                                                                 W,
                                                                 Spheral.Vector2d(xmin.x, xmin.y),
                                                                 Spheral.Vector2d(xmax.x, xmax.y),
                                                                 (nsample[0], nsample[1]))
    else:
        sampledFields = Spheral.sampleMultipleFields2LatticeMash(fieldListSet,
                                                                 r,
                                                                 w,
                                                                 H,
                                                                 mask,
                                                                 W,
                                                                 xmin,
                                                                 xmax,
                                                                 nsample)

    # The above sampling method returns the result distributed over all processors.
    # For the Visit dumping routine below to work though, we need to stitch the
    # fields back together.  Note this won't work if there are too many grid points!
    ntot = nsample[0]*nsample[1]*nsample[2]
    for name, field in (sampledFields[0] +
                        sampledFields[1] +
                        sampledFields[2] +
                        sampledFields[3]):
        if mpi.rank == 0:
            for sendProc in xrange(1, mpi.procs):
                vals = mpi.recv(sendProc)[0]
                print "Received %i values from processor %i" % (len(vals), sendProc)
                field.extend(vals)
        else:
            mpi.send(list(field), 0)
        if mpi.rank == 0:
            print name, len(field), ntot
            assert len(field) == ntot

    # We have to remove the spaces from field names before passing them
    # to visit.
    thpt = ([(name.replace(" ", "_"), val) for (name, val) in sampledFields[0]],
            [(name.replace(" ", "_"), val) for (name, val) in sampledFields[1]],
            [(name.replace(" ", "_"), val) for (name, val) in sampledFields[2]],
            [(name.replace(" ", "_"), val) for (name, val) in sampledFields[3]])
    sampledFields = thpt

    # Make sure the output directory exists.
    import mpi
    import os
    if mpi.rank == 0 and not os.path.exists(baseDirectory):
        try:
            os.makedirs(baseDirectory)
        except:
            raise "Cannot create output directory %s" % baseDirectory
    mpi.barrier()

    # Dump the sucker.
    baseName = "%s-time=%g-cycle=%i" % (baseFileName, currentTime, currentCycle)
    filename = os.path.join(baseDirectory, baseName)
    mastername = os.path.join(baseDirectory, baseFileName + "-samplemaster.visit")
    print "Preparing to write %i fields" % (len(sampledFields[0]) +
                                            len(sampledFields[1]) +
                                            len(sampledFields[2]) +
                                            len(sampledFields[3]))
    if mpi.rank == 0:
        Spheral.writeRectilinearMesh(filename,
                                     binary,
                                     nsample,
                                     coords,
                                     sampledFields)

        # Add this file to the master file.
        mf = open(mastername, "a")
        mf.write("%s.vtk\n" % baseName)
        mf.close()

    mpi.barrier()
    print "Finished: required %0.2f seconds" % (time.clock() - t0)

    return
