#-------------------------------------------------------------------------------
# DumpGzipFileNodeGenerator
#
# Write out the state from a set of NodeLists to a file that is parsable by our
# GzipFileNodeGenerator.
#-------------------------------------------------------------------------------
from math import *
import mpi, gzip
import Spheral

class DumpGzipFileNodeGenerator:

    #---------------------------------------------------------------------------
    # Constructor.
    #---------------------------------------------------------------------------
    def __init__(self,
                 nodeLists,
                 filename,
                 precision = 20,
                 serialize = False,
                 extraFields = [],
                 ):
        self.nodeLists = nodeLists
        self.filename = filename
        self.precision = "%" + "%i.%ie" % (precision + 3, precision)
        self.serialize = serialize
        self.extraFields = extraFields
        self.delimiter = "$"
        self.writing = (mpi.rank == 0 or not serialize)

        # Some mappings to help us write out various data types.
        self._toString = {int                         : self._int2String,
                          bool                        : self._bool2String,
                          float                       : self._float2String,
                          str                         : self._str2String,
                          type(Spheral.Vector1d())    : self._ContainerFloats2String,
                          type(Spheral.Tensor1d())    : self._ContainerFloats2String,
                          type(Spheral.SymTensor1d()) : self._ContainerFloats2String,
                          type(Spheral.Vector2d())    : self._ContainerFloats2String,
                          type(Spheral.Tensor2d())    : self._ContainerFloats2String,
                          type(Spheral.SymTensor2d()) : self._ContainerFloats2String,
                          type(Spheral.Vector3d())    : self._ContainerFloats2String,
                          type(Spheral.Tensor3d())    : self._ContainerFloats2String,
                          type(Spheral.SymTensor3d()) : self._ContainerFloats2String,
                          }

        # Open the file.
        if self.writing:
            self.f = gzip.open(filename, "w")
        else:
            self.f = None

        if isinstance(nodeLists[0], Spheral.NodeList1d):
            SymTensorField = Spheral.SymTensorField1d
        elif isinstance(nodeLists[0], Spheral.NodeList2d):
            SymTensorField = Spheral.SymTensorField2d
        elif isinstance(nodeLists[0], Spheral.NodeList3d):
            SymTensorField = Spheral.SymTensorField3d

        # Walk the NodeLists and write the standard fields.
        for nodes in self.nodeLists:
            self.writeField(nodes.positions(), nodes.name, "positions")
            self.writeField(nodes.mass(), nodes.name, "mass")
            self.writeField(nodes.massDensity(), nodes.name, "density")
            self.writeField(nodes.velocity(), nodes.name, "velocity")
            self.writeField(nodes.specificThermalEnergy(), nodes.name, "specificThermalEnergy")
            Hinv2 = SymTensorField("Hinv2", nodes)
            nodes.Hinverse(Hinv2)
            for i in range(nodes.numInternalNodes):
                thpt = (Hinv2[i]*Hinv2[i]).Symmetric()
                Hinv2[i] = thpt
            self.writeField(Hinv2, nodes.name, "Hinverse2")

        # Add any extra fields requested.
        for field in self.extraFields:
            nodes = field.nodeList()
            self.writeField(field, nodes.name, field.name)

        return
    
    #---------------------------------------------------------------------------
    # Define our internal conversion functions.
    #---------------------------------------------------------------------------
    def _int2String(self, x):
        return str(x)

    def _bool2String(self, x):
        return str(x)

    def _float2String(self, x):
        return self.precision % x

    def _str2String(self, x):
        return x

    def _ContainerFloats2String(self, x):
        result = ""
        for xi in x:
            result += self.precision % xi + " "
        return result

    #---------------------------------------------------------------------------
    # Write a field in our format.
    #---------------------------------------------------------------------------
    def writeField(self, 
                   field,
                   materialName,
                   fieldName):
        print("Writing %s for %s." % (fieldName, materialName))
        vals = list(field.internalValues())
        n = len(vals)
        if self.serialize:
            n = mpi.allreduce(n, mpi.SUM)
        if self.writing:
            self.f.write(materialName + self.delimiter +
                         fieldName)
        if self.serialize:
            for sendProc in range(mpi.procs):
                print("Collecting from ", sendProc)
                otherVals = mpi.bcast(vals, sendProc)
                print("  Received %i values" % len(otherVals))
                if self.writing:
                    for x in otherVals:
                        self.f.write(self.delimiter + self._toString[type(x)](x))
        else:
            if self.writing:
                for x in vals:
                    self.f.write(self.delimiter + self._toString[type(x)](x))
        if self.writing:
            self.f.write("\n")
        return
