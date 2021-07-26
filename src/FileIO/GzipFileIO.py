#-------------------------------------------------------------------------------
# GzipFileIO.
#
# An implementation of FileIO to read/write Spheral++ data in gziped file
# format.  Uses the Python gzip module.
#-------------------------------------------------------------------------------
import struct
import gzip
import pickle
import time
import sys

from SpheralCompiledPackages import *

#-------------------------------------------------------------------------------
# Generic class definition (all dimensions).
#-------------------------------------------------------------------------------
class GzipFileIO(PyFileIO):

    #---------------------------------------------------------------------------
    # Constructor.
    #---------------------------------------------------------------------------
    def __init__(self,
                 fileName,
                 access,
                 precision = 20,
                 readToMemory = True):
        PyFileIO.__init__(self, fileName, access)

        # We enforce the convention that gziped files will have the .gz
        # extension.
        if fileName[-3:] != ".gz":
            fileName += ".gz"
        
        self.precision = precision # "%" + "%i.%ie" % (precision + 3, precision)
        self.readToMemory = readToMemory
        self.open(fileName, access)

        return

    #---------------------------------------------------------------------------
    # Open
    #---------------------------------------------------------------------------
    def open(self,
             fileName,
             access):

        # Open up a file with the requested access pattern.
        mode = None
        if (access == Create or
            access == Write):
            mode = "w"
        elif access == Read:
            mode = "r"
        assert not mode is None
        self.f = gzip.GzipFile(fileName, mode)

        # If requested, read the file to memory.
        if access == Read and self.readToMemory:
            self.lines = {}
            for line in self.f:
                i = line.index(self.terminator())
                assert i < len(line)
                self.lines[line[:i+1]] = line[i+1:-1]

        return

    #---------------------------------------------------------------------------
    # Close
    #---------------------------------------------------------------------------
    def close(self, printTimes = False):
        self.f.close()
        if hasattr(self, "lines"):
            del self.lines
        return

##     #---------------------------------------------------------------------------
##     # Support for writing pickleable Python objects.
##     #---------------------------------------------------------------------------
##     def writeObject(self, object, pathName):
##         self.write(pickle.dumps(object).replace("\n", "\\n"), pathName)
##         return

##     #---------------------------------------------------------------------------
##     # Support for reading pickleable Python objects.
##     #---------------------------------------------------------------------------
##     def readObject(self, pathName):
##         stringvar = self.findPath(pathName).replace("\\n", "\n")
##         return pickle.loads(stringvar)

    #---------------------------------------------------------------------------
    # Support for writing and reading Fields.
    #---------------------------------------------------------------------------
    def writeFieldObject(self, val, pathName):
        self.writeObject(val.string(self.precision), pathName)
        return

    def readFieldObject(self, val, pathName):
        val.string(self.readObject(pathName))
        #print list(val)
        return

    #---------------------------------------------------------------------------
    # Scan the current file and return the line corresponding to the given
    # path.  If it doesn't exist, raise an error.
    #---------------------------------------------------------------------------
    def findPath(self, pathName):
        # Convert the pathName to a string (seems silly, but this hides the
        # details of binary IO).
        #pathstring = self._string2string(pathName) + self.terminator()
        pathstring = str(pathName) + self.terminator()
        npath = len(pathstring)

        # If we read everything to memory, just scan that.
        if self.readToMemory and pathstring in self.lines:
            return self.lines[pathstring]

        else:
            # Go to the beginning of the file.
            self.f.seek(0)

            # Iterate over the lines in the file until we find what we want.
            for line in self.f:
                if line[:npath] == pathstring:
                    return line[npath:-1]

        # Uh oh!  We didn't find the requested path.  Raise an error.
        raise ValueError, "GzipFileIO.findPath : unable to locate %s in %s!" % (pathName,
                                                                                self.fileName)

    #---------------------------------------------------------------------------
    # Standard terminator for encoded strings.
    #---------------------------------------------------------------------------
    def terminator(self):
#         if self.binary:
#             return struct.pack("s", "\0")
#         else:
            return "\0"

    #---------------------------------------------------------------------------
    # We have to actually provide the string write and read methods, since these
    # are used by the writeObject/readObject pickling approach used by the other
    # methods.
    #---------------------------------------------------------------------------
    def write_string(self, val, pathName):
        pathString = str(pathName) + self.terminator()
        valString = str(val) + "\n"
        self.f.write(pathString + valString)
        del valString
        return

    def read_string(self, pathName):
        try:
            result = self.findPath(pathName)
            return result
        except Exception as excp:
            print "WARNING : Unable to restore %s due to exception message: %s" % (pathName, excp)
            pass

    #---------------------------------------------------------------------------
    # pathExists
    #---------------------------------------------------------------------------
    def pathExists(self, pathName):
        try:
            p = self.findPath(pathName)
            return True
        except:
            return False

    #---------------------------------------------------------------------------
    # Use pickling for the majority of the write methods.  Most objects we just
    # convert to strings and pickle that.
    #---------------------------------------------------------------------------
    def write_unsigned_int(self, val, pathName):
        self.writeObject(val, pathName)

    def write_size_t(self, val, pathName):
        self.writeObject(val, pathName)

    def write_int(self, val, pathName):
        self.writeObject(val, pathName)

    def write_bool(self, val, pathName):
        self.writeObject(val, pathName)

    def write_double(self, val, pathName):
        self.writeObject(val, pathName)

    def write_vector_int(self, val, pathName):
        self.writeObject(list(val), pathName)

    def write_vector_double(self, val, pathName):
        self.writeObject(list(val), pathName)

    def write_vector_string(self, val, pathName):
        self.writeObject(list(val), pathName)

    def write_Vector1d(self, val, pathName):
        self.writeObject(val, pathName)

    def write_Tensor1d(self, val, pathName):
        self.writeObject(val, pathName)

    def write_SymTensor1d(self, val, pathName):
        self.writeObject(val, pathName)

    def write_ThirdRankTensor1d(self, val, pathName):
        self.writeObject(val, pathName)

    def write_Vector2d(self, val, pathName):
        self.writeObject(val, pathName)

    def write_Tensor2d(self, val, pathName):
        self.writeObject(val, pathName)

    def write_SymTensor2d(self, val, pathName):
        self.writeObject(val, pathName)

    def write_ThirdRankTensor2d(self, val, pathName):
        self.writeObject(val, pathName)

    def write_Vector3d(self, val, pathName):
        self.writeObject(val, pathName)

    def write_Tensor3d(self, val, pathName):
        self.writeObject(val, pathName)

    def write_SymTensor3d(self, val, pathName):
        self.writeObject(val, pathName)

    def write_ThirdRankTensor3d(self, val, pathName):
        self.writeObject(val, pathName)

    def write_ScalarField1d(self, val, pathName):
        self.writeFieldObject(val, pathName)

    def write_VectorField1d(self, val, pathName):
        self.writeFieldObject(val, pathName)

    def write_TensorField1d(self, val, pathName):
        self.writeFieldObject(val, pathName)

    def write_SymTensorField1d(self, val, pathName):
        self.writeFieldObject(val, pathName)

    def write_ThirdRankTensorField1d(self, val, pathName):
        self.writeFieldObject(val, pathName)

    def write_IntField1d(self, val, pathName):
        self.writeFieldObject(val, pathName)

    def write_UnsignedField1d(self, val, pathName):
        self.writeFieldObject(val, pathName)

    def write_ScalarField2d(self, val, pathName):
        self.writeFieldObject(val, pathName)

    def write_VectorField2d(self, val, pathName):
        self.writeFieldObject(val, pathName)

    def write_TensorField2d(self, val, pathName):
        self.writeFieldObject(val, pathName)

    def write_SymTensorField2d(self, val, pathName):
        self.writeFieldObject(val, pathName)

    def write_ThirdRankTensorField2d(self, val, pathName):
        self.writeFieldObject(val, pathName)

    def write_IntField2d(self, val, pathName):
        self.writeFieldObject(val, pathName)

    def write_UnsignedField2d(self, val, pathName):
        self.writeFieldObject(val, pathName)

    def write_ScalarField3d(self, val, pathName):
        self.writeFieldObject(val, pathName)

    def write_VectorField3d(self, val, pathName):
        self.writeFieldObject(val, pathName)

    def write_TensorField3d(self, val, pathName):
        self.writeFieldObject(val, pathName)

    def write_SymTensorField3d(self, val, pathName):
        self.writeFieldObject(val, pathName)

    def write_ThirdRankTensorField3d(self, val, pathName):
        self.writeFieldObject(val, pathName)

    def write_IntField3d(self, val, pathName):
        self.writeFieldObject(val, pathName)

    def write_UnsignedField3d(self, val, pathName):
        self.writeFieldObject(val, pathName)

    #---------------------------------------------------------------------------
    # We now use unpickling to read objects.
    #---------------------------------------------------------------------------
    def read_unsigned_int(self, pathName):
        return self.readObject(pathName)

    def read_size_t(self, pathName):
        return self.readObject(pathName)

    def read_int(self, pathName):
        return self.readObject(pathName)

    def read_bool(self, pathName):
        return self.readObject(pathName)

    def read_double(self, pathName):
        return self.readObject(pathName)

    def read_vector_int(self, val, pathName):
        self.copyContainer(self.readObject(pathName), val)

    def read_vector_double(self, val, pathName):
        self.copyContainer(self.readObject(pathName), val)

    def read_vector_string(self, val, pathName):
        self.copyContainer(self.readObject(pathName), val)

    def read_Vector1d(self, pathName):
        return self.readObject(pathName)

    def read_Tensor1d(self, pathName):
        return self.readObject(pathName)

    def read_SymTensor1d(self, pathName):
        return self.readObject(pathName)

    def read_ThirdRankTensor1d(self, pathName):
        return self.readObject(pathName)

    def read_Vector2d(self, pathName):
        return self.readObject(pathName)

    def read_Tensor2d(self, pathName):
        return self.readObject(pathName)

    def read_SymTensor2d(self, pathName):
        return self.readObject(pathName)

    def read_ThirdRankTensor2d(self, pathName):
        return self.readObject(pathName)

    def read_Vector3d(self, pathName):
        return self.readObject(pathName)

    def read_Tensor3d(self, pathName):
        return self.readObject(pathName)

    def read_SymTensor3d(self, pathName):
        return self.readObject(pathName)

    def read_ThirdRankTensor3d(self, pathName):
        return self.readObject(pathName)

    def read_ScalarField1d(self, val, pathName):
        self.readFieldObject(val, pathName)

    def read_VectorField1d(self, val, pathName):
        self.readFieldObject(val, pathName)

    def read_TensorField1d(self, val, pathName):
        self.readFieldObject(val, pathName)

    def read_SymTensorField1d(self, val, pathName):
        self.readFieldObject(val, pathName)

    def read_ThirdRankTensorField1d(self, val, pathName):
        self.readFieldObject(val, pathName)

    def read_IntField1d(self, val, pathName):
        self.readFieldObject(val, pathName)

    def read_UnsignedField1d(self, val, pathName):
        self.readFieldObject(val, pathName)

    def read_ScalarField2d(self, val, pathName):
        self.readFieldObject(val, pathName)

    def read_VectorField2d(self, val, pathName):
        self.readFieldObject(val, pathName)

    def read_TensorField2d(self, val, pathName):
        self.readFieldObject(val, pathName)

    def read_SymTensorField2d(self, val, pathName):
        self.readFieldObject(val, pathName)

    def read_ThirdRankTensorField2d(self, val, pathName):
        self.readFieldObject(val, pathName)

    def read_IntField2d(self, val, pathName):
        self.readFieldObject(val, pathName)

    def read_UnsignedField2d(self, val, pathName):
        self.readFieldObject(val, pathName)

    def read_ScalarField3d(self, val, pathName):
        self.readFieldObject(val, pathName)

    def read_VectorField3d(self, val, pathName):
        self.readFieldObject(val, pathName)

    def read_TensorField3d(self, val, pathName):
        self.readFieldObject(val, pathName)

    def read_SymTensorField3d(self, val, pathName):
        self.readFieldObject(val, pathName)

    def read_ThirdRankTensorField3d(self, val, pathName):
        self.readFieldObject(val, pathName)

    def read_IntField3d(self, val, pathName):
        self.readFieldObject(val, pathName)

    def read_UnsignedField3d(self, val, pathName):
        self.readFieldObject(val, pathName)

    #---------------------------------------------------------------------------
    # Copy Vectors.
    #---------------------------------------------------------------------------
    def copyContainer(self, v0, v1):
        for x in v0:
            v1.append(x)
        return
