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
                 precision = 20):
        PyFileIO.__init__(self, fileName, access)
        self.encoding = "utf-8"
        self.terminator = b'\0'
        self.mungenewline = b'<><><>n<><><>'

        # We enforce the convention that gziped files will have the .gz
        # extension.
        if fileName[-3:] != ".gz":
            fileName += ".gz"
        
        self.precision = precision # "%" + "%i.%ie" % (precision + 3, precision)
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
            mode = "wb"
        elif access == Read:
            mode = "rb"
        assert not mode is None
        self.f = gzip.open(fileName, mode=mode)

        # If requested, read the file to memory.
        if access == Read:
            self.lines = {}
            for line in self.f:
                #print("read line: ", line)
                assert line[0:1] == self.terminator
                assert line[-2:-1] == self.terminator
                i = line[1:].index(self.terminator)
                assert i < len(line)
                key = str(line[1:i+1], self.encoding)
                val = line[i+2:-2].replace(self.mungenewline, b'\n')
                self.lines[key] = val
            #print("Final self.lines: \n", self.lines)
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

    # #---------------------------------------------------------------------------
    # # Support for writing and reading Fields.
    # #---------------------------------------------------------------------------
    # def writeFieldObject(self, val, pathName):
    #     self.writeObject(val.string(self.precision), pathName)
    #     return

    # def readFieldObject(self, val, pathName):
    #     val.string(self.readObject(pathName))
    #     #print list(val)
    #     return

    #---------------------------------------------------------------------------
    # Scan the current file and return the line corresponding to the given
    # path.  If it doesn't exist, raise an error.
    #---------------------------------------------------------------------------
    def findPath(self, pathName):
        if not pathName in self.lines:
            raise ValueError("GzipFileIO unknown key: " + pathName)
        #print("findPath: ", pathName, self.lines[pathName])
        return self.lines[pathName]

    #---------------------------------------------------------------------------
    # writeBytes
    #---------------------------------------------------------------------------
    def writeBytes(self, val, pathName):
        #print("writeBytes: ", pathName)
        #print("writeBytes: ", val)
        self.f.write(self.terminator +
                     bytes(pathName, self.encoding) + self.terminator + 
                     val.replace(b'\n', self.mungenewline) + self.terminator + b'\n')
        return

    def readBytes(self, pathName):
        return self.findPath(pathName)

    # #---------------------------------------------------------------------------
    # # pathExists
    # #---------------------------------------------------------------------------
    # def pathExists(self, pathName):
    #     try:
    #         p = self.findPath(pathName)
    #         return True
    #     except:
    #         return False

    # #---------------------------------------------------------------------------
    # # _write: the generic write using pickle
    # #---------------------------------------------------------------------------
    # def _write(self, val, pathName):
    #     self.write_string(pickle.dumps(val).replace(b"\n", b"<<<<<n>>>>>"),
    #                       pathName)

    # #---------------------------------------------------------------------------
    # # _read: the generic read using pickle
    # #---------------------------------------------------------------------------
    # def _read(self, pathName):
    #     stuff = self.readBytes(pathName)
    #     return pickle.loads(stuff.replace(b"<<<<<n>>>>>", b"\n"))

    # #---------------------------------------------------------------------------
    # # Use pickling for the majority of the write methods.  Most objects we just
    # # convert to strings and pickle that.
    # #---------------------------------------------------------------------------
    # def write_unsigned_int(self, val, pathName):
    #     self._write(val, pathName)

    # def write_size_t(self, val, pathName):
    #     self._write(val, pathName)

    # def write_int(self, val, pathName):
    #     self._write(val, pathName)

    # def write_bool(self, val, pathName):
    #     self._write(val, pathName)

    # def write_double(self, val, pathName):
    #     self._write(val, pathName)

    # def write_vector_int(self, val, pathName):
    #     self._write(list(val), pathName)

    # def write_vector_double(self, val, pathName):
    #     self._write(list(val), pathName)

    # def write_vector_string(self, val, pathName):
    #     self._write(list(val), pathName)

    # def write_Vector1d(self, val, pathName):
    #     self._write(val, pathName)

    # def write_Tensor1d(self, val, pathName):
    #     self._write(val, pathName)

    # def write_SymTensor1d(self, val, pathName):
    #     self._write(val, pathName)

    # def write_ThirdRankTensor1d(self, val, pathName):
    #     self._write(val, pathName)

    # def write_Vector2d(self, val, pathName):
    #     self._write(val, pathName)

    # def write_Tensor2d(self, val, pathName):
    #     self._write(val, pathName)

    # def write_SymTensor2d(self, val, pathName):
    #     self._write(val, pathName)

    # def write_ThirdRankTensor2d(self, val, pathName):
    #     self._write(val, pathName)

    # def write_Vector3d(self, val, pathName):
    #     self._write(val, pathName)

    # def write_Tensor3d(self, val, pathName):
    #     self._write(val, pathName)

    # def write_SymTensor3d(self, val, pathName):
    #     self._write(val, pathName)

    # def write_ThirdRankTensor3d(self, val, pathName):
    #     self._write(val, pathName)

    # def write_ScalarField1d(self, val, pathName):
    #     self.writeFieldObject(val, pathName)

    # def write_VectorField1d(self, val, pathName):
    #     self.writeFieldObject(val, pathName)

    # def write_TensorField1d(self, val, pathName):
    #     self.writeFieldObject(val, pathName)

    # def write_SymTensorField1d(self, val, pathName):
    #     self.writeFieldObject(val, pathName)

    # def write_ThirdRankTensorField1d(self, val, pathName):
    #     self.writeFieldObject(val, pathName)

    # def write_IntField1d(self, val, pathName):
    #     self.writeFieldObject(val, pathName)

    # def write_UnsignedField1d(self, val, pathName):
    #     self.writeFieldObject(val, pathName)

    # def write_ScalarField2d(self, val, pathName):
    #     self.writeFieldObject(val, pathName)

    # def write_VectorField2d(self, val, pathName):
    #     self.writeFieldObject(val, pathName)

    # def write_TensorField2d(self, val, pathName):
    #     self.writeFieldObject(val, pathName)

    # def write_SymTensorField2d(self, val, pathName):
    #     self.writeFieldObject(val, pathName)

    # def write_ThirdRankTensorField2d(self, val, pathName):
    #     self.writeFieldObject(val, pathName)

    # def write_IntField2d(self, val, pathName):
    #     self.writeFieldObject(val, pathName)

    # def write_UnsignedField2d(self, val, pathName):
    #     self.writeFieldObject(val, pathName)

    # def write_ScalarField3d(self, val, pathName):
    #     self.writeFieldObject(val, pathName)

    # def write_VectorField3d(self, val, pathName):
    #     self.writeFieldObject(val, pathName)

    # def write_TensorField3d(self, val, pathName):
    #     self.writeFieldObject(val, pathName)

    # def write_SymTensorField3d(self, val, pathName):
    #     self.writeFieldObject(val, pathName)

    # def write_ThirdRankTensorField3d(self, val, pathName):
    #     self.writeFieldObject(val, pathName)

    # def write_IntField3d(self, val, pathName):
    #     self.writeFieldObject(val, pathName)

    # def write_UnsignedField3d(self, val, pathName):
    #     self.writeFieldObject(val, pathName)

    # #---------------------------------------------------------------------------
    # # We now use unpickling to read objects.
    # #---------------------------------------------------------------------------
    # def read_unsigned_int(self, pathName):
    #     return self._read(pathName)

    # def read_size_t(self, pathName):
    #     return self._read(pathName)

    # def read_int(self, pathName):
    #     return self._read(pathName)

    # def read_bool(self, pathName):
    #     return self._read(pathName)

    # def read_double(self, pathName):
    #     return self._read(pathName)

    # def read_vector_int(self, val, pathName):
    #     self.copyContainer(self._read(pathName), val)

    # def read_vector_double(self, val, pathName):
    #     self.copyContainer(self._read(pathName), val)

    # def read_vector_string(self, val, pathName):
    #     self.copyContainer(self._read(pathName), val)

    # def read_Vector1d(self, pathName):
    #     return self._read(pathName)

    # def read_Tensor1d(self, pathName):
    #     return self._read(pathName)

    # def read_SymTensor1d(self, pathName):
    #     return self._read(pathName)

    # def read_ThirdRankTensor1d(self, pathName):
    #     return self._read(pathName)

    # def read_Vector2d(self, pathName):
    #     return self._read(pathName)

    # def read_Tensor2d(self, pathName):
    #     return self._read(pathName)

    # def read_SymTensor2d(self, pathName):
    #     return self._read(pathName)

    # def read_ThirdRankTensor2d(self, pathName):
    #     return self._read(pathName)

    # def read_Vector3d(self, pathName):
    #     return self._read(pathName)

    # def read_Tensor3d(self, pathName):
    #     return self._read(pathName)

    # def read_SymTensor3d(self, pathName):
    #     return self._read(pathName)

    # def read_ThirdRankTensor3d(self, pathName):
    #     return self._read(pathName)

    # def read_ScalarField1d(self, val, pathName):
    #     self.readFieldObject(val, pathName)

    # def read_VectorField1d(self, val, pathName):
    #     self.readFieldObject(val, pathName)

    # def read_TensorField1d(self, val, pathName):
    #     self.readFieldObject(val, pathName)

    # def read_SymTensorField1d(self, val, pathName):
    #     self.readFieldObject(val, pathName)

    # def read_ThirdRankTensorField1d(self, val, pathName):
    #     self.readFieldObject(val, pathName)

    # def read_IntField1d(self, val, pathName):
    #     self.readFieldObject(val, pathName)

    # def read_UnsignedField1d(self, val, pathName):
    #     self.readFieldObject(val, pathName)

    # def read_ScalarField2d(self, val, pathName):
    #     self.readFieldObject(val, pathName)

    # def read_VectorField2d(self, val, pathName):
    #     self.readFieldObject(val, pathName)

    # def read_TensorField2d(self, val, pathName):
    #     self.readFieldObject(val, pathName)

    # def read_SymTensorField2d(self, val, pathName):
    #     self.readFieldObject(val, pathName)

    # def read_ThirdRankTensorField2d(self, val, pathName):
    #     self.readFieldObject(val, pathName)

    # def read_IntField2d(self, val, pathName):
    #     self.readFieldObject(val, pathName)

    # def read_UnsignedField2d(self, val, pathName):
    #     self.readFieldObject(val, pathName)

    # def read_ScalarField3d(self, val, pathName):
    #     self.readFieldObject(val, pathName)

    # def read_VectorField3d(self, val, pathName):
    #     self.readFieldObject(val, pathName)

    # def read_TensorField3d(self, val, pathName):
    #     self.readFieldObject(val, pathName)

    # def read_SymTensorField3d(self, val, pathName):
    #     self.readFieldObject(val, pathName)

    # def read_ThirdRankTensorField3d(self, val, pathName):
    #     self.readFieldObject(val, pathName)

    # def read_IntField3d(self, val, pathName):
    #     self.readFieldObject(val, pathName)

    # def read_UnsignedField3d(self, val, pathName):
    #     self.readFieldObject(val, pathName)

    # #---------------------------------------------------------------------------
    # # Copy Vectors.
    # #---------------------------------------------------------------------------
    # def copyContainer(self, v0, v1):
    #     for x in v0:
    #         v1.append(x)
    #     return
