# ------------------------------------------------------------------------------
# Make modified versions of the FileIO objects extended to handle arbitrary
# python types.
# ------------------------------------------------------------------------------
from FileIO import FlatFileIO as CXXFlatFileIO
from FileIO import FlatFileFormat, AccessType
import pickle, string

class FlatFileIO(CXXFlatFileIO):

    def __init__(self,
                 access,
                 format = FlatFileFormat.ascii):
        CXXFlatFileIO.__init__(self, access, format)
        return

    def writeObject(self, object, pathName):
        self.write(string.replace(pickle.dumps(object), "\n", "\\n"), pathName)
        return

    def readObject(self, pathName):
        stringvar = string.replace(self.readString(pathName), "\\n", "\n")
        return pickle.loads(stringvar)
