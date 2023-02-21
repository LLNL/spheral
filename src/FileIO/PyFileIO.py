#-------------------------------------------------------------------------------
# PyFileIO
# Base class for Python FileIO implementations.  Assumes you just want to
# serialize things using pickle->bytes.  Python descendents need to provide:
#   write_bytes(self, xbytes, path)
#   read_bytes(self, path)
#-------------------------------------------------------------------------------
from SpheralCompiledPackages import *

class PyFileIO(FileIO):

    #---------------------------------------------------------------------------
    # Constructor.
    #---------------------------------------------------------------------------
    def __init__(self, fileName, access):
        FileIO.__init__(self, fileName, access)
        return

    #---------------------------------------------------------------------------
    # write
    #---------------------------------------------------------------------------
    def write(self, x, path):
        self.write_object(x, path)
        return

    #---------------------------------------------------------------------------
    # read
    #---------------------------------------------------------------------------
    def read(self, path):
        return self.read_object(path)
