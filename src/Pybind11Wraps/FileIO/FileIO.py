#-------------------------------------------------------------------------------
# FileIO abstract class
#-------------------------------------------------------------------------------
from PYB11Generator import *
from spheralDimensions import *
dims = spheralDimensions()

class FileIO:

    #...........................................................................
    # Constructors
    def pyinit0(self):
        "Default constructor"

    def pyinit1(self,
                filename = "const std::string",
                access = "AccessType"):
        "Construct with a given file name and access"

    #...........................................................................
    # Abstract interface
    @PYB11pure_virtual
    def open(self,
             fileName = "const std::string",
             access = "AccessType"):
        "Open a file for IO"
        return "void"

    @PYB11pure_virtual
    def close(self):
        "Close the current file we're pointing at"
        return "void"

    #...........................................................................
    # Abstract interface (primitives)
    types = ["unsigned", "int", "bool", "double", "std::string"]
    for ndim in dims:
        types += ["Dim<%i>::Vector" % ndim,
                  "Dim<%i>::Tensor" % ndim,
                  "Dim<%i>::SymTensor" % ndim,
                  "Dim<%i>::ThirdRankTensor" % ndim]

    for T in types:
        exec("""
@PYB11pure_virtual
@PYB11pycppname("write")
def write%(Tmangle)s(self,
    value = "const %(T)s&",
    pathName = "const std::string"):
    "Write %(T)s"
    return "void"

@PYB11pure_virtual
@PYB11pycppname("read")
@PYB11const
def read%(Tmangle)s(self,
    value = "%(T)s&",
    pathName = "const std::string"):
    "Read %(T)s"
    return "void"
""" % {"T"       : T,
       "Tmangle" : T.replace(":", "_").replace("<", "_").replace(">", "_")})

    #...........................................................................
    # Abstract interface (std::vector<primitives>)
    types = ["int", "double", "std::string"]
    for ndim in dims:
        types += ["Dim<%i>::Vector" % ndim,
                  "Dim<%i>::Tensor" % ndim,
                  "Dim<%i>::SymTensor" % ndim,
                  "Dim<%i>::ThirdRankTensor" % ndim]

    for T in types:
        exec("""
@PYB11pure_virtual
@PYB11pycppname("write")
def write%(Tmangle)s(self,
    value = "const %(T)s&",
    pathName = "const std::string"):
    "Write %(T)s"
    return "void"

@PYB11pure_virtual
@PYB11pycppname("read")
@PYB11const
def read%(Tmangle)s(self,
    value = "%(T)s&",
    pathName = "const std::string"):
    "Read %(T)s"
    return "void"
""" % {"T"       : "std::vector<%s>" % T,
       "Tmangle" : ("vector<%s>" % T).replace(":", "_").replace("<", "_").replace(">", "_")})

    #...........................................................................
    # Abstract interface (Field<primitives>)
    for ndim in dims:
        types = ["int",
                 "Dim<%i>::Scalar" % ndim,
                 "Dim<%i>::Vector" % ndim,
                 "Dim<%i>::Tensor" % ndim,
                 "Dim<%i>::SymTensor" % ndim,
                 "Dim<%i>::ThirdRankTensor" % ndim]

        for T in types:
            exec("""
@PYB11pure_virtual
@PYB11pycppname("write")
def writeField%(Tmangle)s(self,
    value = "const Field<Dim<%(ndim)i>, %(T)s>&",
    pathName = "const std::string"):
    "Write Field<Dim<%(ndim)i, %(T)s>"
    return "void"

@PYB11pure_virtual
@PYB11pycppname("read")
@PYB11const
def readField%(Tmangle)s(self,
    value = "Field<Dim<%(ndim)i>, %(T)s>&",
    pathName = "const std::string"):
    "Read %(T)s"
    return "void"
""" % {"ndim" : ndim,
       "T"       : T,
       "Tmangle" : ("Field<%i%s>" % (ndim, T)).replace(":", "_").replace("<", "_").replace(">", "_")})

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def write_unsigned_int(self,
                           value = "const unsigned int",
                           pathName = "const std::string"):
        "Write an unsigned int"
        return "void"

    @PYB11virtual
    def write_int(self,
                  value = "const int",
                  pathName = "const std::string"):
        "Write an int"
        return "void"

    @PYB11virtual
    def write_bool(self,
                   value = "const bool",
                   pathName = "const std::string"):
        "Write a bool"
        return "void"

    @PYB11virtual
    def write_double(self,
                     value = "const double",
                     pathName = "const std::string"):
        "Write a double"
        return "void"

    @PYB11virtual
    def write_string(self,
                     value = "const std::string",
                     pathName = "const std::string"):
        "Write a std::string"
        return "void"

    @PYB11virtual
    @PYB11const
    def read_unsigned_int(self,
                          pathName = "const std::string"):
        "Read an unsigned int"
        return "unsigned int"

    @PYB11virtual
    @PYB11const
    def read_int(self,
                 pathName = "const std::string"):
        "Read an int"
        return "int"

    @PYB11virtual
    @PYB11const
    def read_bool(self,
                  pathName = "const std::string"):
        "Read a bool"
        return "bool"

    @PYB11virtual
    @PYB11const
    def read_double(self,
                    pathName = "const std::string"):
        "Read a double"
        return "double"

    @PYB11virtual
    @PYB11const
    def read_string(self,
                    pathName = "const std::string"):
        "Read a std::string"
        return "std::string"

    #...........................................................................
    # Methods
    for ndim in xrange(1,4):   # all three always required
        exec('''
@PYB11pycppname("write")
def writePlane%(ndim)i(self,
                       value = "const GeomPlane<Dim<%(ndim)i>>&",
                       pathName = "const std::string"):
    "Write a Plane%(ndim)id"
    return "void"

@PYB11pycppname("read")
@PYB11const
def readPlane%(ndim)i(self,
                      value = "GeomPlane<Dim<%(ndim)i>>&",
                      pathName = "const std::string"):
    "Read a Plane%(ndim)id"
    return "void"
''' % {"ndim" : ndim})
