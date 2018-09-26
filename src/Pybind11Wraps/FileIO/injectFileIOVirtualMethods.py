#-------------------------------------------------------------------------------
# Helper to inject common virtual methods for equations of state
#-------------------------------------------------------------------------------
import inspect
from PYB11Generator import *
from spheralDimensions import *
dims = spheralDimensions()

@PYB11ignore
def injectFileIOVirtualMethods(cls):

    @PYB11virtual
    def open(self,
             fileName = "const std::string",
             access = "AccessType"):
        "Open a file for IO"
        return "void"

    @PYB11virtual
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
@PYB11virtual
@PYB11pycppname("write")
def write%(Tmangle)s(self,
    value = "const %(T)s&",
    pathName = "const std::string"):
    "Write %(T)s"
    return "void"

@PYB11virtual
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
@PYB11virtual
@PYB11pycppname("write")
def write%(Tmangle)s(self,
    value = "const %(T)s&",
    pathName = "const std::string"):
    "Write %(T)s"
    return "void"

@PYB11virtual
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
@PYB11virtual
@PYB11pycppname("write")
def writeField%(Tmangle)s(self,
    value = "const Field<Dim<%(ndim)i>, %(T)s>&",
    pathName = "const std::string"):
    "Write Field<Dim<%(ndim)i, %(T)s>"
    return "void"

@PYB11virtual
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

    # Inject em...
    names = [x for x in dir() if inspect.isfunction(eval(x))]
    for name in names:
        exec('cls.%(name)s = eval("%(name)s")' % {"name": name})
    return
