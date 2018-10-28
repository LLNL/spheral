#-------------------------------------------------------------------------------
# FileIO abstract class
#-------------------------------------------------------------------------------
from PYB11Generator import *
from FileIOAbstractMethods import *
from spheralDimensions import *
dims = spheralDimensions()

@PYB11module("SpheralFileIO")
class FileIO:
    "Abstract base class for FileIO objects"

    #...........................................................................
    # Constructors
    def pyinit0(self):
        "Default constructor"

    def pyinit1(self,
                filename = "const std::string",
                access = "AccessType"):
        "Construct with a given file name and access"

    #...........................................................................
    # Abstract methods
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

    @PYB11const
    def splitPathComponents(self, pathName="const std::string&"):
        "A helper function to split a string up into substrings delimited by '/'."
        return "std::vector<std::string>"

    @PYB11const
    def groupName(self, pathName="const std::string"):
        "Return the group (directory) component of a path."
        return "std::string"

    @PYB11const
    def variableName(self, pathName="const std::string"):
        "Return the variable component of a path."
        return "std::string"

    @PYB11implementation("[](FileIO& self, py::handle thing, py::handle path) { self.writeObject(thing.ptr(), path.ptr()); }")
    def writeObject(self,
                    thing = "py::handle",
                    path = "py::handle"):
        "Handle a generic python object through serialization"
        return "void"

    @PYB11returnpolicy("take_ownership")
    @PYB11const
    @PYB11implementation("[](FileIO& self, py::handle path) { return py::handle(self.readObject(path.ptr())); }")
    def readObject(self,
                   path = "py::handle"):
        "Return a generic python object from deserialization."
        return "py::handle"

    #...........................................................................
    # We have to explicitly list the template instatiations for std::vector<>
    # first due to poor overloading selection in pybind11
    @PYB11pycppname("write")
    def writeVecString(self,
                 x = "const std::vector<std::string>&",
                 pathName = "const std::string"):
        "Write a std::vector<std::string>"
        return "void"

    @PYB11const
    @PYB11pycppname("read")
    def readVecString(self,
                x = "std::vector<std::string>&",
                pathName = "const std::string"):
        "Read a std::vector<std::string>"
        return "void"
    
    @PYB11pycppname("write")
    def writeVecInt(self,
                      x = "const std::vector<int>&",
                      pathName = "const std::string"):
        "Write a std::vector<int>"
        return "void"

    @PYB11const
    @PYB11pycppname("read")
    def readVecInt(self,
                   x = "std::vector<int>&",
                   pathName = "const std::string"):
        "Read a std::vector<int>"
        return "void"
    
    @PYB11pycppname("write")
    def writeVecDouble(self,
                      x = "const std::vector<double>&",
                      pathName = "const std::string"):
        "Write a std::vector<double>"
        return "void"

    @PYB11const
    @PYB11pycppname("read")
    def readVecDouble(self,
                   x = "std::vector<double>&",
                   pathName = "const std::string"):
        "Read a std::vector<double>"
        return "void"
    
    @PYB11pycppname("write")
    def writeVecVec1d(self,
                      x = "const std::vector<Dim<1>::Vector>&",
                      pathName = "const std::string"):
        "Write a std::vector<Vector1d>"
        return "void"

    @PYB11const
    @PYB11pycppname("read")
    def readVecVec1d(self,
                     x = "std::vector<Dim<1>::Vector>&",
                     pathName = "const std::string"):
        "Read a std::vector<Vector1d>"
        return "void"
    
    @PYB11pycppname("write")
    def writeVecTen1d(self,
                      x = "const std::vector<Dim<1>::Tensor>&",
                      pathName = "const std::string"):
        "Write a std::vector<Tensor1d>"
        return "void"

    @PYB11const
    @PYB11pycppname("read")
    def readVecTen1d(self,
                     x = "std::vector<Dim<1>::Tensor>&",
                     pathName = "const std::string"):
        "Read a std::vector<Tensor1d>"
        return "void"
    
    @PYB11pycppname("write")
    def writeVecSymTens1d(self,
                          x = "const std::vector<Dim<1>::SymTensor>&",
                          pathName = "const std::string"):
        "Write a std::vector<SymTensor1d>"
        return "void"

    @PYB11const
    @PYB11pycppname("read")
    def readVecSymTen1d(self,
                        x = "std::vector<Dim<1>::SymTensor>&",
                        pathName = "const std::string"):
        "Read a std::vector<SymTensor1d>"
        return "void"
    
    @PYB11pycppname("write")
    def writeVecThirdRankTen1d(self,
                               x = "const std::vector<Dim<1>::ThirdRankTensor>&",
                               pathName = "const std::string"):
        "Write a std::vector<ThirdRankTensor1d>"
        return "void"

    @PYB11const
    @PYB11pycppname("read")
    def readVecThirdRankTen1d(self,
                              x = "std::vector<Dim<1>::ThirdRankTensor>&",
                              pathName = "const std::string"):
        "Read a std::vector<ThirdRankTensor1d>"
        return "void"
    
    @PYB11pycppname("write")
    def writeVecVec2d(self,
                      x = "const std::vector<Dim<2>::Vector>&",
                      pathName = "const std::string"):
        "Write a std::vector<Vector2d>"
        return "void"

    @PYB11const
    @PYB11pycppname("read")
    def readVecVec2d(self,
                     x = "std::vector<Dim<2>::Vector>&",
                     pathName = "const std::string"):
        "Read a std::vector<Vector2d>"
        return "void"
    
    @PYB11pycppname("write")
    def writeVecTen2d(self,
                      x = "const std::vector<Dim<2>::Tensor>&",
                      pathName = "const std::string"):
        "Write a std::vector<Tensor2d>"
        return "void"

    @PYB11const
    @PYB11pycppname("read")
    def readVecTen2d(self,
                     x = "std::vector<Dim<2>::Tensor>&",
                     pathName = "const std::string"):
        "Read a std::vector<Tensor2d>"
        return "void"
    
    @PYB11pycppname("write")
    def writeVecSymTens2d(self,
                          x = "const std::vector<Dim<2>::SymTensor>&",
                          pathName = "const std::string"):
        "Write a std::vector<SymTensor2d>"
        return "void"

    @PYB11const
    @PYB11pycppname("read")
    def readVecSymTen2d(self,
                        x = "std::vector<Dim<2>::SymTensor>&",
                        pathName = "const std::string"):
        "Read a std::vector<SymTensor2d>"
        return "void"
    
    @PYB11pycppname("write")
    def writeVecThirdRankTen2d(self,
                               x = "const std::vector<Dim<2>::ThirdRankTensor>&",
                               pathName = "const std::string"):
        "Write a std::vector<ThirdRankTensor2d>"
        return "void"

    @PYB11const
    @PYB11pycppname("read")
    def readVecThirdRankTen2d(self,
                              x = "std::vector<Dim<2>::ThirdRankTensor>&",
                              pathName = "const std::string"):
        "Read a std::vector<ThirdRankTensor2d>"
        return "void"
    
    @PYB11pycppname("write")
    def writeVecVec3d(self,
                      x = "const std::vector<Dim<3>::Vector>&",
                      pathName = "const std::string"):
        "Write a std::vector<Vector3d>"
        return "void"

    @PYB11const
    @PYB11pycppname("read")
    def readVecVec3d(self,
                     x = "std::vector<Dim<3>::Vector>&",
                     pathName = "const std::string"):
        "Read a std::vector<Vector3d>"
        return "void"
    
    @PYB11pycppname("write")
    def writeVecTen3d(self,
                      x = "const std::vector<Dim<3>::Tensor>&",
                      pathName = "const std::string"):
        "Write a std::vector<Tensor3d>"
        return "void"

    @PYB11const
    @PYB11pycppname("read")
    def readVecTen3d(self,
                     x = "std::vector<Dim<3>::Tensor>&",
                     pathName = "const std::string"):
        "Read a std::vector<Tensor3d>"
        return "void"
    
    @PYB11pycppname("write")
    def writeVecSymTens3d(self,
                          x = "const std::vector<Dim<3>::SymTensor>&",
                          pathName = "const std::string"):
        "Write a std::vector<SymTensor3d>"
        return "void"

    @PYB11const
    @PYB11pycppname("read")
    def readVecSymTen3d(self,
                        x = "std::vector<Dim<3>::SymTensor>&",
                        pathName = "const std::string"):
        "Read a std::vector<SymTensor3d>"
        return "void"
    
    @PYB11pycppname("write")
    def writeVecThirdRankTen3d(self,
                               x = "const std::vector<Dim<3>::ThirdRankTensor>&",
                               pathName = "const std::string"):
        "Write a std::vector<ThirdRankTensor3d>"
        return "void"

    @PYB11const
    @PYB11pycppname("read")
    def readVecThirdRankTen3d(self,
                              x = "std::vector<Dim<3>::ThirdRankTensor>&",
                              pathName = "const std::string"):
        "Read a std::vector<ThirdRankTensor3d>"
        return "void"
    
    #...........................................................................
    # Templates
    @PYB11template("Dimension", "Value")
    @PYB11pycppname("write")
    def writeFieldList(self,
                       fieldList = "const FieldList<%(Dimension)s, %(Value)s>&",
                       pathName = "const std::string"):
        "Write a FieldList<%(Dimension)s, %(Value)s>"
        return "void"

    @PYB11template("Dimension", "Value")
    @PYB11const
    @PYB11pycppname("read")
    def readFieldList(self,
                      fieldList = "FieldList<%(Dimension)s, %(Value)s>&",
                      pathName = "const std::string"):
        "Read a FieldList<%(Dimension)s, %(Value)s>"
        return "void"

    @PYB11template("Dimension", "Value")
    @PYB11pycppname("write")
    def writeFieldVec(self,
                      field = "const Field<%(Dimension)s, std::vector<%(Value)s>>&",
                      pathName = "const std::string"):
        "Write a Field<%(Dimension)s, std::vector<%(Value)s>>"
        return "void"

    @PYB11template("Dimension", "Value")
    @PYB11const
    @PYB11pycppname("read")
    def readFieldVec(self,
                     field = "Field<%(Dimension)s, std::vector<%(Value)s>>&",
                     pathName = "const std::string"):
        "Read a Field<%(Dimension)s, std::vector<%(Value)s>>"
        return "void"

    for ndim in dims:
        for T in ["int",
                  "Dim<%i>::Scalar" % ndim,
                  "Dim<%i>::Vector" % ndim,
                  "Dim<%i>::Tensor" % ndim,
                  "Dim<%i>::SymTensor" % ndim,
                  "Dim<%i>::ThirdRankTensor" % ndim]:
            exec('''
writeFieldList%(Tmangle)s = PYB11TemplateMethod(writeFieldList, template_parameters=("%(Dimension)s", "%(T)s"), pyname="write")
readFieldList%(Tmangle)s =  PYB11TemplateMethod(readFieldList,  template_parameters=("%(Dimension)s", "%(T)s"), pyname="read")

writeFieldVec%(Tmangle)s = PYB11TemplateMethod(writeFieldVec, template_parameters=("%(Dimension)s", "%(T)s"), pyname="write")
readFieldVec%(Tmangle)s =  PYB11TemplateMethod(readFieldVec,  template_parameters=("%(Dimension)s", "%(T)s"), pyname="read")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">",
       "T"         : T,
       "Tmangle"   : ("Field<%i%s>" % (ndim, T)).replace(":", "_").replace("<", "_").replace(">", "_")})

    #...........................................................................
    # Properties
    fileName = PYB11property("const std::string&", "fileName", doc="The current file name")
    access = PYB11property("AccessType", "access", doc="The access type of the currently open file")
    fileOpen = PYB11property("bool", "fileOpen", doc="Is the file currently open?")

#-------------------------------------------------------------------------------
# Inject the abstract interface
#-------------------------------------------------------------------------------
PYB11inject(FileIOAbstractMethods, FileIO, virtual=False, pure_virtual=True)
