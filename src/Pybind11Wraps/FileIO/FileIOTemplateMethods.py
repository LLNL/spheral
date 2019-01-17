#-------------------------------------------------------------------------------
# Define the common template methods of FileIO
#-------------------------------------------------------------------------------
from PYB11Generator import *
from spheralDimensions import *
dims = spheralDimensions()

@PYB11ignore
class FileIOTemplateMethods:

    #...........................................................................
    # std::vector<Value>
    @PYB11template("Value")
    @PYB11cppname("write")
    @PYB11noconvert
    def writeVec(self,
                 x = "const std::vector<%(Value)s>&",
                 pathName = "const std::string"):
        "Write a std::vector<%(Value)s>"
        return "void"

    @PYB11template("Value")
    @PYB11const
    @PYB11cppname("read")
    @PYB11noconvert
    def readVec(self,
                x = "std::vector<%(Value)s>&",
                pathName = "const std::string"):
        "Read a std::vector<%(Value)s>"
        return "void"

    #writeVecInt = PYB11TemplateMethod(writeVec, template_parameters="int", pyname="write")
    #readVecInt  = PYB11TemplateMethod( readVec, template_parameters="int", pyname="read")
    #writeVecDouble = PYB11TemplateMethod(writeVec, template_parameters="double", pyname="write")
    #readVecDouble  = PYB11TemplateMethod( readVec, template_parameters="double", pyname="read")
    #writeVecString = PYB11TemplateMethod(writeVec, template_parameters="std::string", pyname="write")
    #readVecString  = PYB11TemplateMethod( readVec, template_parameters="std::string", pyname="read")

    for ndim in dims:
        for T in ["Dim<%i>::Vector" % ndim,
                  "Dim<%i>::Tensor" % ndim,
                  "Dim<%i>::SymTensor" % ndim,
                  "Dim<%i>::ThirdRankTensor" % ndim]:
            exec('''
writeVec%(Tmangle)s = PYB11TemplateMethod(writeVec, template_parameters="%(T)s", pyname="write")
readVec%(Tmangle)s  = PYB11TemplateMethod( readVec, template_parameters="%(T)s", pyname="read")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">",
       "T"         : T,
       "Tmangle"   : ("Field<%i%s>" % (ndim, T)).replace(":", "_").replace("<", "_").replace(">", "_")})
