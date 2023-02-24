#-------------------------------------------------------------------------------
# Define the common template methods of FileIO
#-------------------------------------------------------------------------------
from PYB11Generator import *
from spheralDimensions import *
dims = spheralDimensions()

@PYB11ignore
class FileIOTemplateMethods:

    #...........................................................................
    # Field<Dim, Value>
    @PYB11template("Dimension", "Value")
    @PYB11noconvert
    @PYB11cppname("write")
    def writeField(self,
                   x = "const Field<%(Dimension)s, %(Value)s>&",
                   path = "const std::string"):
        "Write a Field<%(Dimension)s, %(Value)s"
        return "void"

    @PYB11template("Dimension", "Value")
    @PYB11const
    @PYB11noconvert
    @PYB11cppname("read")
    def readField(self,
                  x = "Field<%(Dimension)s, %(Value)s>&",
                  path = "const std::string"):
        "Read a Field<%(Dimension)s, %(Value)s>"
        return "void"

    #...........................................................................
    # FieldList<Dim, Value>
    @PYB11template("Dimension", "Value")
    @PYB11noconvert
    @PYB11cppname("write")
    def writeFieldList(self,
                       x = "const FieldList<%(Dimension)s, %(Value)s>&",
                       path = "const std::string"):
        "Write a FieldList<%(Dimension)s, %(Value)s"
        return "void"

    @PYB11template("Dimension", "Value")
    @PYB11const
    @PYB11noconvert
    @PYB11cppname("read")
    def readFieldList(self,
                      x = "FieldList<%(Dimension)s, %(Value)s>&",
                      path = "const std::string"):
        "Read a FieldList<%(Dimension)s, %(Value)s>"
        return "void"

    #...........................................................................
    # std::vector<Value>
    @PYB11template("Value")
    @PYB11noconvert
    @PYB11cppname("write")
    def writeVec(self,
                 x = "const std::vector<%(Value)s>&",
                 path = "const std::string"):
        "Write a std::vector<%(Value)s>"
        return "void"

    @PYB11template("Value")
    @PYB11const
    @PYB11noconvert
    @PYB11cppname("read")
    def readVec(self,
                x = "std::vector<%(Value)s>&",
                path = "const std::string"):
        "Read a std::vector<%(Value)s>"
        return "void"

#-------------------------------------------------------------------------------
# Instantiations
#-------------------------------------------------------------------------------
    for ndim in dims:

        # std::vector
        for T in ["Dim<%i>::Vector" % ndim,
                  "Dim<%i>::Tensor" % ndim,
                  "Dim<%i>::SymTensor" % ndim,
                  "Dim<%i>::ThirdRankTensor" % ndim,
                  "Dim<%i>::FacetedVolume" % ndim]:
            exec('''
writeVec{Tmangle} = PYB11TemplateMethod(writeVec, template_parameters="{T}", pyname="write")
readVec{Tmangle}  = PYB11TemplateMethod( readVec, template_parameters="{T}", pyname="read")
'''.format(ndim      = ndim,
           T         = T,
           Tmangle   = ("<%i%s>" % (ndim, T)).replace(":", "_").replace("<", "_").replace(">", "_")))

        # Field/FieldList
        for T in ["int",
                  "unsigned",
                  "Dim<%i>::Scalar" % ndim,
                  "Dim<%i>::Vector" % ndim,
                  "Dim<%i>::Tensor" % ndim,
                  "Dim<%i>::SymTensor" % ndim,
                  "Dim<%i>::ThirdRankTensor" % ndim,
                  "Dim<%i>::FacetedVolume" % ndim]:
            exec('''
writeField{Tmangle} = PYB11TemplateMethod(writeField, template_parameters=("{dim}", "{T}"), pyname="write")
readField{Tmangle} = PYB11TemplateMethod(readField, template_parameters=("{dim}", "{T}"), pyname="read")

writeFieldList{Tmangle} = PYB11TemplateMethod(writeFieldList, template_parameters=("{dim}", "{T}"), pyname="write")
readFieldList{Tmangle} = PYB11TemplateMethod(readFieldList, template_parameters=("{dim}", "{T}"), pyname="read")
'''.format(dim       = "Dim<%i>" % ndim,
           T         = T,
           Tmangle   = ("<%i%s>" % (ndim, T)).replace(":", "_").replace("<", "_").replace(">", "_")))
