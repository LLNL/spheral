from PYB11Generator import *

#-------------------------------------------------------------------------------
# EigenStruct template
#-------------------------------------------------------------------------------
@PYB11template("ndim")
class EigenStruct:
    "Holds the eigenvalues and eigenvectors for a second-rank tensor."

    # Constructors
    def pyinit0(self):
        "Default constructor"

    def pyinit1(self,
                rhs = "const EigenStruct<%(ndim)s>"):
        "Copy constructor"

    # Attributes
    @PYB11readwrite
    def eigenValues(self):
        "The vector of eigenvalues."
        return "Dim<%(ndim)s>::Vector"

    @PYB11readwrite
    def eigenVectors(self):
        "The matrix of eigenvectors as column vectors."
        return "Dim<%(ndim)s>::Vector"

    # String representation
    @PYB11implementation("""
[](const EigenStruct<%(ndim)s>& self) {
  std::ostringstream ss;
  ss << "EigenStruct%(ndim)sd[" << self.eigenValues << ", " << self.eigenVectors << "]";
  return ss.str();
}""")
    def __repr__(self):
        return

#-------------------------------------------------------------------------------
# EigenStruct instantiations.
#-------------------------------------------------------------------------------
EigenStruct1d = PYB11TemplateClass(EigenStruct,
                                   template_parameters = ("1"),
                                   cppname = "EigenStruct<1>",
                                   pyname = "EigenStruct1d")
EigenStruct2d = PYB11TemplateClass(EigenStruct,
                                   template_parameters = ("2"),
                                   cppname = "EigenStruct<2>",
                                   pyname = "EigenStruct2d")
EigenStruct3d = PYB11TemplateClass(EigenStruct,
                                   template_parameters = ("3"),
                                   cppname = "EigenStruct<3>",
                                   pyname = "EigenStruct3d")
