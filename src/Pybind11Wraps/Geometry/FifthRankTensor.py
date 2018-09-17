from PYB11Decorators import *
from PYB11property import *
from PYB11class import *

#-------------------------------------------------------------------------------
# FifthRankTensor template
#-------------------------------------------------------------------------------
@PYB11template("ndim")
class FifthRankTensor:
    "Spheral fifth rank tensor (%(ndim)sx%(ndim)sx%(ndim)sx%(ndim)s) class"

    # Static attributes
    @PYB11static
    @PYB11readonly
    def nDimensions(self):
        "Number of dimensions"

    @PYB11static
    @PYB11readonly
    def numElements(self):
        "Number of elements stored in the type."

    # Constructors
    def pyinit0(self):
        "Default constructor"

    def pyinit1(self,
                rhs = "const Dim<%(ndim)s>::FifthRankTensor"):
        "Copy constructor"

    def pyinit2(self,
                rhs="double"):
        "Construct setting the element values to a constant value."

    # Sequence methods
    @PYB11implementation("[](const Dim<%(ndim)s>::FifthRankTensor& self) { return Dim<%(ndim)s>::FifthRankTensor::numElements; }")
    def __len__(self):
        "The size (number of elements) of the FifthRankTensor."

    @PYB11implementation("[](const Dim<%(ndim)s>::FifthRankTensor &s, size_t i) { if (i >= Dim<%(ndim)s>::FifthRankTensor::numElements) throw py::index_error(); return s[i]; }") 
    def __getitem__(self):
        "Python indexing to get an element."

    @PYB11implementation("[](Dim<%(ndim)s>::FifthRankTensor &s, size_t i, float v) { if (i >= Dim<%(ndim)s>::FifthRankTensor::numElements) throw py::index_error(); s[i] = v; }") 
    def __setitem__(self):
        "Python indexing to set an element."

    @PYB11implementation("[](const Dim<%(ndim)s>::FifthRankTensor &s) { return py::make_iterator(s.begin(), s.end()); }")
    def __iter__(self):
        "Python iteration through a FifthRankTensor."
    @PYB11const
    def __call__(self,
                 i="Dim<%(ndim)s>::FifthRankTensor::size_type", 
                 j="Dim<%(ndim)s>::FifthRankTensor::size_type",
                 k="Dim<%(ndim)s>::FifthRankTensor::size_type",
                 m="Dim<%(ndim)s>::FifthRankTensor::size_type",
                 n="Dim<%(ndim)s>::FifthRankTensor::size_type"):
        "Extract the (i,j,k,m,n) element."
        return "double"

    # String representation
    @PYB11implementation("""
[](const Dim<%(ndim)s>::FifthRankTensor& self) {
  std::string result = "FifthRankTensor" + std::to_string(%(ndim)s) + "d(";
  for (auto val: self) result += (" " + std::to_string(val) + " ");
  result += ")";
  return result;
}""")
    def __repr__(self):
        return

#-------------------------------------------------------------------------------
# FifthRankTensor instantiations.
#-------------------------------------------------------------------------------
FifthRankTensor1d = PYB11TemplateClass(FifthRankTensor,
                                       template_parameters = ("1"),
                                       cppname = "Dim<1>::FifthRankTensor",
                                       pyname = "FifthRankTensor1d")
FifthRankTensor2d = PYB11TemplateClass(FifthRankTensor,
                                       template_parameters = ("2"),
                                       cppname = "Dim<2>::FifthRankTensor",
                                       pyname = "FifthRankTensor2d")
FifthRankTensor3d = PYB11TemplateClass(FifthRankTensor,
                                       template_parameters = ("3"),
                                       cppname = "Dim<3>::FifthRankTensor",
                                       pyname = "FifthRankTensor3d")
