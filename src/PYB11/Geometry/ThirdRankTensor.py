from PYB11Generator import *

#-------------------------------------------------------------------------------
# ThirdRankTensor template
#-------------------------------------------------------------------------------
@PYB11template("ndim")
class ThirdRankTensor:
    "Spheral third rank tensor (%(ndim)sx%(ndim)sx%(ndim)s) class"

    # Static attributes
    nrank = PYB11readonly(static=True, doc="Rank of the tensor", returnpolicy="copy")
    nDimensions = PYB11readonly(static=True, doc="Number of dimensions", returnpolicy="copy")
    numElements = PYB11readonly(static=True, doc="Number of elements stored in the type", returnpolicy="copy")
    zero = PYB11readonly(static=True, doc="The zero value equivalent", returnpolicy="copy")

    # Constructors
    def pyinit0(self):
        "Default constructor"

    def pyinit1(self,
                rhs = "const Dim<%(ndim)s>::ThirdRankTensor"):
        "Copy constructor"

    def pyinit2(self,
                rhs="double"):
        "Construct setting the element values to a constant value."

    # Sequence methods
    @PYB11implementation("[](const Dim<%(ndim)s>::ThirdRankTensor&) { return Dim<%(ndim)s>::ThirdRankTensor::numElements; }")
    def __len__(self):
        "The size (number of elements) of the ThirdRankTensor."

    @PYB11implementation("[](const Dim<%(ndim)s>::ThirdRankTensor &s, size_t i) { if (i >= Dim<%(ndim)s>::ThirdRankTensor::numElements) throw py::index_error(); return s[i]; }") 
    @PYB11returnpolicy("reference_internal")
    def __getitem__(self):
        "Python indexing to get an element."

    @PYB11implementation("[](Dim<%(ndim)s>::ThirdRankTensor &s, size_t i, double v) { if (i >= Dim<%(ndim)s>::ThirdRankTensor::numElements) throw py::index_error(); s[i] = v; }") 
    def __setitem__(self):
        "Python indexing to set an element."

    @PYB11implementation("[](const Dim<%(ndim)s>::ThirdRankTensor &s) { return py::make_iterator(s.begin(), s.end()); }, py::keep_alive<0,1>()")
    def __iter__(self):
        "Python iteration through a ThirdRankTensor."

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def __call__(self,
                 i="Dim<%(ndim)s>::ThirdRankTensor::size_type", 
                 j="Dim<%(ndim)s>::ThirdRankTensor::size_type",
                 k="Dim<%(ndim)s>::ThirdRankTensor::size_type"):
        "Extract the (i,j,k) element."
        return "double"

    @PYB11pycppname("__call__")
    @PYB11implementation("""[](Dim<%(ndim)s>::ThirdRankTensor& self, 
                               Dim<%(ndim)s>::ThirdRankTensor::size_type i,
                               Dim<%(ndim)s>::ThirdRankTensor::size_type j,
                               Dim<%(ndim)s>::ThirdRankTensor::size_type k,
                               double val) { self(i,j,k) = val; }""")
    def assignCall(self,
                   i="Dim<%(ndim)s>::ThirdRankTensor::size_type", 
                   j="Dim<%(ndim)s>::ThirdRankTensor::size_type",
                   k="Dim<%(ndim)s>::ThirdRankTensor::size_type",
                   val="double"):
        return "void"

    # Methods
    def Zero(self):
        "Zero out the elements"
        return "void"

    @PYB11const
    def doubledot(self, rhs="const RankNTensor<%(ndim)s, 3, GeomThirdRankTensor<%(ndim)s>>& rhs"):
        return "double"

    @PYB11const
    def squareElements(self):
        return "const Dim<%(ndim)s>::ThirdRankTensor"

    @PYB11const
    def maxAbsElement(self):
        return "double"

    # Operators
    def __neg__(self):
        return
    def __iadd__(self):
        return
    def __isub__(self):
        return
    def __add__(self):
        return
    def __sub__(self):
        return
    def __imul__(self, rhs="double()"):
        return
    def __idiv__(self, rhs="double()"):
        return
    def __mul__(self, rhs="double()"):
        return
    def __div__(self, rhs="double()"):
        return
                 
    # Comparisons
    def __eq__(self):
        return
    def __ne__(self):
        return
    def __lt__(self):
        return
    def __gt__(self):
        return
    def __le__(self):
        return
    def __ge__(self):
        return

    # String representation
    @PYB11implementation("""
[](const Dim<%(ndim)s>::ThirdRankTensor& self) {
  std::string result = "ThirdRankTensor" + std::to_string(%(ndim)s) + "d(";
  for (auto val: self) result += (" " + std::to_string(val) + " ");
  result += ")";
  return result;
}""")
    def __repr__(self):
        return

#-------------------------------------------------------------------------------
# ThirdRankTensor instantiations.
#-------------------------------------------------------------------------------
ThirdRankTensor1d = PYB11TemplateClass(ThirdRankTensor,
                                       template_parameters = ("1"),
                                       cppname = "Dim<1>::ThirdRankTensor",
                                       pyname = "ThirdRankTensor1d")
ThirdRankTensor2d = PYB11TemplateClass(ThirdRankTensor,
                                       template_parameters = ("2"),
                                       cppname = "Dim<2>::ThirdRankTensor",
                                       pyname = "ThirdRankTensor2d")
ThirdRankTensor3d = PYB11TemplateClass(ThirdRankTensor,
                                       template_parameters = ("3"),
                                       cppname = "Dim<3>::ThirdRankTensor",
                                       pyname = "ThirdRankTensor3d")
