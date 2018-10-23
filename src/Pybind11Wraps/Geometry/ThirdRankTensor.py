from PYB11Generator import *

#-------------------------------------------------------------------------------
# ThirdRankTensor template
#-------------------------------------------------------------------------------
@PYB11template("ndim")
class ThirdRankTensor:
    "Spheral third rank tensor (%(ndim)sx%(ndim)sx%(ndim)s) class"

    # Static attributes
    nDimensions = PYB11readonly(static=True, doc="Number of dimensions")
    numElements = PYB11readonly(static=True, doc="Number of elements stored in the type")
    zero = PYB11readonly(static=True, doc="The zero value equivalent")

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
    @PYB11implementation("[](const Dim<%(ndim)s>::ThirdRankTensor& self) { return Dim<%(ndim)s>::ThirdRankTensor::numElements; }")
    def __len__(self):
        "The size (number of elements) of the ThirdRankTensor."

    @PYB11implementation("[](const Dim<%(ndim)s>::ThirdRankTensor &s, size_t i) { if (i >= Dim<%(ndim)s>::ThirdRankTensor::numElements) throw py::index_error(); return s[i]; }") 
    @PYB11returnpolicy("reference_internal")
    def __getitem__(self):
        "Python indexing to get an element."

    @PYB11implementation("[](Dim<%(ndim)s>::ThirdRankTensor &s, size_t i, float v) { if (i >= Dim<%(ndim)s>::ThirdRankTensor::numElements) throw py::index_error(); s[i] = v; }") 
    def __setitem__(self):
        "Python indexing to set an element."

    @PYB11implementation("[](const Dim<%(ndim)s>::ThirdRankTensor &s) { return py::make_iterator(s.begin(), s.end()); }")
    def __iter__(self):
        "Python iteration through a ThirdRankTensor."
    @PYB11const
    def __call__(self,
                 i="Dim<%(ndim)s>::ThirdRankTensor::size_type", 
                 j="Dim<%(ndim)s>::ThirdRankTensor::size_type",
                 k="Dim<%(ndim)s>::ThirdRankTensor::size_type"):
        "Extract the (i,j,k) element."
        return "double"

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
