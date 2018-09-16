from PYB11Decorators import *
from PYB11property import *
from PYB11class import *

#-------------------------------------------------------------------------------
# Tensor (rank 2) template
#-------------------------------------------------------------------------------
@PYB11template("ndim")
class Tensor:
    "Spheral Geometric Tensor (rank 2: %(ndim)sx%(ndim)s) class"

    # Static attributes
    @PYB11static
    @PYB11readonly
    def nDimensions(self):
        "Number of dimensions"

    @PYB11static
    @PYB11readonly
    def numElements(self):
        "Number of elements stored in the type."

    @PYB11static
    @PYB11readonly
    def zero(self):
        "The zero value equivalent."

    @PYB11static
    @PYB11readonly
    def one(self):
        "The unit value equivalent."

    # Constructors
    def pyinit0(self):
        "Default constructor"

    def pyinit1(self,
                rhs = "const Dim<%(ndim)s>::Tensor"):
        "Copy constructor (tensor)"

    def pyinit2(self,
                rhs = "const Dim<%(ndim)s>::SymTensor"):
        "Copy constructor (symmetric tensor)"

    def pyinit3(self,
                xx="double"):
        "Construct with element values."

    def pyinit3(self,
                xx="double", xy="double",
                yx="double", yy="double"):
        "Construct with element values."

    def pyinit4(self,
                xx="double", xy="double", xz="double",
                yx="double", yy="double", yz="double",
                zx="double", zy="double", zz="double"):
        "Construct with element values."

    # Sequence methods
    @PYB11implementation("[](const Dim<%(ndim)s>::Tensor& self) { return Dim<%(ndim)s>::Tensor::numElements; }")
    def __len__(self):
        "The size (number of elements) of the Tensor."

    @PYB11implementation("[](const Dim<%(ndim)s>::Tensor &s, size_t i) { if (i >= Dim<%(ndim)s>::Tensor::numElements) throw py::index_error(); return s[i]; }") 
    def __getitem__(self):
        "Python indexing to get an element."

    @PYB11implementation("[](Dim<%(ndim)s>::Tensor &s, size_t i, float v) { if (i >= Dim<%(ndim)s>::Tensor::numElements) throw py::index_error(); s[i] = v; }") 
    def __setitem__(self):
        "Python indexing to set an element."

    @PYB11implementation("[](const Dim<%(ndim)s>::Tensor &s) { return py::make_iterator(s.begin(), s.end()); }")
    def __iter__(self):
        "Python iteration through a Tensor."

    @PYB11const
    def __call__(self,
                 row="Dim<%(ndim)s>::Tensor::size_type", 
                 col="Dim<%(ndim)s>::Tensor::size_type"):
        "Extract the (row, column) element."
        return "double"

    # String representation
    @PYB11implementation("""
[](const Dim<%(ndim)s>::Tensor& self) {
  std::string result = "Tensor" + std::to_string(%(ndim)s) + "d(";
  for (auto val: self) result += (" " + std::to_string(val) + " ");
  result += ")";
  return result;
}""")
    def __repr__(self):
        return

    # Operators
    def __neg__(self):
        return
    def __add__(self):
        return
    def __sub__(self):
        return
    def __mul__(self):
        return
    def __iadd__(self):
        return
    def __isub__(self):
        return

    def __add__(self, rhs="Dim<%(ndim)s>::SymTensor()"):
        return
    def __sub__(self, rhs="Dim<%(ndim)s>::SymTensor()"):
        return
    def __mul__(self, rhs="Dim<%(ndim)s>::SymTensor()"):
        return

    def __mul__(self, rhs="float()"):
        return
    def __rmul__(self, rhs="float()"):
        return
    def __div__(self, rhs="float()"):
        return
    def __imul__(self, rhs="float()"):
        return
    def __idiv__(self, rhs="float()"):
        return

    # Comparison
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

    # Methods
    def getRow(self):
        "Extract the indexed row as a Vector%(ndim)sd."
    def getColumn(self):
        "Extract the indexed column as a Vector%(ndim)sd."
    def Zero(self):
        "Zero out the elements."
    def Identity(self):
        "Set equal to the I tensor."
    def Symmetric(self):
        "Return a SymTensor with the symmetric portion of this tensor."
    def SkewSymmetric(self):
        "Return a Tensor with the skew-symmetric portion of this tensor."
    def Transpose(self):
        "Return the transpose of this Tensor."
    def Inverse(self):
        "Return the inverse of the tensor."
    def diagonalElements(self):
        "Return a Vector%(ndim)sd with the diagonal elements of this tensor."
    def Trace(self):
        "Compute the trace (sum of diagonal elements)."
    def Determinant(self):
        "Compute the determinant of the tensor."
    @PYB11const
    def dot(self, rhs="const Dim<%(ndim)s>::Vector&"):
        "Product with a Vector%(ndim)sd."
        return "Dim<%(ndim)s>::Vector"
    @PYB11const
    def doubledot(self, rhs="const Dim<%(ndim)s>::Tensor&"):
        "Double dot contraction with another tensor (returns a scalar)."
        return "double"
    @PYB11const
    def doubledot(self, rhs="const Dim<%(ndim)s>::SymTensor&"):
        "Double dot contraction with a SymTensor (returns a scalar)."
        return "double"
    def selfDoubledot(self):
        "Double dot contraction with ourself (returns a scalar)."
    def square(self):
        "Compute the product with ourself as a matrix product."
    def squareElements(self):
        "Returns the element-wise square of the tensor."
    def eigenValues(self):
        "Return a Vector%(ndim)sd with the eigenvalues of this tensor."
    def rotationalTransform(self):
        "Apply the given rotational transform to the tensor."
    def maxAbsElement(self):
        "Return the maximum of the absolute values of the elements."

    # xx
    @PYB11cppname("xx")
    @PYB11const
    @PYB11ignore
    def getxx(self):
        return "double"

    @PYB11cppname("xx")
    @PYB11ignore
    def setxx(self, val="double"):
        return "void"

    # xy
    @PYB11cppname("xy")
    @PYB11const
    @PYB11ignore
    def getxy(self):
        return "double"

    @PYB11cppname("xy")
    @PYB11ignore
    def setxy(self, val="double"):
        return "void"

    # xz
    @PYB11cppname("xz")
    @PYB11const
    @PYB11ignore
    def getxz(self):
        return "double"

    @PYB11cppname("xz")
    @PYB11ignore
    def setxz(self, val="double"):
        return "void"

    # yx
    @PYB11cppname("yx")
    @PYB11const
    @PYB11ignore
    def getyx(self):
        return "double"

    @PYB11cppname("yx")
    @PYB11ignore
    def setyx(self, val="double"):
        return "void"

    # yy
    @PYB11cppname("yy")
    @PYB11const
    @PYB11ignore
    def getyy(self):
        return "double"

    @PYB11cppname("yy")
    @PYB11ignore
    def setyy(self, val="double"):
        return "void"

    # yz
    @PYB11cppname("yz")
    @PYB11const
    @PYB11ignore
    def getyz(self):
        return "double"

    @PYB11cppname("yz")
    @PYB11ignore
    def setyz(self, val="double"):
        return "void"

    # zx
    @PYB11cppname("zx")
    @PYB11const
    @PYB11ignore
    def getzx(self):
        return "double"

    @PYB11cppname("zx")
    @PYB11ignore
    def setzx(self, val="double"):
        return "void"

    # zy
    @PYB11cppname("zy")
    @PYB11const
    @PYB11ignore
    def getzy(self):
        return "double"

    @PYB11cppname("zy")
    @PYB11ignore
    def setzy(self, val="double"):
        return "void"

    # zz
    @PYB11cppname("zz")
    @PYB11const
    @PYB11ignore
    def getzz(self):
        return "double"

    @PYB11cppname("zz")
    @PYB11ignore
    def setzz(self, val="double"):
        return "void"

    # Properties
    xx = property(getxx, setxx, doc="The xx element.")
    xy = property(getxy, setxy, doc="The xy element.")
    xz = property(getxz, setxz, doc="The xz element.")
    yx = property(getyx, setyx, doc="The yx element.")
    yy = property(getyy, setyy, doc="The yy element.")
    yz = property(getyz, setyz, doc="The yz element.")
    zx = property(getzx, setzx, doc="The zx element.")
    zy = property(getzy, setzy, doc="The zy element.")
    zz = property(getzz, setzz, doc="The zz element.")

#-------------------------------------------------------------------------------
# Tensor instantiations.
#-------------------------------------------------------------------------------
Tensor1d = PYB11TemplateClass(Tensor,
                              template_parameters = ("1"),
                              cppname = "Dim<1>::Tensor",
                              pyname = "Tensor1d")
Tensor2d = PYB11TemplateClass(Tensor,
                              template_parameters = ("2"),
                              cppname = "Dim<2>::Tensor",
                              pyname = "Tensor2d")
Tensor3d = PYB11TemplateClass(Tensor,
                              template_parameters = ("3"),
                              cppname = "Dim<3>::Tensor",
                              pyname = "Tensor3d")
