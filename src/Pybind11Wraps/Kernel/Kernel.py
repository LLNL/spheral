#-------------------------------------------------------------------------------
# Generic Kernel bindings.
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11template("Dimension", "Descendant")
class Kernel:

    typedefs = """
    typedef %(Dimension)s::Vector Vector;
    typedef %(Dimension)s::SymTensor SymTensor;
"""

    # def pyinit(self):
    #     "Default constructor"

    @PYB11const
    def __call__(self,
                 etaMagnitude = "double",
                 H = "const SymTensor&"):
        "Return the kernel value at the given normalized radius"
        return "double"

    #...........................................................................
    @PYB11const
    def kernelValue(self,
                    etaMagnitude = "double",
                    Hdet = "const double&"):
        "Return the value of the kernel"
        return "double"

    @PYB11pycppname("kernelValue")
    @PYB11const
    def kernelValue1(self,
              etaMagnitude = "double",
              H = "const SymTensor&"):
        "Return the value of the kernel"
        return "double"

    @PYB11pycppname("kernelValue")
    @PYB11const
    def kernelValue2(self,
              eta = "const Vector&",
              Hdet = "const double&"):
        "Return the value of the kernel"
        return "double"

    @PYB11pycppname("kernelValue")
    @PYB11const
    def kernelValue3(self,
              eta = "const Vector&",
              H = "const SymTensor&"):
        "Return the value of the kernel"
        return "double"

    #...........................................................................
    @PYB11const
    def gradValue(self,
                    etaMagnitude = "double",
                    Hdet = "const double&"):
        "Return the value of the gradient of the kernel"
        return "double"

    @PYB11pycppname("gradValue")
    @PYB11const
    def gradValue1(self,
              etaMagnitude = "double",
              H = "const SymTensor&"):
        "Return the value of the gradient of the kernel"
        return "double"

    @PYB11pycppname("gradValue")
    @PYB11const
    def gradValue2(self,
              eta = "const Vector&",
              Hdet = "const double&"):
        "Return the value of the gradient of the kernel"
        return "double"

    @PYB11pycppname("gradValue")
    @PYB11const
    def gradValue3(self,
              eta = "const Vector&",
              H = "const SymTensor&"):
        "Return the value of the gradient of the kernel"
        return "double"

    #...........................................................................
    @PYB11const
    def grad2Value(self,
                   etaMagnitude = "double",
                   Hdet = "const double&"):
        "Return the value of the second derivative of the kernel"
        return "double"

    @PYB11pycppname("grad2Value")
    @PYB11const
    def grad2Value1(self,
              etaMagnitude = "double",
              H = "const SymTensor&"):
        "Return the value of the second derivative of the kernel"
        return "double"

    @PYB11pycppname("grad2Value")
    @PYB11const
    def grad2Value2(self,
              eta = "const Vector&",
              Hdet = "const double&"):
        "Return the value of the second derivative of the kernel"
        return "double"

    @PYB11pycppname("grad2Value")
    @PYB11const
    def grad2Value3(self,
              eta = "const Vector&",
              H = "const SymTensor&"):
        "Return the value of the second derivative of the kernel"
        return "double"

    #...........................................................................
    @PYB11pycppname("grad")
    @PYB11const
    def grad_0(self,
               etaMagnitude = "double",
               Hdet = "const double&"):
        "Return the gradient of the kernel"
        return "double"

    @PYB11pycppname("grad")
    @PYB11const
    def grad_1(self,
               etaMagnitude = "double",
               H = "const SymTensor&"):
        "Return the gradient of the kernel"
        return "double"

    @PYB11pycppname("grad")
    @PYB11const
    def grad_2(self,
               eta = "const Vector&",
               Hdet = "const double&"):
        "Return the gradient of the kernel"
        return "double"

    @PYB11pycppname("grad")
    @PYB11const
    def grad_3(self,
               eta = "const Vector&",
               H = "const SymTensor&"):
        "Return the gradient of the kernel"
        return "double"

    #...........................................................................
    @PYB11pycppname("grad2")
    @PYB11const
    def grad2_0(self,
               etaMagnitude = "double",
               Hdet = "const double&"):
        "Return the second derivative of the kernel"
        return "double"

    @PYB11pycppname("grad2")
    @PYB11const
    def grad2_1(self,
               etaMagnitude = "double",
               H = "const SymTensor&"):
        "Return the second derivative of the kernel"
        return "double"

    @PYB11pycppname("grad2")
    @PYB11const
    def grad2_2(self,
               eta = "const Vector&",
               Hdet = "const double&"):
        "Return the second derivative of the kernel"
        return "double"

    @PYB11pycppname("grad2")
    @PYB11const
    def grad2_3(self,
               eta = "const Vector&",
               H = "const SymTensor&"):
        "Return the second derivative of the kernel"
        return "double"

    #...........................................................................
    @PYB11pycppname("gradh")
    @PYB11const
    def gradh_0(self,
                etaMagnitude = "double",
                Hdet = "const double&"):
        "Return the gradient with respect to h of the kernel"
        return "double"

    @PYB11pycppname("gradh")
    @PYB11const
    def gradh_1(self,
               etaMagnitude = "double",
               H = "const SymTensor&"):
        "Return the gradient with respect to h of the kernel"
        return "double"

    @PYB11pycppname("gradh")
    @PYB11const
    def gradh_2(self,
               eta = "const Vector&",
               Hdet = "const double&"):
        "Return the gradient with respect to h of the kernel"
        return "double"

    @PYB11pycppname("gradh")
    @PYB11const
    def gradh_3(self,
               eta = "const Vector&",
               H = "const SymTensor&"):
        "Return the gradient with respect to h of the kernel"
        return "double"

    #...........................................................................
    @PYB11implementation("&KernelPublicist<%(Dimension)s,%(Descendant)s>::setVolumeNormalization")
    @PYB11const
    def setVolumeNormalization(self,
                               val = "double"):
        "Set the volume normalization constant."
        return "void"

    @PYB11implementation("&KernelPublicist<%(Dimension)s,%(Descendant)s>::setKernelExtent")
    @PYB11const
    def setKernelExtent(self,
                        val = "double"):
        "Set the volume normalization constant."
        return "void"

    @PYB11implementation("&KernelPublicist<%(Dimension)s,%(Descendant)s>::setInflectionPoint")
    @PYB11const
    def setInflectionPoint(self,
                           val = "double"):
        "Set the volume normalization constant."
        return "void"

    @PYB11ignore
    @PYB11cppname("volumeNormalization")
    def getvolumeNormalization(self):
        return "double"

    @PYB11ignore
    @PYB11cppname("kernelExtent")
    def getkernelExtent(self):
        return "double"

    @PYB11ignore
    @PYB11cppname("inflectionPoint")
    def getinflectionPoint(self):
        return "double"

    #...........................................................................
    # Properties
    volumeNormalization = property(getvolumeNormalization, doc="The volume normalization to ensure the integral of the kennel is unity.")
    kernelExtent = property(getkernelExtent, doc="The spatial radius of the kernel in dimensionless (x/h) units.")
    inflectionPoint = property(getinflectionPoint, doc="The inflection point radius in dimensionless (x/h) units.")

#-------------------------------------------------------------------------------
# BSpline
#-------------------------------------------------------------------------------
PYB11template("Dimension")
class BSplineKernel(Kernel):

    def pyinit(self):
        "Default constructor"

#-------------------------------------------------------------------------------
# W4Spline
#-------------------------------------------------------------------------------
PYB11template("Dimension")
class W4SplineKernel(Kernel):

    def pyinit(self):
        "Default constructor"

#-------------------------------------------------------------------------------
# Gaussian
#-------------------------------------------------------------------------------
PYB11template("Dimension")
class GaussianKernel(Kernel):

    def pyinit(self,
               extent = "double"):
        "Constructor with the given eta cutoff"

#-------------------------------------------------------------------------------
# SuperGaussian
#-------------------------------------------------------------------------------
PYB11template("Dimension")
class SuperGaussianKernel(Kernel):

    def pyinit(self):
        "Default constructor"

#-------------------------------------------------------------------------------
# PiGaussianSpline
#-------------------------------------------------------------------------------
PYB11template("Dimension")
class PiGaussianKernel(Kernel):

    def pyinit(self):
        "Default constructor"

    def pyinit(self,
               K = "double"):
        "Construct with K"

    @PYB11ignore
    @PYB11const
    def getK(self):
        return "double"

    @PYB11ignore
    def setK(self, val="double"):
        return "void"

    K = property(getK, setK, "Set the K constant")

#-------------------------------------------------------------------------------
# Hat
#-------------------------------------------------------------------------------
PYB11template("Dimension")
class HatKernel(Kernel):

    def pyinit(self,
               eta0 = "double",
               W0 = "double"):
        "Construct with (eta0, K)"

    @PYB11ignore
    @PYB11pycppname("eta0")
    @PYB11const
    def geteta0(self):
        return "double"

    @PYB11ignore
    @PYB11cppname("eta0")
    def seteta0(self, val="double"):
        return "void"
    
    @PYB11ignore
    @PYB11pycppname("W0")
    @PYB11const
    def getW0(self):
        return "double"

    @PYB11ignore
    @PYB11cppname("W0")
    def setW0(self, val="double"):
        return "void"
    
    eta0 = property(geteta0, seteta0, "eta0")
    W0 = property(geteta0, seteta0, "W0")

#-------------------------------------------------------------------------------
# Sinc
#-------------------------------------------------------------------------------
PYB11template("Dimension")
class SincKernel(Kernel):

    def pyinit(self, extent="double"):
        "Construct with the given extent in eta"

#-------------------------------------------------------------------------------
# NSincPolynomial
#-------------------------------------------------------------------------------
PYB11template("Dimension")
class NSincPolynomialKernel(Kernel):

    def pyinit(self, order="int"):
        "Construct with the specified sinc order"

#-------------------------------------------------------------------------------
# NBSpline
#-------------------------------------------------------------------------------
PYB11template("Dimension")
class NBSplineKernel(Kernel):

    def pyinit(self, order="int"):
        "Construct using the given order of b-spline"

    @PYB11const
    def factorial(self, n="int"):
        "Factorial n!"
        return "int"

    @PYB11const
    def binomialCoefficient(self, n="int", m="int"):
        "Compute the binomial coefficient of (n,m)"
        return "int"

    @PYB11const
    def oneSidedPowerFunction(self, s="double", exponent="int"):
        return "double"

    @PYB11ignore
    @PYB11const
    @PYB11pycppname("order")
    def getOrder(self):
        return "int"

    @PYB11ignore
    def setOrder(self, val="int"):
        return "void"

    order = property(getOrder, setOrder, "The order of the kernel")

#-------------------------------------------------------------------------------
# Table
#-------------------------------------------------------------------------------
PYB11template("Dimension")
class TableKernel(Kernel):

    #...........................................................................
    # Constructors
    def pyinit(self,
               kernel = "const BSplineKernel<%(Dimension)s>&",
               numPoints = ("int", "1000"),
               hmult = ("double", "1.0")):
        "Construct with BSpline kernel"

    def pyinita(self,
                kernel = "const W4SplineKernel<%(Dimension)s>&",
                numPoints = ("int", "1000"),
                hmult = ("double", "1.0")):
        "Construct with W4Spline kernel"

    def pyinitb(self,
                kernel = "const GaussianKernel<%(Dimension)s>&",
                numPoints = ("int", "1000"),
                hmult = ("double", "1.0")):
        "Construct with Gaussian kernel"

    def pyinitc(self,
                kernel = "const SuperGaussianKernel<%(Dimension)s>&",
                numPoints = ("int", "1000"),
                hmult = ("double", "1.0")):
        "Construct with SuperGaussian kernel"

    def pyinitd(self,
                kernel = "const PiGaussianKernel<%(Dimension)s>&",
                numPoints = ("int", "1000"),
                hmult = ("double", "1.0")):
        "Construct with PiGaussian kernel"

    def pyinite(self,
                kernel = "const HatKernel<%(Dimension)s>&",
                numPoints = ("int", "1000"),
                hmult = ("double", "1.0")):
        "Construct with Hat kernel"

    def pyinitf(self,
                kernel = "const SincKernel<%(Dimension)s>&",
                numPoints = ("int", "1000"),
                hmult = ("double", "1.0")):
        "Construct with Sinc kernel"

    def pyinitg(self,
                kernel = "const NSincPolynomialKernel<%(Dimension)s>&",
                numPoints = ("int", "1000"),
                hmult = ("double", "1.0")):
        "Construct with NSincPolynomial kernel"

    def pyinith(self,
                kernel = "const QuarticSplineKernel<%(Dimension)s>&",
                numPoints = ("int", "1000"),
                hmult = ("double", "1.0")):
        "Construct with Quartic spline kernel"

    def pyiniti(self,
                kernel = "const QuinticSplineKernel<%(Dimension)s>&",
                numPoints = ("int", "1000"),
                hmult = ("double", "1.0")):
        "Construct with Quintic spline kernel"

    def pyinitj(self,
                kernel = "const NBSplineKernel<%(Dimension)s>&",
                numPoints = ("int", "1000"),
                hmult = ("double", "1.0")):
        "Construct with NBSpline kernel"

    def pyinitk(self,
                kernel = "const WendlandC2Kernel<%(Dimension)s>&",
                numPoints = ("int", "1000"),
                hmult = ("double", "1.0")):
        "Construct with WendlandC2 kernel"

    def pyinitl(self,
                kernel = "const WendlandC4Kernel<%(Dimension)s>&",
                numPoints = ("int", "1000"),
                hmult = ("double", "1.0")):
        "Construct with WendlandC4 kernel"

    def pyinitm(self,
                kernel = "const WendlandC6Kernel<%(Dimension)s>&",
                numPoints = ("int", "1000"),
                hmult = ("double", "1.0")):
        "Construct with WendlandC6 kernel"

    #...........................................................................
    # Methods
    @PYB11const
    def kernelAndGradValues(self,
                            etaMagnitude = "double",
                            Hdet = "double"):
        return "std::pair<double, double>"

    @PYB11pycppname("kernelAndGradValues")
    @PYB11const
    def kernelAndGradValues1(self,
                             etaMagnitude = "const std::vector<double>&",
                             Hdets = "const std::vector<double>&",
                             kernelValues = "std::vector<double>&",
                             gradValues = "std::vector<double>&"):
        return "void"

    @PYB11const
    def equivalentNodesPerSmoothingScale(self,
                                         Wsum = "double"):
        "Compute the nPerh that corresponds to the Wsum value"
        return "double"

    @PYB11const
    def equivalentWsum(self,
                       nPerh = "double"):
        "Compute the Wsum that corresponds to the  nPerh value"
        return "double"

    def f1(self):
        "Look up the f1 (RZ) Garcia-Senz factor"
        return

    def f2(self):
        "Look up the f2 (RZ) Garcia-Senz factor"
        return

    def gradf1(self):
        "Look up the gradf1 (RZ) Garcia-Senz factor"
        return

    def gradf2(self):
        "Look up the gradf2 (RZ) Garcia-Senz factor"
        return

    def f1Andf2(self):
        "Look up the f1 & f2 (RZ) Garcia-Senz factors simultaneously"

    @PYB11const
    def lowerBound(self, etaMagnitude="double"):
        "Find the lower bound index in a lookup table"
        return "int"

    #...........................................................................
    # Properties
    @PYB11pycppname("nperhValues")
    @PYB11const
    def getnperhValues(self):
        return "std::vector<double>"

    @PYB11pycppname("WsumValues")
    @PYB11const
    def getWsumValues(self):
        return "std::vector<double>"

    @PYB11pycppname("numPoints")
    @PYB11const
    def getnumPoints(self):
        return "int"

    @PYB11pycppname("stepSize")
    @PYB11const
    def getstepSize(self):
        return "double"

    @PYB11pycppname("stepSizeInv")
    @PYB11const
    def getstepSizeInv(self):
        return "double"

    nperhValues = property(getnperhValues, doc="The lookup table used for finding nperh")
    WsumValues = property(getWsumValues, doc="The lookup table of Wsum values")
    numPoints = property(getnumPoints, doc="The number of points in the table")
    stepSize = property(getstepSize, doc="The step size in eta for the lookup table")
    stepSizeInv = property(getstepSizeInv, doc="1/(step size) in eta for the lookup table")

#-------------------------------------------------------------------------------
# WendlandC2
#-------------------------------------------------------------------------------
PYB11template("Dimension")
class WendlandC2Kernel(Kernel):

    def pyinit(self):
        "Default constructor"

#-------------------------------------------------------------------------------
# WendlandC4
#-------------------------------------------------------------------------------
PYB11template("Dimension")
class WendlandC4Kernel(Kernel):

    def pyinit(self):
        "Default constructor"

#-------------------------------------------------------------------------------
# WendlandC6
#-------------------------------------------------------------------------------
PYB11template("Dimension")
class WendlandC6Kernel(Kernel):

    def pyinit(self):
        "Default constructor"

#-------------------------------------------------------------------------------
# QuarticSpline
#-------------------------------------------------------------------------------
PYB11template("Dimension")
class QuarticSplineKernel(Kernel):

    def pyinit(self):
        "Default constructor"

#-------------------------------------------------------------------------------
# QuinticSpline
#-------------------------------------------------------------------------------
PYB11template("Dimension")
class QuinticSplineKernel(Kernel):

    def pyinit(self):
        "Default constructor"

#-------------------------------------------------------------------------------
# ExpInv
#-------------------------------------------------------------------------------
PYB11template("Dimension")
class ExpInvKernel(Kernel):

    def pyinit(self):
        "Default constructor"

