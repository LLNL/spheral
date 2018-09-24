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

    def pyinit(self):
        "Default constructor"

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

    @PYB11cppname("kernelValue")
    @PYB11const
    def kernelValue1(self,
              etaMagnitude = "double",
              H = "const SymTensor&"):
        "Return the value of the kernel"
        return "double"

    @PYB11cppname("kernelValue")
    @PYB11const
    def kernelValue2(self,
              eta = "const Vector&",
              Hdet = "const double&"):
        "Return the value of the kernel"
        return "double"

    @PYB11cppname("kernelValue")
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

    @PYB11cppname("gradValue")
    @PYB11const
    def gradValue1(self,
              etaMagnitude = "double",
              H = "const SymTensor&"):
        "Return the value of the gradient of the kernel"
        return "double"

    @PYB11cppname("gradValue")
    @PYB11const
    def gradValue2(self,
              eta = "const Vector&",
              Hdet = "const double&"):
        "Return the value of the gradient of the kernel"
        return "double"

    @PYB11cppname("gradValue")
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

    @PYB11cppname("grad2Value")
    @PYB11const
    def grad2Value1(self,
              etaMagnitude = "double",
              H = "const SymTensor&"):
        "Return the value of the second derivative of the kernel"
        return "double"

    @PYB11cppname("grad2Value")
    @PYB11const
    def grad2Value2(self,
              eta = "const Vector&",
              Hdet = "const double&"):
        "Return the value of the second derivative of the kernel"
        return "double"

    @PYB11cppname("grad2Value")
    @PYB11const
    def grad2Value3(self,
              eta = "const Vector&",
              H = "const SymTensor&"):
        "Return the value of the second derivative of the kernel"
        return "double"

    #...........................................................................
    @PYB11cppname("grad")
    @PYB11const
    def grad_0(self,
               etaMagnitude = "double",
               Hdet = "const double&"):
        "Return the gradient of the kernel"
        return "double"

    @PYB11cppname("grad")
    @PYB11const
    def grad_1(self,
               etaMagnitude = "double",
               H = "const SymTensor&"):
        "Return the gradient of the kernel"
        return "double"

    @PYB11cppname("grad")
    @PYB11const
    def grad_2(self,
               eta = "const Vector&",
               Hdet = "const double&"):
        "Return the gradient of the kernel"
        return "double"

    @PYB11cppname("grad")
    @PYB11const
    def grad_3(self,
               eta = "const Vector&",
               H = "const SymTensor&"):
        "Return the gradient of the kernel"
        return "double"

    #...........................................................................
    @PYB11cppname("grad2")
    @PYB11const
    def grad2_0(self,
               etaMagnitude = "double",
               Hdet = "const double&"):
        "Return the second derivative of the kernel"
        return "double"

    @PYB11cppname("grad2")
    @PYB11const
    def grad2_1(self,
               etaMagnitude = "double",
               H = "const SymTensor&"):
        "Return the second derivative of the kernel"
        return "double"

    @PYB11cppname("grad2")
    @PYB11const
    def grad2_2(self,
               eta = "const Vector&",
               Hdet = "const double&"):
        "Return the second derivative of the kernel"
        return "double"

    @PYB11cppname("grad2")
    @PYB11const
    def grad2_3(self,
               eta = "const Vector&",
               H = "const SymTensor&"):
        "Return the second derivative of the kernel"
        return "double"

    #...........................................................................
    @PYB11cppname("gradh")
    @PYB11const
    def gradh_0(self,
                etaMagnitude = "double",
                Hdet = "const double&"):
        "Return the gradient with respect to h of the kernel"
        return "double"

    @PYB11cppname("gradh")
    @PYB11const
    def gradh_1(self,
               etaMagnitude = "double",
               H = "const SymTensor&"):
        "Return the gradient with respect to h of the kernel"
        return "double"

    @PYB11cppname("gradh")
    @PYB11const
    def gradh_2(self,
               eta = "const Vector&",
               Hdet = "const double&"):
        "Return the gradient with respect to h of the kernel"
        return "double"

    @PYB11cppname("gradh")
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
