from pybindgen import *

from PBGutils import *

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class Kernel:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir, dims):

        self.dims = dims

        # Includes.
        mod.add_include('"%s/KernelTypes.hh"' % srcdir)
    
        # Namespace.
        space = mod.add_cpp_namespace("Spheral")

        # Expose types.
        self.types = ("BSpline", "W4Spline", "Gaussian", "SuperGaussian", "PiGaussian",
                      "Hat", "Sinc", "NSincPolynomial", "NBSpline", "QuarticSpline",
                      "QuinticSpline", "Table", "WendlandC2", "WendlandC4", "WendlandC6", "ExpInv")
        for type in self.types:
            for ndim in self.dims:
                dim = "%id" % ndim
                name = type + "Kernel" + dim
                exec('self.%(name)s = addObject(space, "%(name)s", allow_subclassing=True)' % {"name" : name})
        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        for ndim in self.dims:
            dim = "%id" % ndim

            # Generic Kernel types.
            for type in ("BSpline", "W4Spline", "SuperGaussian", "WendlandC2", "WendlandC4", "WendlandC6", "QuarticSpline", "QuinticSpline", "ExpInv"):
                name = type + "Kernel" + dim
                exec("self.generateDefaultKernelBindings(self.%s, %i)" % (name, ndim))
            # Now some specialized bindings for kernels.
            exec("""
self.generateGaussianKernelBindings(self.GaussianKernel%(dim)s, %(ndim)i)
self.generatePiGaussianKernelBindings(self.PiGaussianKernel%(dim)s, %(ndim)i)
self.generateHatKernelBindings(self.HatKernel%(dim)s, %(ndim)i)
self.generateSincKernelBindings(self.SincKernel%(dim)s, %(ndim)i)
self.generateNSincPolynomialKernelBindings(self.NSincPolynomialKernel%(dim)s, %(ndim)i)
self.generateNBSplineKernelBindings(self.NBSplineKernel%(dim)s, %(ndim)i)
self.generateTableKernelBindings(self.TableKernel%(dim)s, %(ndim)i)
""" % {"dim" : dim, "ndim" : ndim})

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return ["KernelSpace"]

    #---------------------------------------------------------------------------
    # Default (kernels with just the default constructor).
    #---------------------------------------------------------------------------
    def generateDefaultKernelBindings(self, x, ndim):
        x.add_constructor([])
        self.generateGenericKernelBindings(x, ndim)
        return

    #---------------------------------------------------------------------------
    # Gaussian
    #---------------------------------------------------------------------------
    def generateGaussianKernelBindings(self, x, ndim):
        x.add_constructor([param("double", "extent")])
        self.generateGenericKernelBindings(x, ndim)
        return

    #---------------------------------------------------------------------------
    # PiGaussian
    #---------------------------------------------------------------------------
    def generatePiGaussianKernelBindings(self, x, ndim):
        x.add_constructor([])
        x.add_constructor([param("double", "K")])
        x.add_instance_attribute("K", "double", getter="getK", setter="setK")
        self.generateGenericKernelBindings(x, ndim)
        return

    #---------------------------------------------------------------------------
    # Hat
    #---------------------------------------------------------------------------
    def generateHatKernelBindings(self, x, ndim):
        x.add_constructor([param("double", "eta0"), param("double", "W0")])
        x.add_instance_attribute("eta0", "double", getter="eta0", is_const=True)
        x.add_instance_attribute("W0", "double", getter="W0", is_const=True)
        self.generateGenericKernelBindings(x, ndim)
        return

    #---------------------------------------------------------------------------
    # Sinc
    #---------------------------------------------------------------------------
    def generateSincKernelBindings(self, x, ndim):
        x.add_constructor([param("double", "extent")])
        self.generateGenericKernelBindings(x, ndim)
        return

    #---------------------------------------------------------------------------
    # NSincPolynomial
    #---------------------------------------------------------------------------
    def generateNSincPolynomialKernelBindings(self, x, ndim):
        x.add_constructor([param("int", "order")])
        self.generateGenericKernelBindings(x, ndim)
        return

    #---------------------------------------------------------------------------
    # NBSpline
    #---------------------------------------------------------------------------
    def generateNBSplineKernelBindings(self, x, ndim):
        x.add_constructor([param("int", "order")])
        x.add_method("factorial", "int", [param("int", "n")], is_const=True)
        x.add_method("binomialCoefficient", "int", [param("int", "n"), param("int", "m")], is_const=True)
        x.add_method("oneSidedPowerFunction", "double", [param("double", "s"), param("int", "exponent")], is_const=True)
        x.add_instance_attribute("order", "int", getter="order", setter="setOrder")
        self.generateGenericKernelBindings(x, ndim)
        return

    #---------------------------------------------------------------------------
    # TableKernel
    #---------------------------------------------------------------------------
    def generateTableKernelBindings(self, x, ndim):

        bsplinekernel = "Spheral::BSplineKernel%id" % ndim
        w4splinekernel = "Spheral::W4SplineKernel%id" % ndim
        gaussiankernel = "Spheral::GaussianKernel%id" % ndim
        supergaussiankernel = "Spheral::SuperGaussianKernel%id" % ndim
        pigaussiankernel = "Spheral::PiGaussianKernel%id" % ndim
        hatkernel = "Spheral::HatKernel%id" % ndim
        sinckernel = "Spheral::SincKernel%id" % ndim
        nsincpolynomialkernel = "Spheral::NSincPolynomialKernel%id" % ndim
        quarticsplinekernel = "Spheral::QuarticSplineKernel%id" % ndim
        quinticsplinekernel = "Spheral::QuinticSplineKernel%id" % ndim
        nbsplinekernel = "Spheral::NBSplineKernel%id" % ndim
        wendlandc2kernel = "Spheral::WendlandC2Kernel%id" % ndim
        wendlandc4kernel = "Spheral::WendlandC4Kernel%id" % ndim
        wendlandc6kernel = "Spheral::WendlandC6Kernel%id" % ndim

        # Constructors.
        for W in (bsplinekernel, w4splinekernel, gaussiankernel, supergaussiankernel, pigaussiankernel,
                  hatkernel, sinckernel, nsincpolynomialkernel, quarticsplinekernel, quinticsplinekernel, nbsplinekernel, 
                  wendlandc2kernel,wendlandc4kernel,wendlandc6kernel):
            x.add_constructor([constrefparam(W, "kernel"),
                               param("int", "numPoints", default_value="1000"),
                               param("double", "hmult", default_value="1.0")])
            #x.add_method("augment", None, [constrefparam(W, "W")])

        # Methods.
        x.add_method("kernelAndGradValue", "pair_double_double", [param("double", "etaMagnitude"), param("double", "Hdet")], is_const=True)
        x.add_method("kernelAndGradValues", None, [constrefparam("vector_of_double", "etaMagnitudes"),
                                                   constrefparam("vector_of_double", "Hdets"),
                                                   refparam("vector_of_double", "kernelValues"),
                                                   refparam("vector_of_double", "gradValues"),], is_const=True)
        x.add_method("equivalentNodesPerSmoothingScale", "double", [param("double", "Wsum")], is_const=True)
        x.add_method("equivalentWsum", "double", [param("double", "nPerh")], is_const=True)
        x.add_method("f1", "double", [param("double", "etaMagnitude")], is_const=True)
        x.add_method("f2", "double", [param("double", "etaMagnitude")], is_const=True)
        x.add_method("gradf1", "double", [param("double", "etaMagnitude")], is_const=True)
        x.add_method("gradf2", "double", [param("double", "etaMagnitude")], is_const=True)
        x.add_method("f1Andf2", None, [param("double", "etaMagnitude"),
                                       refparam("double", "f1"),
                                       refparam("double", "f2"),
                                       refparam("double", "gradf1"),
                                       refparam("double", "gradf2")], is_const=True)
        x.add_method("lowerBound", "int", [param("double", "etaMagnitude")], is_const=True)
        x.add_method("valid", "bool", [], is_const=True, is_virtual=True)

        # Attributes.
        x.add_instance_attribute("nperhValues", "vector_of_double", getter="nperhValues", is_const=True)
        x.add_instance_attribute("WsumValues", "vector_of_double", getter="WsumValues", is_const=True)
        x.add_instance_attribute("numPoints", "int", getter="numPoints", is_const=True)
        x.add_instance_attribute("stepSize", "double", getter="stepSize", is_const=True)
        x.add_instance_attribute("stepSizeInv", "double", getter="stepSizeInv", is_const=True)

        # Generic methods.
        self.generateGenericKernelBindings(x, ndim)

    #---------------------------------------------------------------------------
    # Add generic Kernel methods.
    #---------------------------------------------------------------------------
    def generateGenericKernelBindings(self, x, ndim):

        # Objects.
        vector = "Vector%id" % ndim
        symtensor = "SymTensor%id" % ndim

        # Methods.
        x.add_method("operator()", "double", [param("double", "etaMagnitude"), param("double", "Hdet")], is_const=True, custom_name = "__call__")
        x.add_method("operator()", "double", [constrefparam(vector, "eta"), param("double", "Hdet")], is_const=True, custom_name="__call__")
        x.add_method("operator()", "double", [param("double", "etaMagnitude"), param(symtensor, "H")], is_const=True, custom_name="__call__")
        x.add_method("operator()", "double", [constrefparam(vector, "eta"), param(symtensor, "H")], is_const=True, custom_name="__call__")

        x.add_method("grad", "double", [param("double", "etaMagnitude"), param("double", "Hdet")], is_const=True)
        x.add_method("grad", "double", [constrefparam(vector, "eta"), param("double", "Hdet")], is_const=True)
        x.add_method("grad", "double", [param("double", "etaMagnitude"), param(symtensor, "H")], is_const=True)
        x.add_method("grad", "double", [constrefparam(vector, "eta"), param(symtensor, "H")], is_const=True)

        x.add_method("grad2", "double", [param("double", "etaMagnitude"), param("double", "Hdet")], is_const=True)
        x.add_method("grad2", "double", [constrefparam(vector, "eta"), param("double", "Hdet")], is_const=True)
        x.add_method("grad2", "double", [param("double", "etaMagnitude"), param(symtensor, "H")], is_const=True)
        x.add_method("grad2", "double", [constrefparam(vector, "eta"), param(symtensor, "H")], is_const=True)

        x.add_method("gradh", "double", [param("double", "etaMagnitude"), param("double", "Hdet")], is_const=True)
        x.add_method("gradh", "double", [constrefparam(vector, "eta"), param("double", "Hdet")], is_const=True)
        x.add_method("gradh", "double", [param("double", "etaMagnitude"), param(symtensor, "H")], is_const=True)
        x.add_method("gradh", "double", [constrefparam(vector, "eta"), param(symtensor, "H")], is_const=True)

        x.add_method("kernelValue", "double", [param("double", "etaMagnitude"), param("double", "etaMagnitude")], is_const=True)
        x.add_method("gradValue", "double", [param("double", "etaMagnitude"), param("double", "etaMagnitude")], is_const=True)
        x.add_method("grad2Value", "double", [param("double", "etaMagnitude"), param("double", "etaMagnitude")], is_const=True)
        x.add_method("gradhValue", "double", [param("double", "etaMagnitude"), param("double", "etaMagnitude")], is_const=True)

        x.add_method("valid", "bool", [], is_const=True, is_virtual=True)

        x.add_method("setVolumeNormalization", None, [param("double", "value")], visibility="protected")
        x.add_method("setKernelExtent", None, [param("double", "value")], visibility="protected")
        x.add_method("setInflectionPoint", None, [param("double", "value")], visibility="protected")

        # Attributes.
        x.add_instance_attribute("volumeNormalization", "double", getter="volumeNormalization", is_const=True)
        x.add_instance_attribute("kernelExtent", "double", getter="kernelExtent", is_const=True)
        x.add_instance_attribute("inflectionPoint", "double", getter="inflectionPoint", is_const=True)

        return

