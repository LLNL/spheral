#-------------------------------------------------------------------------------
# ReproducingKernelMethods
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11template("Dimension")
class ReproducingKernelMethods:
    """Provides the reproducing kernel methods, analogous to the Kernel class for SPH

This is really just a convenient front-end for the methods in RKUtilities"""
    
    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    typedef typename Eigen::SparseMatrix<double> TransformationMatrix;
"""

    def pyinit(self,
               order = "const RKOrder"):
        "Constructor"
        
    @PYB11const
    def transformationMatrix(T = "const Tensor&",
                             needHessian = "const bool"):
        "Compute the transformation matrix to apply to the RK coefficients for the given Tensor"
        return "TransformationMatrix"

    @PYB11const
    def applyTransformation(self,
                            T = "const TransformationMatrix&",
                            corrections = "RKCoefficients<%(Dimension)s>&"):
        "Apply the transformation T to the corrections"
        return "void"

    #..........................................................................
    # Attributes
    order = PYB11property(doc="order to which we are enforcing reproducibility")
    gradCorrectionsSize = PYB11property(doc="The size of the RKCoefficients for corrections + grad")
    hessCorrectionsSize = PYB11property(doc="The size of the RKCoefficients for corrections + grad + hessian")
