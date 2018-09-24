import inspect
from PYB11Generator import *

#-------------------------------------------------------------------------------
# FieldListSet
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
class FieldListSet:

    # Constructors
    def pyinit(self,):
        "Default constructor"

    # Attributes
    @PYB11readwrite
    def ScalarFieldLists(self):
        "The FieldList<Dim, double> set"

    @PYB11readwrite
    def VectorFieldLists(self):
        "The FieldList<Dim, Vector> set"

    @PYB11readwrite
    def TensorFieldLists(self):
        "The FieldList<Dim, Tensor> set"

    @PYB11readwrite
    def SymTensorFieldLists(self):
        "The FieldList<Dim, SymTensor> set"
