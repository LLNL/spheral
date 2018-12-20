import inspect
from PYB11Generator import *

#-------------------------------------------------------------------------------
# FieldListSet
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralField")
class FieldListSet:

    # Constructors
    def pyinit(self,):
        "Default constructor"

    # Attributes
    # ScalarFieldLists = PYB11readwrite(doc="The FieldList<Dim, double> set")
    # VectorFieldLists = PYB11readwrite(doc="The FieldList<Dim, Vector> set")
    # TensorFieldLists = PYB11readwrite(doc="The FieldList<Dim, Tensor> set")
    # SymTensorFieldLists = PYB11readwrite(doc="The FieldList<Dim, SymTensor> set")

    ScalarFieldLists = PYB11property("std::vector<FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>>*",
                                     getterraw = "[](FieldListSet<%(Dimension)s>& self) { std::cerr << "<" << self.ScalarFieldLists.size() << ">" << std::endl; return &(self.ScalarFieldLists); }",
                                     setterraw = "[](FieldListSet<%(Dimension)s>& self, std::vector<FieldList<%(Dimension)s, %(Dimension)s::Scalar>>* x) { self.ScalarFieldLists = *x; }",
                                     returnpolicy = "reference")
    VectorFieldLists = PYB11property("std::vector<FieldList<%(Dimension)s, typename %(Dimension)s::Vector>>*",
                                     getterraw = "[](FieldListSet<%(Dimension)s>& self) { return &(self.VectorFieldLists); }",
                                     setterraw = "[](FieldListSet<%(Dimension)s>& self, std::vector<FieldList<%(Dimension)s, %(Dimension)s::Vector>>* x) { self.VectorFieldLists = *x; }",
                                     returnpolicy = "reference")
    TensorFieldLists = PYB11property("std::vector<FieldList<%(Dimension)s, typename %(Dimension)s::Tensor>>*",
                                     getterraw = "[](FieldListSet<%(Dimension)s>& self) { return &(self.TensorFieldLists); }",
                                     setterraw = "[](FieldListSet<%(Dimension)s>& self, std::vector<FieldList<%(Dimension)s, %(Dimension)s::Tensor>>* x) { self.TensorFieldLists = *x; }",
                                     returnpolicy = "reference")
    SymTensorFieldLists = PYB11property("std::vector<FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>>*",
                                        getterraw = "[](FieldListSet<%(Dimension)s>& self) { return &(self.SymTensorFieldLists); }",
                                        setterraw = "[](FieldListSet<%(Dimension)s>& self, std::vector<FieldList<%(Dimension)s, %(Dimension)s::SymTensor>>* x) { self.SymTensorFieldLists = *x; }",
                                        returnpolicy = "reference")


    # VectorFieldLists = PYB11property("std::vector<FieldList<%(Dimension)s, typename %(Dimension)s::Vector>>&",
    #                                  getterraw = "[](FieldListSet<%(Dimension)s>& self) { return self.VectorFieldLists; }",
    #                                  setterraw = "[](FieldListSet<%(Dimension)s>& self, std::vector<FieldList<%(Dimension)s, %(Dimension)s::Vector>>& x) { self.VectorFieldLists = x; }",
    #                                  returnpolicy = "reference_internal")
    # TensorFieldLists = PYB11property("std::vector<FieldList<%(Dimension)s, typename %(Dimension)s::Tensor>>&",
    #                                  getterraw = "[](FieldListSet<%(Dimension)s>& self) { return self.TensorFieldLists; }",
    #                                  setterraw = "[](FieldListSet<%(Dimension)s>& self, std::vector<FieldList<%(Dimension)s, %(Dimension)s::Tensor>>& x) { self.TensorFieldLists = x; }",
    #                                  returnpolicy = "reference_internal")
    # SymTensorFieldLists = PYB11property("std::vector<FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>>&",
    #                                     getterraw = "[](FieldListSet<%(Dimension)s>& self) { return self.SymTensorFieldLists; }",
    #                                     setterraw = "[](FieldListSet<%(Dimension)s>& self, std::vector<FieldList<%(Dimension)s, %(Dimension)s::SymTensor>>& x) { self.SymTensorFieldLists = x; }",
    #                                     returnpolicy = "reference_internal")
