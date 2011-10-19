from pybindgen import *

import sys
sys.path.append("..")
from PBGutils import *
from CXXTypesModule import generateStdVectorBindings
from ref_return_value import *

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class Field:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod):

        # Includes.
        mod.add_include('"Field/FieldTypes.hh"')
    
        # Namespace.
        Spheral = mod.add_cpp_namespace("Spheral")
        space = Spheral.add_cpp_namespace("FieldSpace")

        # Expose types.
        self.FieldBase1d = addObject(space, "FieldBase1d", allow_subclassing=True)
        self.FieldBase2d = addObject(space, "FieldBase2d", allow_subclassing=True)
        self.FieldBase3d = addObject(space, "FieldBase3d", allow_subclassing=True)

        self.FieldListBase = addObject(space, "FieldListBase", allow_subclassing=True)

        self.FieldStorageType = space.add_enum("FieldStorageType", ["Reference", "Copy"], outer_class=self.FieldListBase)

        self.IntField1d =             addObject(space, "IntField1d", parent=self.FieldBase1d)
        self.ULLField1d =             addObject(space, "ULLField1d", parent=self.FieldBase1d)
        self.ScalarField1d =          addObject(space, "ScalarField1d", parent=self.FieldBase1d)
        self.VectorField1d =          addObject(space, "VectorField1d", parent=self.FieldBase1d)
        self.Vector3dField1d =        addObject(space, "Vector3dField1d", parent=self.FieldBase1d)
        self.TensorField1d =          addObject(space, "TensorField1d", parent=self.FieldBase1d)
        self.SymTensorField1d =       addObject(space, "SymTensorField1d", parent=self.FieldBase1d)
        self.ThirdRankTensorField1d = addObject(space, "ThirdRankTensorField1d", parent=self.FieldBase1d)
        self.VectorDoubleField1d =    addObject(space, "VectorDoubleField1d", parent=self.FieldBase1d)
        self.VectorVectorField1d =    addObject(space, "VectorVectorField1d", parent=self.FieldBase1d)
        self.VectorTensorField1d =    addObject(space, "VectorTensorField1d", parent=self.FieldBase1d)
        self.VectorSymTensorField1d = addObject(space, "VectorSymTensorField1d", parent=self.FieldBase1d)

        self.IntField2d =             addObject(space, "IntField2d", parent=self.FieldBase2d)
        self.ULLField2d =             addObject(space, "ULLField2d", parent=self.FieldBase2d)
        self.ScalarField2d =          addObject(space, "ScalarField2d", parent=self.FieldBase2d)
        self.VectorField2d =          addObject(space, "VectorField2d", parent=self.FieldBase2d)
        self.Vector3dField2d =        addObject(space, "Vector3dField2d", parent=self.FieldBase2d)
        self.TensorField2d =          addObject(space, "TensorField2d", parent=self.FieldBase2d)
        self.SymTensorField2d =       addObject(space, "SymTensorField2d", parent=self.FieldBase2d)
        self.ThirdRankTensorField2d = addObject(space, "ThirdRankTensorField2d", parent=self.FieldBase2d)
        self.VectorDoubleField2d =    addObject(space, "VectorDoubleField2d", parent=self.FieldBase2d)
        self.VectorVectorField2d =    addObject(space, "VectorVectorField2d", parent=self.FieldBase2d)
        self.VectorTensorField2d =    addObject(space, "VectorTensorField2d", parent=self.FieldBase2d)
        self.VectorSymTensorField2d = addObject(space, "VectorSymTensorField2d", parent=self.FieldBase2d)

        self.IntField3d =             addObject(space, "IntField3d", parent=self.FieldBase3d)
        self.ULLField3d =             addObject(space, "ULLField3d", parent=self.FieldBase3d)
        self.ScalarField3d =          addObject(space, "ScalarField3d", parent=self.FieldBase3d)
        self.VectorField3d =          addObject(space, "VectorField3d", parent=self.FieldBase3d)
        self.Vector3dField3d =        addObject(space, "Vector3dField3d", parent=self.FieldBase3d)
        self.TensorField3d =          addObject(space, "TensorField3d", parent=self.FieldBase3d)
        self.SymTensorField3d =       addObject(space, "SymTensorField3d", parent=self.FieldBase3d)
        self.ThirdRankTensorField3d = addObject(space, "ThirdRankTensorField3d", parent=self.FieldBase3d)
        self.VectorDoubleField3d =    addObject(space, "VectorDoubleField3d", parent=self.FieldBase3d)
        self.VectorVectorField3d =    addObject(space, "VectorVectorField3d", parent=self.FieldBase3d)
        self.VectorTensorField3d =    addObject(space, "VectorTensorField3d", parent=self.FieldBase3d)
        self.VectorSymTensorField3d = addObject(space, "VectorSymTensorField3d", parent=self.FieldBase3d)

        self.IntFieldList1d =             addObject(space, "IntFieldList1d", parent=self.FieldListBase)
        self.ULLFieldList1d =             addObject(space, "ULLFieldList1d", parent=self.FieldListBase)
        self.ScalarFieldList1d =          addObject(space, "ScalarFieldList1d", parent=self.FieldListBase)
        self.VectorFieldList1d =          addObject(space, "VectorFieldList1d", parent=self.FieldListBase)
        self.Vector3dFieldList1d =        addObject(space, "Vector3dFieldList1d", parent=self.FieldListBase)
        self.TensorFieldList1d =          addObject(space, "TensorFieldList1d", parent=self.FieldListBase)
        self.SymTensorFieldList1d =       addObject(space, "SymTensorFieldList1d", parent=self.FieldListBase)
        self.ThirdRankTensorFieldList1d = addObject(space, "ThirdRankTensorFieldList1d", parent=self.FieldListBase)
        self.VectorDoubleFieldList1d =    addObject(space, "VectorDoubleFieldList1d", parent=self.FieldListBase)
        self.VectorVectorFieldList1d =    addObject(space, "VectorVectorFieldList1d", parent=self.FieldListBase)
        self.VectorTensorFieldList1d =    addObject(space, "VectorTensorFieldList1d", parent=self.FieldListBase)
        self.VectorSymTensorFieldList1d = addObject(space, "VectorSymTensorFieldList1d", parent=self.FieldListBase)

        self.IntFieldList2d =             addObject(space, "IntFieldList2d", parent=self.FieldListBase)
        self.ULLFieldList2d =             addObject(space, "ULLFieldList2d", parent=self.FieldListBase)
        self.ScalarFieldList2d =          addObject(space, "ScalarFieldList2d", parent=self.FieldListBase)
        self.VectorFieldList2d =          addObject(space, "VectorFieldList2d", parent=self.FieldListBase)
        self.Vector3dFieldList2d =        addObject(space, "Vector3dFieldList2d", parent=self.FieldListBase)
        self.TensorFieldList2d =          addObject(space, "TensorFieldList2d", parent=self.FieldListBase)
        self.SymTensorFieldList2d =       addObject(space, "SymTensorFieldList2d", parent=self.FieldListBase)
        self.ThirdRankTensorFieldList2d = addObject(space, "ThirdRankTensorFieldList2d", parent=self.FieldListBase)
        self.VectorDoubleFieldList2d =    addObject(space, "VectorDoubleFieldList2d", parent=self.FieldListBase)
        self.VectorVectorFieldList2d =    addObject(space, "VectorVectorFieldList2d", parent=self.FieldListBase)
        self.VectorTensorFieldList2d =    addObject(space, "VectorTensorFieldList2d", parent=self.FieldListBase)
        self.VectorSymTensorFieldList2d = addObject(space, "VectorSymTensorFieldList2d", parent=self.FieldListBase)

        self.IntFieldList3d =             addObject(space, "IntFieldList3d", parent=self.FieldListBase)
        self.ULLFieldList3d =             addObject(space, "ULLFieldList3d", parent=self.FieldListBase)
        self.ScalarFieldList3d =          addObject(space, "ScalarFieldList3d", parent=self.FieldListBase)
        self.VectorFieldList3d =          addObject(space, "VectorFieldList3d", parent=self.FieldListBase)
        self.Vector3dFieldList3d =        addObject(space, "Vector3dFieldList3d", parent=self.FieldListBase)
        self.TensorFieldList3d =          addObject(space, "TensorFieldList3d", parent=self.FieldListBase)
        self.SymTensorFieldList3d =       addObject(space, "SymTensorFieldList3d", parent=self.FieldListBase)
        self.ThirdRankTensorFieldList3d = addObject(space, "ThirdRankTensorFieldList3d", parent=self.FieldListBase)
        self.VectorDoubleFieldList3d =    addObject(space, "VectorDoubleFieldList3d", parent=self.FieldListBase)
        self.VectorVectorFieldList3d =    addObject(space, "VectorVectorFieldList3d", parent=self.FieldListBase)
        self.VectorTensorFieldList3d =    addObject(space, "VectorTensorFieldList3d", parent=self.FieldListBase)
        self.VectorSymTensorFieldList3d = addObject(space, "VectorSymTensorFieldList3d", parent=self.FieldListBase)

        self.FieldListSet1d = addObject(space, "FieldListSet1d", allow_subclassing=True)
        self.FieldListSet2d = addObject(space, "FieldListSet2d", allow_subclassing=True)
        self.FieldListSet3d = addObject(space, "FieldListSet3d", allow_subclassing=True)

        # std::vector<Field*>
        for element in ("Int", "Scalar", "Vector", "Tensor", "SymTensor"):
            for dim in ("1d", "2d", "3d"):
                exec("""
self.vector_of_%(element)sFieldPtr%(dim)s = addObject(mod, "vector_of_%(element)sFieldPtr%(dim)s", allow_subclassing=True)
self.vector_of_%(element)sFieldList%(dim)s = addObject(mod, "vector_of_%(element)sFieldList%(dim)s", allow_subclassing=True)
""" % {"element" : element, "dim" : dim})

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):
        self.addFieldBaseMethods(self.FieldBase1d, 1)
        self.addFieldBaseMethods(self.FieldBase2d, 2)
        self.addFieldBaseMethods(self.FieldBase3d, 3)

        self.addFieldMethods(self.IntField1d, "int", "IntField1d", 1, applyNumberMethods=True)
        self.addFieldMethods(self.ULLField1d, "unsigned long long", "ULLField1d", 1, applyNumberMethods=True)
        self.addFieldMethods(self.ScalarField1d, "double", "ScalarField1d", 1, applyNumberMethods=True)
        self.addFieldMethods(self.VectorField1d, "Vector1d", "VectorField1d", 1, applyNumberMethods=True, indexAsPointer=True)
        self.addFieldMethods(self.Vector3dField1d, "Geom3Vector", "Vector3dField1d", 1, applyNumberMethods=True, indexAsPointer=True)
        self.addFieldMethods(self.TensorField1d, "Tensor1d", "TensorField1d", 1, applyNumberMethods=True, indexAsPointer=True)
        self.addFieldMethods(self.SymTensorField1d, "SymTensor1d", "SymTensorField1d", 1, applyNumberMethods=True, indexAsPointer=True)
        self.addFieldMethods(self.ThirdRankTensorField1d, "ThirdRankTensor1d", "ThirdRankTensorField1d", 1, applyNumberMethods=True, indexAsPointer=True)
        self.addFieldMethods(self.VectorDoubleField1d, "vector_of_double", "VectorDoubleField1d", 1, indexAsPointer=True)
        self.addFieldMethods(self.VectorVectorField1d, "vector_of_Vector1d", "VectorVectorField1d", 1, indexAsPointer=True)
        self.addFieldMethods(self.VectorTensorField1d, "vector_of_Tensor1d", "VectorTensorField1d", 1, indexAsPointer=True)
        self.addFieldMethods(self.VectorSymTensorField1d, "vector_of_SymTensor1d", "VectorSymTensorField1d", 1, indexAsPointer=True)

        self.addFieldMethods(self.IntField2d, "int", "IntField2d", 2, applyNumberMethods=True)
        self.addFieldMethods(self.ULLField2d, "unsigned long long", "ULLField2d", 2, applyNumberMethods=True)
        self.addFieldMethods(self.ScalarField2d, "double", "ScalarField2d", 2, applyNumberMethods=True)
        self.addFieldMethods(self.VectorField2d, "Vector2d", "VectorField2d", 2, applyNumberMethods=True, indexAsPointer=True)
        self.addFieldMethods(self.Vector3dField2d, "Geom3Vector", "Vector3dField2d", 2, applyNumberMethods=True, indexAsPointer=True)
        self.addFieldMethods(self.TensorField2d, "Tensor2d", "TensorField2d", 2, applyNumberMethods=True, indexAsPointer=True)
        self.addFieldMethods(self.SymTensorField2d, "SymTensor2d", "SymTensorField2d", 2, applyNumberMethods=True, indexAsPointer=True)
        self.addFieldMethods(self.ThirdRankTensorField2d, "ThirdRankTensor2d", "ThirdRankTensorField2d", 2, applyNumberMethods=True, indexAsPointer=True)
        self.addFieldMethods(self.VectorDoubleField2d, "vector_of_double", "VectorDoubleField2d", 2, indexAsPointer=True)
        self.addFieldMethods(self.VectorVectorField2d, "vector_of_Vector2d", "VectorVectorField2d", 2, indexAsPointer=True)
        self.addFieldMethods(self.VectorTensorField2d, "vector_of_Tensor2d", "VectorTensorField2d", 2, indexAsPointer=True)
        self.addFieldMethods(self.VectorSymTensorField2d, "vector_of_SymTensor2d", "VectorSymTensorField2d", 2, indexAsPointer=True)

        self.addFieldMethods(self.IntField3d, "int", "IntField3d", 3, applyNumberMethods=True)
        self.addFieldMethods(self.ULLField3d, "unsigned long long", "ULLField3d", 3, applyNumberMethods=True)
        self.addFieldMethods(self.ScalarField3d, "double", "ScalarField3d", 3, applyNumberMethods=True)
        self.addFieldMethods(self.VectorField3d, "Vector3d", "VectorField3d", 3, applyNumberMethods=True, indexAsPointer=True)
        self.addFieldMethods(self.Vector3dField3d, "Geom3Vector", "Vector3dField3d", 3, applyNumberMethods=True, indexAsPointer=True)
        self.addFieldMethods(self.TensorField3d, "Tensor3d", "TensorField3d", 3, applyNumberMethods=True, indexAsPointer=True)
        self.addFieldMethods(self.SymTensorField3d, "SymTensor3d", "SymTensorField3d", 3, applyNumberMethods=True, indexAsPointer=True)
        self.addFieldMethods(self.ThirdRankTensorField3d, "ThirdRankTensor3d", "ThirdRankTensorField3d", 3, applyNumberMethods=True, indexAsPointer=True)
        self.addFieldMethods(self.VectorDoubleField3d, "vector_of_double", "VectorDoubleField3d", 3, indexAsPointer=True)
        self.addFieldMethods(self.VectorVectorField3d, "vector_of_Vector3d", "VectorVectorField3d", 3, indexAsPointer=True)
        self.addFieldMethods(self.VectorTensorField3d, "vector_of_Tensor3d", "VectorTensorField3d", 3, indexAsPointer=True)
        self.addFieldMethods(self.VectorSymTensorField3d, "vector_of_SymTensor3d", "VectorSymTensorField3d", 3, indexAsPointer=True)

        self.addFieldListMethods(self.IntFieldList1d, "int", "IntFieldList1d", 1, applyNumberMethods=True)
        self.addFieldListMethods(self.ULLFieldList1d, "unsigned long long", "ULLFieldList1d", 1, applyNumberMethods=True)
        self.addFieldListMethods(self.ScalarFieldList1d, "double", "ScalarFieldList1d", 1, applyNumberMethods=True)
        self.addFieldListMethods(self.VectorFieldList1d, "Vector1d", "VectorFieldList1d", 1, applyNumberMethods=True, indexAsPointer=True)
        self.addFieldListMethods(self.Vector3dFieldList1d, "Geom3Vector", "Vector3dFieldList1d", 1, applyNumberMethods=True, indexAsPointer=True)
        self.addFieldListMethods(self.TensorFieldList1d, "Tensor1d", "TensorFieldList1d", 1, applyNumberMethods=True, indexAsPointer=True)
        self.addFieldListMethods(self.SymTensorFieldList1d, "SymTensor1d", "SymTensorFieldList1d", 1, applyNumberMethods=True, indexAsPointer=True)
        self.addFieldListMethods(self.ThirdRankTensorFieldList1d, "ThirdRankTensor1d", "ThirdRankTensorFieldList1d", 1, applyNumberMethods=True, indexAsPointer=True)
        self.addFieldListMethods(self.VectorDoubleFieldList1d, "vector_of_double", "VectorDoubleFieldList1d", 1, indexAsPointer=True)
        self.addFieldListMethods(self.VectorVectorFieldList1d, "vector_of_Vector1d", "VectorVectorFieldList1d", 1, indexAsPointer=True)
        self.addFieldListMethods(self.VectorTensorFieldList1d, "vector_of_Tensor1d", "VectorTensorFieldList1d", 1, indexAsPointer=True)
        self.addFieldListMethods(self.VectorSymTensorFieldList1d, "vector_of_SymTensor1d", "VectorSymTensorFieldList1d", 1, indexAsPointer=True)

        self.addFieldListMethods(self.IntFieldList2d, "int", "IntFieldList2d", 2, applyNumberMethods=True)
        self.addFieldListMethods(self.ULLFieldList2d, "unsigned long long", "ULLFieldList2d", 2, applyNumberMethods=True)
        self.addFieldListMethods(self.ScalarFieldList2d, "double", "ScalarFieldList2d", 2, applyNumberMethods=True)
        self.addFieldListMethods(self.VectorFieldList2d, "Vector2d", "VectorFieldList2d", 2, applyNumberMethods=True, indexAsPointer=True)
        self.addFieldListMethods(self.Vector3dFieldList2d, "Geom3Vector", "Vector3dFieldList2d", 2, applyNumberMethods=True, indexAsPointer=True)
        self.addFieldListMethods(self.TensorFieldList2d, "Tensor2d", "TensorFieldList2d", 2, applyNumberMethods=True, indexAsPointer=True)
        self.addFieldListMethods(self.SymTensorFieldList2d, "SymTensor2d", "SymTensorFieldList2d", 2, applyNumberMethods=True, indexAsPointer=True)
        self.addFieldListMethods(self.ThirdRankTensorFieldList2d, "ThirdRankTensor2d", "ThirdRankTensorFieldList2d", 2, applyNumberMethods=True, indexAsPointer=True)
        self.addFieldListMethods(self.VectorDoubleFieldList2d, "vector_of_double", "VectorDoubleFieldList2d", 2, indexAsPointer=True)
        self.addFieldListMethods(self.VectorVectorFieldList2d, "vector_of_Vector2d", "VectorVectorFieldList2d", 2, indexAsPointer=True)
        self.addFieldListMethods(self.VectorTensorFieldList2d, "vector_of_Tensor2d", "VectorTensorFieldList2d", 2, indexAsPointer=True)
        self.addFieldListMethods(self.VectorSymTensorFieldList2d, "vector_of_SymTensor2d", "VectorSymTensorFieldList2d", 2, indexAsPointer=True)

        self.addFieldListMethods(self.IntFieldList3d, "int", "IntFieldList3d", 3, applyNumberMethods=True)
        self.addFieldListMethods(self.ULLFieldList3d, "unsigned long long", "ULLFieldList3d", 3, applyNumberMethods=True)
        self.addFieldListMethods(self.ScalarFieldList3d, "double", "ScalarFieldList3d", 3, applyNumberMethods=True)
        self.addFieldListMethods(self.VectorFieldList3d, "Vector3d", "VectorFieldList3d", 3, applyNumberMethods=True, indexAsPointer=True)
        self.addFieldListMethods(self.Vector3dFieldList3d, "Geom3Vector", "Vector3dFieldList3d", 3, applyNumberMethods=True, indexAsPointer=True)
        self.addFieldListMethods(self.TensorFieldList3d, "Tensor3d", "TensorFieldList3d", 3, applyNumberMethods=True, indexAsPointer=True)
        self.addFieldListMethods(self.SymTensorFieldList3d, "SymTensor3d", "SymTensorFieldList3d", 3, applyNumberMethods=True, indexAsPointer=True)
        self.addFieldListMethods(self.ThirdRankTensorFieldList3d, "ThirdRankTensor3d", "ThirdRankTensorFieldList3d", 3, applyNumberMethods=True, indexAsPointer=True)
        self.addFieldListMethods(self.VectorDoubleFieldList3d, "vector_of_double", "VectorDoubleFieldList3d", 3, indexAsPointer=True)
        self.addFieldListMethods(self.VectorVectorFieldList3d, "vector_of_Vector3d", "VectorVectorFieldList3d", 3, indexAsPointer=True)
        self.addFieldListMethods(self.VectorTensorFieldList3d, "vector_of_Tensor3d", "VectorTensorFieldList3d", 3, indexAsPointer=True)
        self.addFieldListMethods(self.VectorSymTensorFieldList3d, "vector_of_SymTensor3d", "VectorSymTensorFieldList3d", 3, indexAsPointer=True)

        self.addFieldListSetMethods(self.FieldListSet1d, 1)
        self.addFieldListSetMethods(self.FieldListSet2d, 2)
        self.addFieldListSetMethods(self.FieldListSet3d, 3)

        # std::vector<Field*> and std::vector<FieldList>
        for element in ["Int", "Scalar", "Vector", "Tensor", "SymTensor"]:
            for dim in ["1d", "2d", "3d"]:
                exec("""
generateStdVectorBindings(self.vector_of_%(element)sFieldPtr%(dim)s, "Spheral::FieldSpace::%(element)sField%(dim)s*", "vector_of_%(element)sFieldPtr%(dim)s", indexAsPointer=True)
generateStdVectorBindings(self.vector_of_%(element)sFieldList%(dim)s, "Spheral::FieldSpace::%(element)sFieldList%(dim)s", "vector_of_%(element)sFieldList%(dim)s", indexAsPointer=True)
""" % {"element" : element, "dim" : dim})

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return ["FieldSpace"]

    #---------------------------------------------------------------------------
    # Add methods (Fields).
    #---------------------------------------------------------------------------
    def addFieldBaseMethods(self, x, ndim):

        # Object names.
        me = "Spheral::FieldSpace::FieldBase%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        fieldbase = "Spheral::FieldSpace::FieldBase%id" % ndim
        intfield = "Spheral::FieldSpace::IntField%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        vectorfield = "Spheral::FieldSpace::VectorField%id" % ndim
        vector3dfield = "Spheral::FieldSpace::Vector3dField%id" % ndim
        tensorfield = "Spheral::FieldSpace::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::FieldSpace::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::FieldSpace::VectorDoubleField%id" % ndim
        vectorvectorfield = "Spheral::FieldSpace::VectorVectorField%id" % ndim
        vectorsymtensorfield = "Spheral::FieldSpace::VectorSymTensorField%id" % ndim
        symtensorfield = "Spheral::FieldSpace::SymTensorField%id" % ndim
        intfieldlist = "Spheral::FieldSpace::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::FieldSpace::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::FieldSpace::VectorFieldList%id" % ndim
        vector3dfieldlist = "Spheral::FieldSpace::Vector3dFieldList%id" % ndim
        tensorfieldlist = "Spheral::FieldSpace::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::FieldSpace::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::FieldSpace::ThirdRankTensorFieldList%id" % ndim
        vectordoublefieldlist = "Spheral::FieldSpace::VectorDoubleFieldList%id" % ndim
        vectorvectorfieldlist = "Spheral::FieldSpace::VectorVectorFieldList%id" % ndim
        vectorsymtensorfieldlist = "Spheral::FieldSpace::VectorSymTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeSpace::NodeList%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        database = "Spheral::DataBaseSpace::DataBase%id" % ndim
        connectivitymap = "Spheral::NeighborSpace::ConnectivityMap%id" % ndim
        key = "pair_NodeList%id_string" % ndim
        vectorkeys = "vector_of_pair_NodeList%id_string" % ndim

        # # Constructors.
        # x.add_constructor([param("std::string", "name")])
        # x.add_constructor([param("std::string", "name"),
        #                    constrefparam(nodelist, "nodeList")])

        # Methods.
        const_ref_return_value(x, me, "%s::nodeList" % me, nodelist, [], "nodeList")
        x.add_method("size", "int", [], is_const=True, is_pure_virtual=True)
        x.add_method("Zero", None, [], is_pure_virtual=True)
        x.add_method("setNodeList", None, [constrefparam(nodelist, "nodeList")], is_pure_virtual=True)
        x.add_method("resizeField", None, [param("int", "size")], is_pure_virtual=True)
        x.add_method("resizeFieldInternal", None, [param("int", "size"),
                                                   param("int", "oldFirstGhostNode")], is_pure_virtual=True)
        x.add_method("resizeFieldGhost", None, [param("int", "size")], is_pure_virtual=True)
        x.add_method("deleteElement", None, [param("int", "nodeID")], is_pure_virtual=True)
        x.add_method("deleteElements", None, [constrefparam("vector_of_int", "nodeIDs")], is_pure_virtual=True)
        x.add_method("packValues", "vector_of_char", [constrefparam("vector_of_int", "nodeIDs")], is_const=True, is_pure_virtual=True)
        x.add_method("unpackValues", None, [param("int", "numElements"),
                                            param("int", "beginInsertionIndex"),
                                            constrefparam("vector_of_char", "buffer")], is_pure_virtual=True)

        # Attributes.
        x.add_instance_attribute("name", "std::string", getter="name", setter="name")

        return

    #---------------------------------------------------------------------------
    # Add methods (Fields).
    #---------------------------------------------------------------------------
    def addFieldMethods(self, x, val, me, ndim, applyNumberMethods=False, indexAsPointer=False):

        # Object names.
        me = "Spheral::FieldSpace::%s" % me
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        fieldbase = "Spheral::FieldSpace::FieldBase%id" % ndim
        intfield = "Spheral::FieldSpace::IntField%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        vectorfield = "Spheral::FieldSpace::VectorField%id" % ndim
        vector3dfield = "Spheral::FieldSpace::Vector3dField%id" % ndim
        tensorfield = "Spheral::FieldSpace::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::FieldSpace::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::FieldSpace::VectorDoubleField%id" % ndim
        vectorvectorfield = "Spheral::FieldSpace::VectorVectorField%id" % ndim
        vectorsymtensorfield = "Spheral::FieldSpace::VectorSymTensorField%id" % ndim
        symtensorfield = "Spheral::FieldSpace::SymTensorField%id" % ndim
        intfieldlist = "Spheral::FieldSpace::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::FieldSpace::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::FieldSpace::VectorFieldList%id" % ndim
        vector3dfieldlist = "Spheral::FieldSpace::Vector3dFieldList%id" % ndim
        tensorfieldlist = "Spheral::FieldSpace::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::FieldSpace::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::FieldSpace::ThirdRankTensorFieldList%id" % ndim
        vectordoublefieldlist = "Spheral::FieldSpace::VectorDoubleFieldList%id" % ndim
        vectorvectorfieldlist = "Spheral::FieldSpace::VectorVectorFieldList%id" % ndim
        vectorsymtensorfieldlist = "Spheral::FieldSpace::VectorSymTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeSpace::NodeList%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        database = "Spheral::DataBaseSpace::DataBase%id" % ndim
        connectivitymap = "Spheral::NeighborSpace::ConnectivityMap%id" % ndim
        key = "pair_NodeList%id_string" % ndim
        vectorkeys = "vector_of_pair_NodeList%id_string" % ndim
        if val == "unsigned long long":
            vector_of_value = "vector_of_ULL"
        else:
            vector_of_value = "vector_of_%s" % val

        # Constructors.
        x.add_constructor([param("std::string", "name")])
        x.add_constructor([param("std::string", "name"),
                           constrefparam(me, "field")])
        x.add_constructor([param("std::string", "name"),
                           constrefparam(nodelist, "nodeList")])
        x.add_constructor([param("std::string", "name"),
                           constrefparam(nodelist, "nodeList"),
                           param(val, "value")])

        # Methods.
        x.add_method("size", "int", [], is_const=True, custom_name="__len__")
        x.add_method("Zero", None, [], is_virtual=True)
        x.add_method("valid", "bool", [], is_const=True)
        x.add_method("internalValues", vector_of_value, [], is_const=True)
        x.add_method("ghostValues", vector_of_value, [], is_const=True)
        x.add_method("allValues", vector_of_value, [], is_const=True)

        # Comparison operators.
        x.add_binary_comparison_operator("==")
        x.add_binary_comparison_operator("!=")
        x.add_binary_comparison_operator("<")
        x.add_binary_comparison_operator(">")
        x.add_binary_comparison_operator("<=")
        x.add_binary_comparison_operator(">=")

        # These methods can only be used with Field values that can be used as numbers.
        if applyNumberMethods:
            x.add_method("applyMin", None, [param(val, "dataMin")])
            x.add_method("applyMax", None, [param(val, "dataMax")])
            x.add_method("sumElements", val, [], is_const=True)
            x.add_method("min", val, [], is_const=True)
            x.add_method("max", val, [], is_const=True)
            x.add_method("localSumElements", val, [], is_const=True)
            x.add_method("localMin", val, [], is_const=True)
            x.add_method("localMax", val, [], is_const=True)
            x.add_method("string", "std::string", [param("int", "precision", default_value="20")], is_const=True)
            x.add_method("string", None, [constrefparam("std::string", "string")])

            # We allow a few special methods for Symmetric tensors.
            if "SymTensor" in val:
                x.add_method("applyScalarMin", None, [param("double", "dataMin")])
                x.add_method("applyScalarMax", None, [param("double", "dataMax")])

        # Math operators.
        # These don't work yet, mostly 'cause the pybindgen implementation assumes that
        # the classes will be default constructable.
#         x.add_binary_numeric_operator("+")
#         x.add_binary_numeric_operator("-")
#         x.add_inplace_numeric_operator("+=")
#         x.add_inplace_numeric_operator("-=")
#         x.add_binary_numeric_operator("+", right_cppclass=val)
#         x.add_binary_numeric_operator("-", right_cppclass=val)
#         x.add_inplace_numeric_operator("+=", right_cppclass=val)
#         x.add_inplace_numeric_operator("-=", right_cppclass=val)

        # Indexing.
        if indexAsPointer:
            x.add_function_as_method("indexFieldAsPointer",
                                     retval(ptr(val), reference_existing_object=True),
                                     [param(me, "self"), param("int", "index")],
                                     template_parameters = [me],
                                     custom_name = "__getitem__")
        else:
            x.add_function_as_method("indexField",
                                     val,
                                     [param(me, "self"), param("int", "index")],
                                     template_parameters = [me],
                                     custom_name = "__getitem__")
        x.add_function_as_method("assignToFieldIndex", None, [param(me, "self"),
                                                              param("int", "index"),
                                                              param(val, "value")],
                                 template_parameters = [me],
                                 custom_name = "__setitem__")

        # The virtual methods from FieldBase.
        x.add_method("size", "int", [], is_const=True, is_virtual=True)
        x.add_method("Zero", None, [], is_virtual=True)
        x.add_method("setNodeList", None, [constrefparam(nodelist, "nodeList")], is_virtual=True)
        x.add_method("resizeField", None, [param("int", "size")], is_virtual=True)
        x.add_method("resizeFieldInternal", None, [param("int", "size"),
                                                   param("int", "oldFirstGhostNode")], is_virtual=True)
        x.add_method("resizeFieldGhost", None, [param("int", "size")], is_virtual=True)
        x.add_method("deleteElement", None, [param("int", "nodeID")], is_virtual=True)
        x.add_method("deleteElements", None, [constrefparam("vector_of_int", "nodeIDs")], is_virtual=True)
        x.add_method("packValues", "vector_of_char", [constrefparam("vector_of_int", "nodeIDs")], is_const=True, is_virtual=True)
        x.add_method("unpackValues", None, [param("int", "numElements"),
                                            param("int", "beginInsertionIndex"),
                                            constrefparam("vector_of_char", "buffer")], is_virtual=True)

        # Attributes.
        x.add_instance_attribute("numElements", "int", is_const=True, getter="numElements")
        x.add_instance_attribute("numInternalElements", "int", is_const=True, getter="numInternalElements")

        return

    #---------------------------------------------------------------------------
    # Add methods (FieldLists).
    #---------------------------------------------------------------------------
    def addFieldListMethods(self, x, val, me, ndim, applyNumberMethods=False, indexAsPointer=False):

        # Object names.
        me = "Spheral::FieldSpace::%s" % me
        field = me.replace("List", "")
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        fieldbase = "Spheral::FieldSpace::FieldBase%id" % ndim
        intfield = "Spheral::FieldSpace::IntField%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        vectorfield = "Spheral::FieldSpace::VectorField%id" % ndim
        vector3dfield = "Spheral::FieldSpace::Vector3dField%id" % ndim
        tensorfield = "Spheral::FieldSpace::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::FieldSpace::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::FieldSpace::VectorDoubleField%id" % ndim
        vectorvectorfield = "Spheral::FieldSpace::VectorVectorField%id" % ndim
        vectorsymtensorfield = "Spheral::FieldSpace::VectorSymTensorField%id" % ndim
        symtensorfield = "Spheral::FieldSpace::SymTensorField%id" % ndim
        intfieldlist = "Spheral::FieldSpace::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::FieldSpace::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::FieldSpace::VectorFieldList%id" % ndim
        vector3dfieldlist = "Spheral::FieldSpace::Vector3dFieldList%id" % ndim
        tensorfieldlist = "Spheral::FieldSpace::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::FieldSpace::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::FieldSpace::ThirdRankTensorFieldList%id" % ndim
        vectordoublefieldlist = "Spheral::FieldSpace::VectorDoubleFieldList%id" % ndim
        vectorvectorfieldlist = "Spheral::FieldSpace::VectorVectorFieldList%id" % ndim
        vectorsymtensorfieldlist = "Spheral::FieldSpace::VectorSymTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeSpace::NodeList%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        database = "Spheral::DataBaseSpace::DataBase%id" % ndim
        connectivitymap = "Spheral::NeighborSpace::ConnectivityMap%id" % ndim
        key = "pair_NodeList%id_string" % ndim
        vectorkeys = "vector_of_pair_NodeList%id_string" % ndim

        # Constructors.
        x.add_constructor([])
        x.add_constructor([param("FieldStorageType", "aStorageType")])

        # Methods.
        x.add_method("copyFields", None, [])
        x.add_method("haveField", "bool", [constrefparam(field, "field")], is_const=True)
        x.add_method("haveNodeList", "bool", [constrefparam(nodelist, "nodeList")], is_const=True)
        x.add_method("assignFields", None, [constrefparam(me, "fieldList")])
        x.add_method("appendField", None, [constrefparam(field, "field")])
        x.add_method("deleteField", None, [constrefparam(field, "field")])
        x.add_method("appendNewField", None, [param("std::string", "name"), constrefparam(nodelist, "nodeList"), param(val, "value")])
        x.add_method("setMasterNodeLists", None, [constrefparam(vector, "r"), constrefparam(symtensor, "H")], is_const=True)
        x.add_method("setMasterNodeLists", None, [constrefparam(vector, "r")], is_const=True)
        x.add_method("setRefineNodeLists", None, [constrefparam(vector, "r"), constrefparam(symtensor, "H")], is_const=True)
        x.add_method("setRefineNodeLists", None, [constrefparam(vector, "r")], is_const=True)
        x.add_method("Zero", None, [])
        x.add_method("size", "int", [], is_const=True, custom_name="__len__")

        # Comparison operators.
        x.add_binary_comparison_operator("==")
        x.add_binary_comparison_operator("!=")
        x.add_binary_comparison_operator("<")
        x.add_binary_comparison_operator(">")
        x.add_binary_comparison_operator("<=")
        x.add_binary_comparison_operator(">=")

        # These methods can only be used with Field values that can be used as numbers.
        if applyNumberMethods:
            x.add_method("applyMin", None, [param(val, "dataMin")])
            x.add_method("applyMax", None, [param(val, "dataMax")])
            x.add_method("sumElements", val, [], is_const=True)
            x.add_method("min", val, [], is_const=True)
            x.add_method("max", val, [], is_const=True)
            x.add_method("localSumElements", val, [], is_const=True)
            x.add_method("localMin", val, [], is_const=True)
            x.add_method("localMax", val, [], is_const=True)

            # We allow a few special methods for Symmetric tensors.
            if "SymTensor" in val:
                x.add_method("applyScalarMin", None, [param("double", "dataMin")])
                x.add_method("applyScalarMax", None, [param("double", "dataMax")])

        # Indexing.
        x.add_function_as_method("indexFieldList",
                                 retval(ptr(field), reference_existing_object=True),
                                 [param(me, "self"), param("int", "index")],
                                 template_parameters = [me],
                                 custom_name = "__getitem__")
        x.add_function_as_method("fieldForNodeList",
                                 retval(ptr(field), reference_existing_object=True),
                                 [param(me, "self"), 
                                  constrefparam(nodelist, "nodeList")],
                                 template_parameters = [dim, me],
                                 custom_name = "fieldForNodeList")
        if indexAsPointer:
            x.add_function_as_method("indexFieldListForValuePointer",
                                     retval(ptr(val), reference_existing_object=True),
                                     [param(me, "self"), param("int", "fieldIndex"), param("int", "nodeIndex")],
                                     template_parameters = [me],
                                     custom_name = "__call__")
        else:
            x.add_function_as_method("indexFieldListForValue",
                                     val,
                                     [param(me, "self"), param("int", "fieldIndex"), param("int", "nodeIndex")],
                                     template_parameters = [me],
                                     custom_name = "__call__")

        # Attributes.
        x.add_instance_attribute("storageType", "FieldStorageType", is_const=True, getter="storageType")
        x.add_instance_attribute("numFields", "int", is_const=True, getter="numFields")
        x.add_instance_attribute("numNodes", "int", is_const=True, getter="numNodes")
        x.add_instance_attribute("numInternalNodes", "int", is_const=True, getter="numInternalNodes")
        x.add_instance_attribute("numGhostNodes", "int", is_const=True, getter="numGhostNodes")

        return

    #---------------------------------------------------------------------------
    # Add methods (FieldListSet).
    #---------------------------------------------------------------------------
    def addFieldListSetMethods(self, x, ndim):

        # Object names.
        vector_of_scalarfieldlist = "vector_of_ScalarFieldList%id" % ndim
        vector_of_vectorfieldlist = "vector_of_VectorFieldList%id" % ndim
        vector_of_tensorfieldlist = "vector_of_TensorFieldList%id" % ndim
        vector_of_symtensorfieldlist = "vector_of_SymTensorFieldList%id" % ndim

        # Constructors.
        x.add_constructor([])

        # Attributes.
        x.add_instance_attribute("ScalarFieldLists", retval(ptr(vector_of_scalarfieldlist), reference_existing_object=True), getter="ScalarFieldListPtrs", is_const=True)
        x.add_instance_attribute("VectorFieldLists", retval(ptr(vector_of_vectorfieldlist), reference_existing_object=True), getter="VectorFieldListPtrs", is_const=True)
        x.add_instance_attribute("TensorFieldLists", retval(ptr(vector_of_tensorfieldlist), reference_existing_object=True), getter="TensorFieldListPtrs", is_const=True)
        x.add_instance_attribute("SymTensorFieldLists", retval(ptr(vector_of_symtensorfieldlist), reference_existing_object=True), getter="SymTensorFieldListPtrs", is_const=True)

        return

