from pybindgen import *

import sys
sys.path.extend(["..", "CXXTypes"])
from PBGutils import *
from CXXTypesModule import generateStdVectorBindings
from ref_return_value import *

#---------------------------------------------------------------------------
# Boundary abstract bindings.
#---------------------------------------------------------------------------
def generateBoundaryVirtualBindings(x, ndim, pureVirtual):

    # Object names.
    vector = "Vector%id" % ndim
    tensor = "Tensor%id" % ndim
    symtensor = "SymTensor%id" % ndim
    intfield = "Spheral::FieldSpace::IntField%id" % ndim
    scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
    vectorfield = "Spheral::FieldSpace::VectorField%id" % ndim
    tensorfield = "Spheral::FieldSpace::TensorField%id" % ndim
    thirdranktensorfield = "Spheral::FieldSpace::ThirdRankTensorField%id" % ndim
    vectordoublefield = "Spheral::FieldSpace::VectorDoubleField%id" % ndim
    symtensorfield = "Spheral::FieldSpace::SymTensorField%id" % ndim
    intfieldlist = "Spheral::FieldSpace::IntFieldList%id" % ndim
    scalarfieldlist = "Spheral::FieldSpace::ScalarFieldList%id" % ndim
    vectorfieldlist = "Spheral::FieldSpace::VectorFieldList%id" % ndim
    tensorfieldlist = "Spheral::FieldSpace::TensorFieldList%id" % ndim
    symtensorfieldlist = "Spheral::FieldSpace::SymTensorFieldList%id" % ndim
    thirdranktensorfieldlist = "Spheral::FieldSpace::ThirdRankTensorFieldList%id" % ndim
    nodelist = "Spheral::NodeSpace::NodeList%id" % ndim
    state = "Spheral::State%id" % ndim
    derivatives = "Spheral::StateDerivatives%id" % ndim
    database = "Spheral::DataBaseSpace::DataBase%id" % ndim
    connectivitymap = "Spheral::NeighborSpace::ConnectivityMap%id" % ndim

    # Virtual methods.
    x.add_method("setGhostNodes", None, [refparam(nodelist, "nodeList")], is_virtual=True, is_pure_virtual=pureVirtual)
    x.add_method("updateGhostNodes", None, [refparam(nodelist, "nodeList")], is_virtual=True, is_pure_virtual=pureVirtual)

    x.add_method("setViolationNodes", None, [refparam(nodelist, "nodeList")], is_virtual=True, is_pure_virtual=pureVirtual)
    x.add_method("updateViolationNodes", None, [refparam(nodelist, "nodeList")], is_virtual=True, is_pure_virtual=pureVirtual)

    x.add_method("applyGhostBoundary", None, [refparam(intfield, "field")], is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
    x.add_method("applyGhostBoundary", None, [refparam(scalarfield, "field")], is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
    x.add_method("applyGhostBoundary", None, [refparam(vectorfield, "field")], is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
    x.add_method("applyGhostBoundary", None, [refparam(tensorfield, "field")], is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
    x.add_method("applyGhostBoundary", None, [refparam(symtensorfield, "field")], is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
    x.add_method("applyGhostBoundary", None, [refparam(thirdranktensorfield, "field")], is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
    x.add_method("applyGhostBoundary", None, [refparam(vectordoublefield, "field")], is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)

    x.add_method("enforceBoundary", None, [refparam(intfield, "field")], is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
    x.add_method("enforceBoundary", None, [refparam(scalarfield, "field")], is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
    x.add_method("enforceBoundary", None, [refparam(vectorfield, "field")], is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
    x.add_method("enforceBoundary", None, [refparam(tensorfield, "field")], is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
    x.add_method("enforceBoundary", None, [refparam(symtensorfield, "field")], is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)
    x.add_method("enforceBoundary", None, [refparam(thirdranktensorfield, "field")], is_const=True, is_virtual=True, is_pure_virtual=pureVirtual)

    return

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class Boundary:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod):

        # Includes.
        mod.add_include('"Boundary/BoundaryTypes.hh"')
        mod.add_include('"Boundary/Boundary.hh"')
    
        # Namespace.
        Spheral = mod.add_cpp_namespace("Spheral")
        space = Spheral.add_cpp_namespace("BoundarySpace")

        # Expose types.
        self.Boundary1d = addObject(space, "Boundary1d", allow_subclassing=True)
        self.Boundary2d = addObject(space, "Boundary2d", allow_subclassing=True)
        self.Boundary3d = addObject(space, "Boundary3d", allow_subclassing=True)

        self.PlanarBoundary1d = addObject(space, "PlanarBoundary1d", parent=self.Boundary1d)
        self.PlanarBoundary2d = addObject(space, "PlanarBoundary2d", parent=self.Boundary2d)
        self.PlanarBoundary3d = addObject(space, "PlanarBoundary3d", parent=self.Boundary3d)

        self.ReflectingBoundary1d = addObject(space, "ReflectingBoundary1d", parent=self.PlanarBoundary1d)
        self.ReflectingBoundary2d = addObject(space, "ReflectingBoundary2d", parent=self.PlanarBoundary2d)
        self.ReflectingBoundary3d = addObject(space, "ReflectingBoundary3d", parent=self.PlanarBoundary3d)

        self.RigidBoundary1d = addObject(space, "RigidBoundary1d", parent=self.PlanarBoundary1d)
        self.RigidBoundary2d = addObject(space, "RigidBoundary2d", parent=self.PlanarBoundary2d)
        self.RigidBoundary3d = addObject(space, "RigidBoundary3d", parent=self.PlanarBoundary3d)

        self.PeriodicBoundary1d = addObject(space, "PeriodicBoundary1d", parent=self.PlanarBoundary1d)
        self.PeriodicBoundary2d = addObject(space, "PeriodicBoundary2d", parent=self.PlanarBoundary2d)
        self.PeriodicBoundary3d = addObject(space, "PeriodicBoundary3d", parent=self.PlanarBoundary3d)

        self.ConstantVelocityBoundary1d = addObject(space, "ConstantVelocityBoundary1d", parent=self.Boundary1d)
        self.ConstantVelocityBoundary2d = addObject(space, "ConstantVelocityBoundary2d", parent=self.Boundary2d)
        self.ConstantVelocityBoundary3d = addObject(space, "ConstantVelocityBoundary3d", parent=self.Boundary3d)

        self.ConstantXVelocityBoundary1d = addObject(space, "ConstantXVelocityBoundary1d", parent=self.ConstantVelocityBoundary1d)
        self.ConstantXVelocityBoundary2d = addObject(space, "ConstantXVelocityBoundary2d", parent=self.ConstantVelocityBoundary2d)
        self.ConstantXVelocityBoundary3d = addObject(space, "ConstantXVelocityBoundary3d", parent=self.ConstantVelocityBoundary3d)

        self.ConstantYVelocityBoundary2d = addObject(space, "ConstantYVelocityBoundary2d", parent=self.ConstantVelocityBoundary2d)
        self.ConstantYVelocityBoundary3d = addObject(space, "ConstantYVelocityBoundary3d", parent=self.ConstantVelocityBoundary3d)

        self.ConstantZVelocityBoundary3d = addObject(space, "ConstantZVelocityBoundary3d", parent=self.ConstantVelocityBoundary3d)

        self.ConstantRVelocityBoundary1d = addObject(space, "ConstantRVelocityBoundary1d", parent=self.ConstantVelocityBoundary1d)
        self.ConstantRVelocityBoundary2d = addObject(space, "ConstantRVelocityBoundary2d", parent=self.ConstantVelocityBoundary2d)
        self.ConstantRVelocityBoundary3d = addObject(space, "ConstantRVelocityBoundary3d", parent=self.ConstantVelocityBoundary3d)

        self.ConstantBoundary1d = addObject(space, "ConstantBoundary1d", parent=self.Boundary1d)
        self.ConstantBoundary2d = addObject(space, "ConstantBoundary2d", parent=self.Boundary2d)
        self.ConstantBoundary3d = addObject(space, "ConstantBoundary3d", parent=self.Boundary3d)

        self.SphericalBoundary = addObject(space, "SphericalBoundary", parent=self.Boundary3d)

        self.CylindricalBoundary = addObject(space, "CylindricalBoundary", parent=self.Boundary3d)

        self.vecBound1d = addObject(mod, "vector_of_Boundary1d", allow_subclassing=True)
        self.vecBound2d = addObject(mod, "vector_of_Boundary2d", allow_subclassing=True)
        self.vecBound3d = addObject(mod, "vector_of_Boundary3d", allow_subclassing=True)

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        self.generateBoundaryBindings(self.Boundary1d, 1)
        self.generateBoundaryBindings(self.Boundary2d, 2)
        self.generateBoundaryBindings(self.Boundary3d, 3)

        self.generatePlanarBoundaryBindings(self.PlanarBoundary1d, 1)
        self.generatePlanarBoundaryBindings(self.PlanarBoundary2d, 2)
        self.generatePlanarBoundaryBindings(self.PlanarBoundary3d, 3)

        self.generateReflectingBoundaryBindings(self.ReflectingBoundary1d, 1)
        self.generateReflectingBoundaryBindings(self.ReflectingBoundary2d, 2)
        self.generateReflectingBoundaryBindings(self.ReflectingBoundary3d, 3)

        self.generateRigidBoundaryBindings(self.RigidBoundary1d, 1)
        self.generateRigidBoundaryBindings(self.RigidBoundary2d, 2)
        self.generateRigidBoundaryBindings(self.RigidBoundary3d, 3)

        self.generatePeriodicBoundaryBindings(self.PeriodicBoundary1d, 1)
        self.generatePeriodicBoundaryBindings(self.PeriodicBoundary2d, 2)
        self.generatePeriodicBoundaryBindings(self.PeriodicBoundary3d, 3)

        self.generateConstantVelocityBoundaryBindings(self.ConstantVelocityBoundary1d, 1)
        self.generateConstantVelocityBoundaryBindings(self.ConstantVelocityBoundary2d, 2)
        self.generateConstantVelocityBoundaryBindings(self.ConstantVelocityBoundary3d, 3)

        self.generateConstantXVelocityBoundaryBindings(self.ConstantXVelocityBoundary1d, 1)
        self.generateConstantXVelocityBoundaryBindings(self.ConstantXVelocityBoundary2d, 2)
        self.generateConstantXVelocityBoundaryBindings(self.ConstantXVelocityBoundary3d, 3)

        self.generateConstantYVelocityBoundaryBindings(self.ConstantYVelocityBoundary2d, 2)
        self.generateConstantYVelocityBoundaryBindings(self.ConstantYVelocityBoundary3d, 3)

        self.generateConstantZVelocityBoundaryBindings(self.ConstantZVelocityBoundary3d, 3)

        self.generateConstantBoundaryBindings(self.ConstantBoundary1d, 1)
        self.generateConstantBoundaryBindings(self.ConstantBoundary2d, 2)
        self.generateConstantBoundaryBindings(self.ConstantBoundary3d, 3)

        self.generateSphericalBoundaryBindings(self.SphericalBoundary)

        self.generateCylindricalBoundaryBindings(self.CylindricalBoundary)

        generateStdVectorBindings(self.vecBound1d, "Spheral::BoundarySpace::Boundary1d*", "vector_of_Boundary1d")
        generateStdVectorBindings(self.vecBound2d, "Spheral::BoundarySpace::Boundary2d*", "vector_of_Boundary2d")
        generateStdVectorBindings(self.vecBound3d, "Spheral::BoundarySpace::Boundary3d*", "vector_of_Boundary3d")
        
        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return ["BoundarySpace"]

    #---------------------------------------------------------------------------
    # Boundary bindings.
    #---------------------------------------------------------------------------
    def generateBoundaryBindings(self, x, ndim):

        # Object names.
        me = "Boundary%id" % ndim
        dim = "Spheral::Dim<%i> " % ndim
        boundarynodes = "Spheral::BoundarySpace::%s::BoundaryNodes" % me
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        intfield = "Spheral::FieldSpace::IntField%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        vectorfield = "Spheral::FieldSpace::VectorField%id" % ndim
        tensorfield = "Spheral::FieldSpace::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::FieldSpace::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::FieldSpace::VectorDoubleField%id" % ndim
        symtensorfield = "Spheral::FieldSpace::SymTensorField%id" % ndim
        intfieldlist = "Spheral::FieldSpace::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::FieldSpace::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::FieldSpace::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::FieldSpace::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::FieldSpace::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::FieldSpace::ThirdRankTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeSpace::NodeList%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        database = "Spheral::DataBaseSpace::DataBase%id" % ndim
        connectivitymap = "Spheral::NeighborSpace::ConnectivityMap%id" % ndim

        # Add the subclass.
        bn = x.add_class("BoundaryNodes", allow_subclassing=True)
        bn.add_function_as_method("controlNodesFromBoundaryNodes",
                                  retval(ptr("vector_of_int"), reference_existing_object=True),
                                  [param(boundarynodes, "self")],
                                  template_parameters = [dim],
                                  custom_name = "controlNodes")
        bn.add_function_as_method("ghostNodesFromBoundaryNodes",
                                  retval(ptr("vector_of_int"), reference_existing_object=True),
                                  [param(boundarynodes, "self")],
                                  template_parameters = [dim],
                                  custom_name = "ghostNodes")
        bn.add_function_as_method("violationNodesFromBoundaryNodes",
                                  retval(ptr("vector_of_int"), reference_existing_object=True),
                                  [param(boundarynodes, "self")],
                                  template_parameters = [dim],
                                  custom_name = "violationNodes")
        
        # Constructors.
        x.add_constructor([])

        # Virtual methods.
        x.add_method("setAllGhostNodes", None, [refparam(database, "dataBase")], is_virtual=True)

        x.add_method("setAllViolationNodes", None, [refparam(database, "dataBase")], is_virtual=True)

        x.add_method("enforceFieldListBoundary", None, [refparam(intfieldlist, "fieldList")], custom_template_method_name="enforceIntFieldListBoundary", is_const=True)
        x.add_method("enforceFieldListBoundary", None, [refparam(scalarfieldlist, "fieldList")], custom_template_method_name="enforceScalarFieldListBoundary", is_const=True)
        x.add_method("enforceFieldListBoundary", None, [refparam(vectorfieldlist, "fieldList")], custom_template_method_name="enforceVectorFieldListBoundary", is_const=True)
        x.add_method("enforceFieldListBoundary", None, [refparam(tensorfieldlist, "fieldList")], custom_template_method_name="enforceTensorFieldListBoundary", is_const=True)
        x.add_method("enforceFieldListBoundary", None, [refparam(symtensorfieldlist, "fieldList")], custom_template_method_name="enforceSymTensorFieldListBoundary", is_const=True)
        x.add_method("enforceFieldListBoundary", None, [refparam(thirdranktensorfieldlist, "fieldList")], custom_template_method_name="enforceThirdRankTensorFieldListBoundary", is_const=True)

        x.add_method("finalizeGhostBoundary", None, [], is_const=True, is_virtual=True)
        x.add_method("clip", None, [refparam(vector, "xmin"), refparam(vector, "xmax")], is_const=True, is_virtual=True)

        # x.add_method("meshGhostNodes", "bool", [], is_virtual=True, is_const=True)

        x.add_method("reset", None, [refparam(database, "dataBase")], is_virtual=True)

        x.add_method("cullGhostNodes", None, [constrefparam(intfieldlist, "flagSet"),
                                              refparam(intfieldlist, "old2newIndexMap"),
                                              refparam("vector_of_int", "numNodesRemoved")], is_virtual=True)

        # Methods.
        x.add_method("haveNodeList", "bool", [constrefparam(nodelist, "nodeList")], is_const=True)
        x.add_method("controlNodes", "vector_of_int", [refparam(nodelist, "nodeList")], is_const=True)
        x.add_method("ghostNodes", "vector_of_int", [refparam(nodelist, "nodeList")], is_const=True)
        x.add_method("violationNodes", "vector_of_int", [refparam(nodelist, "nodeList")], is_const=True)

        x.add_method("applyFieldListGhostBoundary", None, [refparam(scalarfieldlist, "fieldList")], custom_template_method_name="applyScalarFieldListGhostBoundary", is_const=True)
        x.add_method("applyFieldListGhostBoundary", None, [refparam(vectorfieldlist, "fieldList")], custom_template_method_name="applyVectorFieldListGhostBoundary", is_const=True)
        x.add_method("applyFieldListGhostBoundary", None, [refparam(tensorfieldlist, "fieldList")], custom_template_method_name="applyTensorFieldListGhostBoundary", is_const=True)
        x.add_method("applyFieldListGhostBoundary", None, [refparam(symtensorfieldlist, "fieldList")], custom_template_method_name="applySymTensorFieldListGhostBoundary", is_const=True)

        # Protected methods.
        x.add_function_as_method("accessBoundaryNodesFromBoundary",
                                 retval(ptr(boundarynodes), reference_existing_object=True),
                                 [param(me, "self"), param(nodelist, "nodeList")],
                                 template_parameters = [dim],
                                 custom_name = "accessBoundaryNodes")
        x.add_method("addNodeList", None, [refparam(nodelist, "nodeList")], is_virtual=True)

        # Attributes.
        x.add_instance_attribute("numGhostNodes", "int", getter="numGhostNodes", is_const=True)

        # Generate the abstract interface.
        generateBoundaryVirtualBindings(x, ndim, True)

        return

    #---------------------------------------------------------------------------
    # PlanarBoundary bindings.
    #---------------------------------------------------------------------------
    def generatePlanarBoundaryBindings(self, x, ndim):

        # Object names.
        me = "PlanarBoundary%id" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        plane = "Plane%id" % ndim
        intfield = "Spheral::FieldSpace::IntField%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        vectorfield = "Spheral::FieldSpace::VectorField%id" % ndim
        tensorfield = "Spheral::FieldSpace::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::FieldSpace::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::FieldSpace::VectorDoubleField%id" % ndim
        symtensorfield = "Spheral::FieldSpace::SymTensorField%id" % ndim
        intfieldlist = "Spheral::FieldSpace::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::FieldSpace::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::FieldSpace::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::FieldSpace::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::FieldSpace::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::FieldSpace::ThirdRankTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeSpace::NodeList%id" % ndim
        state = "State%id" % ndim
        derivatives = "StateDerivatives%id" % ndim
        database = "Spheral::DataBaseSpace::DataBase%id" % ndim
        connectivitymap = "Spheral::NeighborSpace::ConnectivityMap%id" % ndim
        fileio = "Spheral::FileIOSpace::FileIO"

        # Constructors.
        x.add_constructor([])
        x.add_constructor([refparam(plane, "enterPlane"), refparam(plane, "exitPlane")])

        # Attributes.
        x.add_instance_attribute("enterPlane", plane, getter="enterPlane", setter="setEnterPlane")
        x.add_instance_attribute("exitPlane", plane, getter="exitPlane", setter="setExitPlane")

        # Methods.
        x.add_method("setGhostNodes", None, [refparam(nodelist, "nodeList"), refparam("vector_of_int", "presetControlNodes")])
        x.add_method("mapPosition", vector, [refparam(vector, "position"),
                                             refparam(plane, "enterPlane"),
                                             refparam(plane, "exitPlane")], is_const=True)
        
        # Virtual methods.
        x.add_method("valid", "bool", [], is_const=True, is_virtual=True)
        x.add_method("label", "std::string", [], is_const=True, is_virtual=True)
        x.add_method("dumpState", None, [refparam(fileio, "fileIO"),
                                         refparam("std::string", "pathName")], is_const=True, is_virtual=True)
        x.add_method("restoreState", None, [refparam(fileio, "fileIO"),
                                            refparam("std::string", "pathName")], is_virtual=True)

        # Override the abstract interface.
        x.add_method("setGhostNodes", None, [refparam(nodelist, "nodeList")], is_virtual=True)
        x.add_method("updateGhostNodes", None, [refparam(nodelist, "nodeList")], is_virtual=True)
        x.add_method("setViolationNodes", None, [refparam(nodelist, "nodeList")], is_virtual=True)
        x.add_method("updateViolationNodes", None, [refparam(nodelist, "nodeList")], is_virtual=True)

        return

    #---------------------------------------------------------------------------
    # ReflectingBoundary bindings.
    #---------------------------------------------------------------------------
    def generateReflectingBoundaryBindings(self, x, ndim):

        # Object names.
        me = "ReflectingBoundary%id" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        plane = "Plane%id" % ndim
        intfield = "Spheral::FieldSpace::IntField%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        vectorfield = "Spheral::FieldSpace::VectorField%id" % ndim
        tensorfield = "Spheral::FieldSpace::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::FieldSpace::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::FieldSpace::VectorDoubleField%id" % ndim
        symtensorfield = "Spheral::FieldSpace::SymTensorField%id" % ndim
        intfieldlist = "Spheral::FieldSpace::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::FieldSpace::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::FieldSpace::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::FieldSpace::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::FieldSpace::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::FieldSpace::ThirdRankTensorFieldList%id" % ndim
        nodelist = "NodeList%id" % ndim
        state = "State%id" % ndim
        derivatives = "StateDerivatives%id" % ndim
        database = "Spheral::DataBaseSpace::DataBase%id" % ndim
        connectivitymap = "ConnectivityMap%id" % ndim

        # Constructors.
        x.add_constructor([])
        x.add_constructor([refparam(plane, "plane")])

        # Attributes.
        x.add_instance_attribute("reflectOperator", tensor, getter="reflectOperator", is_const=True)

        # Override the abstract interface.
        x.add_method("applyGhostBoundary", None, [refparam(intfield, "field")], is_const=True, is_virtual=True)
        x.add_method("applyGhostBoundary", None, [refparam(scalarfield, "field")], is_const=True, is_virtual=True)
        x.add_method("applyGhostBoundary", None, [refparam(vectorfield, "field")], is_const=True, is_virtual=True)
        x.add_method("applyGhostBoundary", None, [refparam(tensorfield, "field")], is_const=True, is_virtual=True)
        x.add_method("applyGhostBoundary", None, [refparam(symtensorfield, "field")], is_const=True, is_virtual=True)
        x.add_method("applyGhostBoundary", None, [refparam(thirdranktensorfield, "field")], is_const=True, is_virtual=True)
        x.add_method("applyGhostBoundary", None, [refparam(vectordoublefield, "field")], is_const=True, is_virtual=True)

        x.add_method("enforceBoundary", None, [refparam(intfield, "field")], is_const=True, is_virtual=True)
        x.add_method("enforceBoundary", None, [refparam(scalarfield, "field")], is_const=True, is_virtual=True)
        x.add_method("enforceBoundary", None, [refparam(vectorfield, "field")], is_const=True, is_virtual=True)
        x.add_method("enforceBoundary", None, [refparam(tensorfield, "field")], is_const=True, is_virtual=True)
        x.add_method("enforceBoundary", None, [refparam(symtensorfield, "field")], is_const=True, is_virtual=True)
        x.add_method("enforceBoundary", None, [refparam(thirdranktensorfield, "field")], is_const=True, is_virtual=True)

        return

    #---------------------------------------------------------------------------
    # RigidBoundary bindings.
    #---------------------------------------------------------------------------
    def generateRigidBoundaryBindings(self, x, ndim):

        # Object names.
        me = "RigidBoundary%id" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        plane = "Plane%id" % ndim
        intfield = "Spheral::FieldSpace::IntField%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        vectorfield = "Spheral::FieldSpace::VectorField%id" % ndim
        tensorfield = "Spheral::FieldSpace::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::FieldSpace::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::FieldSpace::VectorDoubleField%id" % ndim
        symtensorfield = "Spheral::FieldSpace::SymTensorField%id" % ndim
        intfieldlist = "Spheral::FieldSpace::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::FieldSpace::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::FieldSpace::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::FieldSpace::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::FieldSpace::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::FieldSpace::ThirdRankTensorFieldList%id" % ndim
        nodelist = "NodeList%id" % ndim
        state = "State%id" % ndim
        derivatives = "StateDerivatives%id" % ndim
        database = "Spheral::DataBaseSpace::DataBase%id" % ndim
        connectivitymap = "ConnectivityMap%id" % ndim

        # Constructors.
        x.add_constructor([])
        x.add_constructor([refparam(plane, "plane")])

        # Attributes.
        x.add_instance_attribute("reflectOperator", tensor, getter="reflectOperator", is_const=True)

        # Override the abstract interface.
        x.add_method("applyGhostBoundary", None, [refparam(intfield, "field")], is_const=True, is_virtual=True)
        x.add_method("applyGhostBoundary", None, [refparam(scalarfield, "field")], is_const=True, is_virtual=True)
        x.add_method("applyGhostBoundary", None, [refparam(vectorfield, "field")], is_const=True, is_virtual=True)
        x.add_method("applyGhostBoundary", None, [refparam(tensorfield, "field")], is_const=True, is_virtual=True)
        x.add_method("applyGhostBoundary", None, [refparam(symtensorfield, "field")], is_const=True, is_virtual=True)
        x.add_method("applyGhostBoundary", None, [refparam(thirdranktensorfield, "field")], is_const=True, is_virtual=True)
        x.add_method("applyGhostBoundary", None, [refparam(vectordoublefield, "field")], is_const=True, is_virtual=True)

        x.add_method("enforceBoundary", None, [refparam(intfield, "field")], is_const=True, is_virtual=True)
        x.add_method("enforceBoundary", None, [refparam(scalarfield, "field")], is_const=True, is_virtual=True)
        x.add_method("enforceBoundary", None, [refparam(vectorfield, "field")], is_const=True, is_virtual=True)
        x.add_method("enforceBoundary", None, [refparam(tensorfield, "field")], is_const=True, is_virtual=True)
        x.add_method("enforceBoundary", None, [refparam(symtensorfield, "field")], is_const=True, is_virtual=True)
        x.add_method("enforceBoundary", None, [refparam(thirdranktensorfield, "field")], is_const=True, is_virtual=True)

        return

    #---------------------------------------------------------------------------
    # PeriodicBoundary bindings.
    #---------------------------------------------------------------------------
    def generatePeriodicBoundaryBindings(self, x, ndim):

        # Object names.
        plane = "Plane%id" % ndim

        # Constructors.
        x.add_constructor([])
        x.add_constructor([refparam(plane, "plane1"), refparam(plane, "plane2")])

        # Generate the abstract interface.
        generateBoundaryVirtualBindings(x, ndim, False)

        return

    #---------------------------------------------------------------------------
    # ConstantVelocityBoundary bindings.
    #---------------------------------------------------------------------------
    def generateConstantVelocityBoundaryBindings(self, x, ndim):

        # Object names.
        me = "Spheral::BoundarySpace::ConstantVelocityBoundary%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        plane = "Plane%id" % ndim
        intfield = "Spheral::FieldSpace::IntField%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        vectorfield = "Spheral::FieldSpace::VectorField%id" % ndim
        tensorfield = "Spheral::FieldSpace::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::FieldSpace::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::FieldSpace::VectorDoubleField%id" % ndim
        symtensorfield = "Spheral::FieldSpace::SymTensorField%id" % ndim
        intfieldlist = "Spheral::FieldSpace::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::FieldSpace::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::FieldSpace::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::FieldSpace::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::FieldSpace::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::FieldSpace::ThirdRankTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeSpace::NodeList%id" % ndim
        state = "State%id" % ndim
        derivatives = "StateDerivatives%id" % ndim
        database = "Spheral::DataBaseSpace::DataBase%id" % ndim
        connectivitymap = "Spheral::NeighborSpace::ConnectivityMap%id" % ndim
        fileio = "Spheral::FileIOSpace::FileIO"

        # Constructors.
        x.add_constructor([refparam(nodelist, "nodeList"), 
                           refparam("vector_of_int", "nodeIndicies")])

        # Attributes.
        #x.add_instance_attribute("nodeList", retval(const_ptr(nodelist), caller_owns_return=True), getter="nodeListPtr", is_const=True)
        x.add_instance_attribute("nodeIndicies", "vector_of_int", getter="nodeIndicies", is_const=True)

        # Methods.
        const_ref_return_value(x, me, "%s::nodeList" % me, nodelist, [], "nodeList")

        # Virtual methods.
        x.add_method("valid", "bool", [], is_const=True, is_virtual=True)
        x.add_method("label", "std::string", [], is_const=True, is_virtual=True)
        x.add_method("dumpState", None, [refparam(fileio, "fileIO"),
                                         refparam("std::string", "pathName")], is_const=True, is_virtual=True)
        x.add_method("restoreState", None, [refparam(fileio, "fileIO"),
                                            refparam("std::string", "pathName")], is_virtual=True)

        # Generate the abstract interface.
        generateBoundaryVirtualBindings(x, ndim, False)

        return

    #---------------------------------------------------------------------------
    # ConstantXVelocityBoundary bindings.
    #---------------------------------------------------------------------------
    def generateConstantXVelocityBoundaryBindings(self, x, ndim):

        # Object names.
        me = "Spheral::BoundarySpace::ConstantXVelocityBoundary%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        plane = "Plane%id" % ndim
        intfield = "Spheral::FieldSpace::IntField%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        vectorfield = "Spheral::FieldSpace::VectorField%id" % ndim
        tensorfield = "Spheral::FieldSpace::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::FieldSpace::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::FieldSpace::VectorDoubleField%id" % ndim
        symtensorfield = "Spheral::FieldSpace::SymTensorField%id" % ndim
        intfieldlist = "Spheral::FieldSpace::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::FieldSpace::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::FieldSpace::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::FieldSpace::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::FieldSpace::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::FieldSpace::ThirdRankTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeSpace::NodeList%id" % ndim
        state = "State%id" % ndim
        derivatives = "StateDerivatives%id" % ndim
        database = "Spheral::DataBaseSpace::DataBase%id" % ndim
        connectivitymap = "ConnectivityMap%id" % ndim

        # Constructors.
        x.add_constructor([refparam(nodelist, "nodeList"), refparam("vector_of_int", "nodeIndicies")])

        return

    #---------------------------------------------------------------------------
    # ConstantYVelocityBoundary bindings.
    #---------------------------------------------------------------------------
    def generateConstantYVelocityBoundaryBindings(self, x, ndim):

        # Object names.
        me = "Spheral::BoundarySpace::ConstantYVelocityBoundary%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        plane = "Plane%id" % ndim
        intfield = "Spheral::FieldSpace::IntField%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        vectorfield = "Spheral::FieldSpace::VectorField%id" % ndim
        tensorfield = "Spheral::FieldSpace::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::FieldSpace::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::FieldSpace::VectorDoubleField%id" % ndim
        symtensorfield = "Spheral::FieldSpace::SymTensorField%id" % ndim
        intfieldlist = "Spheral::FieldSpace::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::FieldSpace::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::FieldSpace::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::FieldSpace::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::FieldSpace::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::FieldSpace::ThirdRankTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeSpace::NodeList%id" % ndim
        state = "State%id" % ndim
        derivatives = "StateDerivatives%id" % ndim
        database = "Spheral::DataBaseSpace::DataBase%id" % ndim
        connectivitymap = "ConnectivityMap%id" % ndim

        # Constructors.
        x.add_constructor([refparam(nodelist, "nodeList"), refparam("vector_of_int", "nodeIndicies")])

        return

    #---------------------------------------------------------------------------
    # ConstantZVelocityBoundary bindings.
    #---------------------------------------------------------------------------
    def generateConstantZVelocityBoundaryBindings(self, x, ndim):

        # Object names.
        me = "Spheral::BoundarySpace::ConstantZVelocityBoundary%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        plane = "Plane%id" % ndim
        intfield = "Spheral::FieldSpace::IntField%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        vectorfield = "Spheral::FieldSpace::VectorField%id" % ndim
        tensorfield = "Spheral::FieldSpace::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::FieldSpace::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::FieldSpace::VectorDoubleField%id" % ndim
        symtensorfield = "Spheral::FieldSpace::SymTensorField%id" % ndim
        intfieldlist = "Spheral::FieldSpace::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::FieldSpace::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::FieldSpace::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::FieldSpace::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::FieldSpace::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::FieldSpace::ThirdRankTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeSpace::NodeList%id" % ndim
        state = "State%id" % ndim
        derivatives = "StateDerivatives%id" % ndim
        database = "Spheral::DataBaseSpace::DataBase%id" % ndim
        connectivitymap = "ConnectivityMap%id" % ndim

        # Constructors.
        x.add_constructor([refparam(nodelist, "nodeList"), refparam("vector_of_int", "nodeIndicies")])

        return

    #---------------------------------------------------------------------------
    # ConstantBoundary bindings.
    #---------------------------------------------------------------------------
    def generateConstantBoundaryBindings(self, x, ndim):

        # Object names.
        me = "Spheral::BoundarySpace::ConstantBoundary%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        plane = "Plane%id" % ndim
        intfield = "Spheral::FieldSpace::IntField%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        vectorfield = "Spheral::FieldSpace::VectorField%id" % ndim
        tensorfield = "Spheral::FieldSpace::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::FieldSpace::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::FieldSpace::VectorDoubleField%id" % ndim
        symtensorfield = "Spheral::FieldSpace::SymTensorField%id" % ndim
        intfieldlist = "Spheral::FieldSpace::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::FieldSpace::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::FieldSpace::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::FieldSpace::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::FieldSpace::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::FieldSpace::ThirdRankTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeSpace::NodeList%id" % ndim
        state = "State%id" % ndim
        derivatives = "StateDerivatives%id" % ndim
        database = "Spheral::DataBaseSpace::DataBase%id" % ndim
        connectivitymap = "ConnectivityMap%id" % ndim

        # Constructors.
        x.add_constructor([refparam(nodelist, "nodeList"), refparam("vector_of_int", "nodeIndicies")])

        # Attributes.
        x.add_instance_attribute("numConstantNodes", "int", getter="numConstantNodes", is_const=True)

        # Methods.
#         x.add_function_as_method("getNodeListPtr", retval(ptr(nodelist), reference_existing_object=True), [param(me, "self")],
#                                  template_parameters = [dim, me],
#                                  custom_name = "nodeList")
        
        # Virtual methods.
        x.add_method("valid", "bool", [], is_const=True, is_virtual=True)

        # Generate the abstract interface.
        generateBoundaryVirtualBindings(x, ndim, False)

        return

    #---------------------------------------------------------------------------
    # SphericalBoundary bindings.
    #---------------------------------------------------------------------------
    def generateSphericalBoundaryBindings(self, x):

        # Object names.
        ndim = 3
        me = "SphericalBoundary"
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        plane = "Plane%id" % ndim
        intfield = "Spheral::FieldSpace::IntField%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        vectorfield = "Spheral::FieldSpace::VectorField%id" % ndim
        tensorfield = "Spheral::FieldSpace::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::FieldSpace::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::FieldSpace::VectorDoubleField%id" % ndim
        symtensorfield = "Spheral::FieldSpace::SymTensorField%id" % ndim
        intfieldlist = "Spheral::FieldSpace::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::FieldSpace::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::FieldSpace::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::FieldSpace::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::FieldSpace::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::FieldSpace::ThirdRankTensorFieldList%id" % ndim
        nodelist = "NodeList%id" % ndim
        state = "State%id" % ndim
        derivatives = "StateDerivatives%id" % ndim
        database = "Spheral::DataBaseSpace::DataBase%id" % ndim
        connectivitymap = "Spheral::NeighborSpace::ConnectivityMap%id" % ndim
        fileio = "Spheral::FileIOSpace::FileIO"

        # Constructors.
        x.add_constructor([refparam(database, "dataBase")])

        # Methods.
        x.add_method("reflectOperator", tensor, [refparam(vector, "r0"), refparam(vector, "r1")], is_const=True)
        
        # Virtual methods.
        x.add_method("label", "std::string", [], is_const=True, is_virtual=True)
        x.add_method("dumpState", None, [refparam(fileio, "fileIO"),
                                         refparam("std::string", "pathName")], is_const=True, is_virtual=True)
        x.add_method("restoreState", None, [refparam(fileio, "fileIO"),
                                            refparam("std::string", "pathName")], is_virtual=True)

        # Generate the abstract interface.
        generateBoundaryVirtualBindings(x, ndim, False)

        return

    #---------------------------------------------------------------------------
    # CylindricalBoundary bindings.
    #---------------------------------------------------------------------------
    def generateCylindricalBoundaryBindings(self, x):

        # Object names.
        ndim = 3
        me = "CylindricalBoundary"
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        plane = "Plane%id" % ndim
        intfield = "Spheral::FieldSpace::IntField%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        vectorfield = "Spheral::FieldSpace::VectorField%id" % ndim
        tensorfield = "Spheral::FieldSpace::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::FieldSpace::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::FieldSpace::VectorDoubleField%id" % ndim
        symtensorfield = "Spheral::FieldSpace::SymTensorField%id" % ndim
        intfieldlist = "Spheral::FieldSpace::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::FieldSpace::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::FieldSpace::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::FieldSpace::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::FieldSpace::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::FieldSpace::ThirdRankTensorFieldList%id" % ndim
        nodelist = "NodeList%id" % ndim
        state = "State%id" % ndim
        derivatives = "StateDerivatives%id" % ndim
        database = "Spheral::DataBaseSpace::DataBase%id" % ndim
        connectivitymap = "Spheral::NeighborSpace::ConnectivityMap%id" % ndim
        fileio = "Spheral::FileIOSpace::FileIO"

        # Constructors.
        x.add_constructor([refparam(database, "dataBase")])

        # Methods.
        x.add_method("reflectOperator", tensor, [refparam(vector, "r0"), refparam(vector, "r1")], is_static=True)
        x.add_method("angularSpacing", "double", [param("double", "ri"),
                                                  param("double", "hzi"),
                                                  param("double", "nodePerh"),
                                                  param("double", "kernelExtent")], is_static=True)
        
        # Virtual methods.
        x.add_method("label", "std::string", [], is_const=True, is_virtual=True)
        x.add_method("dumpState", None, [refparam(fileio, "fileIO"),
                                         refparam("std::string", "pathName")], is_const=True, is_virtual=True)
        x.add_method("restoreState", None, [refparam(fileio, "fileIO"),
                                            refparam("std::string", "pathName")], is_virtual=True)

        # Generate the abstract interface.
        generateBoundaryVirtualBindings(x, ndim, False)

        return
