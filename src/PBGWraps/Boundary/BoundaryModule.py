from pybindgen import *

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
    intfield = "Spheral::IntField%id" % ndim
    scalarfield = "Spheral::ScalarField%id" % ndim
    vectorfield = "Spheral::VectorField%id" % ndim
    tensorfield = "Spheral::TensorField%id" % ndim
    thirdranktensorfield = "Spheral::ThirdRankTensorField%id" % ndim
    vectordoublefield = "Spheral::VectorDoubleField%id" % ndim
    symtensorfield = "Spheral::SymTensorField%id" % ndim
    intfieldlist = "Spheral::IntFieldList%id" % ndim
    scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
    vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
    tensorfieldlist = "Spheral::TensorFieldList%id" % ndim
    symtensorfieldlist = "Spheral::SymTensorFieldList%id" % ndim
    thirdranktensorfieldlist = "Spheral::ThirdRankTensorFieldList%id" % ndim
    nodelist = "Spheral::NodeList%id" % ndim
    state = "Spheral::State%id" % ndim
    derivatives = "Spheral::StateDerivatives%id" % ndim
    database = "Spheral::DataBase%id" % ndim
    connectivitymap = "Spheral::ConnectivityMap%id" % ndim

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
    def __init__(self, mod, srcdir, topsrcdir, dims):

        self.dims = dims

        # Includes.
        mod.add_include('"%s/BoundaryTypes.hh"' % srcdir)
        mod.add_include('"%s/Boundary/Boundary.hh"' % topsrcdir)
    
        # Namespace.
        self.space = mod.add_cpp_namespace("Spheral")

        # Expose types.
        for ndim in self.dims:
            exec("""
self.Boundary%(ndim)id = addObject(self.space, "Boundary%(ndim)id", allow_subclassing=True)
self.PlanarBoundary%(ndim)id = addObject(self.space, "PlanarBoundary%(ndim)id", parent=self.Boundary%(ndim)id)
self.ReflectingBoundary%(ndim)id = addObject(self.space, "ReflectingBoundary%(ndim)id", parent=self.PlanarBoundary%(ndim)id)
self.RigidBoundary%(ndim)id = addObject(self.space, "RigidBoundary%(ndim)id", parent=self.PlanarBoundary%(ndim)id)
self.PeriodicBoundary%(ndim)id = addObject(self.space, "PeriodicBoundary%(ndim)id", parent=self.PlanarBoundary%(ndim)id)
self.ConstantVelocityBoundary%(ndim)id = addObject(self.space, "ConstantVelocityBoundary%(ndim)id", parent=self.Boundary%(ndim)id)
self.ConstantXVelocityBoundary%(ndim)id = addObject(self.space, "ConstantXVelocityBoundary%(ndim)id", parent=self.ConstantVelocityBoundary%(ndim)id)
self.ConstantRVelocityBoundary%(ndim)id = addObject(self.space, "ConstantRVelocityBoundary%(ndim)id", parent=self.ConstantVelocityBoundary%(ndim)id)
self.ConstantBoundary%(ndim)id = addObject(self.space, "ConstantBoundary%(ndim)id", parent=self.Boundary%(ndim)id)
self.CRKSPHVoidBoundary%(ndim)id = addObject(self.space, "CRKSPHVoidBoundary%(ndim)id", parent=self.Boundary%(ndim)id)
self.vecBound%(ndim)id = addObject(mod, "vector_of_Boundary%(ndim)id", allow_subclassing=True)
""" % {"ndim" : ndim})

        if 2 in self.dims:
            self.ConstantYVelocityBoundary2d = addObject(self.space, "ConstantYVelocityBoundary2d", parent=self.ConstantVelocityBoundary2d)

        if 3 in self.dims:
            self.AxisBoundaryRZ = addObject(self.space, "AxisBoundaryRZ", parent=self.ReflectingBoundary2d)
            self.ConstantYVelocityBoundary3d = addObject(self.space, "ConstantYVelocityBoundary3d", parent=self.ConstantVelocityBoundary3d)
            self.ConstantZVelocityBoundary3d = addObject(self.space, "ConstantZVelocityBoundary3d", parent=self.ConstantVelocityBoundary3d)
            self.SphericalBoundary = addObject(self.space, "SphericalBoundary", parent=self.Boundary3d)
            self.CylindricalBoundary = addObject(self.space, "CylindricalBoundary", parent=self.Boundary3d)

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        for ndim in self.dims:
            exec("""
self.generateBoundaryBindings(self.Boundary%(ndim)id, %(ndim)i)
self.generatePlanarBoundaryBindings(self.PlanarBoundary%(ndim)id, %(ndim)i)
self.generateReflectingBoundaryBindings(self.ReflectingBoundary%(ndim)id, %(ndim)i)
self.generateRigidBoundaryBindings(self.RigidBoundary%(ndim)id, %(ndim)i)
self.generatePeriodicBoundaryBindings(self.PeriodicBoundary%(ndim)id, %(ndim)i)
self.generateConstantVelocityBoundaryBindings(self.ConstantVelocityBoundary%(ndim)id, %(ndim)i)
self.generateConstantXVelocityBoundaryBindings(self.ConstantXVelocityBoundary%(ndim)id, %(ndim)i)
self.generateConstantBoundaryBindings(self.ConstantBoundary%(ndim)id, %(ndim)i)
self.generateCRKSPHVoidBoundaryBindings(self.CRKSPHVoidBoundary%(ndim)id, %(ndim)i)
generateStdVectorBindings(self.vecBound%(ndim)id, "Spheral::Boundary%(ndim)id*", "vector_of_Boundary%(ndim)id")
self.space.add_function("dynamicCastBoundary", 
                        retval("Spheral::PlanarBoundary%(ndim)id*", reference_existing_object=True), # caller_owns_return=False),
                        [param("Spheral::Boundary%(ndim)id*", "boundary", transfer_ownership=False)],
                        template_parameters=["Spheral::Boundary%(ndim)id", "Spheral::PlanarBoundary%(ndim)id"],
                        custom_name = "dynamicCastBoundaryToPlanarBoundary%(ndim)id")
""" % {"ndim": ndim})

        if 2 in self.dims:
            self.generateConstantYVelocityBoundaryBindings(self.ConstantYVelocityBoundary2d, 2)

        if 3 in self.dims:
            self.generateAxisBoundaryRZBindings(self.AxisBoundaryRZ)
            self.generateConstantYVelocityBoundaryBindings(self.ConstantYVelocityBoundary3d, 3)
            self.generateConstantZVelocityBoundaryBindings(self.ConstantZVelocityBoundary3d, 3)
            self.generateSphericalBoundaryBindings(self.SphericalBoundary)
            self.generateCylindricalBoundaryBindings(self.CylindricalBoundary)
        
        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return []

    #---------------------------------------------------------------------------
    # Boundary bindings.
    #---------------------------------------------------------------------------
    def generateBoundaryBindings(self, x, ndim):

        # Object names.
        me = "Boundary%id" % ndim
        dim = "Spheral::Dim<%i> " % ndim
        boundarynodes = "Spheral::%s::BoundaryNodes" % me
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        intfield = "Spheral::IntField%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        tensorfield = "Spheral::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::VectorDoubleField%id" % ndim
        vectorvectorfield = "Spheral::VectorVectorField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::ThirdRankTensorFieldList%id" % ndim
        vector_of_Vector = "vector_of_Vector%id" % ndim
        vector_of_Tensor = "vector_of_Tensor%id" % ndim
        vector_of_SymTensor = "vector_of_SymTensor%id" % ndim
        vector_of_ThirdRankTensor = "vector_of_ThirdRankTensor%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        connectivitymap = "Spheral::ConnectivityMap%id" % ndim
        mesh = "Spheral::" + {1 : "LineMesh", 2 : "PolygonalMesh", 3 : "PolyhedralMesh"}[ndim]

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

        x.add_method("enforceFieldListBoundary", None, [refparam(intfieldlist, "fieldList")], is_const=True)
        x.add_method("enforceFieldListBoundary", None, [refparam(scalarfieldlist, "fieldList")], is_const=True)
        x.add_method("enforceFieldListBoundary", None, [refparam(vectorfieldlist, "fieldList")], is_const=True)
        x.add_method("enforceFieldListBoundary", None, [refparam(tensorfieldlist, "fieldList")], is_const=True)
        x.add_method("enforceFieldListBoundary", None, [refparam(symtensorfieldlist, "fieldList")], is_const=True)
        x.add_method("enforceFieldListBoundary", None, [refparam(thirdranktensorfieldlist, "fieldList")], is_const=True)

        x.add_method("initializeProblemStartup", None, [], is_virtual=True)
        x.add_method("finalizeGhostBoundary", None, [], is_const=True, is_virtual=True)
        x.add_method("clip", None, [refparam(vector, "xmin"), refparam(vector, "xmax")], is_const=True, is_virtual=True)

        x.add_method("meshGhostNodes", "bool", [], is_virtual=True, is_const=True)

        x.add_method("reset", None, [refparam(database, "dataBase")], is_virtual=True)

        x.add_method("cullGhostNodes", None, [constrefparam(intfieldlist, "flagSet"),
                                              refparam(intfieldlist, "old2newIndexMap"),
                                              refparam("vector_of_int", "numNodesRemoved")], is_virtual=True)

        x.add_method("enforceBoundary", None, [refparam("vector_of_int", "faceField"), constrefparam(mesh, "mesh")], is_const=True, is_virtual=True)
        x.add_method("enforceBoundary", None, [refparam("vector_of_double", "faceField"), constrefparam(mesh, "mesh")], is_const=True, is_virtual=True)
        x.add_method("enforceBoundary", None, [refparam(vector_of_Vector, "faceField"), constrefparam(mesh, "mesh")], is_const=True, is_virtual=True)
        x.add_method("enforceBoundary", None, [refparam(vector_of_Tensor, "faceField"), constrefparam(mesh, "mesh")], is_const=True, is_virtual=True)
        x.add_method("enforceBoundary", None, [refparam(vector_of_SymTensor, "faceField"), constrefparam(mesh, "mesh")], is_const=True, is_virtual=True)
        x.add_method("enforceBoundary", None, [refparam(vector_of_ThirdRankTensor, "faceField"), constrefparam(mesh, "mesh")], is_const=True, is_virtual=True)

        x.add_method("swapFaceValues", None, [refparam(vectordoublefield, "field"), constrefparam(mesh, "mesh")], is_const=True, is_virtual=True)
        x.add_method("swapFaceValues", None, [refparam(vectorvectorfield, "field"), constrefparam(mesh, "mesh")], is_const=True, is_virtual=True)

        x.add_method("applyGhostBoundary", None, [refparam(vectordoublefield, "field")], is_const=True, is_virtual=True)
        x.add_method("applyGhostBoundary", None, [refparam(vectorvectorfield, "field")], is_const=True, is_virtual=True)

        # Methods.
        x.add_method("haveNodeList", "bool", [constrefparam(nodelist, "nodeList")], is_const=True)
        x.add_method("controlNodes", "vector_of_int", [refparam(nodelist, "nodeList")], is_const=True)
        x.add_method("ghostNodes", "vector_of_int", [refparam(nodelist, "nodeList")], is_const=True)
        x.add_method("violationNodes", "vector_of_int", [refparam(nodelist, "nodeList")], is_const=True)

        x.add_method("applyFieldListGhostBoundary", None, [refparam(scalarfieldlist, "fieldList")], is_const=True)
        x.add_method("applyFieldListGhostBoundary", None, [refparam(vectorfieldlist, "fieldList")], is_const=True)
        x.add_method("applyFieldListGhostBoundary", None, [refparam(tensorfieldlist, "fieldList")], is_const=True)
        x.add_method("applyFieldListGhostBoundary", None, [refparam(symtensorfieldlist, "fieldList")], is_const=True)

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
        intfield = "Spheral::IntField%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        tensorfield = "Spheral::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::VectorDoubleField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::ThirdRankTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        state = "State%id" % ndim
        derivatives = "StateDerivatives%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        connectivitymap = "Spheral::ConnectivityMap%id" % ndim
        fileio = "Spheral::FileIO"
        mesh = "Spheral::" + {1 : "LineMesh", 2 : "PolygonalMesh", 3 : "PolyhedralMesh"}[ndim]

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
        x.add_method("facesOnPlane", "vector_of_unsigned", [constrefparam(mesh, "mesh"),
                                                            constrefparam(plane, "plane"),
                                                            param("double", "tol")], is_const=True);
        
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
        intfield = "Spheral::IntField%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        tensorfield = "Spheral::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::VectorDoubleField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::ThirdRankTensorFieldList%id" % ndim
        nodelist = "NodeList%id" % ndim
        state = "State%id" % ndim
        derivatives = "StateDerivatives%id" % ndim
        database = "Spheral::DataBase%id" % ndim
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
        intfield = "Spheral::IntField%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        tensorfield = "Spheral::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::VectorDoubleField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::ThirdRankTensorFieldList%id" % ndim
        nodelist = "NodeList%id" % ndim
        state = "State%id" % ndim
        derivatives = "StateDerivatives%id" % ndim
        database = "Spheral::DataBase%id" % ndim
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
        me = "Spheral::ConstantVelocityBoundary%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        plane = "Plane%id" % ndim
        intfield = "Spheral::IntField%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        tensorfield = "Spheral::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::VectorDoubleField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::ThirdRankTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        state = "State%id" % ndim
        derivatives = "StateDerivatives%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        connectivitymap = "Spheral::ConnectivityMap%id" % ndim
        fileio = "Spheral::FileIO"

        # Constructors.
        x.add_constructor([refparam(nodelist, "nodeList"), 
                           refparam("vector_of_int", "nodeIndices")])

        # Attributes.
        #x.add_instance_attribute("nodeList", retval(const_ptr(nodelist), caller_owns_return=True), getter="nodeListPtr", is_const=True)
        x.add_instance_attribute("nodeIndices", "vector_of_int", getter="nodeIndices", is_const=True)

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
        me = "Spheral::ConstantXVelocityBoundary%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        plane = "Plane%id" % ndim
        intfield = "Spheral::IntField%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        tensorfield = "Spheral::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::VectorDoubleField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::ThirdRankTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        state = "State%id" % ndim
        derivatives = "StateDerivatives%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        connectivitymap = "ConnectivityMap%id" % ndim

        # Constructors.
        x.add_constructor([refparam(nodelist, "nodeList"), refparam("vector_of_int", "nodeIndices")])

        return

    #---------------------------------------------------------------------------
    # ConstantYVelocityBoundary bindings.
    #---------------------------------------------------------------------------
    def generateConstantYVelocityBoundaryBindings(self, x, ndim):

        # Object names.
        me = "Spheral::ConstantYVelocityBoundary%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        plane = "Plane%id" % ndim
        intfield = "Spheral::IntField%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        tensorfield = "Spheral::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::VectorDoubleField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::ThirdRankTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        state = "State%id" % ndim
        derivatives = "StateDerivatives%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        connectivitymap = "ConnectivityMap%id" % ndim

        # Constructors.
        x.add_constructor([refparam(nodelist, "nodeList"), refparam("vector_of_int", "nodeIndices")])

        return

    #---------------------------------------------------------------------------
    # ConstantZVelocityBoundary bindings.
    #---------------------------------------------------------------------------
    def generateConstantZVelocityBoundaryBindings(self, x, ndim):

        # Object names.
        me = "Spheral::ConstantZVelocityBoundary%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        plane = "Plane%id" % ndim
        intfield = "Spheral::IntField%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        tensorfield = "Spheral::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::VectorDoubleField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::ThirdRankTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        state = "State%id" % ndim
        derivatives = "StateDerivatives%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        connectivitymap = "ConnectivityMap%id" % ndim

        # Constructors.
        x.add_constructor([refparam(nodelist, "nodeList"), refparam("vector_of_int", "nodeIndices")])

        return

    #---------------------------------------------------------------------------
    # ConstantBoundary bindings.
    #---------------------------------------------------------------------------
    def generateConstantBoundaryBindings(self, x, ndim):

        # Object names.
        me = "Spheral::ConstantBoundary%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        plane = "Plane%id" % ndim
        intfield = "Spheral::IntField%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        tensorfield = "Spheral::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::VectorDoubleField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::ThirdRankTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        state = "State%id" % ndim
        derivatives = "StateDerivatives%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        connectivitymap = "ConnectivityMap%id" % ndim
        fileio = "Spheral::FileIO"

        # Constructors.
        x.add_constructor([refparam(nodelist, "nodeList"), 
                           refparam("vector_of_int", "nodeIndices"),
                           refparam(plane, "denialPlane")])

        # Attributes.
        x.add_instance_attribute("numConstantNodes", "int", getter="numConstantNodes", is_const=True)
        x.add_instance_attribute("reflectOperator", tensor, getter="reflectOperator", is_const=True)

        # Methods.
        x.add_method("nodeIndices", "vector_of_int", [], is_const=True)
        
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
    # ConstantBoundary bindings.
    #---------------------------------------------------------------------------
    def generateCRKSPHVoidBoundaryBindings(self, x, ndim):

        # Object names.
        me = "Spheral::CRKSPHVoidBoundary%id" % ndim
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        plane = "Plane%id" % ndim
        intfield = "Spheral::IntField%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        tensorfield = "Spheral::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::VectorDoubleField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::ThirdRankTensorFieldList%id" % ndim
        vectorvectorfieldlist = "Spheral::VectorVectorFieldList%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        state = "State%id" % ndim
        derivatives = "StateDerivatives%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        connectivitymap = "ConnectivityMap%id" % ndim
        fileio = "Spheral::FileIO"

        # Constructors.
        x.add_constructor([constrefparam(intfieldlist, "surfacePoint"), 
                           constrefparam(vectorvectorfieldlist, "etaVoidPoints")])

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
        intfield = "Spheral::IntField%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        tensorfield = "Spheral::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::VectorDoubleField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::ThirdRankTensorFieldList%id" % ndim
        nodelist = "NodeList%id" % ndim
        state = "State%id" % ndim
        derivatives = "StateDerivatives%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        connectivitymap = "Spheral::ConnectivityMap%id" % ndim
        fileio = "Spheral::FileIO"

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
        intfield = "Spheral::IntField%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        tensorfield = "Spheral::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::VectorDoubleField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::ThirdRankTensorFieldList%id" % ndim
        nodelist = "NodeList%id" % ndim
        state = "State%id" % ndim
        derivatives = "StateDerivatives%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        connectivitymap = "Spheral::ConnectivityMap%id" % ndim
        fileio = "Spheral::FileIO"

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

    #---------------------------------------------------------------------------
    # AxisBoundaryRZ bindings.
    #---------------------------------------------------------------------------
    def generateAxisBoundaryRZBindings(self, x):

        # Object names.
        ndim = 2
        me = "AxisBoundaryRZ"
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        plane = "Plane%id" % ndim
        intfield = "Spheral::IntField%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfield = "Spheral::VectorField%id" % ndim
        tensorfield = "Spheral::TensorField%id" % ndim
        thirdranktensorfield = "Spheral::ThirdRankTensorField%id" % ndim
        vectordoublefield = "Spheral::VectorDoubleField%id" % ndim
        symtensorfield = "Spheral::SymTensorField%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::SymTensorFieldList%id" % ndim
        thirdranktensorfieldlist = "Spheral::ThirdRankTensorFieldList%id" % ndim
        nodelist = "Spheral::NodeList%id" % ndim
        state = "State%id" % ndim
        derivatives = "StateDerivatives%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        connectivitymap = "Spheral::ConnectivityMap%id" % ndim
        fileio = "Spheral::FileIO"
        mesh = "Spheral::" + {1 : "LineMesh", 2 : "PolygonalMesh", 3 : "PolyhedralMesh"}[ndim]

        # Constructors.
        x.add_constructor([param("const double", "etamin", default_value="0.1")])

        # Attributes.
        x.add_instance_attribute("etamin", "double", getter="etamin", setter="etamin")

        # Override the abstract interface.
        x.add_method("setViolationNodes", None, [refparam(nodelist, "nodeList")], is_virtual=True)
        x.add_method("updateViolationNodes", None, [refparam(nodelist, "nodeList")], is_virtual=True)

        return

