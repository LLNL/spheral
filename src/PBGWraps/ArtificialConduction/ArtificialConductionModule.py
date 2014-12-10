from pybindgen import *

import sys
sys.path.append("..")
from PBGutils import *
from ref_return_value import *

sys.path.append("../CXXTypes")
from CXXTypesModule import generateStdVectorBindings

sys.path.append("../Physics")
from PhysicsModule import generatePhysicsVirtualBindings

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class ArtificialConduction:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod):

        # Includes.
        mod.add_include('"ArtificialConduction/ArtificialConductionTypes.hh"')
    
        # Namespace.
        Spheral = mod.add_cpp_namespace("Spheral")
        space = Spheral.add_cpp_namespace("PhysicsSpace")

        Physics1d = findObject(space, "Physics1d")
        Physics2d = findObject(space, "Physics2d")
        Physics3d = findObject(space, "Physics3d")

        # Expose types.
        self.ArtificialConduction1d = addObject(space, "ArtificialConduction1d", allow_subclassing=True, parent=Physics1d)
        self.ArtificialConduction2d = addObject(space, "ArtificialConduction2d", allow_subclassing=True, parent=Physics2d)
        self.ArtificialConduction3d = addObject(space, "ArtificialConduction3d", allow_subclassing=True, parent=Physics3d)

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):
        self.generateArtificialConductionBindings(self.ArtificialConduction1d, 1)
        self.generateArtificialConductionBindings(self.ArtificialConduction2d, 2)
        self.generateArtificialConductionBindings(self.ArtificialConduction3d, 3)
        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return ["PhysicsSpace"]

    #---------------------------------------------------------------------------
    # Artificial Conduction
    #---------------------------------------------------------------------------
    def generateArtificialConductionBindings(self, x, ndim):
        
        # Object names.
        me = "Spheral::PhysicsSpace::ArtificialConduction%id" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        scalarfield = "Spheral::FieldSpace::ScalarField%id" % ndim
        vectorfieldlist = "Spheral::FieldSpace::VectorFieldList%id" % ndim
        database = "Spheral::DataBaseSpace::DataBase%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        boundary = "Spheral::BoundarySpace::Boundary%id" % ndim
        vector_of_boundary = "vector_of_Boundary%id" % ndim
        tablekernel = "Spheral::KernelSpace::TableKernel%id" % ndim
        artificialviscosity = "Spheral::ArtificialViscositySpace::ArtificialViscosity%id" % ndim
        
        # Constructor.
        x.add_constructor([constrefparam(tablekernel, "W"),
                           param("double", "arCondAlpha", default_value="0.5")])
                           
        # Methods.
        x.add_method("initializeProblemStartup", None, [refparam(database, "dataBase")], is_virtual=True)
        generatePhysicsVirtualBindings(x,ndim,pureVirtual=False)
    
                           
        return
