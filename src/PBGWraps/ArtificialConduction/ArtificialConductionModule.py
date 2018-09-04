from pybindgen import *

from PBGutils import *
from ref_return_value import *

from CXXTypesModule import generateStdVectorBindings
from PhysicsModule import generatePhysicsVirtualBindings

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class ArtificialConduction:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir, dims):

        self.dims = dims

        # Includes.
        mod.add_include('"%s/ArtificialConductionTypes.hh"' % srcdir)
    
        # Namespace.
        space = mod.add_cpp_namespace("Spheral")

        # Expose types.
        for dim in self.dims:
            exec('''
Physics%(dim)id = findObject(space, "Physics%(dim)id")
self.ArtificialConduction%(dim)id = addObject(space, "ArtificialConduction%(dim)id", allow_subclassing=True, parent=Physics%(dim)id)
''' % {"dim" : dim})

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):
        for dim in self.dims:
            exec('''
self.generateArtificialConductionBindings(self.ArtificialConduction%(dim)id, %(dim)i)
''' % {"dim" : dim})
        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return []

    #---------------------------------------------------------------------------
    # Artificial Conduction
    #---------------------------------------------------------------------------
    def generateArtificialConductionBindings(self, x, ndim):
        
        # Object names.
        me = "Spheral::ArtificialConduction%id" % ndim
        vector = "Vector%id" % ndim
        tensor = "Tensor%id" % ndim
        symtensor = "SymTensor%id" % ndim
        scalarfield = "Spheral::ScalarField%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        database = "Spheral::DataBase%id" % ndim
        state = "Spheral::State%id" % ndim
        derivatives = "Spheral::StateDerivatives%id" % ndim
        boundary = "Spheral::Boundary%id" % ndim
        vector_of_boundary = "vector_of_Boundary%id" % ndim
        tablekernel = "Spheral::TableKernel%id" % ndim
        artificialviscosity = "Spheral::ArtificialViscosity%id" % ndim
        
        # Constructor.
        x.add_constructor([constrefparam(tablekernel, "W"),
                           param("double", "arCondAlpha", default_value="0.5"),
                           param("CRKOrder", "ACcorrectionOrder", default_value="Spheral::CRKOrder::LinearOrder")])
        # Attributes.
        x.add_instance_attribute("ACcorrectionOrder", "CRKOrder", getter="ACcorrectionOrder", setter="ACcorrectionOrder")
                           
        # Methods.
        x.add_method("initializeProblemStartup", None, [refparam(database, "dataBase")], is_virtual=True)
        generatePhysicsVirtualBindings(x,ndim,pureVirtual=False)
    
                           
        return
