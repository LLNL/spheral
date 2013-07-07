from pybindgen import *

import sys
sys.path.append("..")
from PBGutils import *
from ref_return_value import *

sys.path.append("../Physics")
from PhysicsModule import generatePhysicsVirtualBindings

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class SVPH:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod):

        # Includes.
        mod.add_include('"SVPH/SVPHTypes.hh"')
    
        # Namespace.
        Spheral = mod.add_cpp_namespace("Spheral")
        PhysicsSpace = Spheral.add_cpp_namespace("PhysicsSpace")
        self.space = Spheral.add_cpp_namespace("SVPHSpace")
        
        generichydro1d = findObject(PhysicsSpace, "GenericHydro1d")
        generichydro2d = findObject(PhysicsSpace, "GenericHydro2d")
        generichydro3d = findObject(PhysicsSpace, "GenericHydro3d")

        # # Expose types.
        # self.SVPHHydroBase1d = addObject(self.space, "SVPHHydroBase1d", allow_subclassing=True, parent=generichydro1d)
        # self.SVPHHydroBase2d = addObject(self.space, "SVPHHydroBase2d", allow_subclassing=True, parent=generichydro2d)
        # self.SVPHHydroBase3d = addObject(self.space, "SVPHHydroBase3d", allow_subclassing=True, parent=generichydro3d)

        return

    #---------------------------------------------------------------------------
    # Generate bindings.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        # self.generateSVPHHydroBaseBindings(self.SVPHHydroBase1d, 1)
        # self.generateSVPHHydroBaseBindings(self.SVPHHydroBase2d, 2)
        # self.generateSVPHHydroBaseBindings(self.SVPHHydroBase3d, 3)
        
        self.generateDimBindings(1)
        self.generateDimBindings(2)
        self.generateDimBindings(3)

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return ["SVPHSpace"]

    #---------------------------------------------------------------------------
    # Add the types per dimension.
    #---------------------------------------------------------------------------
    def generateDimBindings(self, ndim):

        dim = "Spheral::Dim<%i>" % ndim
        vector = "Spheral::Dim<%i>::Vector" % ndim
        tensor = "Spheral::Dim<%i>::Tensor" % ndim
        symtensor = "Spheral::Dim<%i>::SymTensor" % ndim
        scalarfieldlist = "Spheral::FieldSpace::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::FieldSpace::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::FieldSpace::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::FieldSpace::SymTensorFieldList%id" % ndim
        connectivitymap = "Spheral::NeighborSpace::ConnectivityMap%id" % ndim
        tablekernel = "Spheral::KernelSpace::TableKernel%id" % ndim
        mesh = {1 : "LineMesh", 2 : "PolygonalMesh", 3 : "PolyhedralMesh"}[ndim]

        # SVPH sampling.
        for fl in [scalarfieldlist,
                   vectorfieldlist,
                   tensorfieldlist,
                   symtensorfieldlist]:
            self.space.add_function("sampleFieldListSVPH", 
                                    fl,
                                    [constrefparam(fl, "fieldList"),
                                     constrefparam(vectorfieldlist, "position"),
                                     constrefparam(symtensorfieldlist, "H"),
                                     constrefparam(connectivitymap, "connectivityMap"),
                                     constrefparam(tablekernel, "W"),
                                     constrefparam(mesh, "mesh"),
                                     param("bool", "firstOrderConsistent", default_value="true")],
                                    template_parameters = [dim],
                                    custom_name = "sampleFieldListSVPH%id" % ndim)

        # SVPH gradient.
        for (fl, gfl) in [(scalarfieldlist, vectorfieldlist),
                          (vectorfieldlist, tensorfieldlist)]:
            self.space.add_function("gradientFieldListSVPH", 
                                    gfl,
                                    [constrefparam(fl, "fieldList"),
                                     constrefparam(vectorfieldlist, "position"),
                                     constrefparam(symtensorfieldlist, "H"),
                                     constrefparam(connectivitymap, "connectivityMap"),
                                     constrefparam(tablekernel, "W"),
                                     constrefparam(mesh, "mesh"),
                                     param("bool", "firstOrderConsistent", default_value="true")],
                                    template_parameters = [dim],
                                    custom_name = "gradientFieldListSVPH%id" % ndim)

        return

