from pybindgen import *

import sys
sys.path.append("..")
from PBGutils import *

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class NodeGenerators:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod):
        mod.add_include('"NodeGenerators/generateCylDistributionFromRZ.hh"')
        mod.add_include('"NodeGenerators/fillFacetedVolume.hh"')
        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):
        Spheral = mod.add_cpp_namespace("Spheral")
        Spheral.add_function("generateCylDistributionFromRZ",
                             None,
                             [refparam("vector_of_double", "x"),
                              refparam("vector_of_double", "y"),
                              refparam("vector_of_double", "z"),
                              refparam("vector_of_double", "m"),
                              refparam("vector_of_SymTensor3d", "H"),
                              refparam("vector_of_int", "globalIDs"),
                              refparam("vector_of_vector_of_double", "extraFields"),
                              param("double", "nNoderPerh"),
                              param("double", "kernelExtent"),
                              param("double", "phi"),
                              param("int", "procID"),
                              param("int", "nProcs")],
                             docstring = "Spin the given given 2D information to a 3D distribution.")

        Spheral.add_function("fillFacetedVolume",
                             "vector_of_Vector3d",
                             [constrefparam("Polyhedron", "outerBoundary"),
                              param("unsigned int", "n1d")],
                             docstring = "Return a vector of positions filling the given polyhedron.")

        Spheral.add_function("fillFacetedVolume",
                             "vector_of_Vector3d",
                             [constrefparam("Polyhedron", "innerBoundary"),
                              constrefparam("Polyhedron", "outerBoundary"),
                              param("unsigned int", "n1d")],
                             docstring = "Return a vector of positions filling the volume between an inner and outer bounding polyhedra.")
        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return []

