from pybindgen import *

from PBGutils import *

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class NodeGenerators:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir, dims):

        self.dims = dims

        mod.add_include('"%s/NodeGeneratorsTypes.hh"' % srcdir)
        Spheral = mod.add_cpp_namespace("Spheral")

        # Expose types.
        if 2 in self.dims:
            self.WeightingFunctor2d = addObject(Spheral, "WeightingFunctor2d", allow_subclassing=True)
        if 3 in self.dims:
            self.WeightingFunctor3d = addObject(Spheral, "WeightingFunctor3d", allow_subclassing=True)

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        Spheral = mod.add_cpp_namespace("Spheral")

        if 2 in self.dims:
            self.addWeightingFunctorMethods(self.WeightingFunctor2d, 2)
        if 3 in self.dims:
            self.addWeightingFunctorMethods(self.WeightingFunctor3d, 3)

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
                              param("unsigned int", "n1d"),
                              param("unsigned int", "domain"),
                              param("unsigned int", "numDomains")],
                             docstring = "Return a vector of positions filling the given polyhedron.")

        Spheral.add_function("fillFacetedVolume",
                             "vector_of_Vector3d",
                             [constrefparam("Polyhedron", "innerBoundary"),
                              constrefparam("Polyhedron", "outerBoundary"),
                              param("unsigned int", "n1d"),
                              param("unsigned int", "domain"),
                              param("unsigned int", "numDomains")],
                             docstring = "Return a vector of positions filling the volume between an inner and outer bounding polyhedra.")

        Spheral.add_function("readSiloPolyMesh",
                             None,
                             [param("std::string", "fileName"),
                              param("std::string", "meshName"),
                              refparam("vector_of_Vector3d", "positions"),
                              refparam("vector_of_double", "volumes"),
                              refparam("vector_of_SymTensor3d", "H")],
                             docstring = "Compute stuff useful for a NodeGenerator from a polyhedral mesh in a silo file.")

        subdims = []
        if 2 in self.dims:
            subdims.append(2)
        if 3 in self.dims:
            subdims.append(3)
        for ndim in subdims:
            poly = "Spheral::" + {2 : "Polygon", 3 : "Polyhedron"}[ndim]
            vector = "Vector%id" % ndim
            database = "Spheral::DataBaseSpace::DataBase%id" % ndim
            boundary = "Spheral::BoundarySpace::Boundary%id" % ndim
            vector_of_boundary = "vector_of_Boundary%id" % ndim
            tablekernel = "Spheral::KernelSpace::TableKernel%id" % ndim
            smoothingscalebase = "Spheral::NodeSpace::SmoothingScaleBase%id" % ndim
            weightingfunctor = "Spheral::WeightingFunctor%id" % ndim

            Spheral.add_function("relaxNodeDistribution", None,
                                 [constrefparam(database, "dataBase"),
                                  constrefparam(poly, "boundary"),
                                  constrefparam(vector_of_boundary, "boundaries"),
                                  constrefparam(tablekernel, "W"),
                                  constrefparam(smoothingscalebase, "smoothingScaleMethod"),
                                  constrefparam(weightingfunctor, "weightingFunctor"),# default_value=weightingfunctor+"()"),
                                  constrefparam(weightingfunctor, "massDensityFunctor"),# default_value=weightingfunctor+"()"),
                                  param("double", "targetMass", default_value="0.0"),
                                  param("int", "maxIterations", default_value="100"),
                                  param("double", "tolerance", default_value="1.0e-4")],
                                 docstring = "Iteratively relax a set of nodes within a boundary.")

        return

    #---------------------------------------------------------------------------
    # Add methods to WeightingFunctor.
    #---------------------------------------------------------------------------
    def addWeightingFunctorMethods(self, x, ndim):

        poly = "Spheral::" + {2 : "Polygon", 3 : "Polyhedron"}[ndim]
        #poly = "Spheral::Dim<%i>::FacetedVolume" % ndim
        vector = "Vector%id" % ndim

        x.add_constructor([])
        x.add_method("__call__", "double", [constrefparam(vector, "position"),
                                            constrefparam(poly, "boundary")],
                     is_const=True, is_virtual=True)

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return []

