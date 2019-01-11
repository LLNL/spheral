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

        # centroidalRelaxNodesImpl.
        for ndim in self.dims:
            Spheral.add_function("centroidalRelaxNodesImpl", "unsigned int",
                                 [refparam("DataBase%id" % ndim, "db"),
                                  constrefparam("vector_of_FacetedVolume%id" % ndim, "volumeBoundaries"),
                                  constrefparam("vector_of_vector_of_FacetedVolume%id" % ndim, "holes"),
                                  constrefparam("TableKernel%id" % ndim, "W"),
                                  constrefparam("Spheral::PythonBoundFunctors::SpheralFunctor<Vector%id, double>" % ndim, "rhofunc"),
                                  constrefparam("Spheral::PythonBoundFunctors::SpheralFunctor<Vector%id, Vector%id>" % (ndim, ndim), "gradrhofunc"),
                                  param("bool", "rhoConst"),
                                  param("bool", "useGradRhoFunc"),
                                  refparam("vector_of_Boundary%id" % ndim, "boundaries"),
                                  param("unsigned int", "maxIterations"),
                                  param("double", "fracTol"),
                                  param("CRKOrder", "correctionOrder"),
                                  param("double", "centroidFrac"),
                                  refparam("ScalarFieldList%id" % ndim, "vol"),
                                  refparam("IntFieldList%id" % ndim, "surfacePoint"),
                                  refparam("FacetedVolumeFieldList%id" % ndim, "cells")],
                                 docstring = "Use Lloyds algorithm to relax point positions.")

        subdims = []
        if 2 in self.dims:
            subdims.append(2)
        if 3 in self.dims:
            subdims.append(3)
        for ndim in subdims:
            dim = "Spheral::Dim<%i>" % ndim
            poly = "Spheral::" + {2 : "Polygon", 3 : "Polyhedron"}[ndim]
            vector_of_facetedvolume = "vector_of_FacetedVolume%id" % ndim
            vector = "Vector%id" % ndim
            database = "Spheral::DataBase%id" % ndim
            boundary = "Spheral::Boundary%id" % ndim
            vector_of_boundary = "vector_of_Boundary%id" % ndim
            vector_of_Vector = "vector_of_Vector%id" % ndim
            tablekernel = "Spheral::TableKernel%id" % ndim
            smoothingscalebase = "Spheral::SmoothingScaleBase%id" % ndim
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

            Spheral.add_function("compactFacetedVolumes", "unsigned int",
                                 [refparam(vector_of_facetedvolume, "shapes"),
                                  refparam(vector_of_Vector, "centers"),
                                  refparam("vector_of_int", "flags"),
                                  constrefparam(poly, "surface"),
                                  param("double", "depthmax"),
                                  param("unsigned int", "surfaceIterations"),
                                  param("unsigned int", "maxIterations"),
                                  param("double", "dispfrac"),
                                  param("double", "maxoverlapfrac")],
                                 template_parameters = [dim],
                                 custom_name = "compactFacetedVolumes%id" % ndim,
                                 docstring = "Iteratively compact a set of shapes into surface polygon/polyhedron.")

            Spheral.add_function("chooseRandomNonoverlappingCenter", "unsigned int",
                                 [refparam(vector, "result"),
                                  constrefparam(poly, "trialShape"),
                                  constrefparam(poly, "boundary"),
                                  constrefparam(vector_of_facetedvolume, "existingShapes"),
                                  param("unsigned int", "maxTries")],
                                 template_parameters = [dim],
                                 custom_name = "chooseRandomNonoverlappingCenter%id" % ndim,
                                 docstring = "Search for a non-overlapping position for a shape in a set of shapes inside a surface polygon/polyhedron.")

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

