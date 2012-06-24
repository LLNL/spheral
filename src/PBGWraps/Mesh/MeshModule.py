from pybindgen import *

import sys
sys.path.append("..")
from PBGutils import *
from ref_return_value import *

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class Mesh:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod):

        # Includes.
        mod.add_include('"Mesh/MeshTypes.hh"')
    
        # Namespace.
        Spheral = mod.add_cpp_namespace("Spheral")
        self.space = Spheral.add_cpp_namespace("MeshSpace")

        # Expose types.
        for prefix, dim in (("Line", "1d"),
                            ("Polygonal", "2d"),
                            ("Polyhedral", "3d")):
            exec("""
self.%(prefix)sMesh = addObject(self.space, "%(prefix)sMesh", allow_subclassing=True)
self.%(prefix)sMeshNode = addObject(self.space, "Node", outer_class=self.%(prefix)sMesh, allow_subclassing=True)
self.%(prefix)sMeshEdge = addObject(self.space, "Edge", outer_class=self.%(prefix)sMesh, allow_subclassing=True)
self.%(prefix)sMeshFace = addObject(self.space, "Face", outer_class=self.%(prefix)sMesh, allow_subclassing=True)
self.%(prefix)sMeshZone = addObject(self.space, "Zone", outer_class=self.%(prefix)sMesh, allow_subclassing=True)

self.MeshWall%(dim)s = addObject(self.space, "MeshWall%(dim)s", allow_subclassing=True)
self.PlanarMeshWall%(dim)s = addObject(self.space, "PlanarMeshWall%(dim)s", allow_subclassing=True, parent=self.MeshWall%(dim)s)
self.FacetedMeshWall%(dim)s = addObject(self.space, "FacetedMeshWall%(dim)s", allow_subclassing=True, parent=self.MeshWall%(dim)s)
""" % {"prefix" : prefix, "dim" : dim})

        self.VoroPP2d = addObject(self.space, "VoroPP2d", allow_subclassing=True)
        self.VoroPP3d = addObject(self.space, "VoroPP3d", allow_subclassing=True)

        return

    #---------------------------------------------------------------------------
    # Generate bindings.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        # Add methods to objects.
        for prefix, ndim in (("Line", 1),
                             ("Polygonal", 2),
                             ("Polyhedral", 3)):
            exec("""
self.generateMeshBindings(self.%(prefix)sMesh, "%(prefix)sMesh", %(ndim)i)
self.generateNodeBindings(self.%(prefix)sMeshNode, "%(prefix)sMesh", %(ndim)i)
self.generateEdgeBindings(self.%(prefix)sMeshEdge, "%(prefix)sMesh", %(ndim)i)
self.generateFaceBindings(self.%(prefix)sMeshFace, "%(prefix)sMesh", %(ndim)i)
self.generateZoneBindings(self.%(prefix)sMeshZone, "%(prefix)sMesh", %(ndim)i)

self.generateMeshWallBindings(self.MeshWall%(dim)s, %(ndim)i)
self.generatePlanarMeshWallBindings(self.PlanarMeshWall%(dim)s, %(ndim)i)
self.generateFacetedMeshWallBindings(self.FacetedMeshWall%(dim)s, %(ndim)i)

self.addFunctions("%(prefix)sMesh", %(ndim)i)
""" % {"prefix" : prefix, 
       "ndim" : ndim, 
       "dim" : "%id" % ndim})

        self.generateVoroPPBindings(self.VoroPP2d, "VoroPP2d", 2)
        self.generateVoroPPBindings(self.VoroPP3d, "VoroPP3d", 3)

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return ["MeshSpace"]

    #---------------------------------------------------------------------------
    # Mesh Bindings.
    #---------------------------------------------------------------------------
    def generateMeshBindings(self, x, name, ndim):

        # Object names.
        dim = "Spheral::Dim<%i>" % ndim
        me = "Spheral::MeshSpace::%s" % name
        node = "Spheral::MeshSpace::%s::Node" % name
        edge = "Spheral::MeshSpace::%s::Edge" % name
        face = "Spheral::MeshSpace::%s::Face" % name
        zone = "Spheral::MeshSpace::%s::Zone" % name
        vector = "Spheral::Vector%id" % ndim
        vector_of_vector = "vector_of_Vector%id" % ndim
        nodelist = "Spheral::NodeSpace::NodeList%id" % ndim
        vector_of_nodelist = "vector_of_NodeList%id" % ndim
        facetedvolume = ("Box1d", "Polygon", "Polyhedron")[ndim - 1]

        # Constructors.
        x.add_constructor([])
        x.add_constructor([constrefparam(vector_of_vector, "generators"),
                           constrefparam(vector, "xmin"),
                           constrefparam(vector, "xmax")])
        x.add_constructor([constrefparam(vector_of_vector, "generators"),
                           constrefparam(facetedvolume, "boundary")])
        x.add_constructor([constrefparam(vector_of_vector, "nodePositions"),
                           constrefparam("vector_of_vector_of_unsigned", "edgeNodes"),
                           constrefparam("vector_of_vector_of_unsigned", "faceEdges"),
                           constrefparam("vector_of_vector_of_unsigned", "zoneFaces")])

        # Methods.
        x.add_method("clear", None, [])
        x.add_method("reconstruct", None, [constrefparam(vector_of_vector, "generators"),
                                           constrefparam(facetedvolume, "boundary")])
        x.add_method("removeZonesByMask", None, [constrefparam("vector_of_unsigned", "mask")])
        x.add_method("cleanEdges", None, [param("double", "edgeTol")])
        x.add_method("node", node, [param("unsigned int", "i")], is_const=True)
        x.add_method("edge", edge, [param("unsigned int", "i")], is_const=True)
        x.add_method("face", face, [param("unsigned int", "i")], is_const=True)
        x.add_method("zone", zone, [param("unsigned int", "i")], is_const=True)
        x.add_method("zone", zone, [constrefparam(nodelist, "nodeList"), param("unsigned int", "i")], is_const=True)
        x.add_method("zone", zone, [param("unsigned int", "nodeListi"), param("unsigned int", "i")], is_const=True)
        x.add_method("offset", "unsigned int", [constrefparam(nodelist, "nodeList")], is_const=True)
        x.add_method("offset", "unsigned int", [param("unsigned int", "nodeListi")], is_const=True)
        x.add_method("generateDomainInfo", None, [])
        x.add_method("globalMeshNodeIDs", "vector_of_unsigned", [], is_const=True)
        x.add_method("validDomainInfo", "std::string",
                     [constrefparam(vector, "xmin"), constrefparam(vector, "xmax"), param("bool", "checkUniqueSendProc")], is_const=True)
        x.add_method("storeNodeListOffsets", None,
                     [constrefparam(vector_of_nodelist, "nodeLists"), constrefparam("vector_of_unsigned", "offsets")])
        x.add_method("boundingSurface", facetedvolume, [], is_const=True)

        # Attributes.
        x.add_instance_attribute("nDim", "int", getter="nDim", is_const=True)
        x.add_instance_attribute("numNodes", "unsigned int", getter="numNodes", is_const=True)
        x.add_instance_attribute("numEdges", "unsigned int", getter="numEdges", is_const=True)
        x.add_instance_attribute("numFaces", "unsigned int", getter="numFaces", is_const=True)
        x.add_instance_attribute("numZones", "unsigned int", getter="numZones", is_const=True)
        x.add_instance_attribute("minimumScale", "double", getter="minimumScale", is_const=True)
        x.add_instance_attribute("neighborDomains", "vector_of_unsigned", getter="neighborDomains", is_const=True)
        x.add_instance_attribute("sharedNodes", "vector_of_vector_of_unsigned", getter="sharedNodes", is_const=True)

        x.add_static_attribute("UNSETID", "unsigned int",  is_const=True)
        x.add_static_attribute("minFacesPerZone", "unsigned int",  is_const=True)
        x.add_static_attribute("minEdgesPerZone", "unsigned int",  is_const=True)
        x.add_static_attribute("minNodesPerZone", "unsigned int",  is_const=True)
        x.add_static_attribute("minEdgesPerFace", "unsigned int",  is_const=True)
        x.add_static_attribute("minNodesPerFace", "unsigned int",  is_const=True)

        return

    #---------------------------------------------------------------------------
    # Mesh::Node Bindings.
    #---------------------------------------------------------------------------
    def generateNodeBindings(self, x, meshname, ndim):

        # Object names.
        mesh = "Spheral::MeshSpace::%s" % meshname
        node = "Spheral::MeshSpace::%s::Node" % meshname
        edge = "Spheral::MeshSpace::%s::Edge" % meshname
        face = "Spheral::MeshSpace::%s::Face" % meshname
        zone = "Spheral::MeshSpace::%s::Zone" % meshname
        vector = "Spheral::Vector%id" % ndim
        vector_of_vector = "vector_of_Vector%id" % ndim

        # Constructors.
        x.add_constructor([constrefparam(mesh, "mesh"),
                           param("unsigned int", "ID"),
                           param("vector_of_unsigned", "zoneIDs")])

        # Methods.
        x.add_method("position", vector, [], is_const=True)

        # Attributes.
        x.add_instance_attribute("ID", "unsigned int", getter="ID", is_const=True)
        x.add_instance_attribute("zoneIDs", "vector_of_unsigned", getter="zoneIDs", is_const=True)

        return

    #---------------------------------------------------------------------------
    # Mesh::Edge Bindings.
    #---------------------------------------------------------------------------
    def generateEdgeBindings(self, x, meshname, ndim):

        # Object names.
        mesh = "Spheral::MeshSpace::%s" % meshname
        node = "Spheral::MeshSpace::%s::Node" % meshname
        edge = "Spheral::MeshSpace::%s::Edge" % meshname
        face = "Spheral::MeshSpace::%s::Face" % meshname
        zone = "Spheral::MeshSpace::%s::Zone" % meshname
        vector = "Spheral::Vector%id" % ndim
        vector_of_vector = "vector_of_Vector%id" % ndim

        # Constructors.
        x.add_constructor([constrefparam(mesh, "mesh"),
                           param("unsigned int", "ID"),
                           param("unsigned int", "node1ID"),
                           param("unsigned int", "node2ID")])

        # Methods.
        x.add_method("position", vector, [], is_const=True)
        x.add_method("length", "double", [], is_const=True)

        # Attributes.
        x.add_instance_attribute("ID", "unsigned int", getter="ID", is_const=True)
        x.add_instance_attribute("node1ID", "unsigned int", getter="node1ID", is_const=True)
        x.add_instance_attribute("node2ID", "unsigned int", getter="node2ID", is_const=True)
        x.add_instance_attribute("node1", node, getter="node1", is_const=True)
        x.add_instance_attribute("node2", node, getter="node2", is_const=True)

        return

    #---------------------------------------------------------------------------
    # Mesh::Face Bindings.
    #---------------------------------------------------------------------------
    def generateFaceBindings(self, x, meshname, ndim):

        # Object names.
        mesh = "Spheral::MeshSpace::%s" % meshname
        node = "Spheral::MeshSpace::%s::Node" % meshname
        edge = "Spheral::MeshSpace::%s::Edge" % meshname
        face = "Spheral::MeshSpace::%s::Face" % meshname
        zone = "Spheral::MeshSpace::%s::Zone" % meshname
        vector = "Spheral::Vector%id" % ndim
        vector_of_vector = "vector_of_Vector%id" % ndim

        # Constructors.
        x.add_constructor([constrefparam(mesh, "mesh"),
                           param("unsigned int", "ID"),
                           param("unsigned int", "zone1ID"),
                           param("unsigned int", "zone2ID"),
                           param("vector_of_unsigned", "edgeIDs")])

        # Methods.
        x.add_method("position", vector, [], is_const=True)
        x.add_method("area", "double", [], is_const=True)
        x.add_method("unitNormal", vector, [], is_const=True)
        x.add_method("oppositeZoneID", "unsigned int", [param("unsigned int", "zoneID")], is_const=True)
        x.add_method("compare", "int", [constrefparam(vector, "point"), param("double", "tol", default_value="1.0e-8")], is_const=True)

        # Attributes.
        x.add_instance_attribute("ID", "unsigned int", getter="ID", is_const=True)
        x.add_instance_attribute("numNodes", "unsigned int", getter="numNodes", is_const=True)
        x.add_instance_attribute("numEdges", "unsigned int", getter="numEdges", is_const=True)
        x.add_instance_attribute("zone1ID", "unsigned int", getter="zone1ID", is_const=True)
        x.add_instance_attribute("zone2ID", "unsigned int", getter="zone2ID", is_const=True)
        x.add_instance_attribute("nodeIDs", "vector_of_unsigned", getter="nodeIDs", is_const=True)
        x.add_instance_attribute("edgeIDs", "vector_of_unsigned", getter="edgeIDs", is_const=True)

        return

    #---------------------------------------------------------------------------
    # Mesh::Zone Bindings.
    #---------------------------------------------------------------------------
    def generateZoneBindings(self, x, meshname, ndim):

        # Object names.
        mesh = "Spheral::MeshSpace::%s" % meshname
        node = "Spheral::MeshSpace::%s::Node" % meshname
        edge = "Spheral::MeshSpace::%s::Edge" % meshname
        face = "Spheral::MeshSpace::%s::Face" % meshname
        zone = "Spheral::MeshSpace::%s::Zone" % meshname
        vector = "Spheral::Vector%id" % ndim
        vector_of_vector = "vector_of_Vector%id" % ndim
        convexhull = ("Box1d", "Polygon", "Polyhedron")[ndim - 1]

        # Constructors.
        x.add_constructor([constrefparam(mesh, "mesh"),
                           param("unsigned int", "ID"),
                           param("vector_of_unsigned", "faceIDs")])

        # Methods.
        x.add_method("position", vector, [], is_const=True)
        x.add_method("volume", "double", [], is_const=True)
        x.add_method("convexHull", convexhull, [], is_const=True)

        # Attributes.
        x.add_instance_attribute("ID", "unsigned int", getter="ID", is_const=True)
        x.add_instance_attribute("numNodes", "unsigned int", getter="numNodes", is_const=True)
        x.add_instance_attribute("numEdges", "unsigned int", getter="numEdges", is_const=True)
        x.add_instance_attribute("numFaces", "unsigned int", getter="numFaces", is_const=True)
        x.add_instance_attribute("nodeIDs", "vector_of_unsigned", getter="nodeIDs", is_const=True)
        x.add_instance_attribute("edgeIDs", "vector_of_unsigned", getter="edgeIDs", is_const=True)
        x.add_instance_attribute("faceIDs", "vector_of_unsigned", getter="faceIDs", is_const=True)

        return

    #---------------------------------------------------------------------------
    # MeshWall Bindings.
    #---------------------------------------------------------------------------
    def generateMeshWallBindings(self, x, ndim):

        # Object names.
        dim = "Spheral::Dim<%i>" % ndim
        me = "Spheral::MeshSpace::MeshWall%s" % dim

        # Constructors.
        if ndim == 1:
            x.add_constructor([param("double", "xmin", default_value="-std::numeric_limits<double>::max()"),
                               param("double", "xmax", default_value="std::numeric_limits<double>::max()")])
        else:
            x.add_constructor([])

        return

    #---------------------------------------------------------------------------
    # PlanarMeshWall Bindings.
    #---------------------------------------------------------------------------
    def generatePlanarMeshWallBindings(self, x, ndim):

        # Object names.
        dim = "Spheral::Dim<%i>" % ndim
        me = "Spheral::MeshSpace::PlanarMeshWall%s" % dim
        plane = "Plane%id" % ndim

        # Constructors.
        x.add_constructor([constrefparam(plane, "plane")])
        x.add_constructor([constrefparam(plane, "plane1"),
                           constrefparam(plane, "plane2")])

        return

    #---------------------------------------------------------------------------
    # FacetedMeshWall Bindings.
    #---------------------------------------------------------------------------
    def generateFacetedMeshWallBindings(self, x, ndim):

        # Object names.
        dim = "Spheral::Dim<%i>" % ndim
        me = "Spheral::MeshSpace::FacetedMeshWall%s" % dim
        volume = ("Box1d", "Polygon", "Polyhedron")[ndim - 1]

        # Constructors.
        x.add_constructor([constrefparam(volume, "volume")])

        return

##     #---------------------------------------------------------------------------
##     # VoroPP2d Bindings.
##     #---------------------------------------------------------------------------
##     def generateVoroPP2dBindings(self, x, name, ndim):

##         # Object names.
##         dim = "Spheral::Dim<%i>" % ndim
##         me = "Spheral::MeshSpace::%s" % name
##         vector = "Spheral::Vector%id" % ndim
##         vector_of_vector = "vector_of_Vector%id" % ndim
##         vector_of_vector_of_vector = "vector_of_vector_of_Vector%id" % ndim
##         meshwall = "Spheral::MeshSpace::MeshWall%id" % ndim

##         # Constructors.
##         x.add_constructor([constrefparam(vector_of_vector, "generators"),
##                            constrefparam(vector, "xmin"),
##                            constrefparam(vector, "xmax"),
##                            param("unsigned int", "nx", default_value="20"),
##                            param("unsigned int", "ny", default_value="20")])

##         # Methods.
##         x.add_method("allCells", None, [refparam(vector_of_vector, "vertices"),
##                                         refparam("vector_of_vector_of_unsigned", "cellVertexIndices")], is_const=True)
##         x.add_method("subRegion", None, [constrefparam(vector, "point"),
##                                          param("unsigned int", "i"),
##                                          param("unsigned int", "j")], is_const=True)

##         return

    #---------------------------------------------------------------------------
    # VoroPP Bindings.
    #---------------------------------------------------------------------------
    def generateVoroPPBindings(self, x, name, ndim):

        # Object names.
        dim = "Spheral::Dim<%i>" % ndim
        me = "Spheral::MeshSpace::%s" % name
        vector = "Spheral::Vector%id" % ndim
        vector_of_vector = "vector_of_Vector%id" % ndim
        vector_of_vector_of_vector = "vector_of_vector_of_Vector%id" % ndim
        meshwall = "Spheral::MeshSpace::MeshWall%id" % ndim

        # Constructors.
        x.add_constructor([constrefparam(vector_of_vector, "generators"),
                           constrefparam(vector, "xmin"),
                           constrefparam(vector, "xmax"),
                           param("unsigned int", "nx", default_value="20"),
                           param("unsigned int", "ny", default_value="20"),
                           param("unsigned int", "nz", default_value="20"),
                           param("double", "edgeTol", default_value="1.0e-8")])

        # Methods.
        x.add_method("addBoundary", None, [constrefparam(meshwall, "boundary")])
        x.add_method("subRegion", None, [constrefparam(vector, "point"),
                                         param("unsigned int", "i"),
                                         param("unsigned int", "j"),
                                         param("unsigned int", "k")], is_const=True)

        return

    #---------------------------------------------------------------------------
    # Add standalone functions.
    #---------------------------------------------------------------------------
    def addFunctions(self, meshname, ndim):

        # Object names.
        dim = "Spheral::Dim<%i>" % ndim
        vector = "Spheral::Vector%id" % ndim
        vector_of_vector = "vector_of_Vector%id" % ndim
        vector_of_symtensor = "vector_of_SymTensor%id" % ndim
        vector_of_nodelist = "vector_of_NodeList%id" % ndim
        vector_of_boundary = "vector_of_Boundary%id" % ndim
        nodelist = "Spheral::NodeSpace::NodeList%id" % ndim

        self.space.add_function("computeGenerators",
                                None,
                                [refparam(vector_of_nodelist, "nodeLists"),
                                 refparam(vector_of_boundary, "boundaries"),
                                 constrefparam(vector, "xmin"),
                                 constrefparam(vector, "xmax"),
                                 refparam(vector_of_vector, "positions"),
                                 refparam(vector_of_symtensor, "Hs"),
                                 refparam("vector_of_unsigned", "offsets")],
                                custom_name = "computeGenerators%id" % ndim,
                                template_parameters = [dim],
                                docstring = "Compute the parallel safe generators for a set of NodeLists.")

        self.space.add_function("generateMesh",
                                None,
                                [refparam(vector_of_nodelist, "nodeLists"),
                                 refparam(vector_of_boundary, "boundaries"),
                                 constrefparam(vector, "xmin"),
                                 constrefparam(vector, "xmax"),
                                 param("bool", "generateVoid"),
                                 param("bool", "generateParallelConnectivity"),
                                 param("bool", "removeBoundaryZones"),
                                 param("double", "voidThreshold"),
                                 refparam(meshname, "mesh"),
                                 refparam(nodelist, "voidNodes")],
                                custom_name = "generateMesh%id" % ndim,
                                template_parameters = [dim],
                                docstring = "Generate a mesh for the set of NodeLists.")

        self.space.add_function("hashPositionWrapper",
                                retval(ptr("PyObject"), caller_owns_return=True),
                                [constrefparam(vector, "position"),
                                 constrefparam(vector, "xmin"),
                                 constrefparam(vector, "xmax"),
                                 constrefparam(vector, "boxInv")],
                                template_parameters = [dim],
                                custom_name = "hashPosition")

        return
