from pybindgen import *

from PBGutils import *

#-------------------------------------------------------------------------------
# The class for wrapping this module.
#-------------------------------------------------------------------------------
class Parmetis:

    #---------------------------------------------------------------------------
    # Add the types to the module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir):

        # Includes.
        mod.add_include('"%s/ParmetisTypes.hh"' % srcdir)

        # Namespaces.
        Spheral = mod.add_cpp_namespace("Spheral")
        space = Spheral.add_cpp_namespace("PartitionSpace")
        RedistributeNodes2d = findObject(space, "RedistributeNodes2d")
        RedistributeNodes3d = findObject(space, "RedistributeNodes3d")

        # Expose types.
        self.ParmetisRedistributeNodes2d = addObject(space, "ParmetisRedistributeNodes2d", parent=RedistributeNodes2d, allow_subclassing=True)
        self.ParmetisRedistributeNodes3d = addObject(space, "ParmetisRedistributeNodes3d", parent=RedistributeNodes3d, allow_subclassing=True)

    #---------------------------------------------------------------------------
    # Generate the bindings.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        self.parmetisBindings(self.ParmetisRedistributeNodes2d, 2)
        self.parmetisBindings(self.ParmetisRedistributeNodes3d, 3)

        return

    #---------------------------------------------------------------------------
    # Parmetis
    #---------------------------------------------------------------------------
    def parmetisBindings(self, x, ndim):

        database = "Spheral::DataBaseSpace::DataBase%id" % ndim
        vector_of_boundary = "vector_of_Boundary%id" % ndim
        
        # Constructors.
        x.add_constructor([param("double", "extent")])

        # Methods.
        x.add_method("redistributeNodes", None, [refparam(database, "dataBase"),
                                                 param(vector_of_boundary, "boundaries", default_value="%s()" % vector_of_boundary)],
                     is_virtual = True)
        x.add_method("refineAndRedistributeNodes", None, [refparam(database, "dataBase"),
                                                          param(vector_of_boundary, "boundaries", default_value="%s()" % vector_of_boundary)])

        # Attributes.
        x.add_instance_attribute("normalizedNodeExtent", "double", getter="normalizedNodeExtent", setter="setNormalizedNodeExtent")

        return
