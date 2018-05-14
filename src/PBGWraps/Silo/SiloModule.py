from pybindgen import *

from PBGutils import *
from enumUtilities import *
from ref_return_value import *

from CXXTypesModule import generateStdVectorBindings

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class Silo:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir, dims):

        # Includes.
        mod.add_include('"%s/SiloTypes.hh"' % srcdir)
    
        # Namespace
        self.space = mod.add_cpp_namespace("silo")

        # Expose types.
        self.DBfile = addObject(self.space, "DBfile", allow_subclassing=True)
        self.DBoptlist = addObject(self.space, "DBoptlist_wrapper", custom_name="DBoptlist", allow_subclassing=True)
        self.DBmrgtree = addObject(self.space, "DBmrgtree_wrapper", custom_name="DBmrgtree", allow_subclassing=True)
        self.SiloAttributes = addStructAsEnumDefinition(self.space, "SiloAttributes", "%s/SiloTypes.hh" % srcdir)
        #self.SiloAttributes = addEnumDefinition(self.space, "SiloAttributes", "Silo/SiloTypes.hh")

        self.vector_of_DBoptlist = addObject(mod, "vector_of_DBoptlist", allow_subclassing=True)

        return

    #---------------------------------------------------------------------------
    # Generate bindings.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        self.generateDBoptlistBindings(self.DBoptlist)
        self.generateDBmrgtreeBindings(self.DBmrgtree)
        generateStdVectorBindings(self.vector_of_DBoptlist, "silo::DBoptlist_wrapper*", "vector_of_DBoptlist")

        # DBCreate
        self.space.add_function("DBCreate_wrap",
                                retval("DBfile*", reference_existing_object=True),
                                [param("std::string", "pathName"),
                                 param("int", "mode"),
                                 param("int", "target"),
                                 param("std::string", "fileInfo"),
                                 param("int", "fileType")],
                                custom_name = "DBCreate",
                                docstring = "Create a SILO file.")

        # DBOpen
        self.space.add_function("DBOpen_wrap",
                                retval("DBfile*", reference_existing_object=True),
                                [param("std::string", "pathName"),
                                 param("int", "type"),
                                 param("int", "mode")],
                                custom_name = "DBOpen",
                                docstring = "Open an existing SILO file.")

        # # DBMakeMrgtree
        # self.space.add_function("DBMakeMrgtree_wrap", 
        #                         retval("DBmrgtree*", reference_existing_object=True),
        #                         [param("int", "mesh_type"),
        #                          param("int", "info_bits"),
        #                          param("int", "max_children"),
        #                          refparam("silo::DBoptlist_wrapper", "optlist", default_value="silo::DBoptlist_wrapper(0)")],
        #                         custom_name = "DBMakeMrgtree",
        #                         docstring = "Create a DBmrgtree object.")

        # # DBFreeMrgtree
        # self.space.add_function("DBFreeMrgtree_wrap", 
        #                         retval("DBmrgtree*", reference_existing_object=True),
        #                         [refparam("DBmrgtree", "tree")],
        #                         docstring = "Free a DBmrgtree object.")

        # DBClose
        self.space.add_function("DBClose", "int", [refparam("DBfile", "file")],
                                docstring = "Close a SILO file.")

        # DBMkDir
        self.space.add_function("DBMkDir", "int", [refparam("DBfile", "file"), param("std::string", "dirname")],
                                docstring = "Create a new directory in a Silo file.")

        # DBSetDir
        self.space.add_function("DBSetDir", "int", [refparam("DBfile", "file"), param("std::string", "dirname")],
                                docstring = "Set the current directory in a Silo file.")

        # DBGetDir
        self.space.add_function("DBGetDir", "std::string", [refparam("DBfile", "file")],
                                docstring = "Get the current directory in a Silo file.")

        # DBCpDir
        self.space.add_function("DBCpDir", "int", [refparam("DBfile", "srcFile"), param("std::string", "srcDir"),
                                                   refparam("DBfile", "dstFile"), param("std::string", "dstDir")],
                                docstring = "Create a directory structure from one Silo file to another.")

        # DBPutMultimesh
        self.space.add_function("DBPutMultimesh", "int",
                                [refparam("DBfile", "file"),
                                 param("std::string", "name"),
                                 constrefparam("vector_of_string", "meshNames"),
                                 constrefparam("vector_of_int", "meshTypes"),
                                 refparam("silo::DBoptlist_wrapper", "optlist", default_value="silo::DBoptlist_wrapper(0)")],
                                docstring = "Write a multi-mesh to a SILO file.")

        # DBPutMultimat
        self.space.add_function("DBPutMultimat", "int",
                                [refparam("DBfile", "file"),
                                 param("std::string", "name"),
                                 constrefparam("vector_of_string", "matNames"),
                                 refparam("silo::DBoptlist_wrapper", "optlist", default_value="silo::DBoptlist_wrapper(0)")],
                                docstring = "Write a multi-material object to a SILO file.")

        # DBPutMultivar
        self.space.add_function("DBPutMultivar", "int",
                                [refparam("DBfile", "file"),
                                 param("std::string", "name"),
                                 refparam("vector_of_string", "varNames"),
                                 refparam("vector_of_int", "varTypes"),
                                 refparam("silo::DBoptlist_wrapper", "optlist", default_value="silo::DBoptlist_wrapper(0)")],
                                docstring = "Write a multi-block variable object.")

        # DBPutMaterial
        self.space.add_function("DBPutMaterial", "int",
                                [refparam("DBfile", "file"),
                                 param("std::string", "name"),
                                 param("std::string", "meshName"),
                                 refparam("vector_of_int", "matnos"),
                                 refparam("vector_of_int", "matlist"),
                                 param("vector_of_int", "dims"),
                                 refparam("vector_of_int", "mix_next"),
                                 refparam("vector_of_int", "mix_mat"),
                                 refparam("vector_of_int", "mix_zone"),
                                 refparam("vector_of_double", "mix_vf"),
                                 refparam("silo::DBoptlist_wrapper", "optlist", default_value="silo::DBoptlist_wrapper(0)")],
                                docstring = "Write a Silo material data object.")

        # DBPutUcdmesh
        self.space.add_function("DBPutUcdmesh", "int",
                                [refparam("DBfile", "file"),
                                 param("std::string", "name"),
                                 refparam("vector_of_vector_of_double", "coords"),
                                 param("int", "nzones"),
                                 param("std::string", "zonel_name"),
                                 param("std::string", "facel_name"),
                                 refparam("silo::DBoptlist_wrapper", "optlist", default_value="silo::DBoptlist_wrapper(0)")],
                                docstring = "Write UCD mesh object.")

        # DBPutQuadmesh
        self.space.add_function("DBPutQuadmesh", "int",
                                [refparam("DBfile", "file"),
                                 param("std::string", "name"),
                                 refparam("vector_of_vector_of_double", "coords"),
                                 refparam("silo::DBoptlist_wrapper", "optlist", default_value="silo::DBoptlist_wrapper(0)")],
                                docstring = "Write Quad mesh object.")

        # DBPutDefvars
        self.space.add_function("DBPutDefvars", "int",
                                [refparam("DBfile", "file"),
                                 param("std::string", "name"),
                                 refparam("vector_of_string", "varNames"),
                                 refparam("vector_of_int", "varTypes"),
                                 refparam("vector_of_string", "varDefs"),
                                 refparam("vector_of_DBoptlist", "optlists")],
                                docstring = "Write a collection of variable definitions.")

        # DBPutZonelist2
        self.space.add_function("DBPutZonelist2", "int",
                                [refparam("DBfile", "file"),
                                 param("std::string", "name"),
                                 param("unsigned int", "ndims"),
                                 refparam("vector_of_vector_of_int", "zoneNodes"),
                                 param("unsigned int", "low_offset"),
                                 param("unsigned int", "high_offset"),
                                 refparam("vector_of_int", "shapetype"),
                                 refparam("vector_of_int", "shapesize"),
                                 refparam("vector_of_int", "shapecount"),
                                 refparam("silo::DBoptlist_wrapper", "optlist", default_value="silo::DBoptlist_wrapper(0)")],
                                docstring = "Write a zone list using DBPutZoneList2.")

        # DBPutPHZonelist
        self.space.add_function("DBPutPHZonelist", "int",
                                [refparam("DBfile", "file"),
                                 param("std::string", "name"),
                                 refparam("vector_of_vector_of_int", "faceNodeLists"),
                                 refparam("vector_of_vector_of_int", "zoneFaceLists"),
                                 param("int", "low_offset"),
                                 param("int", "high_offset"),
                                 refparam("silo::DBoptlist_wrapper", "optlist", default_value="silo::DBoptlist_wrapper(0)")],
                                docstring = "Write polyhedral zone list object.")

        # DBPutPointmesh
        self.space.add_function("DBPutPointmesh", "int",
                                [refparam("DBfile", "file"),
                                 param("std::string", "name"),
                                 refparam("vector_of_vector_of_double", "coords"),
                                 refparam("silo::DBoptlist_wrapper", "optlist", default_value="silo::DBoptlist_wrapper(0)")],
                                docstring = "Write Pointmesh object.")

        # DBAddRegion
        self.space.add_function("DBAddRegion", "int",
                                [refparam("DBmrgtree_wrapper", "tree"),
                                 param("std::string", "reg_name"),
                                 param("int", "info_bits"),
                                 param("int", "max_children"),
                                 param("std::string", "maps_name"),
                                 refparam("vector_of_int", "seg_ids"),
                                 refparam("vector_of_int", "seg_lens"),
                                 refparam("vector_of_int", "seg_types"),
                                 refparam("silo::DBoptlist_wrapper", "optlist", default_value="silo::DBoptlist_wrapper(0)")],
                                docstring = "Add a region to a DBMrgtree.")

        # DBSetCwr
        self.space.add_function("DBSetCwr", "int",
                                [refparam("DBmrgtree_wrapper", "tree"),
                                 param("std::string", "path")],
                                docstring = "Set the current working region in tree.")

        # DBGetCwr
        self.space.add_function("DBGetCwr", retval("const char*", caller_owns_return=False),
                                [refparam("DBmrgtree_wrapper", "tree")],
                                docstring = "Get the current working region in tree.")

        # DBPutMrgtree
        self.space.add_function("DBPutMrgtree", "int",
                                [refparam("DBfile", "file"),
                                 param("std::string", "name"),
                                 param("std::string", "mesh_name"),
                                 refparam("DBmrgtree_wrapper", "tree"),
                                 refparam("silo::DBoptlist_wrapper", "optlist", default_value="silo::DBoptlist_wrapper(0)")],
                                docstring = "Write a DBmrgtree object to a file.")

        # Bind templated types for a variety of Silo data types.
        for type in ("int", "float", "double", "char"):
            ext = "_" + type.replace(" ", "_")
        
            # DBWrite
            self.space.add_function("DBWrite", "int",
                                    [refparam("DBfile", "file"),
                                     param("std::string", "varname"),
                                     param(type, "var")],
                                    template_parameters = [type],
                                    custom_name = "DBWrite",
                                    docstring = "Write a(n) %s to a SILO file." % type)

            # DBWrite
            self.space.add_function("DBWrite_vector", "int",
                                    [refparam("DBfile", "file"),
                                     param("std::string", "varname"),
                                     refparam("vector_of_%s" % type, "var")],
                                    template_parameters = [type],
                                    custom_name = "DBWrite",
                                    docstring = "Write a std::vector<%s> to a SILO file." % type)
            # DBWrite
            self.space.add_function("DBWrite_vector_of_vector", "int",
                                    [refparam("DBfile", "file"),
                                     param("std::string", "varname"),
                                     refparam("vector_of_vector_of_%s" % type, "var")],
                                    template_parameters = [type],
                                    custom_name = "DBWrite",
                                    docstring = "Write a std::vector<std::vector<%s>> to a SILO file." % type)
            # DBPutCompoundarray
            self.space.add_function("DBPutCompoundarray", type,
                                    [refparam("DBfile", "file"),
                                     param("std::string", "name"),
                                     constrefparam("vector_of_string", "elemNames"),
                                     constrefparam("vector_of_vector_of_%s" % type, "values"),
                                     refparam("silo::DBoptlist_wrapper", "optlist", default_value="silo::DBoptlist_wrapper(0)")],
                                    template_parameters = [type],
                                    custom_name = "DBPutCompoundarray",
                                    docstring = "Write a compound array of %s." % type)

            # DBReadVar
            self.space.add_function("DBReadVar", type,
                                    [refparam("DBfile", "file"),
                                     param("std::string", "varname")],
                                    template_parameters = [type],
                                    custom_name = "DBReadVar_%s" % type,
                                    docstring = "Read a(n) %s from a SILO file." % type)

            # DBPutUcdvar1
            self.space.add_function("DBPutUcdvar1", "int",
                                    [refparam("DBfile", "file"),
                                     param("std::string", "name"),
                                     param("std::string", "meshName"),
                                     refparam("vector_of_%s" % type, "values"),
                                     refparam("vector_of_%s" % type, "mixValues"),
                                     param("int", "centering"),
                                     refparam("silo::DBoptlist_wrapper", "optlist", default_value="silo::DBoptlist_wrapper(0)")],
                                    template_parameters = [type],
                                    custom_name = "DBPutUcdvar1",
                                    docstring = "Write a Ucd array of type %s." % type)

            # DBPutQuadvar1
            self.space.add_function("DBPutQuadvar1", "int",
                                    [refparam("DBfile", "file"),
                                     param("std::string", "name"),
                                     param("std::string", "meshName"),
                                     refparam("vector_of_%s" % type, "values"),
                                     refparam("vector_of_%s" % type, "mixValues"),
                                     param("int", "centering"),
                                     refparam("vector_of_int", "vardims"),
                                     refparam("silo::DBoptlist_wrapper", "optlist", default_value="silo::DBoptlist_wrapper(0)")],
                                    template_parameters = [type],
                                    custom_name = "DBPutQuadvar1",
                                    docstring = "Write a Quad array of type %s." % type)

            # DBPutPointvar1
            self.space.add_function("DBPutPointvar1", "int",
                                    [refparam("DBfile", "file"),
                                     param("std::string", "name"),
                                     param("std::string", "meshName"),
                                     refparam("vector_of_%s" % type, "values"),
                                     refparam("silo::DBoptlist_wrapper", "optlist", default_value="silo::DBoptlist_wrapper(0)")],
                                    template_parameters = [type],
                                    custom_name = "DBPutPointvar1",
                                    docstring = "Write a Pointmesh array of type %s." % type)

        # Support for writing compound types.
        for type in ("Vector2d", "Vector3d", "Tensor2d", "Tensor3d", "SymTensor2d", "SymTensor3d"):

            # DBPutUcdvar
            self.space.add_function("DBPutUcdvar", "int",
                                    [refparam("DBfile", "file"),
                                     param("std::string", "name"),
                                     param("std::string", "meshName"),
                                     refparam("vector_of_%s" % type, "values"),
                                     refparam("vector_of_%s" % type, "mixValues"),
                                     param("int", "centering"),
                                     refparam("silo::DBoptlist_wrapper", "optlist", default_value="silo::DBoptlist_wrapper(0)")],
                                    template_parameters = [type],
                                    custom_name = "DBPutUcdvar",
                                    docstring = "Write a Ucd array of type %s." % type)

            # DBPutQuadvar
            self.space.add_function("DBPutQuadvar", "int",
                                    [refparam("DBfile", "file"),
                                     param("std::string", "name"),
                                     param("std::string", "meshName"),
                                     refparam("vector_of_%s" % type, "values"),
                                     refparam("vector_of_%s" % type, "mixValues"),
                                     param("int", "centering"),
                                     refparam("vector_of_int", "vardims"),
                                     refparam("silo::DBoptlist_wrapper", "optlist", default_value="silo::DBoptlist_wrapper(0)")],
                                    template_parameters = [type],
                                    custom_name = "DBPutQuadvar",
                                    docstring = "Write a Quad array of type %s." % type)

            # DBPutPointvar
            self.space.add_function("DBPutPointvar", "int",
                                    [refparam("DBfile", "file"),
                                     param("std::string", "name"),
                                     param("std::string", "meshName"),
                                     refparam("vector_of_%s" % type, "values"),
                                     refparam("silo::DBoptlist_wrapper", "optlist", default_value="silo::DBoptlist_wrapper(0)")],
                                    template_parameters = [type],
                                    custom_name = "DBPutPointvar",
                                    docstring = "Write a Pointmesh array of type %s." % type)

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return ["silo"]

    #---------------------------------------------------------------------------
    # DBoptlist Bindings.
    #---------------------------------------------------------------------------
    def generateDBoptlistBindings(self, x):

        # Constructors.
        x.add_constructor([param("int", "maxopts", default_value="1024")])

        # Methods.
        for ValueType in ("int", "double", "std::string"):
            ValueBase = ValueType.split(":")[-1]
            x.add_method("addOption", "int", [param("int", "option"), param(ValueType, "value")],
                         template_parameters = [ValueType],
                         custom_name = "addOption")
            x.add_method("addOption", "int", [param("int", "option"), 
                                              param("int", "option_size"),
                                              constrefparam("vector_of_%s" % ValueBase, "value")],
                         template_parameters = [ValueType],
                         custom_name = "addOption")

            x.add_method("getOption", ValueType, [param("int", "option")],
                         template_parameters = [ValueType],
                         custom_name = "getOption%s" % ValueBase.upper())
            x.add_method("getOption", "vector_of_%s" % ValueBase, [param("int", "option"),
                                                                   param("int", "option_size")],
                         template_parameters = [ValueType],
                         custom_name = "getOption%s" % ValueBase.upper())

        return

    #---------------------------------------------------------------------------
    # DBmrgtree Bindings.
    #---------------------------------------------------------------------------
    def generateDBmrgtreeBindings(self, x):

        # Constructors.
        x.add_constructor([param("int", "mesh_type", default_value="DB_POINTMESH"),
                           param("int", "info_bits", default_value="0"),
                           param("int", "max_children", default_value="1024"),
                           param("silo::DBoptlist_wrapper", "optlist", default_value="silo::DBoptlist_wrapper(0)")])

        # Attributes.
        x.add_instance_attribute("name", "std::string", getter="name", setter="name")
        x.add_instance_attribute("src_mesh_name", "std::string", getter="src_mesh_name", setter="src_mesh_name")
        x.add_instance_attribute("src_mesh_type", "int", getter="src_mesh_type", setter="src_mesh_type")
        x.add_instance_attribute("type_info_bits", "int", getter="type_info_bits", setter="type_info_bits")
        x.add_instance_attribute("num_nodes", "int", getter="num_nodes", setter="num_nodes")
        #x.add_instance_attribute("mrgvar_onames", retval("char**", caller_owns_return=False))
        #x.add_instance_attribute("mrgvar_rnames", retval("char**", caller_owns_return=False))

        return
