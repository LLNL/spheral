"""
Provides wrappers for the Silo library.
"""

from PYB11Generator import *

includes = ['"Geometry/Dimension.hh"',
            '"Geometry/GeomPlane.hh"',
            '"SiloWrappers.hh"']

namespaces = ["silo",
              "Spheral"]

#-------------------------------------------------------------------------------
class DBfile:
    "Opaque object for silo file struct"

#-------------------------------------------------------------------------------
@PYB11cppname("DBoptlist_wrapper")
class DBoptlist:
    "The silo optlist collection."

    # Constructor
    def pyinit(self,
               maxopts = ("int", "1024")):
        "Construct with the given number of avilable option slots."

    #...........................................................................
    # addOption
    # @PYB11template("T")
    # @PYB11implementation("""[](DBoptlist_wrapper& self,
    #                            int option,
    #                            %(T)s& value) {
    #                                return DBAddOption(self.mOptlistPtr, option, static_cast<void*>(&value));
    #                            }""")
    # def addOption(self, option="int", value="%(T)s&"):
    #     return "int"

    # @PYB11template("T")
    # @PYB11implementation("""[](DBoptlist_wrapper& self,
    #                            int option,
    #                            int option_size,
    #                            std::vector<%(T)s>& value) {
    #                                auto value_size = value.size();
    #                                auto result = DBAddOption(self.mOptlistPtr, option_size, static_cast<void*>(&value_size));
    #                                if (result != 0) return result;
    #                                return DBAddOption(self.mOptlistPtr, option, static_cast<void*>(&(value.front())));
    #                            }""")
    # def addOptionVec(self, option="int", value="std::vector<%(T)s>&"):
    #     return "int"

    # addOptionInt =       PYB11TemplateMethod(addOption,    template_parameters="int",    pyname="addOption")
    # addOptionDouble =    PYB11TemplateMethod(addOption,    template_parameters="double", pyname="addOption")
    # addOptionVecInt =    PYB11TemplateMethod(addOptionVec, template_parameters="int",    pyname="addOption")
    # addOptionVecDouble = PYB11TemplateMethod(addOptionVec, template_parameters="double", pyname="addOption")

    # @PYB11pyname("addOption")
    # @PYB11implementation("""[](DBoptlist_wrapper& self,
    #                            int option,
    #                            std::string& value) {
    #                                return DBAddOption(self.mOptlistPtr, option, const_cast<char*>(value.c_str()));
    #                            }""")
    # def addOptionString(self, option="int", value="std::string&"):
    #     return "int"

    # @PYB11pyname("addOption")
    # @PYB11implementation("""[](DBoptlist_wrapper& self,
    #                            int option,
    #                            int option_size,
    #                            std::vector<std::string>& value) {
    #                                auto value_size = value.size();
    #                                auto result = DBAddOption(self.mOptlistPtr, option_size, static_cast<void*>(&value_size));
    #                                if (result != 0) return result;
    #                                char** charArray = new char*[value.size()];
    #                                for (auto k = 0; k < value.size(); ++k) {
    #                                  charArray[k] = new char[value[k].size() + 1];
    #                                  strcpy(charArray[k], value[k].c_str());
    #                                }
    #                                return DBAddOption(self.mOptlistPtr, option, static_cast<void*>(charArray));
    #                            }""")
    # def addOptionVecString(self, option="int", value="std::vector<std::string>&"):
    #     return "int"

    for (T, Label) in (("int", "Int"),
                       ("double", "Double"),
                       ("std::string", "String")):
        exec('''
@PYB11pycppname("addOption")
def addOption%(Label)s(self,
                       option = "int",
                       value = "const %(T)s&"):
    return "int"

@PYB11pycppname("getOption")
def getOption%(Label)s(self,
                       option = "int"):
            return "%(T)s"

@PYB11pycppname("addOption")
def addOptionVec%(Label)s(self,
                          option = "int",
                          option_size = "int",
                          value = "const std::vector<%(T)s>&"):
    return "int"

@PYB11pycppname("getOption")
def getOptionVec%(Label)s(self,
                          option = "int",
                          option_size = "int"):
    return "std::vector<%(T)s>"
''' % {"T" : T,
       "Label" : Label})

#-------------------------------------------------------------------------------
@PYB11cppname("DBmrgtree_wrapper")
class DBmrgtree:
    "Another silo thingus."

    # Constructor
    def pyinit(self,
               mesh_type = ("int", "DB_POINTMESH"),
               info_bits = ("int", "0"),
               max_children = ("int", "1024"),
               optlist = ("silo::DBoptlist_wrapper", "silo::DBoptlist_wrapper(1024)")):
        "Constructor"

    # Methods for defining properties
    @PYB11const
    @PYB11cppname("name")
    @PYB11ignore
    def getname(self):
        return "std::string"

    @PYB11cppname("name")
    @PYB11ignore
    def setname(self, val="std::string"):
        return "void"

    @PYB11const
    @PYB11cppname("src_mesh_name")
    @PYB11ignore
    def getsrc_mesh_name(self):
        return "std::string"

    @PYB11cppname("src_mesh_name")
    @PYB11ignore
    def setsrc_mesh_name(self, val="std::string"):
        return "void"

    @PYB11const
    @PYB11cppname("src_mesh_type")
    @PYB11ignore
    def getsrc_mesh_type(self):
        return "int"

    @PYB11cppname("src_mesh_type")
    @PYB11ignore
    def setsrc_mesh_type(self, val="int"):
        return "void"

    @PYB11const
    @PYB11cppname("type_info_bits")
    @PYB11ignore
    def gettype_info_bits(self):
        return "int"

    @PYB11cppname("type_info_bits")
    @PYB11ignore
    def settype_info_bits(self, val="int"):
        return "void"

    @PYB11const
    @PYB11cppname("num_nodes")
    @PYB11ignore
    def getnum_nodes(self):
        return "int"

    @PYB11cppname("num_nodes")
    @PYB11ignore
    def setnum_nodes(self, val="int"):
        return "void"

    # Properties
    name = property(getname, setname)
    src_mesh_name = property(getsrc_mesh_name, setsrc_mesh_name)
    src_mesh_type = property(getsrc_mesh_type, setsrc_mesh_type)
    type_info_bits = property(gettype_info_bits, settype_info_bits)
    num_nodes = property(getnum_nodes, setnum_nodes)

#-------------------------------------------------------------------------------
# STL types
# vector_of_DBoptlist = PYB11_bind_vector("silo::DBoptlist_wrapper", opaque=True)

#-------------------------------------------------------------------------------
# Module methods
@PYB11returnpolicy("reference")
@PYB11cppname("DBCreate_wrap")
def DBCreate(pathName = "std::string",
             mode = "int",
             target = "int",
             fileInfo = "std::string",
             fileType = "int"):
    "Create a silo database"
    return "DBfile*"

@PYB11returnpolicy("reference")
@PYB11cppname("DBOpen_wrap")
def DBOpen(name = "std::string",
           type = "int",
           mode = "int"):
    "Open a silo database"
    return "DBfile*"

@PYB11namespace("silo")
def DBClose():
    "Close a silo database"

@PYB11namespace("silo")
def DBMkDir():
    "Make a new path in a silo database"

@PYB11namespace("silo")
def DBSetDir():
    "Go to the given path in a silo database"

@PYB11namespace("silo")
def DBGetDir():
    "Get the current path in the silo database."

@PYB11namespace("silo")
def DBCpDir():
    "Copy a directory to another."

@PYB11namespace("silo")
def DBPutMultimesh():
    "Write a multimesh"

@PYB11namespace("silo")
def DBPutMultimat():
    "Write a multimat"

@PYB11namespace("silo")
def DBPutMultivar():
    "Write a multivar"

@PYB11namespace("silo")
def DBPutMaterial():
    "Write a material"

@PYB11namespace("silo")
def DBPutUcdmesh():
    "Write a UCD mesh"

@PYB11namespace("silo")
def DBPutQuadmesh():
    "Write a quad mesh"

@PYB11namespace("silo")
def DBPutDefvars():
    "Write var definitions"

@PYB11namespace("silo")
def DBPutZonelist2():
    "Write a zonelist"

@PYB11namespace("silo")
def DBPutPHZonelist():
    "Write a polyhedral zonelist"

@PYB11namespace("silo")
def DBPutPointmesh():
    "Write a point mesh"

@PYB11namespace("silo")
def DBAddRegion():
    "Add a region to a silo data base"

@PYB11namespace("silo")
def DBSetCwr():
    "set cwr?"

@PYB11namespace("silo")
def DBGetCwr():
    "Get the cwr"

@PYB11namespace("silo")
def DBPutMrgtree():
    "Write an mrg tree"

#-------------------------------------------------------------------------------
@PYB11template("T")
@PYB11namespace("silo")
def DBWrite():
    "Write a %(T)s to a silo database."

@PYB11template("T")
@PYB11namespace("silo")
def DBWrite_vector():
    "Write a std::vector<%(T)s> to a silo database."

@PYB11template("T")
@PYB11namespace("silo")
def DBWrite_vector_of_vector():
    "Write a std::vector<std::vector<%(T)s>> to a silo database."

@PYB11template("T")
@PYB11namespace("silo")
def DBPutCompoundarray():
    "Write a compound array of %(T)s to a silo database."

@PYB11template("T")
@PYB11namespace("silo")
def DBReadVar():
    "Read a %(T)s from a silo database."

@PYB11template("T")
@PYB11namespace("silo")
def DBPutUcdvar1():
    "Write a UCD scalar variable of %(T)s to a silo database."

@PYB11template("T")
@PYB11namespace("silo")
def DBPutQuadvar1():
    "Write a quad mesh scalar variable of %(T)s to a silo database."

@PYB11template("T")
@PYB11namespace("silo")
def DBPutPointvar1():
    "Write a point mesh scalar variable of %(T)s to a silo database."

@PYB11template("T")
@PYB11namespace("silo")
def DBPutUcdvar():
    "Write a UCD mesh variable of %(T)s to a silo database."

@PYB11template("T")
@PYB11namespace("silo")
def DBPutQuadvar():
    "Write a quad mesh variable of %(T)s to a silo database."

@PYB11template("T")
@PYB11namespace("silo")
def DBPutPointvar():
    "Write a point mesh variable of %(T)s to a silo database."

for d in ("int", "double"):
    exec('''
DBWrite_%(d)s = PYB11TemplateFunction(DBWrite, ("%(d)s",), pyname="DBWrite")
DBWrite_vector_%(d)s = PYB11TemplateFunction(DBWrite_vector, ("%(d)s",), pyname="DBWrite")
DBWrite_vector_of_vector_%(d)s = PYB11TemplateFunction(DBWrite_vector_of_vector, ("%(d)s",), pyname="DBWrite")
DBPutCompoundarray_%(d)s = PYB11TemplateFunction(DBPutCompoundarray, ("%(d)s", ), pyname="DBPutCompoundarray")
DBReadVar_%(d)s = PYB11TemplateFunction(DBReadVar, ("%(d)s",), pyname="DBReadVar")
DBPutUcdvar1_%(d)s = PYB11TemplateFunction(DBPutUcdvar1, ("%(d)s",), pyname="DBPutUcdvar1")
DBPutQuadvar1_%(d)s = PYB11TemplateFunction(DBPutQuadvar1, ("%(d)s",), pyname="DBPutUcdvar1")
DBPutPointvar1_%(d)s = PYB11TemplateFunction(DBPutPointvar1, ("%(d)s",), pyname="DBPutPointvar1")
''' % {"d" : d})

for d in ("Vector", "Tensor", "SymTensor"):
    for ndim in (2, 3):
        exec('''
DBPutUcdvar_%(ndim)i = PYB11TemplateFunction(DBPutUcdvar, "Dim<%(ndim)i>::%(d)s", pyname="DBPutUcdvar")
DBPutQuadvar_%(ndim)i = PYB11TemplateFunction(DBPutQuadvar, "Dim<%(ndim)i>::%(d)s", pyname="DBPutQuadvar")
DBPutPointvar_%(ndim)i = PYB11TemplateFunction(DBPutPointvar, "Dim<%(ndim)i>::%(d)s", pyname="DBPutPointvar")
''') % {"ndim" : ndim,
        "d" : d}

#-------------------------------------------------------------------------------
# Taken from the silo.h file, expose the #define variables as module attributes.
DB_ZONETYPE_POLYHEDRON = PYB11attr()
DB_ZONETYPE_TET = PYB11attr()
DB_ZONETYPE_PYRAMID = PYB11attr()
DB_ZONETYPE_PRISM = PYB11attr()
DB_ZONETYPE_HEX = PYB11attr()

DB_NETCDF = PYB11attr()
DB_PDB = PYB11attr()
DB_TAURUS = PYB11attr()
DB_UNKNOWN = PYB11attr()
DB_DEBUG = PYB11attr()
DB_HDF5X = PYB11attr()
DB_PDBP = PYB11attr()

DB_HDF5_SEC2_OBSOLETE = PYB11attr()
DB_HDF5_STDIO_OBSOLETE = PYB11attr()
DB_HDF5_CORE_OBSOLETE = PYB11attr()
DB_HDF5_MPIO_OBSOLETE = PYB11attr()
DB_HDF5_MPIOP_OBSOLETE = PYB11attr()

DB_H5VFD_DEFAULT = PYB11attr()
DB_H5VFD_SEC2 = PYB11attr()
DB_H5VFD_STDIO = PYB11attr()
DB_H5VFD_CORE = PYB11attr()
DB_H5VFD_LOG = PYB11attr()
DB_H5VFD_SPLIT = PYB11attr()
DB_H5VFD_DIRECT = PYB11attr()
DB_H5VFD_FAMILY = PYB11attr()
DB_H5VFD_MPIO = PYB11attr()
DB_H5VFD_MPIP = PYB11attr()
DB_H5VFD_SILO = PYB11attr()

DB_FILE_OPTS_H5_DEFAULT_DEFAULT = PYB11attr()
DB_FILE_OPTS_H5_DEFAULT_SEC2 = PYB11attr()
DB_FILE_OPTS_H5_DEFAULT_STDIO = PYB11attr()
DB_FILE_OPTS_H5_DEFAULT_CORE = PYB11attr()
DB_FILE_OPTS_H5_DEFAULT_LOG = PYB11attr()
DB_FILE_OPTS_H5_DEFAULT_SPLIT = PYB11attr()
DB_FILE_OPTS_H5_DEFAULT_DIRECT = PYB11attr()
DB_FILE_OPTS_H5_DEFAULT_FAMILY = PYB11attr()
DB_FILE_OPTS_H5_DEFAULT_MPIO = PYB11attr()
DB_FILE_OPTS_H5_DEFAULT_MPIP = PYB11attr()
DB_FILE_OPTS_H5_DEFAULT_SILO = PYB11attr()
DB_FILE_OPTS_H5_DEFAULT_SILO = PYB11attr()

DB_HDF5 = PYB11attr()
DB_HDF5_SEC2 = PYB11attr()
DB_HDF5_STDIO = PYB11attr()
DB_HDF5_CORE = PYB11attr()
DB_HDF5_LOG = PYB11attr()
DB_HDF5_SPLIT = PYB11attr()
DB_HDF5_DIRECT = PYB11attr()
DB_HDF5_FAMILY = PYB11attr()
DB_HDF5_MPIO = PYB11attr()
DB_HDF5_MPIOP = PYB11attr()
DB_HDF5_MPIP = PYB11attr()
DB_HDF5_SILO = PYB11attr()

DB_NFILES = PYB11attr()
DB_NFILTERS = PYB11attr()

DBAll = PYB11attr()
DBNone = PYB11attr()
DBCalc = PYB11attr()
DBMatMatnos = PYB11attr()
DBMatMatlist = PYB11attr()
DBMatMixList = PYB11attr()
DBCurveArrays = PYB11attr()
DBPMCoords = PYB11attr()
DBPVData = PYB11attr()
DBQMCoords = PYB11attr()
DBQVData = PYB11attr()
DBUMCoords = PYB11attr()
DBUMFacelist = PYB11attr()
DBUMZonelist = PYB11attr()
DBUVData = PYB11attr()
DBFacelistInfo = PYB11attr()
DBZonelistInfo = PYB11attr()
DBMatMatnames = PYB11attr()
DBUMGlobNodeNo = PYB11attr()
DBZonelistGlobZoneNo = PYB11attr()
DBMatMatcolors = PYB11attr()
DBCSGMBoundaryInfo = PYB11attr()
DBCSGMZonelist = PYB11attr()
DBCSGMBoundaryNames = PYB11attr()
DBCSGVData = PYB11attr()
DBCSGZonelistZoneNames = PYB11attr()
DBCSGZonelistRegNames = PYB11attr()
DBMMADJNodelists = PYB11attr()
DBMMADJZonelists = PYB11attr()
DBPMGlobNodeNo = PYB11attr()

DBObjectType = PYB11enum(("DB_INVALID_OBJECT",
                          "DB_QUADRECT",
                          "DB_QUADCURV",
                          "DB_QUADMESH",
                          "DB_QUADVAR",
                          "DB_UCDMESH",
                          "DB_UCDVAR",
                          "DB_MULTIMESH",
                          "DB_MULTIVAR",
                          "DB_MULTIMAT",
                          "DB_MULTIMATSPECIES",
                          "DB_MULTIBLOCKMESH",
                          "DB_MULTIBLOCKVAR",
                          "DB_MULTIMESHADJ",
                          "DB_MATERIAL",
                          "DB_MATSPECIES",
                          "DB_FACELIST",
                          "DB_ZONELIST",
                          "DB_EDGELIST",
                          "DB_PHZONELIST",
                          "DB_CSGZONELIST",
                          "DB_CSGMESH",
                          "DB_CSGVAR",
                          "DB_CURVE",
                          "DB_DEFVARS",
                          "DB_POINTMESH",
                          "DB_POINTVAR",
                          "DB_ARRAY",
                          "DB_DIR",
                          "DB_VARIABLE",
                          "DB_MRGTREE",
                          "DB_GROUPELMAP",
                          "DB_MRGVAR",
                          "DB_USERDEF"),
                         export_values = True)

DBdatatype = PYB11enum(("DB_INT",
                        "DB_SHORT",
                        "DB_LONG",
                        "DB_FLOAT",
                        "DB_DOUBLE",
                        "DB_CHAR",
                        "DB_LONG_LONG",
                        "DB_NOTYPE"),
                       export_values = True)

DB_CLOBBER = PYB11attr()
DB_NOCLOBBER = PYB11attr()

DB_READ = PYB11attr()
DB_APPEND = PYB11attr()

DB_LOCAL = PYB11attr()
DB_SUN3 = PYB11attr()
DB_SUN4 = PYB11attr()
DB_SGI = PYB11attr()
DB_RS6000 = PYB11attr()
DB_CRAY = PYB11attr()
DB_INTEL = PYB11attr()

DBOPT_FIRST = PYB11attr()
DBOPT_ALIGN = PYB11attr()
DBOPT_COORDSYS = PYB11attr()
DBOPT_CYCLE = PYB11attr()
DBOPT_FACETYPE = PYB11attr()
DBOPT_HI_OFFSET = PYB11attr()
DBOPT_LO_OFFSET = PYB11attr()
DBOPT_LABEL = PYB11attr()
DBOPT_XLABEL = PYB11attr()
DBOPT_YLABEL = PYB11attr()
DBOPT_ZLABEL = PYB11attr()
DBOPT_MAJORORDER = PYB11attr()
DBOPT_NSPACE = PYB11attr()
DBOPT_ORIGIN = PYB11attr()
DBOPT_PLANAR = PYB11attr()
DBOPT_TIME = PYB11attr()
DBOPT_UNITS = PYB11attr()
DBOPT_XUNITS = PYB11attr()
DBOPT_YUNITS = PYB11attr()
DBOPT_ZUNITS = PYB11attr()
DBOPT_DTIME = PYB11attr()
DBOPT_USESPECMF = PYB11attr()
DBOPT_XVARNAME = PYB11attr()
DBOPT_YVARNAME = PYB11attr()
DBOPT_ZVARNAME = PYB11attr()
DBOPT_ASCII_LABEL = PYB11attr()
DBOPT_MATNOS = PYB11attr()
DBOPT_NMATNOS = PYB11attr()
DBOPT_MATNAME = PYB11attr()
DBOPT_NMAT = PYB11attr()
DBOPT_NMATSPEC = PYB11attr()
DBOPT_BASEINDEX = PYB11attr()
DBOPT_ZONENUM = PYB11attr()
DBOPT_NODENUM = PYB11attr()
DBOPT_BLOCKORIGIN = PYB11attr()
DBOPT_GROUPNUM = PYB11attr()
DBOPT_GROUPORIGIN = PYB11attr()
DBOPT_NGROUPS = PYB11attr()
DBOPT_MATNAMES = PYB11attr()
DBOPT_EXTENTS_SIZE = PYB11attr()
DBOPT_EXTENTS = PYB11attr()
DBOPT_MATCOUNTS = PYB11attr()
DBOPT_MATLISTS = PYB11attr()
DBOPT_MIXLENS = PYB11attr()
DBOPT_ZONECOUNTS = PYB11attr()
DBOPT_HAS_EXTERNAL_ZONES = PYB11attr()
DBOPT_PHZONELIST = PYB11attr()
DBOPT_MATCOLORS = PYB11attr()
DBOPT_BNDNAMES = PYB11attr()
DBOPT_REGNAMES = PYB11attr()
DBOPT_ZONENAMES = PYB11attr()
DBOPT_HIDE_FROM_GUI = PYB11attr()
DBOPT_TOPO_DIM = PYB11attr()
DBOPT_REFERENCE = PYB11attr()
DBOPT_GROUPINGS_SIZE = PYB11attr()
DBOPT_GROUPINGS = PYB11attr()
DBOPT_GROUPINGNAMES = PYB11attr()
DBOPT_ALLOWMAT0 = PYB11attr()
DBOPT_MRGTREE_NAME = PYB11attr()
DBOPT_REGION_PNAMES = PYB11attr()
DBOPT_TENSOR_RANK = PYB11attr()
DBOPT_MMESH_NAME = PYB11attr()
DBOPT_TV_CONNECTIVITY = PYB11attr()
DBOPT_DISJOINT_MODE = PYB11attr()
DBOPT_MRGV_ONAMES = PYB11attr()
DBOPT_MRGV_RNAMES = PYB11attr()
DBOPT_SPECNAMES = PYB11attr()
DBOPT_SPECCOLORS = PYB11attr()
DBOPT_LLONGNZNUM = PYB11attr()
DBOPT_CONSERVED = PYB11attr()
DBOPT_EXTENSIVE = PYB11attr()
DBOPT_MB_FILE_NS = PYB11attr()
DBOPT_MB_BLOCK_NS = PYB11attr()
DBOPT_MB_BLOCK_TYPE = PYB11attr()
DBOPT_MB_EMPTY_LIST = PYB11attr()
DBOPT_MB_EMPTY_COUNT = PYB11attr()
DBOPT_LAST = PYB11attr()
  
DBOPT_H5_FIRST = PYB11attr()
DBOPT_H5_VFD = PYB11attr()
DBOPT_H5_RAW_FILE_OPTS = PYB11attr()
DBOPT_H5_RAW_EXTENSION = PYB11attr()
DBOPT_H5_META_FILE_OPTS = PYB11attr()
DBOPT_H5_META_EXTENSION = PYB11attr()
DBOPT_H5_CORE_ALLOC_INC = PYB11attr()
DBOPT_H5_CORE_NO_BACK_STORE = PYB11attr()
DBOPT_H5_META_BLOCK_SIZE = PYB11attr()
DBOPT_H5_SMALL_RAW_SIZE = PYB11attr()
DBOPT_H5_ALIGN_MIN = PYB11attr()
DBOPT_H5_ALIGN_VAL = PYB11attr()
DBOPT_H5_DIRECT_MEM_ALIGN = PYB11attr()
DBOPT_H5_DIRECT_BLOCK_SIZE = PYB11attr()
DBOPT_H5_DIRECT_BUF_SIZE = PYB11attr()
DBOPT_H5_LOG_NAME = PYB11attr()
DBOPT_H5_LOG_BUF_SIZE = PYB11attr()
DBOPT_H5_MPIO_COMM = PYB11attr()
DBOPT_H5_MPIO_INFO = PYB11attr()
DBOPT_H5_MPIP_NO_GPFS_HINTS = PYB11attr()
DBOPT_H5_SIEVE_BUF_SIZE = PYB11attr()
DBOPT_H5_CACHE_NELMTS = PYB11attr()
DBOPT_H5_CACHE_NBYTES = PYB11attr()
DBOPT_H5_CACHE_POLICY = PYB11attr()
DBOPT_H5_FAM_SIZE = PYB11attr()
DBOPT_H5_FAM_FILE_OPTS = PYB11attr()
DBOPT_H5_USER_DRIVER_ID = PYB11attr()
DBOPT_H5_USER_DRIVER_INFO = PYB11attr()
DBOPT_H5_SILO_BLOCK_SIZE = PYB11attr()
DBOPT_H5_SILO_BLOCK_COUNT = PYB11attr()
DBOPT_H5_SILO_LOG_STATS = PYB11attr()
DBOPT_H5_SILO_USE_DIRECT = PYB11attr()
DBOPT_H5_LAST = PYB11attr()

DB_TOP = PYB11attr()
DB_NONE = PYB11attr()
DB_ALL = PYB11attr()
DB_ABORT = PYB11attr()
DB_SUSPEND = PYB11attr()
DB_RESUME = PYB11attr()
DB_ALL_AND_DRVR = PYB11attr()

DB_ROWMAJOR = PYB11attr()
DB_COLMAJOR = PYB11attr()

DB_NOTCENT = PYB11attr()
DB_NODECENT = PYB11attr()
DB_ZONECENT = PYB11attr()
DB_FACECENT = PYB11attr()
DB_BNDCENT = PYB11attr()
DB_EDGECENT = PYB11attr()
DB_BLOCKCENT = PYB11attr()

DB_CARTESIAN = PYB11attr()
DB_CYLINDRICAL = PYB11attr()
DB_SPHERICAL = PYB11attr()
DB_NUMERICAL = PYB11attr()
DB_OTHER = PYB11attr()

DB_RECTILINEAR = PYB11attr()
DB_CURVILINEAR = PYB11attr()

DB_AREA = PYB11attr()
DB_VOLUME = PYB11attr()

DB_ON = PYB11attr()
DB_OFF = PYB11attr()

DB_ABUTTING = PYB11attr()
DB_FLOATING = PYB11attr()

DB_VARTYPE_SCALAR = PYB11attr()
DB_VARTYPE_VECTOR = PYB11attr()
DB_VARTYPE_TENSOR = PYB11attr()
DB_VARTYPE_SYMTENSOR = PYB11attr()
DB_VARTYPE_ARRAY = PYB11attr()
DB_VARTYPE_MATERIAL = PYB11attr()
DB_VARTYPE_SPECIES = PYB11attr()
DB_VARTYPE_LABEL = PYB11attr()

DBCSG_QUADRIC_G = PYB11attr()
DBCSG_SPHERE_PR = PYB11attr()
DBCSG_ELLIPSOID_PRRR = PYB11attr()
DBCSG_PLANE_G = PYB11attr()
DBCSG_PLANE_X = PYB11attr()
DBCSG_PLANE_Y = PYB11attr()
DBCSG_PLANE_Z = PYB11attr()
DBCSG_PLANE_PN = PYB11attr()
DBCSG_PLANE_PPP = PYB11attr()
DBCSG_CYLINDER_PNLR = PYB11attr()
DBCSG_CYLINDER_PPR = PYB11attr()
DBCSG_BOX_XYZXYZ = PYB11attr()
DBCSG_CONE_PNLA = PYB11attr()
DBCSG_CONE_PPA = PYB11attr()
DBCSG_POLYHEDRON_KF = PYB11attr()
DBCSG_HEX_6F = PYB11attr()
DBCSG_TET_4F = PYB11attr()
DBCSG_PYRAMID_5F = PYB11attr()
DBCSG_PRISM_5F = PYB11attr()

DBCSG_QUADRATIC_G = PYB11attr()
DBCSG_CIRCLE_PR = PYB11attr()
DBCSG_ELLIPSE_PRR = PYB11attr()
DBCSG_LINE_G = PYB11attr()
DBCSG_LINE_X = PYB11attr()
DBCSG_LINE_Y = PYB11attr()
DBCSG_LINE_PN = PYB11attr()
DBCSG_LINE_PP = PYB11attr()
DBCSG_BOX_XYXY = PYB11attr()
DBCSG_ANGLE_PNLA = PYB11attr()
DBCSG_ANGLE_PPA = PYB11attr()
DBCSG_POLYGON_KP = PYB11attr()
DBCSG_TRI_3P = PYB11attr()
DBCSG_QUAD_4P = PYB11attr()

DBCSG_INNER = PYB11attr()
DBCSG_OUTER = PYB11attr()
DBCSG_ON = PYB11attr()
DBCSG_UNION = PYB11attr()
DBCSG_INTERSECT = PYB11attr()
DBCSG_DIFF = PYB11attr()
DBCSG_COMPLIMENT = PYB11attr()
DBCSG_XFORM = PYB11attr()
DBCSG_SWEEP = PYB11attr()

DB_PREORDER = PYB11attr()
DB_POSTORDER = PYB11attr()
DB_FROMCWR = PYB11attr()

DB_ZONETYPE_BEAM = PYB11attr()

DB_ZONETYPE_POLYGON = PYB11attr()
DB_ZONETYPE_TRIANGLE = PYB11attr()
DB_ZONETYPE_QUAD = PYB11attr()
