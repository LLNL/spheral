"""
Provides wrappers for the Silo library.
"""

from PYB11Generator import *

includes = ['"Geometry/Dimension.hh"',
            '"SiloWrappers.hh"']

#-------------------------------------------------------------------------------
class DBfile:
    "Opaque object for silo file struct"

#-------------------------------------------------------------------------------
@PYB11namespace("silo")
@PYB11cppname("DBoptlist_wrapper")
class DBoptlist:
    "The silo optlist collection."

    # Constructor
    def pyinit(self,
               maxopts = ("int", "1024")):
        "Construct with the given number of avilable option slots."

    def addOption(self,
                  option = "int",
                  value = "const int&"):
        return "int"

    @PYB11pycppname("addOption")
    def addOption1(self,
                   option = "int",
                   value = "const double&"):
        return "int"

    @PYB11pycppname("addOption")
    def addOption2(self,
                   option = "int",
                   value = "const std::string&"):
        return "int"

    @PYB11pycppname("addOption")
    def addOption3(self,
                   option = "int",
                   value = "const std::vector<int>&"):
        return "int"

    @PYB11pycppname("addOption")
    def addOption4(self,
                   option = "int",
                   value = "const std::vector<double>&"):
        return "int"

    @PYB11pycppname("addOption")
    def addOption5(self,
                   option = "int",
                   value = "const std::vector<std::string>&"):
        return "int"

    def getOption(self,
                  option = "int"):
        return "int"

    @PYB11pycppname("getOption")
    def getOption1(self,
                   option = "int"):
        return "double"

    @PYB11pycppname("getOption")
    def getOption2(self,
                   option = "int"):
        return "std::string"

    @PYB11pycppname("getOption")
    def getOption3(self,
                   option = "int"):
        return "std::vector<int>"

    @PYB11pycppname("getOption")
    def getOption4(self,
                   option = "int"):
        return "std::vector<double>"

    @PYB11pycppname("getOption")
    def getOption5(self,
                   option = "int"):
        return "std::vector<std::string>"

#-------------------------------------------------------------------------------
@PYB11namespace("silo")
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
vector_of_DBoptlist = PYB11_bind_vector("silo::DBoptlist_wrapper")

#-------------------------------------------------------------------------------
# Module methods
@PYB11cppname("silo::DBCreate_wrap")
def DBCreate():
    "Create a silo database"

@PYB11cppname("silo::DBOpen_wrap")
def DBOpen():
    "Open a silo database"

@PYB11namespace("silo")
def DBClose():
    "Close a silo database"

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

@PYB11template("T")
def DBWrite():
    "Write a %(T)s to a silo database."

DBWrite_int = PYB11TemplateFunction(DBWrite, ("int",), cppname="DBWrite", pyname="DBWrite")

#-------------------------------------------------------------------------------
# Taken from the silo.h file, expose the #define variables as module attributes.
DB_ZONETYPE_POLYHEDRON = PYB11attr("DB_ZONETYPE_POLYHEDRON")
DB_ZONETYPE_TET = PYB11attr("DB_ZONETYPE_TET")
DB_ZONETYPE_PYRAMID = PYB11attr("DB_ZONETYPE_PYRAMID")
DB_ZONETYPE_PRISM = PYB11attr("DB_ZONETYPE_PRISM")
DB_ZONETYPE_HEX = PYB11attr("DB_ZONETYPE_HEX")

DB_NETCDF = PYB11attr("DB_NETCDF")
DB_PDB = PYB11attr("DB_PDB")
DB_TAURUS = PYB11attr("DB_TAURUS")
DB_UNKNOWN = PYB11attr("DB_UNKNOWN")
DB_DEBUG = PYB11attr("DB_DEBUG")
DB_HDF5X = PYB11attr("DB_HDF5X")
DB_PDBP = PYB11attr("DB_PDBP")

DB_HDF5_SEC2_OBSOLETE = PYB11attr("DB_HDF5_SEC2_OBSOLETE")
DB_HDF5_STDIO_OBSOLETE = PYB11attr("DB_HDF5_STDIO_OBSOLETE")
DB_HDF5_CORE_OBSOLETE = PYB11attr("DB_HDF5_CORE_OBSOLETE")
DB_HDF5_MPIO_OBSOLETE = PYB11attr("DB_HDF5_MPIO_OBSOLETE")
DB_HDF5_MPIOP_OBSOLETE = PYB11attr("DB_HDF5_MPIOP_OBSOLETE")

DB_H5VFD_DEFAULT = PYB11attr("DB_H5VFD_DEFAULT")
DB_H5VFD_SEC2 = PYB11attr("DB_H5VFD_SEC2")
DB_H5VFD_STDIO = PYB11attr("DB_H5VFD_STDIO")
DB_H5VFD_CORE = PYB11attr("DB_H5VFD_CORE")
DB_H5VFD_LOG = PYB11attr("DB_H5VFD_LOG")
DB_H5VFD_SPLIT = PYB11attr("DB_H5VFD_SPLIT")
DB_H5VFD_DIRECT = PYB11attr("DB_H5VFD_DIRECT")
DB_H5VFD_FAMILY = PYB11attr("DB_H5VFD_FAMILY")
DB_H5VFD_MPIO = PYB11attr("DB_H5VFD_MPIO")
DB_H5VFD_MPIP = PYB11attr("DB_H5VFD_MPIP")
DB_H5VFD_SILO = PYB11attr("DB_H5VFD_SILO")

DB_FILE_OPTS_H5_DEFAULT_DEFAULT = PYB11attr("DB_FILE_OPTS_H5_DEFAULT_DEFAULT")
DB_FILE_OPTS_H5_DEFAULT_SEC2 = PYB11attr("DB_FILE_OPTS_H5_DEFAULT_SEC2")
DB_FILE_OPTS_H5_DEFAULT_STDIO = PYB11attr("DB_FILE_OPTS_H5_DEFAULT_STDIO")
DB_FILE_OPTS_H5_DEFAULT_CORE = PYB11attr("DB_FILE_OPTS_H5_DEFAULT_CORE")
DB_FILE_OPTS_H5_DEFAULT_LOG = PYB11attr("DB_FILE_OPTS_H5_DEFAULT_LOG")
DB_FILE_OPTS_H5_DEFAULT_SPLIT = PYB11attr("DB_FILE_OPTS_H5_DEFAULT_SPLIT")
DB_FILE_OPTS_H5_DEFAULT_DIRECT = PYB11attr("DB_FILE_OPTS_H5_DEFAULT_DIRECT")
DB_FILE_OPTS_H5_DEFAULT_FAMILY = PYB11attr("DB_FILE_OPTS_H5_DEFAULT_FAMILY")
DB_FILE_OPTS_H5_DEFAULT_MPIO = PYB11attr("DB_FILE_OPTS_H5_DEFAULT_MPIO")
DB_FILE_OPTS_H5_DEFAULT_MPIP = PYB11attr("DB_FILE_OPTS_H5_DEFAULT_MPIP")
DB_FILE_OPTS_H5_DEFAULT_SILO = PYB11attr("DB_FILE_OPTS_H5_DEFAULT_SILO")
DB_FILE_OPTS_H5_DEFAULT_SILO = PYB11attr("DB_FILE_OPTS_H5_DEFAULT_SILO")

DB_HDF5 = PYB11attr("DB_HDF5")
DB_HDF5_SEC2 = PYB11attr("DB_HDF5_SEC2")
DB_HDF5_STDIO = PYB11attr("DB_HDF5_STDIO")
DB_HDF5_CORE = PYB11attr("DB_HDF5_CORE")
DB_HDF5_LOG = PYB11attr("DB_HDF5_LOG")
DB_HDF5_SPLIT = PYB11attr("DB_HDF5_SPLIT")
DB_HDF5_DIRECT = PYB11attr("DB_HDF5_DIRECT")
DB_HDF5_FAMILY = PYB11attr("DB_HDF5_FAMILY")
DB_HDF5_MPIO = PYB11attr("DB_HDF5_MPIO")
DB_HDF5_MPIOP = PYB11attr("DB_HDF5_MPIOP")
DB_HDF5_MPIP = PYB11attr("DB_HDF5_MPIP")
DB_HDF5_SILO = PYB11attr("DB_HDF5_SILO")

DB_NFILES = PYB11attr("DB_NFILES")
DB_NFILTERS = PYB11attr("DB_NFILTERS")

DBAll = PYB11attr("DBAll")
DBNone = PYB11attr("DBNone")
DBCalc = PYB11attr("DBCalc")
DBMatMatnos = PYB11attr("DBMatMatnos")
DBMatMatlist = PYB11attr("DBMatMatlist")
DBMatMixList = PYB11attr("DBMatMixList")
DBCurveArrays = PYB11attr("DBCurveArrays")
DBPMCoords = PYB11attr("DBPMCoords")
DBPVData = PYB11attr("DBPVData")
DBQMCoords = PYB11attr("DBQMCoords")
DBQVData = PYB11attr("DBQVData")
DBUMCoords = PYB11attr("DBUMCoords")
DBUMFacelist = PYB11attr("DBUMFacelist")
DBUMZonelist = PYB11attr("DBUMZonelist")
DBUVData = PYB11attr("DBUVData")
DBFacelistInfo = PYB11attr("DBFacelistInfo")
DBZonelistInfo = PYB11attr("DBZonelistInfo")
DBMatMatnames = PYB11attr("DBMatMatnames")
DBUMGlobNodeNo = PYB11attr("DBUMGlobNodeNo")
DBZonelistGlobZoneNo = PYB11attr("DBZonelistGlobZoneNo")
DBMatMatcolors = PYB11attr("DBMatMatcolors")
DBCSGMBoundaryInfo = PYB11attr("DBCSGMBoundaryInfo")
DBCSGMZonelist = PYB11attr("DBCSGMZonelist")
DBCSGMBoundaryNames = PYB11attr("DBCSGMBoundaryNames")
DBCSGVData = PYB11attr("DBCSGVData")
DBCSGZonelistZoneNames = PYB11attr("DBCSGZonelistZoneNames")
DBCSGZonelistRegNames = PYB11attr("DBCSGZonelistRegNames")
DBMMADJNodelists = PYB11attr("DBMMADJNodelists")
DBMMADJZonelists = PYB11attr("DBMMADJZonelists")
DBPMGlobNodeNo = PYB11attr("DBPMGlobNodeNo")

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

DB_CLOBBER = PYB11attr("DB_CLOBBER")
DB_NOCLOBBER = PYB11attr("DB_NOCLOBBER")

DB_READ = PYB11attr("DB_READ")
DB_APPEND = PYB11attr("DB_APPEND")

DB_LOCAL = PYB11attr("DB_LOCAL")
DB_SUN3 = PYB11attr("DB_SUN3")
DB_SUN4 = PYB11attr("DB_SUN4")
DB_SGI = PYB11attr("DB_SGI")
DB_RS6000 = PYB11attr("DB_RS6000")
DB_CRAY = PYB11attr("DB_CRAY")
DB_INTEL = PYB11attr("DB_INTEL")

DBOPT_FIRST = PYB11attr("DBOPT_FIRST")
DBOPT_ALIGN = PYB11attr("DBOPT_ALIGN")
DBOPT_COORDSYS = PYB11attr("DBOPT_COORDSYS")
DBOPT_CYCLE = PYB11attr("DBOPT_CYCLE")
DBOPT_FACETYPE = PYB11attr("DBOPT_FACETYPE")
DBOPT_HI_OFFSET = PYB11attr("DBOPT_HI_OFFSET")
DBOPT_LO_OFFSET = PYB11attr("DBOPT_LO_OFFSET")
DBOPT_LABEL = PYB11attr("DBOPT_LABEL")
DBOPT_XLABEL = PYB11attr("DBOPT_XLABEL")
DBOPT_YLABEL = PYB11attr("DBOPT_YLABEL")
DBOPT_ZLABEL = PYB11attr("DBOPT_ZLABEL")
DBOPT_MAJORORDER = PYB11attr("DBOPT_MAJORORDER")
DBOPT_NSPACE = PYB11attr("DBOPT_NSPACE")
DBOPT_ORIGIN = PYB11attr("DBOPT_ORIGIN")
DBOPT_PLANAR = PYB11attr("DBOPT_PLANAR")
DBOPT_TIME = PYB11attr("DBOPT_TIME")
DBOPT_UNITS = PYB11attr("DBOPT_UNITS")
DBOPT_XUNITS = PYB11attr("DBOPT_XUNITS")
DBOPT_YUNITS = PYB11attr("DBOPT_YUNITS")
DBOPT_ZUNITS = PYB11attr("DBOPT_ZUNITS")
DBOPT_DTIME = PYB11attr("DBOPT_DTIME")
DBOPT_USESPECMF = PYB11attr("DBOPT_USESPECMF")
DBOPT_XVARNAME = PYB11attr("DBOPT_XVARNAME")
DBOPT_YVARNAME = PYB11attr("DBOPT_YVARNAME")
DBOPT_ZVARNAME = PYB11attr("DBOPT_ZVARNAME")
DBOPT_ASCII_LABEL = PYB11attr("DBOPT_ASCII_LABEL")
DBOPT_MATNOS = PYB11attr("DBOPT_MATNOS")
DBOPT_NMATNOS = PYB11attr("DBOPT_NMATNOS")
DBOPT_MATNAME = PYB11attr("DBOPT_MATNAME")
DBOPT_NMAT = PYB11attr("DBOPT_NMAT")
DBOPT_NMATSPEC = PYB11attr("DBOPT_NMATSPEC")
DBOPT_BASEINDEX = PYB11attr("DBOPT_BASEINDEX")
DBOPT_ZONENUM = PYB11attr("DBOPT_ZONENUM")
DBOPT_NODENUM = PYB11attr("DBOPT_NODENUM")
DBOPT_BLOCKORIGIN = PYB11attr("DBOPT_BLOCKORIGIN")
DBOPT_GROUPNUM = PYB11attr("DBOPT_GROUPNUM")
DBOPT_GROUPORIGIN = PYB11attr("DBOPT_GROUPORIGIN")
DBOPT_NGROUPS = PYB11attr("DBOPT_NGROUPS")
DBOPT_MATNAMES = PYB11attr("DBOPT_MATNAMES")
DBOPT_EXTENTS_SIZE = PYB11attr("DBOPT_EXTENTS_SIZE")
DBOPT_EXTENTS = PYB11attr("DBOPT_EXTENTS")
DBOPT_MATCOUNTS = PYB11attr("DBOPT_MATCOUNTS")
DBOPT_MATLISTS = PYB11attr("DBOPT_MATLISTS")
DBOPT_MIXLENS = PYB11attr("DBOPT_MIXLENS")
DBOPT_ZONECOUNTS = PYB11attr("DBOPT_ZONECOUNTS")
DBOPT_HAS_EXTERNAL_ZONES = PYB11attr("DBOPT_HAS_EXTERNAL_ZONES")
DBOPT_PHZONELIST = PYB11attr("DBOPT_PHZONELIST")
DBOPT_MATCOLORS = PYB11attr("DBOPT_MATCOLORS")
DBOPT_BNDNAMES = PYB11attr("DBOPT_BNDNAMES")
DBOPT_REGNAMES = PYB11attr("DBOPT_REGNAMES")
DBOPT_ZONENAMES = PYB11attr("DBOPT_ZONENAMES")
DBOPT_HIDE_FROM_GUI = PYB11attr("DBOPT_HIDE_FROM_GUI")
DBOPT_TOPO_DIM = PYB11attr("DBOPT_TOPO_DIM")
DBOPT_REFERENCE = PYB11attr("DBOPT_REFERENCE")
DBOPT_GROUPINGS_SIZE = PYB11attr("DBOPT_GROUPINGS_SIZE")
DBOPT_GROUPINGS = PYB11attr("DBOPT_GROUPINGS")
DBOPT_GROUPINGNAMES = PYB11attr("DBOPT_GROUPINGNAMES")
DBOPT_ALLOWMAT0 = PYB11attr("DBOPT_ALLOWMAT0")
DBOPT_MRGTREE_NAME = PYB11attr("DBOPT_MRGTREE_NAME")
DBOPT_REGION_PNAMES = PYB11attr("DBOPT_REGION_PNAMES")
DBOPT_TENSOR_RANK = PYB11attr("DBOPT_TENSOR_RANK")
DBOPT_MMESH_NAME = PYB11attr("DBOPT_MMESH_NAME")
DBOPT_TV_CONNECTIVITY = PYB11attr("DBOPT_TV_CONNECTIVITY")
DBOPT_DISJOINT_MODE = PYB11attr("DBOPT_DISJOINT_MODE")
DBOPT_MRGV_ONAMES = PYB11attr("DBOPT_MRGV_ONAMES")
DBOPT_MRGV_RNAMES = PYB11attr("DBOPT_MRGV_RNAMES")
DBOPT_SPECNAMES = PYB11attr("DBOPT_SPECNAMES")
DBOPT_SPECCOLORS = PYB11attr("DBOPT_SPECCOLORS")
DBOPT_LLONGNZNUM = PYB11attr("DBOPT_LLONGNZNUM")
DBOPT_CONSERVED = PYB11attr("DBOPT_CONSERVED")
DBOPT_EXTENSIVE = PYB11attr("DBOPT_EXTENSIVE")
DBOPT_MB_FILE_NS = PYB11attr("DBOPT_MB_FILE_NS")
DBOPT_MB_BLOCK_NS = PYB11attr("DBOPT_MB_BLOCK_NS")
DBOPT_MB_BLOCK_TYPE = PYB11attr("DBOPT_MB_BLOCK_TYPE")
DBOPT_MB_EMPTY_LIST = PYB11attr("DBOPT_MB_EMPTY_LIST")
DBOPT_MB_EMPTY_COUNT = PYB11attr("DBOPT_MB_EMPTY_COUNT")
DBOPT_LAST = PYB11attr("DBOPT_LAST")
  
DBOPT_H5_FIRST = PYB11attr("DBOPT_H5_FIRST")
DBOPT_H5_VFD = PYB11attr("DBOPT_H5_VFD")
DBOPT_H5_RAW_FILE_OPTS = PYB11attr("DBOPT_H5_RAW_FILE_OPTS")
DBOPT_H5_RAW_EXTENSION = PYB11attr("DBOPT_H5_RAW_EXTENSION")
DBOPT_H5_META_FILE_OPTS = PYB11attr("DBOPT_H5_META_FILE_OPTS")
DBOPT_H5_META_EXTENSION = PYB11attr("DBOPT_H5_META_EXTENSION")
DBOPT_H5_CORE_ALLOC_INC = PYB11attr("DBOPT_H5_CORE_ALLOC_INC")
DBOPT_H5_CORE_NO_BACK_STORE = PYB11attr("DBOPT_H5_CORE_NO_BACK_STORE")
DBOPT_H5_META_BLOCK_SIZE = PYB11attr("DBOPT_H5_META_BLOCK_SIZE")
DBOPT_H5_SMALL_RAW_SIZE = PYB11attr("DBOPT_H5_SMALL_RAW_SIZE")
DBOPT_H5_ALIGN_MIN = PYB11attr("DBOPT_H5_ALIGN_MIN")
DBOPT_H5_ALIGN_VAL = PYB11attr("DBOPT_H5_ALIGN_VAL")
DBOPT_H5_DIRECT_MEM_ALIGN = PYB11attr("DBOPT_H5_DIRECT_MEM_ALIGN")
DBOPT_H5_DIRECT_BLOCK_SIZE = PYB11attr("DBOPT_H5_DIRECT_BLOCK_SIZE")
DBOPT_H5_DIRECT_BUF_SIZE = PYB11attr("DBOPT_H5_DIRECT_BUF_SIZE")
DBOPT_H5_LOG_NAME = PYB11attr("DBOPT_H5_LOG_NAME")
DBOPT_H5_LOG_BUF_SIZE = PYB11attr("DBOPT_H5_LOG_BUF_SIZE")
DBOPT_H5_MPIO_COMM = PYB11attr("DBOPT_H5_MPIO_COMM")
DBOPT_H5_MPIO_INFO = PYB11attr("DBOPT_H5_MPIO_INFO")
DBOPT_H5_MPIP_NO_GPFS_HINTS = PYB11attr("DBOPT_H5_MPIP_NO_GPFS_HINTS")
DBOPT_H5_SIEVE_BUF_SIZE = PYB11attr("DBOPT_H5_SIEVE_BUF_SIZE")
DBOPT_H5_CACHE_NELMTS = PYB11attr("DBOPT_H5_CACHE_NELMTS")
DBOPT_H5_CACHE_NBYTES = PYB11attr("DBOPT_H5_CACHE_NBYTES")
DBOPT_H5_CACHE_POLICY = PYB11attr("DBOPT_H5_CACHE_POLICY")
DBOPT_H5_FAM_SIZE = PYB11attr("DBOPT_H5_FAM_SIZE")
DBOPT_H5_FAM_FILE_OPTS = PYB11attr("DBOPT_H5_FAM_FILE_OPTS")
DBOPT_H5_USER_DRIVER_ID = PYB11attr("DBOPT_H5_USER_DRIVER_ID")
DBOPT_H5_USER_DRIVER_INFO = PYB11attr("DBOPT_H5_USER_DRIVER_INFO")
DBOPT_H5_SILO_BLOCK_SIZE = PYB11attr("DBOPT_H5_SILO_BLOCK_SIZE")
DBOPT_H5_SILO_BLOCK_COUNT = PYB11attr("DBOPT_H5_SILO_BLOCK_COUNT")
DBOPT_H5_SILO_LOG_STATS = PYB11attr("DBOPT_H5_SILO_LOG_STATS")
DBOPT_H5_SILO_USE_DIRECT = PYB11attr("DBOPT_H5_SILO_USE_DIRECT")
DBOPT_H5_LAST = PYB11attr("DBOPT_H5_LAST")

DB_TOP = PYB11attr("DB_TOP")
DB_NONE = PYB11attr("DB_NONE")
DB_ALL = PYB11attr("DB_ALL")
DB_ABORT = PYB11attr("DB_ABORT")
DB_SUSPEND = PYB11attr("DB_SUSPEND")
DB_RESUME = PYB11attr("DB_RESUME")
DB_ALL_AND_DRVR = PYB11attr("DB_ALL_AND_DRVR")

DB_ROWMAJOR = PYB11attr("DB_ROWMAJOR")
DB_COLMAJOR = PYB11attr("DB_COLMAJOR")

DB_NOTCENT = PYB11attr("DB_NOTCENT")
DB_NODECENT = PYB11attr("DB_NODECENT")
DB_ZONECENT = PYB11attr("DB_ZONECENT")
DB_FACECENT = PYB11attr("DB_FACECENT")
DB_BNDCENT = PYB11attr("DB_BNDCENT")
DB_EDGECENT = PYB11attr("DB_EDGECENT")
DB_BLOCKCENT = PYB11attr("DB_BLOCKCENT")

DB_CARTESIAN = PYB11attr("DB_CARTESIAN")
DB_CYLINDRICAL = PYB11attr("DB_CYLINDRICAL")
DB_SPHERICAL = PYB11attr("DB_SPHERICAL")
DB_NUMERICAL = PYB11attr("DB_NUMERICAL")
DB_OTHER = PYB11attr("DB_OTHER")

DB_RECTILINEAR = PYB11attr("DB_RECTILINEAR")
DB_CURVILINEAR = PYB11attr("DB_CURVILINEAR")

DB_AREA = PYB11attr("DB_AREA")
DB_VOLUME = PYB11attr("DB_VOLUME")

DB_ON = PYB11attr("DB_ON")
DB_OFF = PYB11attr("DB_OFF")

DB_ABUTTING = PYB11attr("DB_ABUTTING")
DB_FLOATING = PYB11attr("DB_FLOATING")

DB_VARTYPE_SCALAR = PYB11attr("DB_VARTYPE_SCALAR")
DB_VARTYPE_VECTOR = PYB11attr("DB_VARTYPE_VECTOR")
DB_VARTYPE_TENSOR = PYB11attr("DB_VARTYPE_TENSOR")
DB_VARTYPE_SYMTENSOR = PYB11attr("DB_VARTYPE_SYMTENSOR")
DB_VARTYPE_ARRAY = PYB11attr("DB_VARTYPE_ARRAY")
DB_VARTYPE_MATERIAL = PYB11attr("DB_VARTYPE_MATERIAL")
DB_VARTYPE_SPECIES = PYB11attr("DB_VARTYPE_SPECIES")
DB_VARTYPE_LABEL = PYB11attr("DB_VARTYPE_LABEL")

DBCSG_QUADRIC_G = PYB11attr("DBCSG_QUADRIC_G")
DBCSG_SPHERE_PR = PYB11attr("DBCSG_SPHERE_PR")
DBCSG_ELLIPSOID_PRRR = PYB11attr("DBCSG_ELLIPSOID_PRRR")
DBCSG_PLANE_G = PYB11attr("DBCSG_PLANE_G")
DBCSG_PLANE_X = PYB11attr("DBCSG_PLANE_X")
DBCSG_PLANE_Y = PYB11attr("DBCSG_PLANE_Y")
DBCSG_PLANE_Z = PYB11attr("DBCSG_PLANE_Z")
DBCSG_PLANE_PN = PYB11attr("DBCSG_PLANE_PN")
DBCSG_PLANE_PPP = PYB11attr("DBCSG_PLANE_PPP")
DBCSG_CYLINDER_PNLR = PYB11attr("DBCSG_CYLINDER_PNLR")
DBCSG_CYLINDER_PPR = PYB11attr("DBCSG_CYLINDER_PPR")
DBCSG_BOX_XYZXYZ = PYB11attr("DBCSG_BOX_XYZXYZ")
DBCSG_CONE_PNLA = PYB11attr("DBCSG_CONE_PNLA")
DBCSG_CONE_PPA = PYB11attr("DBCSG_CONE_PPA")
DBCSG_POLYHEDRON_KF = PYB11attr("DBCSG_POLYHEDRON_KF")
DBCSG_HEX_6F = PYB11attr("DBCSG_HEX_6F")
DBCSG_TET_4F = PYB11attr("DBCSG_TET_4F")
DBCSG_PYRAMID_5F = PYB11attr("DBCSG_PYRAMID_5F")
DBCSG_PRISM_5F = PYB11attr("DBCSG_PRISM_5F")

DBCSG_QUADRATIC_G = PYB11attr("DBCSG_QUADRATIC_G")
DBCSG_CIRCLE_PR = PYB11attr("DBCSG_CIRCLE_PR")
DBCSG_ELLIPSE_PRR = PYB11attr("DBCSG_ELLIPSE_PRR")
DBCSG_LINE_G = PYB11attr("DBCSG_LINE_G")
DBCSG_LINE_X = PYB11attr("DBCSG_LINE_X")
DBCSG_LINE_Y = PYB11attr("DBCSG_LINE_Y")
DBCSG_LINE_PN = PYB11attr("DBCSG_LINE_PN")
DBCSG_LINE_PP = PYB11attr("DBCSG_LINE_PP")
DBCSG_BOX_XYXY = PYB11attr("DBCSG_BOX_XYXY")
DBCSG_ANGLE_PNLA = PYB11attr("DBCSG_ANGLE_PNLA")
DBCSG_ANGLE_PPA = PYB11attr("DBCSG_ANGLE_PPA")
DBCSG_POLYGON_KP = PYB11attr("DBCSG_POLYGON_KP")
DBCSG_TRI_3P = PYB11attr("DBCSG_TRI_3P")
DBCSG_QUAD_4P = PYB11attr("DBCSG_QUAD_4P")

DBCSG_INNER = PYB11attr("DBCSG_INNER")
DBCSG_OUTER = PYB11attr("DBCSG_OUTER")
DBCSG_ON = PYB11attr("DBCSG_ON")
DBCSG_UNION = PYB11attr("DBCSG_UNION")
DBCSG_INTERSECT = PYB11attr("DBCSG_INTERSECT")
DBCSG_DIFF = PYB11attr("DBCSG_DIFF")
DBCSG_COMPLIMENT = PYB11attr("DBCSG_COMPLIMENT")
DBCSG_XFORM = PYB11attr("DBCSG_XFORM")
DBCSG_SWEEP = PYB11attr("DBCSG_SWEEP")

DB_PREORDER = PYB11attr("DB_PREORDER")
DB_POSTORDER = PYB11attr("DB_POSTORDER")
DB_FROMCWR = PYB11attr("DB_FROMCWR")

DB_ZONETYPE_BEAM = PYB11attr("DB_ZONETYPE_BEAM")

DB_ZONETYPE_POLYGON = PYB11attr("DB_ZONETYPE_POLYGON")
DB_ZONETYPE_TRIANGLE = PYB11attr("DB_ZONETYPE_TRIANGLE")
DB_ZONETYPE_QUAD = PYB11attr("DB_ZONETYPE_QUAD")
