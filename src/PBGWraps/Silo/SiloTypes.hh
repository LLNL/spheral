#ifndef __PBGWRAPS_SILOTYPES__
#define __PBGWRAPS_SILOTYPES__

#include <stdio.h>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include "boost/shared_ptr.hpp"

#include "silo.h"

#include "Geometry/Dimension.hh"
#include "Utilities/DBC.hh"
#include "PBGWraps/CXXTypes/CXXTypes.hh"

namespace silo {

//------------------------------------------------------------------------------
// Trait class to help mapping Spheral types to silo.
//------------------------------------------------------------------------------
template<typename T>
struct Spheral2Silo;

template<>
struct Spheral2Silo<double> {
  typedef double Value;
  static unsigned numElements() { return 1; }
  static void copyElement(const Value& x, double** silovars, const unsigned i) {
    silovars[0][i] = x;
  }
};

template<>
struct Spheral2Silo<Spheral::Dim<2>::Vector> {
  typedef Spheral::Dim<2>::Vector Value;
  static unsigned numElements() { return 2; }
  static void copyElement(const Value& x, double** silovars, const unsigned i) {
    silovars[0][i] = x.x();
    silovars[1][i] = x.y();
  }
};

template<>
struct Spheral2Silo<Spheral::Dim<2>::Tensor> {
  typedef Spheral::Dim<2>::Tensor Value;
  static unsigned numElements() { return 4; }
  static void copyElement(const Value& x, double** silovars, const unsigned i) {
    silovars[0][i] = x.xx(); silovars[1][i] = x.xy();
    silovars[2][i] = x.yx(); silovars[3][i] = x.yy();
  }
};

template<>
struct Spheral2Silo<Spheral::Dim<2>::SymTensor> {
  typedef Spheral::Dim<2>::SymTensor Value;
  static unsigned numElements() { return 4; }
  static void copyElement(const Value& x, double** silovars, const unsigned i) {
    silovars[0][i] = x.xx(); silovars[1][i] = x.xy();
    silovars[2][i] = x.yx(); silovars[3][i] = x.yy();
  }
};

template<>
struct Spheral2Silo<Spheral::Dim<3>::Vector> {
  typedef Spheral::Dim<3>::Vector Value;
  static unsigned numElements() { return 3; }
  static void copyElement(const Value& x, double** silovars, const unsigned i) {
    silovars[0][i] = x.x();
    silovars[1][i] = x.y();
    silovars[2][i] = x.z();
  }
};

template<>
struct Spheral2Silo<Spheral::Dim<3>::Tensor> {
  typedef Spheral::Dim<3>::Tensor Value;
  static unsigned numElements() { return 9; }
  static void copyElement(const Value& x, double** silovars, const unsigned i) {
    silovars[0][i] = x.xx(); silovars[1][i] = x.xy(); silovars[2][i] = x.xz();
    silovars[3][i] = x.yx(); silovars[4][i] = x.yy(); silovars[5][i] = x.yz();
    silovars[6][i] = x.zx(); silovars[7][i] = x.zy(); silovars[8][i] = x.zz();
  }
};

template<>
struct Spheral2Silo<Spheral::Dim<3>::SymTensor> {
  typedef Spheral::Dim<3>::SymTensor Value;
  static unsigned numElements() { return 9; }
  static void copyElement(const Value& x, double** silovars, const unsigned i) {
    silovars[0][i] = x.xx(); silovars[1][i] = x.xy(); silovars[2][i] = x.xz();
    silovars[3][i] = x.yx(); silovars[4][i] = x.yy(); silovars[5][i] = x.yz();
    silovars[6][i] = x.zx(); silovars[7][i] = x.zy(); silovars[8][i] = x.zz();
  }
};

//------------------------------------------------------------------------------
// Names
//------------------------------------------------------------------------------
typedef ::DBfile DBfile;
typedef ::DBoptlist DBoptlist;
using ::DBMakeOptlist;
using ::DBClearOptlist;
using ::DBClearOption;
using std::string;

//------------------------------------------------------------------------------
// Convert a std::string -> char*
//------------------------------------------------------------------------------
struct ConvertStringToCharStar {
  char* operator()(const std::string& x) {
    return const_cast<char*>(x.c_str());
  }
};

//------------------------------------------------------------------------------
// An struct to help exposing the many silo attributes.
//------------------------------------------------------------------------------
enum SiloAttributes {
  _DB_ZONETYPE_POLYHEDRON = DB_ZONETYPE_POLYHEDRON,
  _DB_ZONETYPE_TET = DB_ZONETYPE_TET,
  _DB_ZONETYPE_PYRAMID = DB_ZONETYPE_PYRAMID,
  _DB_ZONETYPE_PRISM = DB_ZONETYPE_PRISM,
  _DB_ZONETYPE_HEX = DB_ZONETYPE_HEX,

  _DB_NETCDF = DB_NETCDF,
  _DB_PDB = DB_PDB,
  _DB_TAURUS = DB_TAURUS,
  _DB_UNKNOWN = DB_UNKNOWN,
  _DB_DEBUG = DB_DEBUG,
  _DB_HDF5X = DB_HDF5X,
  _DB_PDBP = DB_PDBP,

  _DB_HDF5_SEC2_OBSOLETE = DB_HDF5_SEC2_OBSOLETE,
  _DB_HDF5_STDIO_OBSOLETE = DB_HDF5_STDIO_OBSOLETE,
  _DB_HDF5_CORE_OBSOLETE = DB_HDF5_CORE_OBSOLETE,
  _DB_HDF5_MPIO_OBSOLETE = DB_HDF5_MPIO_OBSOLETE,
  _DB_HDF5_MPIOP_OBSOLETE = DB_HDF5_MPIOP_OBSOLETE,

  _DB_H5VFD_DEFAULT = DB_H5VFD_DEFAULT,
  _DB_H5VFD_SEC2 = DB_H5VFD_SEC2,
  _DB_H5VFD_STDIO = DB_H5VFD_STDIO,
  _DB_H5VFD_CORE = DB_H5VFD_CORE,
  _DB_H5VFD_LOG = DB_H5VFD_LOG,
  _DB_H5VFD_SPLIT = DB_H5VFD_SPLIT,
  _DB_H5VFD_DIRECT = DB_H5VFD_DIRECT,
  _DB_H5VFD_FAMILY = DB_H5VFD_FAMILY,
  _DB_H5VFD_MPIO = DB_H5VFD_MPIO,
  _DB_H5VFD_MPIP = DB_H5VFD_MPIP,
  _DB_H5VFD_SILO = DB_H5VFD_SILO,

  _DB_FILE_OPTS_H5_DEFAULT_DEFAULT = DB_FILE_OPTS_H5_DEFAULT_DEFAULT,
  _DB_FILE_OPTS_H5_DEFAULT_SEC2    = DB_FILE_OPTS_H5_DEFAULT_SEC2,
  _DB_FILE_OPTS_H5_DEFAULT_STDIO   = DB_FILE_OPTS_H5_DEFAULT_STDIO,
  _DB_FILE_OPTS_H5_DEFAULT_CORE    = DB_FILE_OPTS_H5_DEFAULT_CORE,
  _DB_FILE_OPTS_H5_DEFAULT_LOG     = DB_FILE_OPTS_H5_DEFAULT_LOG,
  _DB_FILE_OPTS_H5_DEFAULT_SPLIT   = DB_FILE_OPTS_H5_DEFAULT_SPLIT,
  _DB_FILE_OPTS_H5_DEFAULT_DIRECT  = DB_FILE_OPTS_H5_DEFAULT_DIRECT,
  _DB_FILE_OPTS_H5_DEFAULT_FAMILY  = DB_FILE_OPTS_H5_DEFAULT_FAMILY,
  _DB_FILE_OPTS_H5_DEFAULT_MPIO    = DB_FILE_OPTS_H5_DEFAULT_MPIO,
  _DB_FILE_OPTS_H5_DEFAULT_MPIP    = DB_FILE_OPTS_H5_DEFAULT_MPIP,
  _DB_FILE_OPTS_H5_DEFAULT_SILO    = DB_FILE_OPTS_H5_DEFAULT_SILO,
  _DB_FILE_OPTS_LAST               = DB_FILE_OPTS_H5_DEFAULT_SILO,

  _DB_HDF5 = DB_HDF5,
  _DB_HDF5_SEC2 = DB_HDF5_SEC2,
  _DB_HDF5_STDIO = DB_HDF5_STDIO,
  _DB_HDF5_CORE = DB_HDF5_CORE,
  _DB_HDF5_LOG = DB_HDF5_LOG,
  _DB_HDF5_SPLIT = DB_HDF5_SPLIT,
  _DB_HDF5_DIRECT = DB_HDF5_DIRECT,
  _DB_HDF5_FAMILY = DB_HDF5_FAMILY,
  _DB_HDF5_MPIO = DB_HDF5_MPIO,
  _DB_HDF5_MPIOP = DB_HDF5_MPIOP,
  _DB_HDF5_MPIP = DB_HDF5_MPIP,
  _DB_HDF5_SILO = DB_HDF5_SILO,

  _DB_NFILES = DB_NFILES,
  _DB_NFILTERS = DB_NFILTERS,

  _DBAll = DBAll,
  _DBNone = DBNone,
  _DBCalc = DBCalc,
  _DBMatMatnos = DBMatMatnos,
  _DBMatMatlist = DBMatMatlist,
  _DBMatMixList = DBMatMixList,
  _DBCurveArrays = DBCurveArrays,
  _DBPMCoords = DBPMCoords,
  _DBPVData = DBPVData,
  _DBQMCoords = DBQMCoords,
  _DBQVData = DBQVData,
  _DBUMCoords = DBUMCoords,
  _DBUMFacelist = DBUMFacelist,
  _DBUMZonelist = DBUMZonelist,
  _DBUVData = DBUVData,
  _DBFacelistInfo = DBFacelistInfo,
  _DBZonelistInfo = DBZonelistInfo,
  _DBMatMatnames = DBMatMatnames,
  _DBUMGlobNodeNo = DBUMGlobNodeNo,
  _DBZonelistGlobZoneNo = DBZonelistGlobZoneNo,
  _DBMatMatcolors = DBMatMatcolors,
  _DBCSGMBoundaryInfo = DBCSGMBoundaryInfo,
  _DBCSGMZonelist = DBCSGMZonelist,
  _DBCSGMBoundaryNames = DBCSGMBoundaryNames,
  _DBCSGVData = DBCSGVData,
  _DBCSGZonelistZoneNames = DBCSGZonelistZoneNames,
  _DBCSGZonelistRegNames = DBCSGZonelistRegNames,
  _DBMMADJNodelists = DBMMADJNodelists,
  _DBMMADJZonelists = DBMMADJZonelists,
  _DBPMGlobNodeNo = DBPMGlobNodeNo,

  _DB_INVALID_OBJECT = DB_INVALID_OBJECT,
  _DB_QUADRECT = DB_QUADRECT,
  _DB_QUADCURV = DB_QUADCURV,
  _DB_QUADMESH = DB_QUADMESH,
  _DB_QUADVAR = DB_QUADVAR,
  _DB_UCDMESH = DB_UCDMESH,
  _DB_UCDVAR = DB_UCDVAR,
  _DB_MULTIMESH = DB_MULTIMESH,
  _DB_MULTIVAR = DB_MULTIVAR,
  _DB_MULTIMAT = DB_MULTIMAT,
  _DB_MULTIMATSPECIES = DB_MULTIMATSPECIES,
  _DB_MULTIBLOCKMESH = DB_MULTIBLOCKMESH,
  _DB_MULTIBLOCKVAR = DB_MULTIBLOCKVAR,
  _DB_MULTIMESHADJ = DB_MULTIMESHADJ,
  _DB_MATERIAL = DB_MATERIAL,
  _DB_MATSPECIES = DB_MATSPECIES,
  _DB_FACELIST = DB_FACELIST,
  _DB_ZONELIST = DB_ZONELIST,
  _DB_EDGELIST = DB_EDGELIST,
  _DB_PHZONELIST = DB_PHZONELIST,
  _DB_CSGZONELIST = DB_CSGZONELIST,
  _DB_CSGMESH = DB_CSGMESH,
  _DB_CSGVAR = DB_CSGVAR,
  _DB_CURVE = DB_CURVE,
  _DB_DEFVARS = DB_DEFVARS,
  _DB_POINTMESH = DB_POINTMESH,
  _DB_POINTVAR = DB_POINTVAR,
  _DB_ARRAY = DB_ARRAY,
  _DB_DIR = DB_DIR,
  _DB_VARIABLE = DB_VARIABLE,
  _DB_MRGTREE = DB_MRGTREE,
  _DB_GROUPELMAP = DB_GROUPELMAP,
  _DB_MRGVAR = DB_MRGVAR,
  _DB_USERDEF = DB_USERDEF,

  _DB_INT = DB_INT,
  _DB_SHORT = DB_SHORT,
  _DB_LONG = DB_LONG,
  _DB_FLOAT = DB_FLOAT,
  _DB_DOUBLE = DB_DOUBLE,
  _DB_CHAR = DB_CHAR,
  _DB_LONG_LONG = DB_LONG_LONG,
  _DB_NOTYPE = DB_NOTYPE,

  _DB_CLOBBER = DB_CLOBBER,
  _DB_NOCLOBBER = DB_NOCLOBBER,

  _DB_READ = DB_READ,
  _DB_APPEND = DB_APPEND,

  _DB_LOCAL = DB_LOCAL,
  _DB_SUN3 = DB_SUN3,
  _DB_SUN4 = DB_SUN4,
  _DB_SGI = DB_SGI,
  _DB_RS6000 = DB_RS6000,
  _DB_CRAY = DB_CRAY,
  _DB_INTEL = DB_INTEL,

  _DBOPT_FIRST = DBOPT_FIRST,
  _DBOPT_ALIGN = DBOPT_ALIGN,
  _DBOPT_COORDSYS = DBOPT_COORDSYS,
  _DBOPT_CYCLE = DBOPT_CYCLE,
  _DBOPT_FACETYPE = DBOPT_FACETYPE,
  _DBOPT_HI_OFFSET = DBOPT_HI_OFFSET,
  _DBOPT_LO_OFFSET = DBOPT_LO_OFFSET,
  _DBOPT_LABEL = DBOPT_LABEL,
  _DBOPT_XLABEL = DBOPT_XLABEL,
  _DBOPT_YLABEL = DBOPT_YLABEL,
  _DBOPT_ZLABEL = DBOPT_ZLABEL,
  _DBOPT_MAJORORDER = DBOPT_MAJORORDER,
  _DBOPT_NSPACE = DBOPT_NSPACE,
  _DBOPT_ORIGIN = DBOPT_ORIGIN,
  _DBOPT_PLANAR = DBOPT_PLANAR,
  _DBOPT_TIME = DBOPT_TIME,
  _DBOPT_UNITS = DBOPT_UNITS,
  _DBOPT_XUNITS = DBOPT_XUNITS,
  _DBOPT_YUNITS = DBOPT_YUNITS,
  _DBOPT_ZUNITS = DBOPT_ZUNITS,
  _DBOPT_DTIME = DBOPT_DTIME,
  _DBOPT_USESPECMF = DBOPT_USESPECMF,
  _DBOPT_XVARNAME = DBOPT_XVARNAME,
  _DBOPT_YVARNAME = DBOPT_YVARNAME,
  _DBOPT_ZVARNAME = DBOPT_ZVARNAME,
  _DBOPT_ASCII_LABEL = DBOPT_ASCII_LABEL,
  _DBOPT_MATNOS = DBOPT_MATNOS,
  _DBOPT_NMATNOS = DBOPT_NMATNOS,
  _DBOPT_MATNAME = DBOPT_MATNAME,
  _DBOPT_NMAT = DBOPT_NMAT,
  _DBOPT_NMATSPEC = DBOPT_NMATSPEC,
  _DBOPT_BASEINDEX = DBOPT_BASEINDEX,
  _DBOPT_ZONENUM = DBOPT_ZONENUM,
  _DBOPT_NODENUM = DBOPT_NODENUM,
  _DBOPT_BLOCKORIGIN = DBOPT_BLOCKORIGIN,
  _DBOPT_GROUPNUM = DBOPT_GROUPNUM,
  _DBOPT_GROUPORIGIN = DBOPT_GROUPORIGIN,
  _DBOPT_NGROUPS = DBOPT_NGROUPS,
  _DBOPT_MATNAMES = DBOPT_MATNAMES,
  _DBOPT_EXTENTS_SIZE = DBOPT_EXTENTS_SIZE,
  _DBOPT_EXTENTS = DBOPT_EXTENTS,
  _DBOPT_MATCOUNTS = DBOPT_MATCOUNTS,
  _DBOPT_MATLISTS = DBOPT_MATLISTS,
  _DBOPT_MIXLENS = DBOPT_MIXLENS,
  _DBOPT_ZONECOUNTS = DBOPT_ZONECOUNTS,
  _DBOPT_HAS_EXTERNAL_ZONES = DBOPT_HAS_EXTERNAL_ZONES,
  _DBOPT_PHZONELIST = DBOPT_PHZONELIST,
  _DBOPT_MATCOLORS = DBOPT_MATCOLORS,
  _DBOPT_BNDNAMES = DBOPT_BNDNAMES,
  _DBOPT_REGNAMES = DBOPT_REGNAMES,
  _DBOPT_ZONENAMES = DBOPT_ZONENAMES,
  _DBOPT_HIDE_FROM_GUI = DBOPT_HIDE_FROM_GUI,
  _DBOPT_TOPO_DIM = DBOPT_TOPO_DIM,
  _DBOPT_REFERENCE = DBOPT_REFERENCE,
  _DBOPT_GROUPINGS_SIZE = DBOPT_GROUPINGS_SIZE,
  _DBOPT_GROUPINGS = DBOPT_GROUPINGS,
  _DBOPT_GROUPINGNAMES = DBOPT_GROUPINGNAMES,
  _DBOPT_ALLOWMAT0 = DBOPT_ALLOWMAT0,
  _DBOPT_MRGTREE_NAME = DBOPT_MRGTREE_NAME,
  _DBOPT_REGION_PNAMES = DBOPT_REGION_PNAMES,
  _DBOPT_TENSOR_RANK = DBOPT_TENSOR_RANK,
  _DBOPT_MMESH_NAME = DBOPT_MMESH_NAME,
  _DBOPT_TV_CONNECTIVITY = DBOPT_TV_CONNECTIVITY,
  _DBOPT_DISJOINT_MODE = DBOPT_DISJOINT_MODE,
  _DBOPT_MRGV_ONAMES = DBOPT_MRGV_ONAMES,
  _DBOPT_MRGV_RNAMES = DBOPT_MRGV_RNAMES,
  _DBOPT_SPECNAMES = DBOPT_SPECNAMES,
  _DBOPT_SPECCOLORS = DBOPT_SPECCOLORS,
  _DBOPT_LLONGNZNUM = DBOPT_LLONGNZNUM,
  _DBOPT_CONSERVED = DBOPT_CONSERVED,
  _DBOPT_EXTENSIVE = DBOPT_EXTENSIVE,
  _DBOPT_MB_FILE_NS = DBOPT_MB_FILE_NS,
  _DBOPT_MB_BLOCK_NS = DBOPT_MB_BLOCK_NS,
  _DBOPT_MB_BLOCK_TYPE = DBOPT_MB_BLOCK_TYPE,
  _DBOPT_MB_EMPTY_LIST = DBOPT_MB_EMPTY_LIST,
  _DBOPT_MB_EMPTY_COUNT = DBOPT_MB_EMPTY_COUNT,
  _DBOPT_LAST = DBOPT_LAST,
  
  _DBOPT_H5_FIRST = DBOPT_H5_FIRST,
  _DBOPT_H5_VFD = DBOPT_H5_VFD,
  _DBOPT_H5_RAW_FILE_OPTS = DBOPT_H5_RAW_FILE_OPTS,
  _DBOPT_H5_RAW_EXTENSION = DBOPT_H5_RAW_EXTENSION,
  _DBOPT_H5_META_FILE_OPTS = DBOPT_H5_META_FILE_OPTS,
  _DBOPT_H5_META_EXTENSION = DBOPT_H5_META_EXTENSION,
  _DBOPT_H5_CORE_ALLOC_INC = DBOPT_H5_CORE_ALLOC_INC,
  _DBOPT_H5_CORE_NO_BACK_STORE = DBOPT_H5_CORE_NO_BACK_STORE,
  _DBOPT_H5_META_BLOCK_SIZE = DBOPT_H5_META_BLOCK_SIZE,
  _DBOPT_H5_SMALL_RAW_SIZE = DBOPT_H5_SMALL_RAW_SIZE,
  _DBOPT_H5_ALIGN_MIN = DBOPT_H5_ALIGN_MIN,
  _DBOPT_H5_ALIGN_VAL = DBOPT_H5_ALIGN_VAL,
  _DBOPT_H5_DIRECT_MEM_ALIGN = DBOPT_H5_DIRECT_MEM_ALIGN,
  _DBOPT_H5_DIRECT_BLOCK_SIZE = DBOPT_H5_DIRECT_BLOCK_SIZE,
  _DBOPT_H5_DIRECT_BUF_SIZE = DBOPT_H5_DIRECT_BUF_SIZE,
  _DBOPT_H5_LOG_NAME = DBOPT_H5_LOG_NAME,
  _DBOPT_H5_LOG_BUF_SIZE = DBOPT_H5_LOG_BUF_SIZE,
  _DBOPT_H5_MPIO_COMM = DBOPT_H5_MPIO_COMM,
  _DBOPT_H5_MPIO_INFO = DBOPT_H5_MPIO_INFO,
  _DBOPT_H5_MPIP_NO_GPFS_HINTS = DBOPT_H5_MPIP_NO_GPFS_HINTS,
  _DBOPT_H5_SIEVE_BUF_SIZE = DBOPT_H5_SIEVE_BUF_SIZE,
  _DBOPT_H5_CACHE_NELMTS = DBOPT_H5_CACHE_NELMTS,
  _DBOPT_H5_CACHE_NBYTES = DBOPT_H5_CACHE_NBYTES,
  _DBOPT_H5_CACHE_POLICY = DBOPT_H5_CACHE_POLICY,
  _DBOPT_H5_FAM_SIZE = DBOPT_H5_FAM_SIZE,
  _DBOPT_H5_FAM_FILE_OPTS = DBOPT_H5_FAM_FILE_OPTS,
  _DBOPT_H5_USER_DRIVER_ID = DBOPT_H5_USER_DRIVER_ID,
  _DBOPT_H5_USER_DRIVER_INFO = DBOPT_H5_USER_DRIVER_INFO,
  _DBOPT_H5_SILO_BLOCK_SIZE = DBOPT_H5_SILO_BLOCK_SIZE,
  _DBOPT_H5_SILO_BLOCK_COUNT = DBOPT_H5_SILO_BLOCK_COUNT,
  _DBOPT_H5_SILO_LOG_STATS = DBOPT_H5_SILO_LOG_STATS,
  _DBOPT_H5_SILO_USE_DIRECT = DBOPT_H5_SILO_USE_DIRECT,
  _DBOPT_H5_LAST = DBOPT_H5_LAST,

  _DB_TOP = DB_TOP,
  _DB_NONE = DB_NONE,
  _DB_ALL = DB_ALL,
  _DB_ABORT = DB_ABORT,
  _DB_SUSPEND = DB_SUSPEND,
  _DB_RESUME = DB_RESUME,
  _DB_ALL_AND_DRVR = DB_ALL_AND_DRVR,

  _DB_ROWMAJOR = DB_ROWMAJOR,
  _DB_COLMAJOR = DB_COLMAJOR,

  _DB_NOTCENT = DB_NOTCENT,
  _DB_NODECENT = DB_NODECENT,
  _DB_ZONECENT = DB_ZONECENT,
  _DB_FACECENT = DB_FACECENT,
  _DB_BNDCENT = DB_BNDCENT,
  _DB_EDGECENT = DB_EDGECENT,
  _DB_BLOCKCENT = DB_BLOCKCENT,

  _DB_CARTESIAN = DB_CARTESIAN,
  _DB_CYLINDRICAL = DB_CYLINDRICAL,
  _DB_SPHERICAL = DB_SPHERICAL,
  _DB_NUMERICAL = DB_NUMERICAL,
  _DB_OTHER = DB_OTHER,

  _DB_RECTILINEAR = DB_RECTILINEAR,
  _DB_CURVILINEAR = DB_CURVILINEAR,

  _DB_AREA = DB_AREA,
  _DB_VOLUME = DB_VOLUME,

  _DB_ON = DB_ON,
  _DB_OFF = DB_OFF,

  _DB_ABUTTING = DB_ABUTTING,
  _DB_FLOATING = DB_FLOATING,

  _DB_VARTYPE_SCALAR = DB_VARTYPE_SCALAR,
  _DB_VARTYPE_VECTOR = DB_VARTYPE_VECTOR,
  _DB_VARTYPE_TENSOR = DB_VARTYPE_TENSOR,
  _DB_VARTYPE_SYMTENSOR = DB_VARTYPE_SYMTENSOR,
  _DB_VARTYPE_ARRAY = DB_VARTYPE_ARRAY,
  _DB_VARTYPE_MATERIAL = DB_VARTYPE_MATERIAL,
  _DB_VARTYPE_SPECIES = DB_VARTYPE_SPECIES,
  _DB_VARTYPE_LABEL = DB_VARTYPE_LABEL,

  _DBCSG_QUADRIC_G = DBCSG_QUADRIC_G,
  _DBCSG_SPHERE_PR = DBCSG_SPHERE_PR,
  _DBCSG_ELLIPSOID_PRRR = DBCSG_ELLIPSOID_PRRR,
  _DBCSG_PLANE_G = DBCSG_PLANE_G,
  _DBCSG_PLANE_X = DBCSG_PLANE_X,
  _DBCSG_PLANE_Y = DBCSG_PLANE_Y,
  _DBCSG_PLANE_Z = DBCSG_PLANE_Z,
  _DBCSG_PLANE_PN = DBCSG_PLANE_PN,
  _DBCSG_PLANE_PPP = DBCSG_PLANE_PPP,
  _DBCSG_CYLINDER_PNLR = DBCSG_CYLINDER_PNLR,
  _DBCSG_CYLINDER_PPR = DBCSG_CYLINDER_PPR,
  _DBCSG_BOX_XYZXYZ = DBCSG_BOX_XYZXYZ,
  _DBCSG_CONE_PNLA = DBCSG_CONE_PNLA,
  _DBCSG_CONE_PPA = DBCSG_CONE_PPA,
  _DBCSG_POLYHEDRON_KF = DBCSG_POLYHEDRON_KF,
  _DBCSG_HEX_6F = DBCSG_HEX_6F,
  _DBCSG_TET_4F = DBCSG_TET_4F,
  _DBCSG_PYRAMID_5F = DBCSG_PYRAMID_5F,
  _DBCSG_PRISM_5F = DBCSG_PRISM_5F,

  _DBCSG_QUADRATIC_G = DBCSG_QUADRATIC_G,
  _DBCSG_CIRCLE_PR = DBCSG_CIRCLE_PR,
  _DBCSG_ELLIPSE_PRR = DBCSG_ELLIPSE_PRR,
  _DBCSG_LINE_G = DBCSG_LINE_G,
  _DBCSG_LINE_X = DBCSG_LINE_X,
  _DBCSG_LINE_Y = DBCSG_LINE_Y,
  _DBCSG_LINE_PN = DBCSG_LINE_PN,
  _DBCSG_LINE_PP = DBCSG_LINE_PP,
  _DBCSG_BOX_XYXY = DBCSG_BOX_XYXY,
  _DBCSG_ANGLE_PNLA = DBCSG_ANGLE_PNLA,
  _DBCSG_ANGLE_PPA = DBCSG_ANGLE_PPA,
  _DBCSG_POLYGON_KP = DBCSG_POLYGON_KP,
  _DBCSG_TRI_3P = DBCSG_TRI_3P,
  _DBCSG_QUAD_4P = DBCSG_QUAD_4P,

  _DBCSG_INNER = DBCSG_INNER,
  _DBCSG_OUTER = DBCSG_OUTER,
  _DBCSG_ON = DBCSG_ON,
  _DBCSG_UNION = DBCSG_UNION,
  _DBCSG_INTERSECT = DBCSG_INTERSECT,
  _DBCSG_DIFF = DBCSG_DIFF,
  _DBCSG_COMPLIMENT = DBCSG_COMPLIMENT,
  _DBCSG_XFORM = DBCSG_XFORM,
  _DBCSG_SWEEP = DBCSG_SWEEP,

  _DB_PREORDER = DB_PREORDER,
  _DB_POSTORDER = DB_POSTORDER,
  _DB_FROMCWR = DB_FROMCWR,

  _DB_ZONETYPE_BEAM = DB_ZONETYPE_BEAM,

  _DB_ZONETYPE_POLYGON = DB_ZONETYPE_POLYGON,
  _DB_ZONETYPE_TRIANGLE = DB_ZONETYPE_TRIANGLE,
  _DB_ZONETYPE_QUAD = DB_ZONETYPE_QUAD,
};

//------------------------------------------------------------------------------
// A trait class for for mapping types -> silo types.
//------------------------------------------------------------------------------
template<typename T> struct SiloTraits {};

template<>
struct SiloTraits<int> {
  static std::vector<int> dims() { return std::vector<int>(1, 1); }
  static int datatype() { return DB_INT; }
};

template<>
struct SiloTraits<short> {
  static std::vector<int> dims() { return std::vector<int>(1, 1); }
  static int datatype() { return DB_SHORT; }
};

template<>
struct SiloTraits<long> {
  static std::vector<int> dims() { return std::vector<int>(1, 1); }
  static int datatype() { return DB_LONG; }
};

template<>
struct SiloTraits<float> {
  static std::vector<int> dims() { return std::vector<int>(1, 1); }
  static int datatype() { return DB_FLOAT; }
};

template<>
struct SiloTraits<double> {
  static std::vector<int> dims() { return std::vector<int>(1, 1); }
  static int datatype() { return DB_DOUBLE; }
};

template<>
struct SiloTraits<char> {
  static std::vector<int> dims() { return std::vector<int>(1, 1); }
  static int datatype() { return DB_CHAR; }
};

template<>
struct SiloTraits<long long> {
  static std::vector<int> dims() { return std::vector<int>(1, 1); }
  static int datatype() { return DB_LONG_LONG; }
};

//------------------------------------------------------------------------------
// Wrapper class to handle the memory managemnt necessary with DBoptlist.
//------------------------------------------------------------------------------
struct DBoptlist_wrapper {
  DBoptlist* mOptlistPtr;
  std::vector<boost::shared_ptr<void> > mCache;

  // Constructors.
  DBoptlist_wrapper(const int maxopts=1024):
    mOptlistPtr(DBMakeOptlist(maxopts)),
    mCache() {}

  // Destructor.
  ~DBoptlist_wrapper() {
    VERIFY(DBFreeOptlist(mOptlistPtr) == 0);
  }

  // Generic functor definitions for adding and getting options.
  template<typename Value>
  struct AddOptionFunctor {
    int writeValue(DBoptlist_wrapper& optlist_wrapper,
                   const int option,
                   const Value& value) {
      boost::shared_ptr<void> voidValue(new Value(value));
      optlist_wrapper.mCache.push_back(voidValue);
      return DBAddOption(optlist_wrapper.mOptlistPtr, option, voidValue.get());
    }
    int writeVector(DBoptlist_wrapper& optlist_wrapper,
                    const int option,
                    const int option_size,
                    const std::vector<Value>& value) {
      DBoptlist_wrapper::AddOptionFunctor<int>().writeValue(optlist_wrapper, option_size, value.size());
      boost::shared_ptr<void> voidValue(new std::vector<Value>(value));
      optlist_wrapper.mCache.push_back(voidValue);
      Value* frontPtr = &(((std::vector<Value>*) voidValue.get())->front());
      return DBAddOption(optlist_wrapper.mOptlistPtr, option, frontPtr);
    }
  };

  template<typename Value>
  struct GetOptionFunctor {
    Value readValue(DBoptlist_wrapper& optlist_wrapper,
                    const int option) {
      return *((Value*) DBGetOption(optlist_wrapper.mOptlistPtr, option));
    }
    std::vector<Value> readVector(DBoptlist_wrapper& optlist_wrapper,
                                  const int option,
                                  const int option_size) {
      const unsigned vecsize = DBoptlist_wrapper::GetOptionFunctor<int>().readValue(optlist_wrapper, option_size);
      Value* frontPtr = (Value*) DBGetOption(optlist_wrapper.mOptlistPtr, option);
      return std::vector<Value>(frontPtr, frontPtr + vecsize);
    }
  };

  // Function definitions that use the functors.
  template<typename Value>
  int addOption(const int option,
                const Value& value) {
    return DBoptlist_wrapper::AddOptionFunctor<Value>().writeValue(*this, option, value);
  }

  template<typename Value>
  int addOption(const int option,
                const int option_size,
                const std::vector<Value>& value) {
    return DBoptlist_wrapper::AddOptionFunctor<Value>().writeVector(*this, option, option_size, value);
  }

  template<typename Value>
  Value getOption(const int option) {
    return DBoptlist_wrapper::GetOptionFunctor<Value>().readValue(*this, option);
  }

  template<typename Value>
  std::vector<Value> getOption(const int option,
                               const int option_size) {
    return DBoptlist_wrapper::GetOptionFunctor<Value>().readVector(*this, option, option_size);
  }
};
                              
//..............................................................................
// std::string specializations.
//..............................................................................
template<>
struct
DBoptlist_wrapper::AddOptionFunctor<std::string> {
  int
  writeValue(DBoptlist_wrapper& optlist_wrapper,
             const int option,
             const std::string& value) {
    boost::shared_ptr<void> voidValue(new std::string(value));
    optlist_wrapper.mCache.push_back(voidValue);
    return DBAddOption(optlist_wrapper.mOptlistPtr, option, (char*) ((std::string*) voidValue.get())->c_str());
  }
  int
  writeVector(DBoptlist_wrapper& optlist_wrapper,
              const int option,
              const int option_size,
              const std::vector<std::string>& value) {
    VERIFY(optlist_wrapper.addOption<int>(option_size, value.size()) == 0);
    boost::shared_ptr<void> voidCopy(new std::vector<std::string>(value));
    boost::shared_ptr<void> voidValue(new std::vector<char*>());
    std::vector<std::string>& stringvec = *((std::vector<string>*) voidCopy.get());
    std::vector<char*>& charvec = *((std::vector<char*>*) (voidValue.get()));
    for (unsigned k = 0; k != value.size(); ++k) charvec.push_back(const_cast<char*>(stringvec[k].c_str()));
    VERIFY(charvec.size() == value.size());
    optlist_wrapper.mCache.push_back(voidValue);
    return DBAddOption(optlist_wrapper.mOptlistPtr, option, (char**) &charvec.front());
  }
};

template<>
struct
DBoptlist_wrapper::GetOptionFunctor<std::string> {
  std::string
  readValue(DBoptlist_wrapper& optlist_wrapper,
            const int option) {
    char* result = (char*) DBGetOption(optlist_wrapper.mOptlistPtr, option);
    return std::string(result);
  }
  std::vector<std::string>
  readVector(DBoptlist_wrapper& optlist_wrapper,
             const int option,
             const int option_size) {
    const int resultsize = optlist_wrapper.getOption<int>(option_size);
    VERIFY(resultsize > 0);
    char** chararray = (char**) DBGetOption(optlist_wrapper.mOptlistPtr, option);
    std::vector<std::string> result(chararray, chararray + resultsize);
    // for (unsigned k = 0;  k != resultsize; ++k) result.push_back(std::string(chararray[k]));
    VERIFY(result.size() == resultsize);
    return result;
  }
};

//------------------------------------------------------------------------------
// DBCreate
//------------------------------------------------------------------------------
inline
DBfile*
DBCreate_wrap(std::string pathName,
              int mode,
              int target,
              std::string fileInfo,
              int fileType) {
  return DBCreate(pathName.c_str(), mode, target, fileInfo.c_str(), fileType);
}

//------------------------------------------------------------------------------
// DBClose
//------------------------------------------------------------------------------
inline
int
DBClose(DBfile& file) {
  return DBClose(&file);
}

//------------------------------------------------------------------------------
// DBWrite
//------------------------------------------------------------------------------
template<typename T>
inline
int
DBWrite(DBfile& file,
        std::string varname,
        T& var) {
  return DBWrite(&file, varname.c_str(), (void*) &var, 
                 &(SiloTraits<T>::dims()).front(),
                 SiloTraits<T>::dims().size(),
                 SiloTraits<T>::datatype());
}

//------------------------------------------------------------------------------
// DBReadVar
//------------------------------------------------------------------------------
template<typename T>
inline
T
DBReadVar(DBfile& file,
          std::string varname) {
  T result;
  DBReadVar(&file, varname.c_str(), (void*) &result);
  return result;
}

//------------------------------------------------------------------------------
// DBPutMultimesh
//------------------------------------------------------------------------------
inline
int
DBPutMultimesh(DBfile& file,
               std::string name,
               std::vector<std::string>& meshNames,
               std::vector<int>& meshTypes,
               DBoptlist_wrapper& optlist) {

  // Pre-conditions.
  VERIFY2(meshNames.size() == meshTypes.size(), "meshNames and meshTypes must be same length:  " << meshNames.size() << " != " << meshTypes.size());

  // Convert names to char*.
  std::vector<char*> meshNames1;
  std::transform(meshNames.begin(), meshNames.end(), std::back_inserter(meshNames1),
                 ConvertStringToCharStar());
  CHECK(meshNames1.size() == meshNames.size());

  // Do the deed.
  const int result = DBPutMultimesh(&file, 
                                    name.c_str(), 
                                    meshNames.size(),
                                    &meshNames1.front(), 
                                    &meshTypes.front(),
                                    optlist.mOptlistPtr);

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// DBPutMultimat
//------------------------------------------------------------------------------
inline
int
DBPutMultimat(DBfile& file,
              std::string name,
              std::vector<std::string>& matNames,
              DBoptlist_wrapper& optlist) {

  // Convert names to char*.
  std::vector<char*> matNames1;
  std::transform(matNames.begin(), matNames.end(), std::back_inserter(matNames1),
                 ConvertStringToCharStar());
  CHECK(matNames1.size() == matNames.size());

  // Do the deed.
  const int result = DBPutMultimat(&file, 
                                   name.c_str(), 
                                   matNames.size(),
                                   &matNames1.front(), 
                                   optlist.mOptlistPtr);

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// DBPutCompoundarray
//------------------------------------------------------------------------------
template<typename T>
inline
int
DBPutCompoundarray(DBfile& file,
                   std::string name,
                   std::vector<std::string>& elemNames,
                   std::vector<std::vector<T> >& values,
                   DBoptlist_wrapper& optlist) {

  // Preconditions.
  VERIFY(elemNames.size() == values.size());

  // Convert names to char*.
  std::vector<char*> elemNames1;
  std::transform(elemNames.begin(), elemNames.end(), std::back_inserter(elemNames1),
                 ConvertStringToCharStar());
  CHECK(elemNames1.size() == elemNames.size());

  // Read the sizes of each array.
  std::vector<int> elemLengths;
  elemLengths.reserve(values.size());
  for (unsigned k = 0; k != values.size(); ++k) elemLengths.push_back(values[k].size());
  const unsigned numValues = std::accumulate(elemLengths.begin(), elemLengths.end(), 0);

  // Flatten the values to a single arrray.
  vector<T> flatValues;
  flatValues.reserve(numValues);
  for (unsigned k = 0; k != values.size(); ++k) {
    for (unsigned j = 0; j != values[k].size(); ++j) {
      flatValues.push_back(values[k][j]);
    }
  }
  CHECK(flatValues.size() == numValues);

  // Do the deed.
  CHECK(elemNames.size() == elemLengths.size());
  const int result = DBPutCompoundarray(&file, 
                                        name.c_str(), 
                                        &elemNames1.front(),
                                        &elemLengths.front(),
                                        elemNames1.size(),
                                        &flatValues.front(),
                                        numValues,
                                        SiloTraits<T>::datatype(),
                                        optlist.mOptlistPtr);

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// DBPutMultivar
//------------------------------------------------------------------------------
inline
int
DBPutMultivar(DBfile& file,
              std::string name,
              std::vector<std::string>& varNames,
              std::vector<int>& varTypes,
              DBoptlist_wrapper& optlist) {

  // Preconditions.
  const unsigned numVars = varNames.size();
  VERIFY(varTypes.size() == numVars);

  // Convert names to char*.
  std::vector<char*> varNames1;
  std::transform(varNames.begin(), varNames.end(), std::back_inserter(varNames1),
                 ConvertStringToCharStar());
  CHECK(varNames1.size() == varNames.size());

  // Do the deed.
  const int result = DBPutMultivar(&file, 
                                   const_cast<char*>(name.c_str()), 
                                   numVars,
                                   &varNames1.front(),
                                   &varTypes.front(),
                                   optlist.mOptlistPtr);

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// DBPutMaterial
//------------------------------------------------------------------------------
inline
int
DBPutMaterial(DBfile& file,
              std::string name,
              std::string meshName,
              std::vector<int>& matnos,
              std::vector<int>& matlist,
              std::vector<int>& mix_next,
              std::vector<int>& mix_mat,
              std::vector<int>& mix_zone,
              std::vector<double>& mix_vf,
              DBoptlist_wrapper& optlist) {

  // Preconditions.
  const unsigned nmat = matnos.size();
  const unsigned numMix = mix_next.size();
  VERIFY(mix_mat.size() == numMix);
  VERIFY(mix_zone.size() == numMix);
  VERIFY(mix_vf.size() == numMix);

  // Dimensionality of the matlist.
  // For now we only support 1-D lists (I don't understand what silo does with 
  // greater dimensionality?)
  vector<int> dims(1, matlist.size());

  // Do the deed.
  const int result = DBPutMaterial(&file, 
                                   name.c_str(),
                                   meshName.c_str(),
                                   nmat,
                                   &matnos.front(),
                                   &matlist.front(),
                                   &dims.front(),
                                   dims.size(),
                                   &mix_next.front(),
                                   &mix_mat.front(),
                                   &mix_zone.front(),
                                   &mix_vf.front(),
                                   numMix,
                                   SiloTraits<double>::datatype(),
                                   optlist.mOptlistPtr);

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// DBPutUcdmesh
//------------------------------------------------------------------------------
inline
int
DBPutUcdmesh(DBfile& file,
             std::string name,
             std::vector<std::vector<double> >& coords,
             int nzones,
             std::string zonel_name,
             std::string facel_name,
             DBoptlist_wrapper& optlist) {

  // Preconditions.
  const unsigned ndims = coords.size();
  VERIFY(ndims == 2 or ndims == 3);
  const unsigned nnodes = coords[0].size();
  for (unsigned idim = 0; idim != ndims; ++idim) VERIFY(coords[idim].size() == nnodes);
  VERIFY(nzones > 0);

  // We need the C-stylish pointers to the coordinates.
  double** coordPtrs = new double*[ndims];
  for (unsigned k = 0; k != ndims; ++k) {
    coordPtrs[k] = new double[nnodes];
    std::copy(coords[k].begin(), coords[k].end(), coordPtrs[k]);
  }

  // Convert strings to char*.
  char* zonel_name1 = (zonel_name == "NULL") ? NULL : const_cast<char*>(zonel_name.c_str());
  char* facel_name1 = (facel_name == "NULL") ? NULL : const_cast<char*>(facel_name.c_str());

  // Do the deed.
  const int result = DBPutUcdmesh(&file, 
                                  name.c_str(),
                                  ndims,
                                  NULL,
                                  coordPtrs,
                                  nnodes,
                                  nzones,
                                  zonel_name1,
                                  facel_name1,
                                  SiloTraits<double>::datatype(),
                                  optlist.mOptlistPtr);

  // That's it.
  for (unsigned k = 0; k != ndims; ++k) delete[] coordPtrs[k];
  delete[] coordPtrs;
  return result;
}

//------------------------------------------------------------------------------
// DBPutDefvars
//------------------------------------------------------------------------------
inline
int
DBPutDefvars(DBfile& file,
             std::string name,
             std::vector<std::string>& varNames,
             std::vector<int>& varTypes,
             std::vector<std::string>& varDefs,
             std::vector<DBoptlist_wrapper*>& optlists) {

  // Preconditions.
  const unsigned ndefs = varNames.size();
  VERIFY(varNames.size() == ndefs and
         varTypes.size() == ndefs and
         varDefs.size() == ndefs and
         optlists.size() == ndefs);

  // Convert names to char*'s
  std::vector<char*> names, defns;
  std::transform(varNames.begin(), varNames.end(), std::back_inserter(names), ConvertStringToCharStar());
  std::transform(varDefs.begin(), varDefs.end(), std::back_inserter(defns), ConvertStringToCharStar());
  CHECK(names.size() == varNames.size());
  CHECK(defns.size() == varDefs.size());
  
  // Copy the optlists to an array of pointers.
  std::vector<DBoptlist*> optlistptrs;
  for (std::vector<DBoptlist_wrapper*>::iterator itr = optlists.begin();
       itr != optlists.end();
       ++itr) optlistptrs.push_back((*itr)->mOptlistPtr);
  CHECK(optlistptrs.size() == ndefs);

  return DBPutDefvars(&file,
                      name.c_str(),
                      ndefs,
                      &names.front(),
                      &varTypes.front(),
                      &defns.front(),
                      &optlistptrs.front());
}

//------------------------------------------------------------------------------
// DBPutUcdvar
// We assume here that the underlying element type is double.
//------------------------------------------------------------------------------
template<typename T>
inline
int
DBPutUcdvar(DBfile& file,
            std::string name,
            std::string meshName,
            std::vector<T>& values,
            std::vector<T>& mixValues,
            int centering,
            DBoptlist_wrapper& optlist) {

  // Preconditions.
  VERIFY(centering == _DB_NODECENT or
         centering == _DB_EDGECENT or
         centering == _DB_FACECENT or
         centering == _DB_ZONECENT);

  unsigned i, j;

  // Build the sub-variable names.
  int nvars = Spheral2Silo<T>::numElements();
  int nels = values.size();
  int mixlen = mixValues.size();
  vector<char*> varnames;
  for (i = 0; i != nvars; ++i) {
    varnames.push_back(const_cast<char*>((name + "_").c_str()));
    sprintf(varnames.back(), "%i", i);
  }

  // Build the sub-variables.
  double** vars = new double*[nvars];
  double** mixvars = new double*[nvars];
  for (i = 0; i != nvars; ++i) {
    vars[i] = new double[nels];
    mixvars[i] = new double[mixlen];
  }
  for (j = 0; j != nels; ++j) Spheral2Silo<T>::copyElement(values[j], vars, j);
  for (j = 0; j != mixlen; ++j) Spheral2Silo<T>::copyElement(mixValues[j], mixvars, j);

  const int result = DBPutUcdvar(&file,
                                 name.c_str(),
                                 meshName.c_str(),
                                 nvars,
                                 &varnames.front(),
                                 (void*) vars,
                                 nels,
                                 (void*) mixvars,
                                 mixlen,
                                 SiloTraits<double>::datatype(),
                                 centering,
                                 optlist.mOptlistPtr);

  // That's it.
  for (i = 0; i != nvars; ++i) {
    delete[] vars[i];
    delete[] mixvars[i];
  }
  delete[] vars, mixvars;
  return result;
}

//------------------------------------------------------------------------------
// DBPutUcdvar1
//------------------------------------------------------------------------------
template<typename T>
inline
int
DBPutUcdvar1(DBfile& file,
             std::string name,
             std::string meshName,
             std::vector<T>& values,
             std::vector<T>& mixValues,
             int centering,
             DBoptlist_wrapper& optlist) {

  // Preconditions.
  VERIFY(centering == _DB_NODECENT or
         centering == _DB_EDGECENT or
         centering == _DB_FACECENT or
         centering == _DB_ZONECENT);

  return DBPutUcdvar1(&file,
                      name.c_str(),
                      meshName.c_str(),
                      (void*) &(*values.begin()),
                      values.size(),
                      (void*) &(*mixValues.begin()),
                      mixValues.size(),
                      SiloTraits<T>::datatype(),
                      centering,
                      optlist.mOptlistPtr);
}

//------------------------------------------------------------------------------
// DBPutZonelist2
//------------------------------------------------------------------------------
inline
int
DBPutZonelist2(DBfile& file,
               std::string name,
               unsigned ndims,
               std::vector<std::vector<int> >& zoneNodes,
               unsigned low_offset,
               unsigned high_offset,
               std::vector<int>& shapetype,
               std::vector<int>& shapesize,
               std::vector<int>& shapecount,
               DBoptlist_wrapper& optlist) {

  // Preconditions.
  const unsigned nzones = zoneNodes.size();
  const unsigned nshapes = shapetype.size();
  VERIFY(shapetype.size() <= nzones);
  VERIFY(shapesize.size() <= nzones);
  VERIFY(shapecount.size() <= nzones);
  VERIFY(shapesize.size() == nshapes);
  VERIFY(shapecount.size() == nshapes);

  // Construct the flat array of zone nodes.
  vector<int> nodelist;
  for (unsigned k = 0; k != zoneNodes.size(); ++k) {
    //nodelist.push_back(zoneNodes[k].size());
    std::copy(zoneNodes[k].begin(), zoneNodes[k].end(), std::back_inserter(nodelist));
  }

  // Do the deed.
  const int result = DBPutZonelist2(&file, 
                                    name.c_str(),
                                    nzones,
                                    ndims,
                                    &nodelist.front(),
                                    nodelist.size(),
                                    0,
                                    low_offset,
                                    high_offset,
                                    &shapetype.front(),
                                    &shapesize.front(),
                                    &shapecount.front(),
                                    nshapes,
                                    optlist.mOptlistPtr);

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// DBPutPHZonelist
//------------------------------------------------------------------------------
inline
int
DBPutPHZonelist(DBfile& file,
                std::string name,
                std::vector<std::vector<int> >& faceNodeLists,
                std::vector<std::vector<int> >& zoneFaceLists,
                unsigned low_offset,
                unsigned high_offset,
                DBoptlist_wrapper& optlist) {

  // Preconditions.
  const unsigned nfaces = faceNodeLists.size();
  const unsigned nzones = zoneFaceLists.size();

  // Construct the flat arrays of face-node info and zone-face info.
  vector<int> nodecnts, nodelist, facecnts, facelist;
  for (unsigned k = 0; k != nfaces; ++k) {
    const std::vector<int>& faceNodes = faceNodeLists[k];
    nodecnts.push_back(faceNodes.size());
    std::copy(faceNodes.begin(), faceNodes.end(), std::back_inserter(nodelist));
  }
  for (unsigned k = 0; k != nzones; ++k) {
    const std::vector<int>& zoneFaces = zoneFaceLists[k];
    facecnts.push_back(zoneFaces.size());
    std::copy(zoneFaces.begin(), zoneFaces.end(), std::back_inserter(facelist));
  }
  CHECK(nodecnts.size() == nfaces);
  CHECK(facecnts.size() == nzones);

  // Do the deed.
  const int result = DBPutPHZonelist(&file, 
                                     name.c_str(),
                                     nfaces,
                                     &nodecnts.front(),
                                     nodelist.size(),
                                     &nodelist.front(),
                                     NULL,
                                     nzones,
                                     &facecnts.front(),
                                     facelist.size(),
                                     &facelist.front(),
                                     0,
                                     low_offset,
                                     high_offset,
                                     optlist.mOptlistPtr);

  // That's it.
  return result;
}

}

//------------------------------------------------------------------------------
// Typedef for vectors of optlists.
//------------------------------------------------------------------------------
typedef std::vector<silo::DBoptlist_wrapper*> vector_of_DBoptlist;

#endif
