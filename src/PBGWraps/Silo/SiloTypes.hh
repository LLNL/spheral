#ifndef __PBGWRAPS_SILOTYPES__
#define __PBGWRAPS_SILOTYPES__

#include "Geometry/Dimension.hh"
#include "Utilities/DBC.hh"
#include "PBGWraps/CXXTypes/CXXTypes.hh"

#include "silo.h"

#include <stdio.h>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <memory>

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
typedef ::DBmrgtree DBmrgtree;
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
struct SiloAttributes {
  const static long _DB_ZONETYPE_POLYHEDRON = DB_ZONETYPE_POLYHEDRON;
  const static long _DB_ZONETYPE_TET = DB_ZONETYPE_TET;
  const static long _DB_ZONETYPE_PYRAMID = DB_ZONETYPE_PYRAMID;
  const static long _DB_ZONETYPE_PRISM = DB_ZONETYPE_PRISM;
  const static long _DB_ZONETYPE_HEX = DB_ZONETYPE_HEX;

  const static long _DB_NETCDF = DB_NETCDF;
  const static long _DB_PDB = DB_PDB;
  const static long _DB_TAURUS = DB_TAURUS;
  const static long _DB_UNKNOWN = DB_UNKNOWN;
  const static long _DB_DEBUG = DB_DEBUG;
  const static long _DB_HDF5X = DB_HDF5X;
  const static long _DB_PDBP = DB_PDBP;

  const static long _DB_HDF5_SEC2_OBSOLETE = DB_HDF5_SEC2_OBSOLETE;
  const static long _DB_HDF5_STDIO_OBSOLETE = DB_HDF5_STDIO_OBSOLETE;
  const static long _DB_HDF5_CORE_OBSOLETE = DB_HDF5_CORE_OBSOLETE;
  const static long _DB_HDF5_MPIO_OBSOLETE = DB_HDF5_MPIO_OBSOLETE;
  const static long _DB_HDF5_MPIOP_OBSOLETE = DB_HDF5_MPIOP_OBSOLETE;

  const static long _DB_H5VFD_DEFAULT = DB_H5VFD_DEFAULT;
  const static long _DB_H5VFD_SEC2 = DB_H5VFD_SEC2;
  const static long _DB_H5VFD_STDIO = DB_H5VFD_STDIO;
  const static long _DB_H5VFD_CORE = DB_H5VFD_CORE;
  const static long _DB_H5VFD_LOG = DB_H5VFD_LOG;
  const static long _DB_H5VFD_SPLIT = DB_H5VFD_SPLIT;
  const static long _DB_H5VFD_DIRECT = DB_H5VFD_DIRECT;
  const static long _DB_H5VFD_FAMILY = DB_H5VFD_FAMILY;
  const static long _DB_H5VFD_MPIO = DB_H5VFD_MPIO;
  const static long _DB_H5VFD_MPIP = DB_H5VFD_MPIP;
  const static long _DB_H5VFD_SILO = DB_H5VFD_SILO;

  const static long _DB_FILE_OPTS_H5_DEFAULT_DEFAULT = DB_FILE_OPTS_H5_DEFAULT_DEFAULT;
  const static long _DB_FILE_OPTS_H5_DEFAULT_SEC2    = DB_FILE_OPTS_H5_DEFAULT_SEC2;
  const static long _DB_FILE_OPTS_H5_DEFAULT_STDIO   = DB_FILE_OPTS_H5_DEFAULT_STDIO;
  const static long _DB_FILE_OPTS_H5_DEFAULT_CORE    = DB_FILE_OPTS_H5_DEFAULT_CORE;
  const static long _DB_FILE_OPTS_H5_DEFAULT_LOG     = DB_FILE_OPTS_H5_DEFAULT_LOG;
  const static long _DB_FILE_OPTS_H5_DEFAULT_SPLIT   = DB_FILE_OPTS_H5_DEFAULT_SPLIT;
  const static long _DB_FILE_OPTS_H5_DEFAULT_DIRECT  = DB_FILE_OPTS_H5_DEFAULT_DIRECT;
  const static long _DB_FILE_OPTS_H5_DEFAULT_FAMILY  = DB_FILE_OPTS_H5_DEFAULT_FAMILY;
  const static long _DB_FILE_OPTS_H5_DEFAULT_MPIO    = DB_FILE_OPTS_H5_DEFAULT_MPIO;
  const static long _DB_FILE_OPTS_H5_DEFAULT_MPIP    = DB_FILE_OPTS_H5_DEFAULT_MPIP;
  const static long _DB_FILE_OPTS_H5_DEFAULT_SILO    = DB_FILE_OPTS_H5_DEFAULT_SILO;
  const static long _DB_FILE_OPTS_LAST               = DB_FILE_OPTS_H5_DEFAULT_SILO;

  const static long _DB_HDF5 = DB_HDF5;
  const static long _DB_HDF5_SEC2 = DB_HDF5_SEC2;
  const static long _DB_HDF5_STDIO = DB_HDF5_STDIO;
  const static long _DB_HDF5_CORE = DB_HDF5_CORE;
  const static long _DB_HDF5_LOG = DB_HDF5_LOG;
  const static long _DB_HDF5_SPLIT = DB_HDF5_SPLIT;
  const static long _DB_HDF5_DIRECT = DB_HDF5_DIRECT;
  const static long _DB_HDF5_FAMILY = DB_HDF5_FAMILY;
  const static long _DB_HDF5_MPIO = DB_HDF5_MPIO;
  const static long _DB_HDF5_MPIOP = DB_HDF5_MPIOP;
  const static long _DB_HDF5_MPIP = DB_HDF5_MPIP;
  const static long _DB_HDF5_SILO = DB_HDF5_SILO;

  const static long _DB_NFILES = DB_NFILES;
  const static long _DB_NFILTERS = DB_NFILTERS;

  const static long _DBAll = DBAll;
  const static long _DBNone = DBNone;
  const static long _DBCalc = DBCalc;
  const static long _DBMatMatnos = DBMatMatnos;
  const static long _DBMatMatlist = DBMatMatlist;
  const static long _DBMatMixList = DBMatMixList;
  const static long _DBCurveArrays = DBCurveArrays;
  const static long _DBPMCoords = DBPMCoords;
  const static long _DBPVData = DBPVData;
  const static long _DBQMCoords = DBQMCoords;
  const static long _DBQVData = DBQVData;
  const static long _DBUMCoords = DBUMCoords;
  const static long _DBUMFacelist = DBUMFacelist;
  const static long _DBUMZonelist = DBUMZonelist;
  const static long _DBUVData = DBUVData;
  const static long _DBFacelistInfo = DBFacelistInfo;
  const static long _DBZonelistInfo = DBZonelistInfo;
  const static long _DBMatMatnames = DBMatMatnames;
  const static long _DBUMGlobNodeNo = DBUMGlobNodeNo;
  const static long _DBZonelistGlobZoneNo = DBZonelistGlobZoneNo;
  const static long _DBMatMatcolors = DBMatMatcolors;
  const static long _DBCSGMBoundaryInfo = DBCSGMBoundaryInfo;
  const static long _DBCSGMZonelist = DBCSGMZonelist;
  const static long _DBCSGMBoundaryNames = DBCSGMBoundaryNames;
  const static long _DBCSGVData = DBCSGVData;
  const static long _DBCSGZonelistZoneNames = DBCSGZonelistZoneNames;
  const static long _DBCSGZonelistRegNames = DBCSGZonelistRegNames;
  const static long _DBMMADJNodelists = DBMMADJNodelists;
  const static long _DBMMADJZonelists = DBMMADJZonelists;
  const static long _DBPMGlobNodeNo = DBPMGlobNodeNo;

  const static long _DB_INVALID_OBJECT = DB_INVALID_OBJECT;
  const static long _DB_QUADRECT = DB_QUADRECT;
  const static long _DB_QUADCURV = DB_QUADCURV;
  const static long _DB_QUADMESH = DB_QUADMESH;
  const static long _DB_QUADVAR = DB_QUADVAR;
  const static long _DB_UCDMESH = DB_UCDMESH;
  const static long _DB_UCDVAR = DB_UCDVAR;
  const static long _DB_MULTIMESH = DB_MULTIMESH;
  const static long _DB_MULTIVAR = DB_MULTIVAR;
  const static long _DB_MULTIMAT = DB_MULTIMAT;
  const static long _DB_MULTIMATSPECIES = DB_MULTIMATSPECIES;
  const static long _DB_MULTIBLOCKMESH = DB_MULTIBLOCKMESH;
  const static long _DB_MULTIBLOCKVAR = DB_MULTIBLOCKVAR;
  const static long _DB_MULTIMESHADJ = DB_MULTIMESHADJ;
  const static long _DB_MATERIAL = DB_MATERIAL;
  const static long _DB_MATSPECIES = DB_MATSPECIES;
  const static long _DB_FACELIST = DB_FACELIST;
  const static long _DB_ZONELIST = DB_ZONELIST;
  const static long _DB_EDGELIST = DB_EDGELIST;
  const static long _DB_PHZONELIST = DB_PHZONELIST;
  const static long _DB_CSGZONELIST = DB_CSGZONELIST;
  const static long _DB_CSGMESH = DB_CSGMESH;
  const static long _DB_CSGVAR = DB_CSGVAR;
  const static long _DB_CURVE = DB_CURVE;
  const static long _DB_DEFVARS = DB_DEFVARS;
  const static long _DB_POINTMESH = DB_POINTMESH;
  const static long _DB_POINTVAR = DB_POINTVAR;
  const static long _DB_ARRAY = DB_ARRAY;
  const static long _DB_DIR = DB_DIR;
  const static long _DB_VARIABLE = DB_VARIABLE;
  const static long _DB_MRGTREE = DB_MRGTREE;
  const static long _DB_GROUPELMAP = DB_GROUPELMAP;
  const static long _DB_MRGVAR = DB_MRGVAR;
  const static long _DB_USERDEF = DB_USERDEF;

  const static long _DB_INT = DB_INT;
  const static long _DB_SHORT = DB_SHORT;
  const static long _DB_LONG = DB_LONG;
  const static long _DB_FLOAT = DB_FLOAT;
  const static long _DB_DOUBLE = DB_DOUBLE;
  const static long _DB_CHAR = DB_CHAR;
  const static long _DB_LONG_LONG = DB_LONG_LONG;
  const static long _DB_NOTYPE = DB_NOTYPE;

  const static long _DB_CLOBBER = DB_CLOBBER;
  const static long _DB_NOCLOBBER = DB_NOCLOBBER;

  const static long _DB_READ = DB_READ;
  const static long _DB_APPEND = DB_APPEND;

  const static long _DB_LOCAL = DB_LOCAL;
  const static long _DB_SUN3 = DB_SUN3;
  const static long _DB_SUN4 = DB_SUN4;
  const static long _DB_SGI = DB_SGI;
  const static long _DB_RS6000 = DB_RS6000;
  const static long _DB_CRAY = DB_CRAY;
  const static long _DB_INTEL = DB_INTEL;

  const static long _DBOPT_FIRST = DBOPT_FIRST;
  const static long _DBOPT_ALIGN = DBOPT_ALIGN;
  const static long _DBOPT_COORDSYS = DBOPT_COORDSYS;
  const static long _DBOPT_CYCLE = DBOPT_CYCLE;
  const static long _DBOPT_FACETYPE = DBOPT_FACETYPE;
  const static long _DBOPT_HI_OFFSET = DBOPT_HI_OFFSET;
  const static long _DBOPT_LO_OFFSET = DBOPT_LO_OFFSET;
  const static long _DBOPT_LABEL = DBOPT_LABEL;
  const static long _DBOPT_XLABEL = DBOPT_XLABEL;
  const static long _DBOPT_YLABEL = DBOPT_YLABEL;
  const static long _DBOPT_ZLABEL = DBOPT_ZLABEL;
  const static long _DBOPT_MAJORORDER = DBOPT_MAJORORDER;
  const static long _DBOPT_NSPACE = DBOPT_NSPACE;
  const static long _DBOPT_ORIGIN = DBOPT_ORIGIN;
  const static long _DBOPT_PLANAR = DBOPT_PLANAR;
  const static long _DBOPT_TIME = DBOPT_TIME;
  const static long _DBOPT_UNITS = DBOPT_UNITS;
  const static long _DBOPT_XUNITS = DBOPT_XUNITS;
  const static long _DBOPT_YUNITS = DBOPT_YUNITS;
  const static long _DBOPT_ZUNITS = DBOPT_ZUNITS;
  const static long _DBOPT_DTIME = DBOPT_DTIME;
  const static long _DBOPT_USESPECMF = DBOPT_USESPECMF;
  const static long _DBOPT_XVARNAME = DBOPT_XVARNAME;
  const static long _DBOPT_YVARNAME = DBOPT_YVARNAME;
  const static long _DBOPT_ZVARNAME = DBOPT_ZVARNAME;
  const static long _DBOPT_ASCII_LABEL = DBOPT_ASCII_LABEL;
  const static long _DBOPT_MATNOS = DBOPT_MATNOS;
  const static long _DBOPT_NMATNOS = DBOPT_NMATNOS;
  const static long _DBOPT_MATNAME = DBOPT_MATNAME;
  const static long _DBOPT_NMAT = DBOPT_NMAT;
  const static long _DBOPT_NMATSPEC = DBOPT_NMATSPEC;
  const static long _DBOPT_BASEINDEX = DBOPT_BASEINDEX;
  const static long _DBOPT_ZONENUM = DBOPT_ZONENUM;
  const static long _DBOPT_NODENUM = DBOPT_NODENUM;
  const static long _DBOPT_BLOCKORIGIN = DBOPT_BLOCKORIGIN;
  const static long _DBOPT_GROUPNUM = DBOPT_GROUPNUM;
  const static long _DBOPT_GROUPORIGIN = DBOPT_GROUPORIGIN;
  const static long _DBOPT_NGROUPS = DBOPT_NGROUPS;
  const static long _DBOPT_MATNAMES = DBOPT_MATNAMES;
  const static long _DBOPT_EXTENTS_SIZE = DBOPT_EXTENTS_SIZE;
  const static long _DBOPT_EXTENTS = DBOPT_EXTENTS;
  const static long _DBOPT_MATCOUNTS = DBOPT_MATCOUNTS;
  const static long _DBOPT_MATLISTS = DBOPT_MATLISTS;
  const static long _DBOPT_MIXLENS = DBOPT_MIXLENS;
  const static long _DBOPT_ZONECOUNTS = DBOPT_ZONECOUNTS;
  const static long _DBOPT_HAS_EXTERNAL_ZONES = DBOPT_HAS_EXTERNAL_ZONES;
  const static long _DBOPT_PHZONELIST = DBOPT_PHZONELIST;
  const static long _DBOPT_MATCOLORS = DBOPT_MATCOLORS;
  const static long _DBOPT_BNDNAMES = DBOPT_BNDNAMES;
  const static long _DBOPT_REGNAMES = DBOPT_REGNAMES;
  const static long _DBOPT_ZONENAMES = DBOPT_ZONENAMES;
  const static long _DBOPT_HIDE_FROM_GUI = DBOPT_HIDE_FROM_GUI;
  const static long _DBOPT_TOPO_DIM = DBOPT_TOPO_DIM;
  const static long _DBOPT_REFERENCE = DBOPT_REFERENCE;
  const static long _DBOPT_GROUPINGS_SIZE = DBOPT_GROUPINGS_SIZE;
  const static long _DBOPT_GROUPINGS = DBOPT_GROUPINGS;
  const static long _DBOPT_GROUPINGNAMES = DBOPT_GROUPINGNAMES;
  const static long _DBOPT_ALLOWMAT0 = DBOPT_ALLOWMAT0;
  const static long _DBOPT_MRGTREE_NAME = DBOPT_MRGTREE_NAME;
  const static long _DBOPT_REGION_PNAMES = DBOPT_REGION_PNAMES;
  const static long _DBOPT_TENSOR_RANK = DBOPT_TENSOR_RANK;
  const static long _DBOPT_MMESH_NAME = DBOPT_MMESH_NAME;
  const static long _DBOPT_TV_CONNECTIVITY = DBOPT_TV_CONNECTIVITY;
  const static long _DBOPT_DISJOINT_MODE = DBOPT_DISJOINT_MODE;
  const static long _DBOPT_MRGV_ONAMES = DBOPT_MRGV_ONAMES;
  const static long _DBOPT_MRGV_RNAMES = DBOPT_MRGV_RNAMES;
  const static long _DBOPT_SPECNAMES = DBOPT_SPECNAMES;
  const static long _DBOPT_SPECCOLORS = DBOPT_SPECCOLORS;
  const static long _DBOPT_LLONGNZNUM = DBOPT_LLONGNZNUM;
  const static long _DBOPT_CONSERVED = DBOPT_CONSERVED;
  const static long _DBOPT_EXTENSIVE = DBOPT_EXTENSIVE;
  const static long _DBOPT_MB_FILE_NS = DBOPT_MB_FILE_NS;
  const static long _DBOPT_MB_BLOCK_NS = DBOPT_MB_BLOCK_NS;
  const static long _DBOPT_MB_BLOCK_TYPE = DBOPT_MB_BLOCK_TYPE;
  const static long _DBOPT_MB_EMPTY_LIST = DBOPT_MB_EMPTY_LIST;
  const static long _DBOPT_MB_EMPTY_COUNT = DBOPT_MB_EMPTY_COUNT;
  const static long _DBOPT_LAST = DBOPT_LAST;
  
  const static long _DBOPT_H5_FIRST = DBOPT_H5_FIRST;
  const static long _DBOPT_H5_VFD = DBOPT_H5_VFD;
  const static long _DBOPT_H5_RAW_FILE_OPTS = DBOPT_H5_RAW_FILE_OPTS;
  const static long _DBOPT_H5_RAW_EXTENSION = DBOPT_H5_RAW_EXTENSION;
  const static long _DBOPT_H5_META_FILE_OPTS = DBOPT_H5_META_FILE_OPTS;
  const static long _DBOPT_H5_META_EXTENSION = DBOPT_H5_META_EXTENSION;
  const static long _DBOPT_H5_CORE_ALLOC_INC = DBOPT_H5_CORE_ALLOC_INC;
  const static long _DBOPT_H5_CORE_NO_BACK_STORE = DBOPT_H5_CORE_NO_BACK_STORE;
  const static long _DBOPT_H5_META_BLOCK_SIZE = DBOPT_H5_META_BLOCK_SIZE;
  const static long _DBOPT_H5_SMALL_RAW_SIZE = DBOPT_H5_SMALL_RAW_SIZE;
  const static long _DBOPT_H5_ALIGN_MIN = DBOPT_H5_ALIGN_MIN;
  const static long _DBOPT_H5_ALIGN_VAL = DBOPT_H5_ALIGN_VAL;
  const static long _DBOPT_H5_DIRECT_MEM_ALIGN = DBOPT_H5_DIRECT_MEM_ALIGN;
  const static long _DBOPT_H5_DIRECT_BLOCK_SIZE = DBOPT_H5_DIRECT_BLOCK_SIZE;
  const static long _DBOPT_H5_DIRECT_BUF_SIZE = DBOPT_H5_DIRECT_BUF_SIZE;
  const static long _DBOPT_H5_LOG_NAME = DBOPT_H5_LOG_NAME;
  const static long _DBOPT_H5_LOG_BUF_SIZE = DBOPT_H5_LOG_BUF_SIZE;
  const static long _DBOPT_H5_MPIO_COMM = DBOPT_H5_MPIO_COMM;
  const static long _DBOPT_H5_MPIO_INFO = DBOPT_H5_MPIO_INFO;
  const static long _DBOPT_H5_MPIP_NO_GPFS_HINTS = DBOPT_H5_MPIP_NO_GPFS_HINTS;
  const static long _DBOPT_H5_SIEVE_BUF_SIZE = DBOPT_H5_SIEVE_BUF_SIZE;
  const static long _DBOPT_H5_CACHE_NELMTS = DBOPT_H5_CACHE_NELMTS;
  const static long _DBOPT_H5_CACHE_NBYTES = DBOPT_H5_CACHE_NBYTES;
  const static long _DBOPT_H5_CACHE_POLICY = DBOPT_H5_CACHE_POLICY;
  const static long _DBOPT_H5_FAM_SIZE = DBOPT_H5_FAM_SIZE;
  const static long _DBOPT_H5_FAM_FILE_OPTS = DBOPT_H5_FAM_FILE_OPTS;
  const static long _DBOPT_H5_USER_DRIVER_ID = DBOPT_H5_USER_DRIVER_ID;
  const static long _DBOPT_H5_USER_DRIVER_INFO = DBOPT_H5_USER_DRIVER_INFO;
  const static long _DBOPT_H5_SILO_BLOCK_SIZE = DBOPT_H5_SILO_BLOCK_SIZE;
  const static long _DBOPT_H5_SILO_BLOCK_COUNT = DBOPT_H5_SILO_BLOCK_COUNT;
  const static long _DBOPT_H5_SILO_LOG_STATS = DBOPT_H5_SILO_LOG_STATS;
  const static long _DBOPT_H5_SILO_USE_DIRECT = DBOPT_H5_SILO_USE_DIRECT;
  const static long _DBOPT_H5_LAST = DBOPT_H5_LAST;

  const static long _DB_TOP = DB_TOP;
  const static long _DB_NONE = DB_NONE;
  const static long _DB_ALL = DB_ALL;
  const static long _DB_ABORT = DB_ABORT;
  const static long _DB_SUSPEND = DB_SUSPEND;
  const static long _DB_RESUME = DB_RESUME;
  const static long _DB_ALL_AND_DRVR = DB_ALL_AND_DRVR;

  const static long _DB_ROWMAJOR = DB_ROWMAJOR;
  const static long _DB_COLMAJOR = DB_COLMAJOR;

  const static long _DB_NOTCENT = DB_NOTCENT;
  const static long _DB_NODECENT = DB_NODECENT;
  const static long _DB_ZONECENT = DB_ZONECENT;
  const static long _DB_FACECENT = DB_FACECENT;
  const static long _DB_BNDCENT = DB_BNDCENT;
  const static long _DB_EDGECENT = DB_EDGECENT;
  const static long _DB_BLOCKCENT = DB_BLOCKCENT;

  const static long _DB_CARTESIAN = DB_CARTESIAN;
  const static long _DB_CYLINDRICAL = DB_CYLINDRICAL;
  const static long _DB_SPHERICAL = DB_SPHERICAL;
  const static long _DB_NUMERICAL = DB_NUMERICAL;
  const static long _DB_OTHER = DB_OTHER;

  const static long _DB_RECTILINEAR = DB_RECTILINEAR;
  const static long _DB_CURVILINEAR = DB_CURVILINEAR;

  const static long _DB_AREA = DB_AREA;
  const static long _DB_VOLUME = DB_VOLUME;

  const static long _DB_ON = DB_ON;
  const static long _DB_OFF = DB_OFF;

  const static long _DB_ABUTTING = DB_ABUTTING;
  const static long _DB_FLOATING = DB_FLOATING;

  const static long _DB_VARTYPE_SCALAR = DB_VARTYPE_SCALAR;
  const static long _DB_VARTYPE_VECTOR = DB_VARTYPE_VECTOR;
  const static long _DB_VARTYPE_TENSOR = DB_VARTYPE_TENSOR;
  const static long _DB_VARTYPE_SYMTENSOR = DB_VARTYPE_SYMTENSOR;
  const static long _DB_VARTYPE_ARRAY = DB_VARTYPE_ARRAY;
  const static long _DB_VARTYPE_MATERIAL = DB_VARTYPE_MATERIAL;
  const static long _DB_VARTYPE_SPECIES = DB_VARTYPE_SPECIES;
  const static long _DB_VARTYPE_LABEL = DB_VARTYPE_LABEL;

  const static long _DBCSG_QUADRIC_G = DBCSG_QUADRIC_G;
  const static long _DBCSG_SPHERE_PR = DBCSG_SPHERE_PR;
  const static long _DBCSG_ELLIPSOID_PRRR = DBCSG_ELLIPSOID_PRRR;
  const static long _DBCSG_PLANE_G = DBCSG_PLANE_G;
  const static long _DBCSG_PLANE_X = DBCSG_PLANE_X;
  const static long _DBCSG_PLANE_Y = DBCSG_PLANE_Y;
  const static long _DBCSG_PLANE_Z = DBCSG_PLANE_Z;
  const static long _DBCSG_PLANE_PN = DBCSG_PLANE_PN;
  const static long _DBCSG_PLANE_PPP = DBCSG_PLANE_PPP;
  const static long _DBCSG_CYLINDER_PNLR = DBCSG_CYLINDER_PNLR;
  const static long _DBCSG_CYLINDER_PPR = DBCSG_CYLINDER_PPR;
  const static long _DBCSG_BOX_XYZXYZ = DBCSG_BOX_XYZXYZ;
  const static long _DBCSG_CONE_PNLA = DBCSG_CONE_PNLA;
  const static long _DBCSG_CONE_PPA = DBCSG_CONE_PPA;
  const static long _DBCSG_POLYHEDRON_KF = DBCSG_POLYHEDRON_KF;
  const static long _DBCSG_HEX_6F = DBCSG_HEX_6F;
  const static long _DBCSG_TET_4F = DBCSG_TET_4F;
  const static long _DBCSG_PYRAMID_5F = DBCSG_PYRAMID_5F;
  const static long _DBCSG_PRISM_5F = DBCSG_PRISM_5F;

  const static long _DBCSG_QUADRATIC_G = DBCSG_QUADRATIC_G;
  const static long _DBCSG_CIRCLE_PR = DBCSG_CIRCLE_PR;
  const static long _DBCSG_ELLIPSE_PRR = DBCSG_ELLIPSE_PRR;
  const static long _DBCSG_LINE_G = DBCSG_LINE_G;
  const static long _DBCSG_LINE_X = DBCSG_LINE_X;
  const static long _DBCSG_LINE_Y = DBCSG_LINE_Y;
  const static long _DBCSG_LINE_PN = DBCSG_LINE_PN;
  const static long _DBCSG_LINE_PP = DBCSG_LINE_PP;
  const static long _DBCSG_BOX_XYXY = DBCSG_BOX_XYXY;
  const static long _DBCSG_ANGLE_PNLA = DBCSG_ANGLE_PNLA;
  const static long _DBCSG_ANGLE_PPA = DBCSG_ANGLE_PPA;
  const static long _DBCSG_POLYGON_KP = DBCSG_POLYGON_KP;
  const static long _DBCSG_TRI_3P = DBCSG_TRI_3P;
  const static long _DBCSG_QUAD_4P = DBCSG_QUAD_4P;

  const static long _DBCSG_INNER = DBCSG_INNER;
  const static long _DBCSG_OUTER = DBCSG_OUTER;
  const static long _DBCSG_ON = DBCSG_ON;
  const static long _DBCSG_UNION = DBCSG_UNION;
  const static long _DBCSG_INTERSECT = DBCSG_INTERSECT;
  const static long _DBCSG_DIFF = DBCSG_DIFF;
  const static long _DBCSG_COMPLIMENT = DBCSG_COMPLIMENT;
  const static long _DBCSG_XFORM = DBCSG_XFORM;
  const static long _DBCSG_SWEEP = DBCSG_SWEEP;

  const static long _DB_PREORDER = DB_PREORDER;
  const static long _DB_POSTORDER = DB_POSTORDER;
  const static long _DB_FROMCWR = DB_FROMCWR;

  const static long _DB_ZONETYPE_BEAM = DB_ZONETYPE_BEAM;

  const static long _DB_ZONETYPE_POLYGON = DB_ZONETYPE_POLYGON;
  const static long _DB_ZONETYPE_TRIANGLE = DB_ZONETYPE_TRIANGLE;
  const static long _DB_ZONETYPE_QUAD = DB_ZONETYPE_QUAD;
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
  std::vector<std::shared_ptr<void> > mCache;

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
      std::shared_ptr<void> voidValue(new Value(value));
      optlist_wrapper.mCache.push_back(voidValue);
      return DBAddOption(optlist_wrapper.mOptlistPtr, option, voidValue.get());
    }
    int writeVector(DBoptlist_wrapper& optlist_wrapper,
                    const int option,
                    const int option_size,
                    const std::vector<Value>& value) {
      DBoptlist_wrapper::AddOptionFunctor<int>().writeValue(optlist_wrapper, option_size, value.size());
      std::shared_ptr<void> voidValue(new std::vector<Value>(value));
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
                              
//------------------------------------------------------------------------------
// Wrapper class to handle the memory managemnt necessary with DBmrgtree
//------------------------------------------------------------------------------
struct DBmrgtree_wrapper {
  DBmrgtree* mDBmrgtree;

  // Constructors.
  DBmrgtree_wrapper(int mesh_type,
                    int info_bits,
                    int max_children,
                    DBoptlist_wrapper optlist):
    mDBmrgtree(DBMakeMrgtree(mesh_type, info_bits, max_children, optlist.mOptlistPtr)) {}

  // Destructor.
  ~DBmrgtree_wrapper() {
    DBFreeMrgtree(mDBmrgtree);
  }

  // name
  std::string name() const { return (mDBmrgtree->name != NULL ? 
                                     std::string(mDBmrgtree->name) :
                                     std::string()); }
  void name(std::string val) { 
    mDBmrgtree->name = new char[val.length() + 1];
    strcpy(mDBmrgtree->name, val.c_str());
  }

  // src_mesh_name
  std::string src_mesh_name() const { return (mDBmrgtree->src_mesh_name != NULL ? 
                                              std::string(mDBmrgtree->src_mesh_name) :
                                              std::string()); }
  void src_mesh_name(std::string val) { 
    mDBmrgtree->src_mesh_name = new char[val.length() + 1];
    strcpy(mDBmrgtree->src_mesh_name, val.c_str());
  }

  // src_mesh_type
  int src_mesh_type() const { return mDBmrgtree->src_mesh_type; }
  void src_mesh_type(int val) { mDBmrgtree->src_mesh_type = val; }

  // type_info_bits
  int type_info_bits() const { return mDBmrgtree->type_info_bits; }
  void type_info_bits(int val) { mDBmrgtree->type_info_bits = val; }

  // num_nodes
  int num_nodes() const { return mDBmrgtree->num_nodes; }
  void num_nodes(int val) { mDBmrgtree->num_nodes = val; }
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
    std::shared_ptr<void> voidValue(new std::string(value));
    optlist_wrapper.mCache.push_back(voidValue);
    return DBAddOption(optlist_wrapper.mOptlistPtr, option, (char*) ((std::string*) voidValue.get())->c_str());
  }
  int
  writeVector(DBoptlist_wrapper& optlist_wrapper,
              const int option,
              const int option_size,
              const std::vector<std::string>& value) {
    VERIFY(optlist_wrapper.addOption<int>(option_size, value.size()) == 0);
    std::shared_ptr<void> voidValue(new char*[value.size()]);
    char** charArray = (char**) voidValue.get();
    for (auto k = 0; k < value.size(); ++k) {
      charArray[k] = new char[value[k].size() + 1];
      strcpy(charArray[k], value[k].c_str());
    }
    optlist_wrapper.mCache.push_back(voidValue);
    return DBAddOption(optlist_wrapper.mOptlistPtr, option, charArray);
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
// DBOpen
//------------------------------------------------------------------------------
inline
DBfile*
DBOpen_wrap(std::string pathName,
              int type,
              int mode) {
  return DBOpen(pathName.c_str(), type, mode);
}

//------------------------------------------------------------------------------
// DBMakeMrgtree
//------------------------------------------------------------------------------
// inline
// DBmrgtree*
// DBMakeMrgtree_wrap(int mesh_type,
//                    int info_bits,
//                    int max_children,
//                    DBoptlist_wrapper& optlist) {
//   return DBMakeMrgtree(mesh_type,
//                        info_bits,
//                        max_children,
//                        optlist.mOptlistPtr);
// }

//------------------------------------------------------------------------------
// DBFreeMrgtree
//------------------------------------------------------------------------------
// inline
// void
// DBFreeMrgtree_wrap(DBmrgtree& tree) {
//   DBFreeMrgtree(&tree);
// }

//------------------------------------------------------------------------------
// DBClose
//------------------------------------------------------------------------------
inline
int
DBClose(DBfile& file) {
  return DBClose(&file);
}

//------------------------------------------------------------------------------
// DBMkDir
//------------------------------------------------------------------------------
inline
int
DBMkDir(DBfile& file,
        std::string dirname) {
  return DBMkDir(&file, dirname.c_str());
}

//------------------------------------------------------------------------------
// DBSetDir
//------------------------------------------------------------------------------
inline
int
DBSetDir(DBfile& file,
         std::string dirname) {
  return DBSetDir(&file, dirname.c_str());
}

//------------------------------------------------------------------------------
// DBGetDir
//------------------------------------------------------------------------------
inline
std::string
DBGetDir(DBfile& file) {
  char result[256];
  auto valid = DBGetDir(&file, result);
  VERIFY2(valid == 0, "Silo ERROR: unable to fetch directory name.");
  return std::string(result);
}

//------------------------------------------------------------------------------
// DBCpDir
//------------------------------------------------------------------------------
inline
int
DBCpDir(DBfile& srcFile,
        std::string srcDir,
        DBfile& dstFile,
        std::string dstDir) {
  return DBCpDir(&srcFile, srcDir.c_str(), &dstFile, dstDir.c_str());
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
  return DBWrite(&file,
                 varname.c_str(),
                 (void*) &var, 
                 &(SiloTraits<T>::dims()).front(),
                 SiloTraits<T>::dims().size(),
                 SiloTraits<T>::datatype());
}

//------------------------------------------------------------------------------
// DBWrite vector<T>
//------------------------------------------------------------------------------
template<typename T>
inline
int
DBWrite_vector(DBfile& file,
               std::string varname,
               std::vector<T>& var) {
  auto dims = std::vector<int>(1, var.size());
  return DBWrite(&file,
                 varname.c_str(),
                 (void*) &var.front(), 
                 &dims.front(),
                 1,
                 SiloTraits<T>::datatype());
}

//------------------------------------------------------------------------------
// DBWrite vector<vector<T>>
//------------------------------------------------------------------------------
template<typename T>
inline
int
DBWrite_vector_of_vector(DBfile& file,
                         std::string varname,
                         std::vector<std::vector<T>>& var) {
  auto ndims = var.size();
  auto dims = std::vector<int>(ndims);
  std::vector<T> varlinear;
  for (auto i = 0; i < ndims; ++i) {
    dims[i] = var[i].size();
    auto istart = varlinear.size();
    varlinear.resize(varlinear.size() + dims[i]);
    std::copy(var[i].begin(), var[i].end(), varlinear.begin() + istart);
  }
  return DBWrite(&file,
                 varname.c_str(),
                 (void*) &varlinear.front(),
                 &dims.front(),
                 ndims,
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
  std::vector<T> flatValues;
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
              std::vector<int>  dims,
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

  // If dims is empty, set it as a 1D list based on matlist.
  if (dims.empty()) dims = std::vector<int>(1, matlist.size());

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
// DBPutQuadmesh
// Note we assume just the unique (x,y,z) coordinates are provided, but we
// replicate them here for the silo file writing.
//------------------------------------------------------------------------------
inline
int
DBPutQuadmesh(DBfile& file,
              std::string name,
              std::vector<std::vector<double> >& coords,
              DBoptlist_wrapper& optlist) {

  // Preconditions.
  const auto ndims = coords.size();
  VERIFY(ndims == 2 or ndims == 3);

  // Number of nodes in each dimension.
  auto nxnodes = coords[0].size();
  auto nxynodes = nxnodes*coords[1].size();
  std::vector<int> meshdims(ndims);
  auto nnodes = 1;
  for (auto k = 0; k < ndims; ++k) {
    meshdims[k] = coords[k].size();
    nnodes *= coords[k].size();
  }

  // We need the C-stylish pointers to the coordinates.
  // This is where we flesh out to the nnodes number of values too.
  double** coordPtrs = new double*[ndims];
  for (auto k = 0; k < ndims; ++k) coordPtrs[k] = new double[nnodes];
  for (auto inode = 0; inode < nnodes; ++inode) {
    const size_t index[3] = {inode % nxnodes,
                             (inode % nxynodes) / nxnodes,
                             inode / nxynodes};
    for (auto k = 0; k < ndims; ++k) coordPtrs[k][inode] = coords[k][index[k]];
  }

  // Do the deed.
  const int result = DBPutQuadmesh(&file,                            // dbfile
                                   name.c_str(),                     // name
                                   NULL,                             // coordnames
                                   coordPtrs,                        // coords
                                   &meshdims[0],                     // dims
                                   ndims,                            // ndims
                                   SiloTraits<double>::datatype(),   // datatype
                                   DB_NONCOLLINEAR,                  // coordtype
                                   optlist.mOptlistPtr);             // optlist

  // That's it.
  for (auto k = 0; k < ndims; ++k) delete[] coordPtrs[k];
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
  VERIFY(centering == DB_NODECENT or
         centering == DB_EDGECENT or
         centering == DB_FACECENT or
         centering == DB_ZONECENT);

  unsigned i, j;

  // Build the sub-variable names.
  int nvars = Spheral2Silo<T>::numElements();
  int nels = values.size();
  int mixlen = mixValues.size();
  std::vector<char*> varnames;
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
  delete[] vars;
  delete[] mixvars;
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
  VERIFY(centering == DB_NODECENT or
         centering == DB_EDGECENT or
         centering == DB_FACECENT or
         centering == DB_ZONECENT);

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
// DBPutQuadvar
// We assume here that the underlying element type is double.
//------------------------------------------------------------------------------
template<typename T>
inline
int
DBPutQuadvar(DBfile& file,
             std::string name,
             std::string meshName,
             std::vector<T>& values,
             std::vector<T>& mixValues,
             int centering,
             std::vector<int>& vardims,
             DBoptlist_wrapper& optlist) {

  // Preconditions.
  VERIFY(centering == DB_NODECENT or
         centering == DB_ZONECENT);
  auto ndims = vardims.size();
  VERIFY(ndims == 1 or ndims == 2 or ndims == 3);

  // Build the sub-variable names.
  auto nvars = Spheral2Silo<T>::numElements();
  auto nels = values.size();
  auto mixlen = mixValues.size();
  std::vector<char*> varnames;
  for (auto i = 0; i != nvars; ++i) {
    varnames.push_back(const_cast<char*>((name + "_").c_str()));
    sprintf(varnames.back(), "%i", i);
  }

  // Build the sub-variables.
  double** vars = new double*[nvars];
  double** mixvars = new double*[nvars];
  for (auto i = 0; i != nvars; ++i) {
    vars[i] = new double[nels];
    mixvars[i] = new double[mixlen];
  }
  for (auto j = 0; j != nels; ++j) Spheral2Silo<T>::copyElement(values[j], vars, j);
  for (auto j = 0; j != mixlen; ++j) Spheral2Silo<T>::copyElement(mixValues[j], mixvars, j);

  const auto result = DBPutQuadvar(&file,
                                   name.c_str(),
                                   meshName.c_str(),
                                   nvars,
                                   &varnames.front(),
                                   (void*) vars,
                                   &vardims.front(),
                                   ndims,
                                   (void*) mixvars,
                                   mixlen,
                                   SiloTraits<double>::datatype(),
                                   centering,
                                   optlist.mOptlistPtr);

  // That's it.
  for (auto i = 0; i != nvars; ++i) {
    delete[] vars[i];
    delete[] mixvars[i];
  }
  delete[] vars;
  delete[] mixvars;
  return result;
}

//------------------------------------------------------------------------------
// DBPutQuadvar1
//------------------------------------------------------------------------------
template<typename T>
inline
int
DBPutQuadvar1(DBfile& file,
              std::string name,
              std::string meshName,
              std::vector<T>& values,
              std::vector<T>& mixValues,
              int centering,
              std::vector<int>& vardims,
              DBoptlist_wrapper& optlist) {

  // Preconditions.
  VERIFY(centering == DB_NODECENT or
         centering == DB_EDGECENT or
         centering == DB_FACECENT or
         centering == DB_ZONECENT);
  const auto ndims = vardims.size();
  VERIFY(ndims == 1 or ndims == 2 or ndims == 3);

  return DBPutQuadvar1(&file,
                       name.c_str(),
                       meshName.c_str(),
                       (void*) &values.front(),
                       &vardims.front(),
                       ndims,
                       (void*) &mixValues.front(),
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
  std::vector<int> nodelist;
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
  std::vector<int> nodecnts, nodelist, facecnts, facelist;
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

//------------------------------------------------------------------------------
// DBPutPointmesh
//------------------------------------------------------------------------------
inline
int
DBPutPointmesh(DBfile& file,
               std::string name,
               std::vector<std::vector<double> >& coords,
               DBoptlist_wrapper& optlist) {

  // Preconditions.
  const unsigned ndims = coords.size();
  VERIFY(ndims == 2 or ndims == 3);
  const unsigned npoints = coords[0].size();
  for (unsigned idim = 0; idim != ndims; ++idim) VERIFY(coords[idim].size() == npoints);

  // We need the C-stylish pointers to the coordinates.
  double** coordPtrs = new double*[ndims];
  for (unsigned k = 0; k != ndims; ++k) {
    coordPtrs[k] = new double[npoints];
    std::copy(coords[k].begin(), coords[k].end(), coordPtrs[k]);
  }

  // Do the deed.
  const int result = DBPutPointmesh(&file, 
                                    name.c_str(),
                                    ndims,
                                    coordPtrs,
                                    npoints,
                                    SiloTraits<double>::datatype(),
                                    optlist.mOptlistPtr);

  // That's it.
  for (unsigned k = 0; k != ndims; ++k) delete[] coordPtrs[k];
  delete[] coordPtrs;
  return result;
}

//------------------------------------------------------------------------------
// DBPutPointvar
// We assume here that the underlying element type is double.
//------------------------------------------------------------------------------
template<typename T>
inline
int
DBPutPointvar(DBfile& file,
              std::string name,
              std::string meshName,
              std::vector<T>& values,
              DBoptlist_wrapper& optlist) {

  unsigned i, j;

  // Build the sub-variable names.
  int nvars = Spheral2Silo<T>::numElements();
  int nels = values.size();

  // Build the sub-variables.
  double** vars = new double*[nvars];
  for (i = 0; i != nvars; ++i) {
    vars[i] = new double[nels];
  }
  for (j = 0; j != nels; ++j) Spheral2Silo<T>::copyElement(values[j], vars, j);

  const int result = DBPutPointvar(&file,
                                   name.c_str(),
                                   meshName.c_str(),
                                   nvars,
                                   (void*) vars,
                                   nels,
                                   SiloTraits<double>::datatype(),
                                   optlist.mOptlistPtr);

  // That's it.
  for (i = 0; i != nvars; ++i) {
    delete[] vars[i];
  }
  delete[] vars;
  return result;
}

//------------------------------------------------------------------------------
// DBPutPointvar1
//------------------------------------------------------------------------------
template<typename T>
inline
int
DBPutPointvar1(DBfile& file,
               std::string name,
               std::string meshName,
               std::vector<T>& values,
               DBoptlist_wrapper& optlist) {

  return DBPutPointvar1(&file,
                        name.c_str(),
                        meshName.c_str(),
                        (void*) &(*values.begin()),
                        values.size(),
                        SiloTraits<T>::datatype(),
                        optlist.mOptlistPtr);
}

//------------------------------------------------------------------------------
// DBAddRegion
//------------------------------------------------------------------------------
inline
int
DBAddRegion(DBmrgtree_wrapper& tree,
            std::string reg_name,
            int info_bits,
            int max_children,
            std::string maps_name,
            std::vector<int>& seg_ids,
            std::vector<int>& seg_lens,
            std::vector<int>& seg_types,
            DBoptlist_wrapper& optlist) {
  int nsegs = seg_ids.size();
  VERIFY(seg_lens.size() == nsegs);
  VERIFY(seg_types.size() == nsegs);
  return DBAddRegion(tree.mDBmrgtree,
                     reg_name.c_str(),
                     info_bits,
                     max_children,
                     (maps_name.size() > 0 ? maps_name.c_str() : NULL),
                     nsegs,
                     (nsegs > 0 ? &(*seg_ids.begin()) : NULL),
                     (nsegs > 0 ? &(*seg_lens.begin()) : NULL),
                     (nsegs > 0 ? &(*seg_types.begin()) : NULL),
                     optlist.mOptlistPtr);
}

//------------------------------------------------------------------------------
// DBSetCwr
//------------------------------------------------------------------------------
inline
int
DBSetCwr(DBmrgtree_wrapper& tree,
         std::string path) {
  return DBSetCwr(tree.mDBmrgtree,
                  path.c_str());
}

//------------------------------------------------------------------------------
// DBGetCwr
//------------------------------------------------------------------------------
inline
const char*
DBGetCwr(DBmrgtree_wrapper& tree) {
  return DBGetCwr(tree.mDBmrgtree);
}

//------------------------------------------------------------------------------
// DBPutMrgtree
//------------------------------------------------------------------------------
inline
int
DBPutMrgtree(DBfile& file,
             std::string name,
             std::string mesh_name,
             DBmrgtree_wrapper& tree,
             DBoptlist_wrapper& optlist) {
  return DBPutMrgtree(&file,
                      name.c_str(),
                      mesh_name.c_str(),
                      tree.mDBmrgtree,
                      optlist.mOptlistPtr);
}

}

//------------------------------------------------------------------------------
// Typedef for vectors of optlists.
//------------------------------------------------------------------------------
typedef std::vector<silo::DBoptlist_wrapper*> vector_of_DBoptlist;

#endif
