// Put Python includes first to avoid compile warnings about redefining _POSIX_C_SOURCE
#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/operators.h"

#include "Geometry/Dimension.hh"
#include "SiloWrappers.hh"

namespace py = pybind11;
using namespace pybind11::literals;

namespace {

//------------------------------------------------------------------------------
// Add the type dependent methods to the silo module (single values).
//------------------------------------------------------------------------------
template<typename T>
void
bindSiloTypeMethods(py::module& m) {
  m.def("DBWrite", &silo::DBWrite<T>);
  m.def("DBPutCompountarray", &silo::DBPutCompoundarray<T>);
  m.def("DBReadVar", &silo::DBReadVar<T>);
  m.def("DBPutUcdvar1", &silo::DBPutUcdvar1<T>);
  m.def("DBPutPointvar1", &silo::DBPutPointvar1<T>);
}

//------------------------------------------------------------------------------
// Add the type dependent methods to the silo module (compund types).
//------------------------------------------------------------------------------
template<typename T>
void
bindSiloCompoundTypeMethods(py::module& m) {
  m.def("DBPutUcdvar", &silo::DBPutUcdvar<T>);
  m.def("DBPutPointvar", &silo::DBPutPointvar<T>);
}

}

//------------------------------------------------------------------------------
// Make the module
//------------------------------------------------------------------------------
PYBIND11_PLUGIN(SpheralSilo) {
  py::module m("SpheralSilo", "Spheral module supporting Silo bindings.");

  // Silo structs
  py::class_<DBfile>(m, "DBfile");

  py::class_<silo::DBoptlist_wrapper>(m ,"DBoptlist")
    .def(py::init<int>(), "maxopts"_a=1024)

    .def("addOption", [](silo::DBoptlist_wrapper& self, int option, int value){ return self.addOption(option, value); })
    .def("addOption", [](silo::DBoptlist_wrapper& self, int option, double value){ return self.addOption(option, value); })
    .def("addOption", [](silo::DBoptlist_wrapper& self, int option, std::string value){ return self.addOption(option, value); })

    .def("addOption", [](silo::DBoptlist_wrapper& self, int option_size, int option, const std::vector<int>& value){ return self.addOption(option, value); })
    .def("addOption", [](silo::DBoptlist_wrapper& self, int option_size, int option, const std::vector<double>& value){ return self.addOption(option, value); })
    .def("addOption", [](silo::DBoptlist_wrapper& self, int option_size, int option, const std::vector<std::string>& value){ return self.addOption(option, value); })

    .def("getOption", [](silo::DBoptlist_wrapper& self, int option){ return self.getOption<int>(option); })
    .def("getOption", [](silo::DBoptlist_wrapper& self, int option){ return self.getOption<double>(option); })
    .def("getOption", [](silo::DBoptlist_wrapper& self, int option){ return self.getOption<std::string>(option); })

    .def("getOption", [](silo::DBoptlist_wrapper& self, int option){ return self.getOption<std::vector<int>>(option); })
    .def("getOption", [](silo::DBoptlist_wrapper& self, int option){ return self.getOption<std::vector<double>>(option); })
    .def("getOption", [](silo::DBoptlist_wrapper& self, int option){ return self.getOption<std::vector<std::string>>(option); })
    ;

  py::class_<silo::DBmrgtree_wrapper>(m, "DBmrgtree")
    .def(py::init<int, int, int, silo::DBoptlist_wrapper>(), "mesh_type"_a="DB_POINTMESH", "info_bits"_a=0, "max_children"_a=1024, "optlist"_a=silo::DBoptlist_wrapper(1024))
    .def_property("name",
                  (std::string (silo::DBmrgtree_wrapper::*)() const) &silo::DBmrgtree_wrapper::name,
                  (void (silo::DBmrgtree_wrapper::*)(std::string)) &silo::DBmrgtree_wrapper::name)
    .def_property("src_mesh_name",
                  (std::string (silo::DBmrgtree_wrapper::*)() const) &silo::DBmrgtree_wrapper::src_mesh_name,
                  (void (silo::DBmrgtree_wrapper::*)(std::string)) &silo::DBmrgtree_wrapper::src_mesh_name)
    .def_property("src_mesh_type",
                  (int (silo::DBmrgtree_wrapper::*)() const) &silo::DBmrgtree_wrapper::src_mesh_type,
                  (void (silo::DBmrgtree_wrapper::*)(int)) &silo::DBmrgtree_wrapper::src_mesh_type)
    .def_property("type_info_bits",
                  (int (silo::DBmrgtree_wrapper::*)() const) &silo::DBmrgtree_wrapper::type_info_bits,
                  (void (silo::DBmrgtree_wrapper::*)(int)) &silo::DBmrgtree_wrapper::type_info_bits)
    .def_property("num_nodes",
                  (int (silo::DBmrgtree_wrapper::*)() const) &silo::DBmrgtree_wrapper::num_nodes,
                  (void (silo::DBmrgtree_wrapper::*)(int)) &silo::DBmrgtree_wrapper::num_nodes)
    ;

  // std::vector<DBoptlist>
  py::bind_vector<std::vector<silo::DBoptlist_wrapper>>(m, "vector_of_DBoptlist");

  // Methods
  m.def("DBCreate", &silo::DBCreate_wrap);
  m.def("DBOpen", &silo::DBOpen_wrap);
  m.def("DBClose", &DBClose);
  m.def("DBPutMultimesh", &silo::DBPutMultimesh);
  m.def("DBPutMultimat", &silo::DBPutMultimat);
  m.def("DBPutMultivar", &silo::DBPutMultivar);
  m.def("DBPutMaterial", &silo::DBPutMaterial);
  m.def("DBPutUcdmesh", &silo::DBPutUcdmesh);
  m.def("DBPutDefvars", &silo::DBPutDefvars);
  m.def("DBPutZonelist2", &silo::DBPutZonelist2);
  m.def("DBPutPHZonelist", &silo::DBPutPHZonelist);
  m.def("DBPutPointmesh", &silo::DBPutPointmesh);
  m.def("DBAddRegion", &silo::DBAddRegion);
  m.def("DBSetCwr", &silo::DBSetCwr);
  m.def("DBGetCwr", &silo::DBGetCwr);
  m.def("DBPutMrgtree", &silo::DBPutMrgtree);
  bindSiloTypeMethods<int>(m);
  bindSiloTypeMethods<float>(m);
  bindSiloTypeMethods<double>(m);
  bindSiloTypeMethods<char>(m);
  bindSiloCompoundTypeMethods<Spheral::Dim<2>::Vector>(m);
  bindSiloCompoundTypeMethods<Spheral::Dim<2>::Tensor>(m);
  bindSiloCompoundTypeMethods<Spheral::Dim<2>::SymTensor>(m);
  bindSiloCompoundTypeMethods<Spheral::Dim<3>::Vector>(m);
  bindSiloCompoundTypeMethods<Spheral::Dim<3>::Tensor>(m);
  bindSiloCompoundTypeMethods<Spheral::Dim<3>::SymTensor>(m);

  //............................................................................
  // Taken from the silo.h file, expose the #define variables as module attributes.
  m.attr("DB_ZONETYPE_POLYHEDRON") = DB_ZONETYPE_POLYHEDRON;
  m.attr("DB_ZONETYPE_TET") = DB_ZONETYPE_TET;
  m.attr("DB_ZONETYPE_PYRAMID") = DB_ZONETYPE_PYRAMID;
  m.attr("DB_ZONETYPE_PRISM") = DB_ZONETYPE_PRISM;
  m.attr("DB_ZONETYPE_HEX") = DB_ZONETYPE_HEX;

  m.attr("DB_NETCDF") = DB_NETCDF;
  m.attr("DB_PDB") = DB_PDB;
  m.attr("DB_TAURUS") = DB_TAURUS;
  m.attr("DB_UNKNOWN") = DB_UNKNOWN;
  m.attr("DB_DEBUG") = DB_DEBUG;
  m.attr("DB_HDF5X") = DB_HDF5X;
  m.attr("DB_PDBP") = DB_PDBP;

  m.attr("DB_HDF5_SEC2_OBSOLETE") = DB_HDF5_SEC2_OBSOLETE;
  m.attr("DB_HDF5_STDIO_OBSOLETE") = DB_HDF5_STDIO_OBSOLETE;
  m.attr("DB_HDF5_CORE_OBSOLETE") = DB_HDF5_CORE_OBSOLETE;
  m.attr("DB_HDF5_MPIO_OBSOLETE") = DB_HDF5_MPIO_OBSOLETE;
  m.attr("DB_HDF5_MPIOP_OBSOLETE") = DB_HDF5_MPIOP_OBSOLETE;

  m.attr("DB_H5VFD_DEFAULT") = DB_H5VFD_DEFAULT;
  m.attr("DB_H5VFD_SEC2") = DB_H5VFD_SEC2;
  m.attr("DB_H5VFD_STDIO") = DB_H5VFD_STDIO;
  m.attr("DB_H5VFD_CORE") = DB_H5VFD_CORE;
  m.attr("DB_H5VFD_LOG") = DB_H5VFD_LOG;
  m.attr("DB_H5VFD_SPLIT") = DB_H5VFD_SPLIT;
  m.attr("DB_H5VFD_DIRECT") = DB_H5VFD_DIRECT;
  m.attr("DB_H5VFD_FAMILY") = DB_H5VFD_FAMILY;
  m.attr("DB_H5VFD_MPIO") = DB_H5VFD_MPIO;
  m.attr("DB_H5VFD_MPIP") = DB_H5VFD_MPIP;
  m.attr("DB_H5VFD_SILO") = DB_H5VFD_SILO;

  m.attr("DB_FILE_OPTS_H5_DEFAULT_DEFAULT") = DB_FILE_OPTS_H5_DEFAULT_DEFAULT;
  m.attr("DB_FILE_OPTS_H5_DEFAULT_SEC2") = DB_FILE_OPTS_H5_DEFAULT_SEC2;
  m.attr("DB_FILE_OPTS_H5_DEFAULT_STDIO") = DB_FILE_OPTS_H5_DEFAULT_STDIO;
  m.attr("DB_FILE_OPTS_H5_DEFAULT_CORE") = DB_FILE_OPTS_H5_DEFAULT_CORE;
  m.attr("DB_FILE_OPTS_H5_DEFAULT_LOG") = DB_FILE_OPTS_H5_DEFAULT_LOG;
  m.attr("DB_FILE_OPTS_H5_DEFAULT_SPLIT") = DB_FILE_OPTS_H5_DEFAULT_SPLIT;
  m.attr("DB_FILE_OPTS_H5_DEFAULT_DIRECT") = DB_FILE_OPTS_H5_DEFAULT_DIRECT;
  m.attr("DB_FILE_OPTS_H5_DEFAULT_FAMILY") = DB_FILE_OPTS_H5_DEFAULT_FAMILY;
  m.attr("DB_FILE_OPTS_H5_DEFAULT_MPIO") = DB_FILE_OPTS_H5_DEFAULT_MPIO;
  m.attr("DB_FILE_OPTS_H5_DEFAULT_MPIP") = DB_FILE_OPTS_H5_DEFAULT_MPIP;
  m.attr("DB_FILE_OPTS_H5_DEFAULT_SILO") = DB_FILE_OPTS_H5_DEFAULT_SILO;
  m.attr("DB_FILE_OPTS_H5_DEFAULT_SILO") = DB_FILE_OPTS_H5_DEFAULT_SILO;

  m.attr("DB_HDF5") = DB_HDF5;
  m.attr("DB_HDF5_SEC2") = DB_HDF5_SEC2;
  m.attr("DB_HDF5_STDIO") = DB_HDF5_STDIO;
  m.attr("DB_HDF5_CORE") = DB_HDF5_CORE;
  m.attr("DB_HDF5_LOG") = DB_HDF5_LOG;
  m.attr("DB_HDF5_SPLIT") = DB_HDF5_SPLIT;
  m.attr("DB_HDF5_DIRECT") = DB_HDF5_DIRECT;
  m.attr("DB_HDF5_FAMILY") = DB_HDF5_FAMILY;
  m.attr("DB_HDF5_MPIO") = DB_HDF5_MPIO;
  m.attr("DB_HDF5_MPIOP") = DB_HDF5_MPIOP;
  m.attr("DB_HDF5_MPIP") = DB_HDF5_MPIP;
  m.attr("DB_HDF5_SILO") = DB_HDF5_SILO;

  m.attr("DB_NFILES") = DB_NFILES;
  m.attr("DB_NFILTERS") = DB_NFILTERS;

  m.attr("DBAll") = DBAll;
  m.attr("DBNone") = DBNone;
  m.attr("DBCalc") = DBCalc;
  m.attr("DBMatMatnos") = DBMatMatnos;
  m.attr("DBMatMatlist") = DBMatMatlist;
  m.attr("DBMatMixList") = DBMatMixList;
  m.attr("DBCurveArrays") = DBCurveArrays;
  m.attr("DBPMCoords") = DBPMCoords;
  m.attr("DBPVData") = DBPVData;
  m.attr("DBQMCoords") = DBQMCoords;
  m.attr("DBQVData") = DBQVData;
  m.attr("DBUMCoords") = DBUMCoords;
  m.attr("DBUMFacelist") = DBUMFacelist;
  m.attr("DBUMZonelist") = DBUMZonelist;
  m.attr("DBUVData") = DBUVData;
  m.attr("DBFacelistInfo") = DBFacelistInfo;
  m.attr("DBZonelistInfo") = DBZonelistInfo;
  m.attr("DBMatMatnames") = DBMatMatnames;
  m.attr("DBUMGlobNodeNo") = DBUMGlobNodeNo;
  m.attr("DBZonelistGlobZoneNo") = DBZonelistGlobZoneNo;
  m.attr("DBMatMatcolors") = DBMatMatcolors;
  m.attr("DBCSGMBoundaryInfo") = DBCSGMBoundaryInfo;
  m.attr("DBCSGMZonelist") = DBCSGMZonelist;
  m.attr("DBCSGMBoundaryNames") = DBCSGMBoundaryNames;
  m.attr("DBCSGVData") = DBCSGVData;
  m.attr("DBCSGZonelistZoneNames") = DBCSGZonelistZoneNames;
  m.attr("DBCSGZonelistRegNames") = DBCSGZonelistRegNames;
  m.attr("DBMMADJNodelists") = DBMMADJNodelists;
  m.attr("DBMMADJZonelists") = DBMMADJZonelists;
  m.attr("DBPMGlobNodeNo") = DBPMGlobNodeNo;

  m.attr("DB_INVALID_OBJECT") = DB_INVALID_OBJECT;
  m.attr("DB_QUADRECT") = DB_QUADRECT;
  m.attr("DB_QUADCURV") = DB_QUADCURV;
  m.attr("DB_QUADMESH") = DB_QUADMESH;
  m.attr("DB_QUADVAR") = DB_QUADVAR;
  m.attr("DB_UCDMESH") = DB_UCDMESH;
  m.attr("DB_UCDVAR") = DB_UCDVAR;
  m.attr("DB_MULTIMESH") = DB_MULTIMESH;
  m.attr("DB_MULTIVAR") = DB_MULTIVAR;
  m.attr("DB_MULTIMAT") = DB_MULTIMAT;
  m.attr("DB_MULTIMATSPECIES") = DB_MULTIMATSPECIES;
  m.attr("DB_MULTIBLOCKMESH") = DB_MULTIBLOCKMESH;
  m.attr("DB_MULTIBLOCKVAR") = DB_MULTIBLOCKVAR;
  m.attr("DB_MULTIMESHADJ") = DB_MULTIMESHADJ;
  m.attr("DB_MATERIAL") = DB_MATERIAL;
  m.attr("DB_MATSPECIES") = DB_MATSPECIES;
  m.attr("DB_FACELIST") = DB_FACELIST;
  m.attr("DB_ZONELIST") = DB_ZONELIST;
  m.attr("DB_EDGELIST") = DB_EDGELIST;
  m.attr("DB_PHZONELIST") = DB_PHZONELIST;
  m.attr("DB_CSGZONELIST") = DB_CSGZONELIST;
  m.attr("DB_CSGMESH") = DB_CSGMESH;
  m.attr("DB_CSGVAR") = DB_CSGVAR;
  m.attr("DB_CURVE") = DB_CURVE;
  m.attr("DB_DEFVARS") = DB_DEFVARS;
  m.attr("DB_POINTMESH") = DB_POINTMESH;
  m.attr("DB_POINTVAR") = DB_POINTVAR;
  m.attr("DB_ARRAY") = DB_ARRAY;
  m.attr("DB_DIR") = DB_DIR;
  m.attr("DB_VARIABLE") = DB_VARIABLE;
  m.attr("DB_MRGTREE") = DB_MRGTREE;
  m.attr("DB_GROUPELMAP") = DB_GROUPELMAP;
  m.attr("DB_MRGVAR") = DB_MRGVAR;
  m.attr("DB_USERDEF") = DB_USERDEF;

  m.attr("DB_INT") = DB_INT;
  m.attr("DB_SHORT") = DB_SHORT;
  m.attr("DB_LONG") = DB_LONG;
  m.attr("DB_FLOAT") = DB_FLOAT;
  m.attr("DB_DOUBLE") = DB_DOUBLE;
  m.attr("DB_CHAR") = DB_CHAR;
  m.attr("DB_LONG_LONG") = DB_LONG_LONG;
  m.attr("DB_NOTYPE") = DB_NOTYPE;

  m.attr("DB_CLOBBER") = DB_CLOBBER;
  m.attr("DB_NOCLOBBER") = DB_NOCLOBBER;

  m.attr("DB_READ") = DB_READ;
  m.attr("DB_APPEND") = DB_APPEND;

  m.attr("DB_LOCAL") = DB_LOCAL;
  m.attr("DB_SUN3") = DB_SUN3;
  m.attr("DB_SUN4") = DB_SUN4;
  m.attr("DB_SGI") = DB_SGI;
  m.attr("DB_RS6000") = DB_RS6000;
  m.attr("DB_CRAY") = DB_CRAY;
  m.attr("DB_INTEL") = DB_INTEL;

  m.attr("DBOPT_FIRST") = DBOPT_FIRST;
  m.attr("DBOPT_ALIGN") = DBOPT_ALIGN;
  m.attr("DBOPT_COORDSYS") = DBOPT_COORDSYS;
  m.attr("DBOPT_CYCLE") = DBOPT_CYCLE;
  m.attr("DBOPT_FACETYPE") = DBOPT_FACETYPE;
  m.attr("DBOPT_HI_OFFSET") = DBOPT_HI_OFFSET;
  m.attr("DBOPT_LO_OFFSET") = DBOPT_LO_OFFSET;
  m.attr("DBOPT_LABEL") = DBOPT_LABEL;
  m.attr("DBOPT_XLABEL") = DBOPT_XLABEL;
  m.attr("DBOPT_YLABEL") = DBOPT_YLABEL;
  m.attr("DBOPT_ZLABEL") = DBOPT_ZLABEL;
  m.attr("DBOPT_MAJORORDER") = DBOPT_MAJORORDER;
  m.attr("DBOPT_NSPACE") = DBOPT_NSPACE;
  m.attr("DBOPT_ORIGIN") = DBOPT_ORIGIN;
  m.attr("DBOPT_PLANAR") = DBOPT_PLANAR;
  m.attr("DBOPT_TIME") = DBOPT_TIME;
  m.attr("DBOPT_UNITS") = DBOPT_UNITS;
  m.attr("DBOPT_XUNITS") = DBOPT_XUNITS;
  m.attr("DBOPT_YUNITS") = DBOPT_YUNITS;
  m.attr("DBOPT_ZUNITS") = DBOPT_ZUNITS;
  m.attr("DBOPT_DTIME") = DBOPT_DTIME;
  m.attr("DBOPT_USESPECMF") = DBOPT_USESPECMF;
  m.attr("DBOPT_XVARNAME") = DBOPT_XVARNAME;
  m.attr("DBOPT_YVARNAME") = DBOPT_YVARNAME;
  m.attr("DBOPT_ZVARNAME") = DBOPT_ZVARNAME;
  m.attr("DBOPT_ASCII_LABEL") = DBOPT_ASCII_LABEL;
  m.attr("DBOPT_MATNOS") = DBOPT_MATNOS;
  m.attr("DBOPT_NMATNOS") = DBOPT_NMATNOS;
  m.attr("DBOPT_MATNAME") = DBOPT_MATNAME;
  m.attr("DBOPT_NMAT") = DBOPT_NMAT;
  m.attr("DBOPT_NMATSPEC") = DBOPT_NMATSPEC;
  m.attr("DBOPT_BASEINDEX") = DBOPT_BASEINDEX;
  m.attr("DBOPT_ZONENUM") = DBOPT_ZONENUM;
  m.attr("DBOPT_NODENUM") = DBOPT_NODENUM;
  m.attr("DBOPT_BLOCKORIGIN") = DBOPT_BLOCKORIGIN;
  m.attr("DBOPT_GROUPNUM") = DBOPT_GROUPNUM;
  m.attr("DBOPT_GROUPORIGIN") = DBOPT_GROUPORIGIN;
  m.attr("DBOPT_NGROUPS") = DBOPT_NGROUPS;
  m.attr("DBOPT_MATNAMES") = DBOPT_MATNAMES;
  m.attr("DBOPT_EXTENTS_SIZE") = DBOPT_EXTENTS_SIZE;
  m.attr("DBOPT_EXTENTS") = DBOPT_EXTENTS;
  m.attr("DBOPT_MATCOUNTS") = DBOPT_MATCOUNTS;
  m.attr("DBOPT_MATLISTS") = DBOPT_MATLISTS;
  m.attr("DBOPT_MIXLENS") = DBOPT_MIXLENS;
  m.attr("DBOPT_ZONECOUNTS") = DBOPT_ZONECOUNTS;
  m.attr("DBOPT_HAS_EXTERNAL_ZONES") = DBOPT_HAS_EXTERNAL_ZONES;
  m.attr("DBOPT_PHZONELIST") = DBOPT_PHZONELIST;
  m.attr("DBOPT_MATCOLORS") = DBOPT_MATCOLORS;
  m.attr("DBOPT_BNDNAMES") = DBOPT_BNDNAMES;
  m.attr("DBOPT_REGNAMES") = DBOPT_REGNAMES;
  m.attr("DBOPT_ZONENAMES") = DBOPT_ZONENAMES;
  m.attr("DBOPT_HIDE_FROM_GUI") = DBOPT_HIDE_FROM_GUI;
  m.attr("DBOPT_TOPO_DIM") = DBOPT_TOPO_DIM;
  m.attr("DBOPT_REFERENCE") = DBOPT_REFERENCE;
  m.attr("DBOPT_GROUPINGS_SIZE") = DBOPT_GROUPINGS_SIZE;
  m.attr("DBOPT_GROUPINGS") = DBOPT_GROUPINGS;
  m.attr("DBOPT_GROUPINGNAMES") = DBOPT_GROUPINGNAMES;
  m.attr("DBOPT_ALLOWMAT0") = DBOPT_ALLOWMAT0;
  m.attr("DBOPT_MRGTREE_NAME") = DBOPT_MRGTREE_NAME;
  m.attr("DBOPT_REGION_PNAMES") = DBOPT_REGION_PNAMES;
  m.attr("DBOPT_TENSOR_RANK") = DBOPT_TENSOR_RANK;
  m.attr("DBOPT_MMESH_NAME") = DBOPT_MMESH_NAME;
  m.attr("DBOPT_TV_CONNECTIVITY") = DBOPT_TV_CONNECTIVITY;
  m.attr("DBOPT_DISJOINT_MODE") = DBOPT_DISJOINT_MODE;
  m.attr("DBOPT_MRGV_ONAMES") = DBOPT_MRGV_ONAMES;
  m.attr("DBOPT_MRGV_RNAMES") = DBOPT_MRGV_RNAMES;
  m.attr("DBOPT_SPECNAMES") = DBOPT_SPECNAMES;
  m.attr("DBOPT_SPECCOLORS") = DBOPT_SPECCOLORS;
  m.attr("DBOPT_LLONGNZNUM") = DBOPT_LLONGNZNUM;
  m.attr("DBOPT_CONSERVED") = DBOPT_CONSERVED;
  m.attr("DBOPT_EXTENSIVE") = DBOPT_EXTENSIVE;
  m.attr("DBOPT_MB_FILE_NS") = DBOPT_MB_FILE_NS;
  m.attr("DBOPT_MB_BLOCK_NS") = DBOPT_MB_BLOCK_NS;
  m.attr("DBOPT_MB_BLOCK_TYPE") = DBOPT_MB_BLOCK_TYPE;
  m.attr("DBOPT_MB_EMPTY_LIST") = DBOPT_MB_EMPTY_LIST;
  m.attr("DBOPT_MB_EMPTY_COUNT") = DBOPT_MB_EMPTY_COUNT;
  m.attr("DBOPT_LAST") = DBOPT_LAST;
  
  m.attr("DBOPT_H5_FIRST") = DBOPT_H5_FIRST;
  m.attr("DBOPT_H5_VFD") = DBOPT_H5_VFD;
  m.attr("DBOPT_H5_RAW_FILE_OPTS") = DBOPT_H5_RAW_FILE_OPTS;
  m.attr("DBOPT_H5_RAW_EXTENSION") = DBOPT_H5_RAW_EXTENSION;
  m.attr("DBOPT_H5_META_FILE_OPTS") = DBOPT_H5_META_FILE_OPTS;
  m.attr("DBOPT_H5_META_EXTENSION") = DBOPT_H5_META_EXTENSION;
  m.attr("DBOPT_H5_CORE_ALLOC_INC") = DBOPT_H5_CORE_ALLOC_INC;
  m.attr("DBOPT_H5_CORE_NO_BACK_STORE") = DBOPT_H5_CORE_NO_BACK_STORE;
  m.attr("DBOPT_H5_META_BLOCK_SIZE") = DBOPT_H5_META_BLOCK_SIZE;
  m.attr("DBOPT_H5_SMALL_RAW_SIZE") = DBOPT_H5_SMALL_RAW_SIZE;
  m.attr("DBOPT_H5_ALIGN_MIN") = DBOPT_H5_ALIGN_MIN;
  m.attr("DBOPT_H5_ALIGN_VAL") = DBOPT_H5_ALIGN_VAL;
  m.attr("DBOPT_H5_DIRECT_MEM_ALIGN") = DBOPT_H5_DIRECT_MEM_ALIGN;
  m.attr("DBOPT_H5_DIRECT_BLOCK_SIZE") = DBOPT_H5_DIRECT_BLOCK_SIZE;
  m.attr("DBOPT_H5_DIRECT_BUF_SIZE") = DBOPT_H5_DIRECT_BUF_SIZE;
  m.attr("DBOPT_H5_LOG_NAME") = DBOPT_H5_LOG_NAME;
  m.attr("DBOPT_H5_LOG_BUF_SIZE") = DBOPT_H5_LOG_BUF_SIZE;
  m.attr("DBOPT_H5_MPIO_COMM") = DBOPT_H5_MPIO_COMM;
  m.attr("DBOPT_H5_MPIO_INFO") = DBOPT_H5_MPIO_INFO;
  m.attr("DBOPT_H5_MPIP_NO_GPFS_HINTS") = DBOPT_H5_MPIP_NO_GPFS_HINTS;
  m.attr("DBOPT_H5_SIEVE_BUF_SIZE") = DBOPT_H5_SIEVE_BUF_SIZE;
  m.attr("DBOPT_H5_CACHE_NELMTS") = DBOPT_H5_CACHE_NELMTS;
  m.attr("DBOPT_H5_CACHE_NBYTES") = DBOPT_H5_CACHE_NBYTES;
  m.attr("DBOPT_H5_CACHE_POLICY") = DBOPT_H5_CACHE_POLICY;
  m.attr("DBOPT_H5_FAM_SIZE") = DBOPT_H5_FAM_SIZE;
  m.attr("DBOPT_H5_FAM_FILE_OPTS") = DBOPT_H5_FAM_FILE_OPTS;
  m.attr("DBOPT_H5_USER_DRIVER_ID") = DBOPT_H5_USER_DRIVER_ID;
  m.attr("DBOPT_H5_USER_DRIVER_INFO") = DBOPT_H5_USER_DRIVER_INFO;
  m.attr("DBOPT_H5_SILO_BLOCK_SIZE") = DBOPT_H5_SILO_BLOCK_SIZE;
  m.attr("DBOPT_H5_SILO_BLOCK_COUNT") = DBOPT_H5_SILO_BLOCK_COUNT;
  m.attr("DBOPT_H5_SILO_LOG_STATS") = DBOPT_H5_SILO_LOG_STATS;
  m.attr("DBOPT_H5_SILO_USE_DIRECT") = DBOPT_H5_SILO_USE_DIRECT;
  m.attr("DBOPT_H5_LAST") = DBOPT_H5_LAST;

  m.attr("DB_TOP") = DB_TOP;
  m.attr("DB_NONE") = DB_NONE;
  m.attr("DB_ALL") = DB_ALL;
  m.attr("DB_ABORT") = DB_ABORT;
  m.attr("DB_SUSPEND") = DB_SUSPEND;
  m.attr("DB_RESUME") = DB_RESUME;
  m.attr("DB_ALL_AND_DRVR") = DB_ALL_AND_DRVR;

  m.attr("DB_ROWMAJOR") = DB_ROWMAJOR;
  m.attr("DB_COLMAJOR") = DB_COLMAJOR;

  m.attr("DB_NOTCENT") = DB_NOTCENT;
  m.attr("DB_NODECENT") = DB_NODECENT;
  m.attr("DB_ZONECENT") = DB_ZONECENT;
  m.attr("DB_FACECENT") = DB_FACECENT;
  m.attr("DB_BNDCENT") = DB_BNDCENT;
  m.attr("DB_EDGECENT") = DB_EDGECENT;
  m.attr("DB_BLOCKCENT") = DB_BLOCKCENT;

  m.attr("DB_CARTESIAN") = DB_CARTESIAN;
  m.attr("DB_CYLINDRICAL") = DB_CYLINDRICAL;
  m.attr("DB_SPHERICAL") = DB_SPHERICAL;
  m.attr("DB_NUMERICAL") = DB_NUMERICAL;
  m.attr("DB_OTHER") = DB_OTHER;

  m.attr("DB_RECTILINEAR") = DB_RECTILINEAR;
  m.attr("DB_CURVILINEAR") = DB_CURVILINEAR;

  m.attr("DB_AREA") = DB_AREA;
  m.attr("DB_VOLUME") = DB_VOLUME;

  m.attr("DB_ON") = DB_ON;
  m.attr("DB_OFF") = DB_OFF;

  m.attr("DB_ABUTTING") = DB_ABUTTING;
  m.attr("DB_FLOATING") = DB_FLOATING;

  m.attr("DB_VARTYPE_SCALAR") = DB_VARTYPE_SCALAR;
  m.attr("DB_VARTYPE_VECTOR") = DB_VARTYPE_VECTOR;
  m.attr("DB_VARTYPE_TENSOR") = DB_VARTYPE_TENSOR;
  m.attr("DB_VARTYPE_SYMTENSOR") = DB_VARTYPE_SYMTENSOR;
  m.attr("DB_VARTYPE_ARRAY") = DB_VARTYPE_ARRAY;
  m.attr("DB_VARTYPE_MATERIAL") = DB_VARTYPE_MATERIAL;
  m.attr("DB_VARTYPE_SPECIES") = DB_VARTYPE_SPECIES;
  m.attr("DB_VARTYPE_LABEL") = DB_VARTYPE_LABEL;

  m.attr("DBCSG_QUADRIC_G") = DBCSG_QUADRIC_G;
  m.attr("DBCSG_SPHERE_PR") = DBCSG_SPHERE_PR;
  m.attr("DBCSG_ELLIPSOID_PRRR") = DBCSG_ELLIPSOID_PRRR;
  m.attr("DBCSG_PLANE_G") = DBCSG_PLANE_G;
  m.attr("DBCSG_PLANE_X") = DBCSG_PLANE_X;
  m.attr("DBCSG_PLANE_Y") = DBCSG_PLANE_Y;
  m.attr("DBCSG_PLANE_Z") = DBCSG_PLANE_Z;
  m.attr("DBCSG_PLANE_PN") = DBCSG_PLANE_PN;
  m.attr("DBCSG_PLANE_PPP") = DBCSG_PLANE_PPP;
  m.attr("DBCSG_CYLINDER_PNLR") = DBCSG_CYLINDER_PNLR;
  m.attr("DBCSG_CYLINDER_PPR") = DBCSG_CYLINDER_PPR;
  m.attr("DBCSG_BOX_XYZXYZ") = DBCSG_BOX_XYZXYZ;
  m.attr("DBCSG_CONE_PNLA") = DBCSG_CONE_PNLA;
  m.attr("DBCSG_CONE_PPA") = DBCSG_CONE_PPA;
  m.attr("DBCSG_POLYHEDRON_KF") = DBCSG_POLYHEDRON_KF;
  m.attr("DBCSG_HEX_6F") = DBCSG_HEX_6F;
  m.attr("DBCSG_TET_4F") = DBCSG_TET_4F;
  m.attr("DBCSG_PYRAMID_5F") = DBCSG_PYRAMID_5F;
  m.attr("DBCSG_PRISM_5F") = DBCSG_PRISM_5F;

  m.attr("DBCSG_QUADRATIC_G") = DBCSG_QUADRATIC_G;
  m.attr("DBCSG_CIRCLE_PR") = DBCSG_CIRCLE_PR;
  m.attr("DBCSG_ELLIPSE_PRR") = DBCSG_ELLIPSE_PRR;
  m.attr("DBCSG_LINE_G") = DBCSG_LINE_G;
  m.attr("DBCSG_LINE_X") = DBCSG_LINE_X;
  m.attr("DBCSG_LINE_Y") = DBCSG_LINE_Y;
  m.attr("DBCSG_LINE_PN") = DBCSG_LINE_PN;
  m.attr("DBCSG_LINE_PP") = DBCSG_LINE_PP;
  m.attr("DBCSG_BOX_XYXY") = DBCSG_BOX_XYXY;
  m.attr("DBCSG_ANGLE_PNLA") = DBCSG_ANGLE_PNLA;
  m.attr("DBCSG_ANGLE_PPA") = DBCSG_ANGLE_PPA;
  m.attr("DBCSG_POLYGON_KP") = DBCSG_POLYGON_KP;
  m.attr("DBCSG_TRI_3P") = DBCSG_TRI_3P;
  m.attr("DBCSG_QUAD_4P") = DBCSG_QUAD_4P;

  m.attr("DBCSG_INNER") = DBCSG_INNER;
  m.attr("DBCSG_OUTER") = DBCSG_OUTER;
  m.attr("DBCSG_ON") = DBCSG_ON;
  m.attr("DBCSG_UNION") = DBCSG_UNION;
  m.attr("DBCSG_INTERSECT") = DBCSG_INTERSECT;
  m.attr("DBCSG_DIFF") = DBCSG_DIFF;
  m.attr("DBCSG_COMPLIMENT") = DBCSG_COMPLIMENT;
  m.attr("DBCSG_XFORM") = DBCSG_XFORM;
  m.attr("DBCSG_SWEEP") = DBCSG_SWEEP;

  m.attr("DB_PREORDER") = DB_PREORDER;
  m.attr("DB_POSTORDER") = DB_POSTORDER;
  m.attr("DB_FROMCWR") = DB_FROMCWR;

  m.attr("DB_ZONETYPE_BEAM") = DB_ZONETYPE_BEAM;

  m.attr("DB_ZONETYPE_POLYGON") = DB_ZONETYPE_POLYGON;
  m.attr("DB_ZONETYPE_TRIANGLE") = DB_ZONETYPE_TRIANGLE;
  m.attr("DB_ZONETYPE_QUAD") = DB_ZONETYPE_QUAD;

  return m.ptr();
}
