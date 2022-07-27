include_directories(.)

set(DEM_inst 
    DEMBase
    ReplacePairFieldList
    IncrementPairFieldList
    ReplaceAndIncrementPairFieldList
    LinearSpringDEM)

set(DEM_sources
    DEMFieldNames.cc
    ContactStorageLocation.cc
    DEMDimension.cc)

instantiate(DEM_inst DEM_sources)

set(DEM_headers 
    ReplacePairFieldList.hh
    ReplaceAndIncrementPairFieldList.hh
    IncrementPairFieldList.hh
    LinearSpringDEM.hh
    DEMBase.hh
    DEMDimension.hh
    DEMFieldNames.hh
    ContactStorageLocation.hh)

spheral_install_python_files(DEM.py)

spheral_add_obj_library(DEM)

