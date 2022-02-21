include_directories(.)

set(DEM_inst 
    DEMBase
    HerzianDEM
    LinearSpringDEM)

set(DEM_sources
    DEMFieldNames.cc
    DEMDimension.cc)

instantiate(DEM_inst DEM_sources)

set(DEM_headers 
    HerzianDEM.hh
    LinearSpringDEM.hh
    DEMBase.hh
    DEMDimension.hh
    DEMFieldNames.hh)

spheral_install_python_files(DEM.py)

spheral_add_cxx_library(DEM)

