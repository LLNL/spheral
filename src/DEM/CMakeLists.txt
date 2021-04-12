include_directories(.)

set(DEM_inst 
    DEMBase
    ContactModelBase)

set(DEM_sources)

instantiate(DEM_inst DEM_sources)

set(DEM_headers 
    DEMBase.hh
    ContactModelBase.hh)

spheral_install_python_files(DEM.py)

spheral_add_cxx_library(DEM)

