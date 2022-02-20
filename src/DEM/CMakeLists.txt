include_directories(.)

set(DEM_inst 
    DEMBase
    ContactModelBase
    DampedLinearSpring
    computeParticleRadius)

set(DEM_sources
    DEMFieldNames.cc
    DEMDimension.cc)

instantiate(DEM_inst DEM_sources)

set(DEM_headers 
    DEMBase.hh
    ContactModelBase.hh
    DampedLinearSpring.hh
    DEMDimension.hh
    computeParticleRadius.hh
    DEMFieldNames.hh)

spheral_install_python_files(DEM.py)

spheral_add_cxx_library(DEM)

