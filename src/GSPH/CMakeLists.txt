include_directories(.)
set(GSPH_inst
    GSPHHydroBase
    GSPHSpecificThermalEnergyPolicy
    )


set(GSPH_sources )

instantiate(GSPH_inst GSPH_sources)

set(GSPH_headers
    GSPHHydroBase.hh
    GSPHSpecificThermalEnergyPolicy.hh
    )

spheral_install_python_files(
    GSPHHydros.py
    )

spheral_add_cxx_library(GSPH)

