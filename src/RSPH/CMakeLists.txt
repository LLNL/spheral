include_directories(.)
set(RSPH_inst
    RSPHHydroBase
    FSISpecificThermalEnergyPolicy
    )


set(RSPH_sources )

instantiate(RSPH_inst RSPH_sources)

set(RSPH_headers
    RSPHHydroBase.hh
    FSISpecificThermalEnergyPolicy.hh
    )

spheral_install_python_files(
    RSPHHydros.py
    )

spheral_add_cxx_library(RSPH)

