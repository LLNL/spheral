include_directories(.)
set(FSISPH_inst
    SolidFSISPHHydroBase
    FSISPHHydroBase
    computeSurfaceNormals
    computeFSISPHSumMassDensity
    FSISpecificThermalEnergyPolicy
    )


set(FSISPH_sources )

if(ENABLE_2D)
  list(APPEND FSISPH_sources
       FSISPHHydroBaseRZ.cc
       SolidFSISPHHydroBaseRZ.cc
      )
endif()


instantiate(FSISPH_inst FSISPH_sources)

set(FSISPH_headers
    SolidFSISPHHydroBase.hh
    SolidFSISPHHydroBaseRZ.hh
    FSISPHHydroBase.hh
    FSISPHHydroBaseRZ.hh
    computeSurfaceNormals.hh
    computeFSISPHSumMassDensity.hh
    FSISpecificThermalEnergyPolicy.hh
    )

spheral_install_python_files(
  FSISPHHydros.py
  )

spheral_add_cxx_library(FSISPH)

