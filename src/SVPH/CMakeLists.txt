include_directories(.)
set(SVPH_inst
    computeSVPHCorrections
    computeSVPHCorrectionsOnFaces
    computeSumVoronoiCellMassDensityFromFaces
    sampleFieldListSVPH
    gradientFieldListSVPH
    SVPHFacetedHydroBase
    SVPHCorrectionsPolicy
    SVPHMassDensityPolicy
    CompatibleFaceSpecificThermalEnergyPolicy
    MeshIdealHPolicy
    SpecificThermalEnergyVolumePolicy
    CellPressurePolicy
   )


set(SVPH_sources )

instantiate(SVPH_inst SVPH_sources)

set(SVPH_headers
    CellPressurePolicy.hh
    CompatibleFaceSpecificThermalEnergyPolicy.hh
    computeSumVoronoiCellMassDensityFromFaces.hh
    computeSVPHCorrections.hh
    computeSVPHCorrectionsOnFaces.hh
    gradientFieldListSVPH.hh
    MeshIdealHPolicy.hh
    sampleFieldListSVPH.hh
    SpecificThermalEnergyVolumePolicy.hh
    SVPHCorrectionsPolicy.hh
    SVPHFacetedHydroBase.hh
    SVPHFacetedHydroBaseInline.hh
    SVPHHydroBase.hh
    SVPHHydroBaseInline.hh
    SVPHMassDensityPolicy.hh
    )

set(SVPH_python_targets
  SVPHHydros.py
  )

spheral_add_cxx_library(SVPH)

install(TARGETS              Spheral_SVPH
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${SVPH_headers}
        DESTINATION include/SVPH
        )

if(NOT ENABLE_CXXONLY)
  install(FILES       ${SVPH_python_targets}
          DESTINATION .
          )
endif()

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_SVPH)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
