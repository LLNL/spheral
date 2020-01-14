include_directories(.)
set(SPH_inst
    SPHHydroBase
    PSPHHydroBase
    SolidSPHHydroBase
    computeSPHOmegaGradhCorrection
    computeSPHSumMassDensity
    correctSPHSumMassDensity
    computePSPHCorrections
    computeSumVoronoiCellMassDensity
   )


set(SPH_sources )

if(ENABLE_2D)
  list(APPEND SPH_sources
       SPHHydroBaseGSRZ.cc
       SPHHydroBaseRZ.cc
       SolidSPHHydroBaseRZ.cc
      )
endif()

if(ENABLE_3D)
  list(APPEND SPH_sources
       computeHydrostaticEquilibriumPressure.cc
      )
endif()

instantiate(SPH_inst SPH_sources)

set(SPH_headers
    OmegaGradhPolicy.hh
    SPHHydroBase.hh
    SPHHydroBaseGSRZ.hh
    SPHHydroBaseInline.hh
    SPHHydroBaseRZ.hh
    SolidSPHHydroBase.hh
    SolidSPHHydroBaseInline.hh
    SolidSPHHydroBaseRZ.hh
    #SolidSPHHydroBaseRZInline.hh
    computeHydrostaticEquilibriumPressure.hh
    computeSPHOmegaGradhCorrection.hh
    computeSPHSumMassDensity.hh
    computeSumVoronoiCellMassDensity.hh
    correctSPHSumMassDensity.hh
    )

set(SPH_python_targets
    SPHHydros.py
    PSPHHydros.py
    SPHUtilities.py
    FacetedSurfaceASPHHydro.py
    )

spheral_add_cxx_library(SPH)

install(TARGETS              Spheral_SPH
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${SPH_headers}
        DESTINATION include/SPH
        )

if(NOT ENABLE_CXXONLY)
  install(FILES       ${SPH_python_targets}
          DESTINATION .
          )
endif()

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_SPH)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
