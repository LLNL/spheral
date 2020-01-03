include_directories(.)
set(CRKSPH_inst
    computeCRKSPHSumMassDensity
    computeSolidCRKSPHSumMassDensity
    computeCRKSPHMoments
    computeCRKSPHCorrections
    CRKSPHHydroBase
    interpolateCRKSPH
    gradientCRKSPH
    detectSurface
    SolidCRKSPHHydroBase
    editMultimaterialSurfaceTopology
    zerothOrderSurfaceCorrections
   )


set(CRKSPH_sources
    centerOfMass.cc
   )

if(ENABLE_2D)
  list(APPEND CRKSPH_sources
       CRKSPHHydroBaseRZ.cc
       SolidCRKSPHHydroBaseRZ.cc
      )
endif()

if(ENABLE_3D)
  list(APPEND CRKSPH_sources
      )
endif()

instantiate(CRKSPH_inst CRKSPH_sources)

set(CRKSPH_headers
    CRKSPHHydroBase.hh
    CRKSPHHydroBaseInline.hh
    CRKSPHHydroBaseRZ.hh
    CRKSPHUtilities.hh
    CRKSPHUtilitiesInline.hh
    SolidCRKSPHHydroBase.hh
    SolidCRKSPHHydroBaseInline.hh
    SolidCRKSPHHydroBaseRZ.hh
    SurfaceNodeCoupling.hh
    centerOfMass.hh
    computeCRKSPHCorrections.hh
    computeCRKSPHEvaluation.hh
    computeCRKSPHIntegral.hh
    computeCRKSPHMoments.hh
    computeCRKSPHSumMassDensity.hh
    computeHullSumMassDensity.hh
    computeSolidCRKSPHSumMassDensity.hh
    computeVoronoiCentroids.hh
    detectSurface.hh
    gradientCRKSPH.hh
    interpolateCRKSPH.hh
    volumeSpacing.hh
    zerothOrderSurfaceCorrections.hh
    )

set(CRKSPH_python_targets
  CRKSPHHydros.py
  )

blt_add_library(NAME        Spheral_CRKSPH
                SOURCES     ${CRKSPH_sources}
                HEADERS     ${CRKSPH_headers}
                DEPENDS_ON  ${spheral_depends}
                SHARED      TRUE
                )

install(TARGETS              Spheral_CRKSPH
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${CRKSPH_headers}
        DESTINATION include/CRKSPH
        )

if(NOT ENABLE_CXXONLY)
  install(FILES       ${CRKSPH_python_targets}
          DESTINATION .
          )
endif()
    
list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_CRKSPH)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
