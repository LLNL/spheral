set(CRKSPH_inst
    computeCRKSPHSumMassDensity
    computeSolidCRKSPHSumMassDensity
    CRKSPHHydroBase
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

instantiate(CRKSPH_inst CRKSPH_sources)

set(CRKSPH_headers
    CRKSPHHydroBase.hh
    CRKSPHHydroBaseInline.hh
    CRKSPHHydroBaseRZ.hh
    CRKSPHVariant.hh
    SolidCRKSPHHydroBase.hh
    SolidCRKSPHHydroBaseInline.hh
    SolidCRKSPHHydroBaseRZ.hh
    centerOfMass.hh
    computeCRKSPHEvaluation.hh
    computeCRKSPHIntegral.hh
    computeCRKSPHSumMassDensity.hh
    computeHullSumMassDensity.hh
    computeSolidCRKSPHSumMassDensity.hh
    detectSurface.hh
    editMultimaterialSurfaceTopology.hh
    volumeSpacing.hh
    zerothOrderSurfaceCorrections.hh
    )

blt_add_library(NAME        Spheral_CRKSPH
                SOURCES     ${CRKSPH_sources}
                HEADERS     ${CRKSPH_headers}
                DEPENDS_ON  ${spheral_depends}
                )

install(TARGETS              Spheral_CRKSPH
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${CRKSPH_headers}
        DESTINATION include/CRKSPH
        )
