set(SPH_inst
    PSPHHydroBase
    SPHHydroBase
    SolidSPHHydroBase
    computePSPHCorrections
    computeSPHOmegaGradhCorrection
    computeSPHSumMassDensity
    computeSumVoronoiCellMassDensity
    correctSPHSumMassDensity
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
    DamagedNodeCoupling.hh
    DamagedNodeCouplingWithFrags.hh
    NodeCoupling.hh
    OmegaGradhPolicy.hh
    PSPHHydroBase.hh
    PSPHHydroBaseInline.hh
    SPHHydroBase.hh
    SPHHydroBaseGSRZ.hh
    SPHHydroBaseInline.hh
    SPHHydroBaseRZ.hh
    SolidSPHHydroBase.hh
    SolidSPHHydroBaseInline.hh
    SolidSPHHydroBaseRZ.hh
    SolidSPHHydroBaseRZInline.hh
    computeHydrostaticEquilibriumPressure.hh
    computePSPHCorrections.hh
    computeSPHOmegaGradhCorrection.hh
    computeSPHSumMassDensity.hh
    computeSumVoronoiCellMassDensity.hh
    correctSPHSumMassDensity.hh
    )

blt_add_library(NAME        Spheral_SPH
                SOURCES     ${SPH_sources}
                HEADERS     ${SPH_headers}
                DEPENDS_ON  ${spheral_depends}
                )

install(TARGETS              Spheral_SPH
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${SPH_headers}
        DESTINATION include/SPH
        )
