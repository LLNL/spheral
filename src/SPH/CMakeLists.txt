set(SPH_sources2d
    PSPHHydroBaseInst.cc.2.cc
    SPHHydroBaseInst.cc.2.cc
    SPHHydroBaseGSRZ.cc
    SPHHydroBaseRZ.cc
    SolidSPHHydroBaseInst.cc.2.cc
    SolidSPHHydroBaseRZ.cc
    computePSPHCorrectionsInst.cc.2.cc
    computeSPHOmegaGradhCorrectionInst.cc.2.cc
    computeSPHSumMassDensityInst.cc.2.cc
    computeSumVoronoiCellMassDensityInst.cc.2.cc
    correctSPHSumMassDensityInst.cc.2.cc
    )

set(SPH_sources3d
    PSPHHydroBaseInst.cc.3.cc
    SPHHydroBaseInst.cc.3.cc
    SolidSPHHydroBaseInst.cc.3.cc
    computeHydrostaticEquilibriumPressure.cc
    computePSPHCorrectionsInst.cc.3.cc
    computeSPHOmegaGradhCorrectionInst.cc.3.cc
    computeSPHSumMassDensityInst.cc.3.cc
    computeSumVoronoiCellMassDensityInst.cc.3.cc
    correctSPHSumMassDensityInst.cc.3.cc
    )

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

INCLUDE_DIRECTORIES(../)

set(spheral_depends BOOST Eigen qhull silo polytope mpi)

blt_add_library(NAME        Spheral_SPH
                SOURCES     ${SPH_sources2d} ${SPH_sources3d}
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
