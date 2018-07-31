set(Strength_sources2d
    BulkModulusPolicyInst.cc.2.cc
    DeviatoricStressPolicyInst.cc.2.cc
    PlasticStrainPolicyInst.cc.2.cc
    RZDeviatoricStressPolicy.cc
    RZPlasticStrainPolicy.cc
    ShearModulusPolicyInst.cc.2.cc
    StrengthSoundSpeedPolicyInst.cc.2.cc
    YieldStrengthPolicyInst.cc.2.cc
    )

set(Strength_sources3d
    BulkModulusPolicyInst.cc.3.cc
    DeviatoricStressPolicyInst.cc.3.cc
    PlasticStrainPolicyInst.cc.3.cc
    ShearModulusPolicyInst.cc.3.cc
    StrengthSoundSpeedPolicyInst.cc.3.cc
    YieldStrengthPolicyInst.cc.3.cc
    )

set(Strength_sources
    SolidFieldNames.cc
    )

set(Strength_headers
    BulkModulusPolicy.hh
    DeviatoricStressPolicy.hh
    PlasticStrainPolicy.hh
    RZDeviatoricStressPolicy.hh
    RZPlasticStrainPolicy.hh
    ShearModulusPolicy.hh
    SolidFieldNames.hh
    StrengthSoundSpeedPolicy.hh
    YieldStrengthPolicy.hh
    effectiveKernelVolume.hh
    )

INCLUDE_DIRECTORIES(../)

set(spheral_depends BOOST Eigen qhull silo polytope mpi)

blt_add_library(NAME        Spheral_Strength
                SOURCES     ${Strength_sources} ${Strength_sources2d} ${Strength_sources3d}
                HEADERS     ${Strength_headers}
                DEPENDS_ON  ${spheral_depends}
                )

install(TARGETS              Spheral_Strength
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Strength_headers}
        DESTINATION include/Strength
        )
