set(ArtificialViscosity_sources2d
    ArtificialViscosityInst.cc.2.cc
    CRKSPHMonaghanGingoldViscosityInst.cc.2.cc
    CullenDehnenViscosityInst.cc.2.cc
    FiniteVolumeViscosityInst.cc.2.cc
    IncrementCullenMultipliersInst.cc.2.cc
    MonaghanGingoldViscosityInst.cc.2.cc
    MonaghanGingoldViscosityGSRZ.cc
    MorrisMonaghanReducingViscosityInst.cc.2.cc
    TensorCRKSPHViscosityInst.cc.2.cc
    TensorMonaghanGingoldViscosityInst.cc.2.cc
    TensorSVPHViscosityInst.cc.2.cc
    VonNeumanViscosityInst.cc.2.cc
    )

set(ArtificialViscosity_sources3d
    ArtificialViscosityInst.cc.3.cc
    CRKSPHMonaghanGingoldViscosityInst.cc.3.cc
    CullenDehnenViscosityInst.cc.3.cc
    FiniteVolumeViscosityInst.cc.3.cc
    IncrementCullenMultipliersInst.cc.3.cc
    MonaghanGingoldViscosityInst.cc.3.cc
    MorrisMonaghanReducingViscosityInst.cc.3.cc
    TensorCRKSPHViscosityInst.cc.3.cc
    TensorMonaghanGingoldViscosityInst.cc.3.cc
    TensorSVPHViscosityInst.cc.3.cc
    VonNeumanViscosityInst.cc.3.cc
    )

set(ArtificialViscosity_headers
    ArtificialViscosity.hh
    ArtificialViscosityInline.hh
    ArtificialViscosityList.hh
    CRKSPHMonaghanGingoldViscosity.hh
    CRKSPHMonaghanGingoldViscosityRZ.hh
    CheapVonNeumanViscosity.hh
    CullenDehnenViscosity.hh
    FiniteVolumeViscosity.hh
    GradPressureViscosity.hh
    IncrementCullenMultipliers.hh
    MASHVonNeumanViscosity.hh
    MonaghanGingoldKurapatenkoViscosity.hh
    MonaghanGingoldSumViscosity.hh
    MonaghanGingoldViscosity.hh
    MonaghanGingoldViscosityGSRZ.hh
    MonaghanGingoldViscosityRZ.hh
    MorrisMonaghanReducingViscosity.hh
    NohViscosity.hh
    RadialViscosity.hh
    TensorCRKSPHViscosity.hh
    TensorMonaghanGingoldViscosity.hh
    TensorSVPHViscosity.hh
    VonNeumanViscosity.hh
    )

INCLUDE_DIRECTORIES(../)

set(spheral_depends BOOST Eigen qhull silo polytope mpi)

blt_add_library(NAME        Spheral_ArtificialViscosity
                SOURCES     ${ArtificialViscosity_sources2d} ${ArtificialViscosity_sources3d}
                HEADERS     ${ArtificialViscosity_headers}
                DEPENDS_ON  ${spheral_depends}
                )

install(TARGETS              Spheral_ArtificialViscosity
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${ArtificialViscosity_headers}
        DESTINATION include/ArtificialViscosity
        )
