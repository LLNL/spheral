
set(ArtificialViscosity_inst
    ArtificialViscosity
    CRKSPHMonaghanGingoldViscosity
    CullenDehnenViscosity
    FiniteVolumeViscosity
    IncrementCullenMultipliers
    MonaghanGingoldViscosity
    MorrisMonaghanReducingViscosity
    TensorCRKSPHViscosity
    TensorMonaghanGingoldViscosity
    TensorSVPHViscosity
    VonNeumanViscosity
   )


set(ArtificialViscosity_sources )

if(ENABLE_2D)
  list(APPEND ArtificialViscosity_sources
       MonaghanGingoldViscosityGSRZ.cc
      )
endif()

instantiate(ArtificialViscosity_inst ArtificialViscosity_sources)

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

blt_add_library(NAME        Spheral_ArtificialViscosity
                SOURCES     ${ArtificialViscosity_sources}
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
