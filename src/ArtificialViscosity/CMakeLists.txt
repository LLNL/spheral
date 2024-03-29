include_directories(.)
set(ArtificialViscosity_inst
    ArtificialViscosity
    LimitedMonaghanGingoldViscosity
    CullenDehnenViscosity
    FiniteVolumeViscosity
    IncrementCullenMultipliers
    MonaghanGingoldViscosity
    MorrisMonaghanReducingViscosity
    TensorSVPHViscosity
    TensorCRKSPHViscosity
    TensorMonaghanGingoldViscosity
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
    LimitedMonaghanGingoldViscosity.hh
    LimitedMonaghanGingoldViscosityRZ.hh
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

spheral_add_obj_library(ArtificialViscosity SPHERAL_OBJ_LIBS)

