set(ArtificialViscosity_inst
  ArtificialViscosity
  CRKSPHMonaghanGingoldViscosity
  CullenDehnenViscosity
  FiniteVolumeViscosity
  IncrementCullenMultipliers
  #MonaghanGingoldViscosity
  MorrisMonaghanReducingViscosity
  TensorCRKSPHViscosity
  TensorMonaghanGingoldViscosity
  TensorSVPHViscosity
  VonNeumanViscosity
  )

set(ArtificialViscosity_sources
  #ArtificialViscosityList.cc
  #CheapVonNeumanViscosity.cc
  #CRKSPHMonaghanGingoldViscosityRZ.cc
  #GradPressureViscosity.cc
  #MASHVonNeumanViscosity.cc
  #MGKViscosity.cc
  #MonaghanGingoldKurapatenkoViscosity.cc
  #MonaghanGingoldSumViscosity.cc
  #MonaghanGingoldViscosityGSRZ.cc
  #MonaghanGingoldViscosityRZ.cc
  #NohViscosity.cc
  #RadialViscosity.cc
  )
instantiate(ArtificialViscosity_inst ArtificialViscosity_sources)

set(ArtificialViscosity_headers
  ArtificialViscosity.hh
  ArtificialViscosityInline.hh
  ArtificialViscosityList.hh
  #CheapVonNeumanViscosity.hh
  CRKSPHMonaghanGingoldViscosity.hh
  #CRKSPHMonaghanGingoldViscosityRZ.hh
  CullenDehnenViscosity.hh
  FiniteVolumeViscosity.hh
  #GradPressureViscosity.hh
  IncrementCullenMultipliers.hh
  #MASHVonNeumanViscosity.hh
  #MonaghanGingoldKurapatenkoViscosity.hh
  #MonaghanGingoldSumViscosity.hh
  #MonaghanGingoldViscosity.hh
  #MonaghanGingoldViscosityGSRZ.hh
  #MonaghanGingoldViscosityRZ.hh
  MorrisMonaghanReducingViscosity.hh
  #NohViscosity.hh
  #RadialViscosity.hh
  TensorCRKSPHViscosity.hh
  TensorMonaghanGingoldViscosity.hh
  TensorSVPHViscosity.hh
  VonNeumanViscosity.hh
  )

set(ArtificialViscosity_includes
  ${SPHERAL_SOURCE_DIR}
  ${SPHERAL_SOURCE_DIR}/ArtificialViscosity
  ${EIGEN_LIB_LOCATION}
  )


blt_add_library(NAME        Spheral_ArtificialViscosity
                HEADERS     ${ArtificialViscosity_headers}
                SOURCES     ${ArtificialViscosity_sources}
                DEPENDS_ON  ${spheral_depends}
                INCLUDES    ${ArtificialViscosity_includes}
                )

install(TARGETS             Spheral_ArtificialViscosity
        EXPORT              ${PROJECT_NAME}-targets
        DESTINATION         lib
        INCLUDES DESTINATION include
        )

install(FILES       ${ArtificialViscosity_headers}
        DESTINATION include/ArtificialViscosity
        )

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_ArtificialViscosity)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
