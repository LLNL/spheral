
set(ArtificialViscosity_inst
    ArtificialViscosity
    CRKSPHMonaghanGingoldViscosity
    MonaghanGingoldViscosity
    TensorCRKSPHViscosity
    TensorMonaghanGingoldViscosity
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
    MonaghanGingoldViscosity.hh
    MonaghanGingoldViscosityGSRZ.hh
    MonaghanGingoldViscosityRZ.hh
    TensorCRKSPHViscosity.hh
    TensorMonaghanGingoldViscosity.hh
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
