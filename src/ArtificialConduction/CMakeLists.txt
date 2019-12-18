include_directories(.)
set(ArtificialConduction_inst
    ArtificialConduction
    ArtificialConductionPolicy
    )

set(ArtificialConduction_sources )
instantiate(ArtificialConduction_inst ArtificialConduction_sources)

set(ArtificialConduction_headers
    ArtificialConduction.hh
    ArtificialConductionPolicy.hh
    ) 

set(ArtificialConduction_includes
  ${SPHERAL_SOURCE_DIR}
  ${EIGEN_LIB_LOCATION}
  )


blt_add_library(NAME        Spheral_ArtificialConduction
                HEADERS     ${ArtificialConduction_headers}
                SOURCES     ${ArtificialConduction_sources}
                DEPENDS_ON  ${spheral_depends}
                INCLUDES    ${ArtificialConduction_includes}
                SHARED      TRUE
                )

install(TARGETS             Spheral_ArtificialConduction
        EXPORT              ${PROJECT_NAME}-targets
        DESTINATION         lib
        INCLUDES DESTINATION include
        )

install(FILES       ${ArtificialConduction_headers}
        DESTINATION include/ArtificialConduction
        )

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_ArtificialConduction)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
