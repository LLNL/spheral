include_directories(.)
set(Damage_inst
    computeFragmentField
   )


set(Damage_sources )

instantiate(Damage_inst Damage_sources)

set(Damage_headers
    computeFragmentField.hh
    )

blt_add_library(NAME        Spheral_Damage
                SOURCES     ${Damage_sources}
                HEADERS     ${Damage_headers}
                DEPENDS_ON  ${spheral_depends}
                SHARED      TRUE
                )

install(TARGETS              Spheral_Damage
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Damage_headers}
        DESTINATION include/Damage
        )

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_Damage)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
