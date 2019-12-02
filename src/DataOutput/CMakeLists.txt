include_directories(.)
set(DataOutput_sources
    RestartRegistrar.cc
    )

set(DataOutput_headers
    Restart.hh
    RestartHandle.hh
    RestartInline.hh
    RestartRegistrar.hh
    RestartRegistrarInline.hh
    registerWithRestart.hh
    )

blt_add_library(NAME        Spheral_DataOutput
                SOURCES     ${DataOutput_sources}
                HEADERS     ${DataOutput_headers}
                DEPENDS_ON  ${spheral_depends}
                SHARED      TRUE
                )

install(TARGETS              Spheral_DataOutput
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${DataOutput_headers}
        DESTINATION include/DataOutput
        )

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_DataOutput)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
