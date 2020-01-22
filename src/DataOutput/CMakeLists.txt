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
                )

install(TARGETS              Spheral_DataOutput
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${DataOutput_headers}
        DESTINATION include/DataOutput
        )
