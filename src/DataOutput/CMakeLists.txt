include_directories(.)
set(DataOutput_sources
    RestartRegistrar.cc
    #RestartableObject.cc
    )

set(DataOutput_headers
    Restart.hh
    RestartHandle.hh
    RestartInline.hh
    RestartRegistrar.hh
    RestartRegistrarInline.hh
    registerWithRestart.hh
    )

spheral_add_cxx_library(DataOutput)

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
