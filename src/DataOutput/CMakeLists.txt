set(DataOutput_inst
    )

set(DataOutput_sources
  RestartableObject.cc
  RestartRegistrar.cc
  )
instantiate(DataOutput_inst DataOutput_sources)

set(DataOutput_headers
  registerWithRestart.hh
  Restart.hh
  RestartableObject.hh
  RestartHandle.hh
  RestartInline.hh
  RestartRegistrar.hh
  RestartRegistrarInline.hh
  )

set(DataOutput_includes
  ${SPHERAL_SOURCE_DIR}
  ${SPHERAL_SOURCE_DIR}/DataOutput
  ${EIGEN_LIB_LOCATION}
  )


blt_add_library(NAME        Spheral_DataOutput
                HEADERS     ${DataOutput_headers}
                SOURCES     ${DataOutput_sources}
                DEPENDS_ON  ${spheral_depends}
                INCLUDES    ${DataOutput_includes}
                )

install(TARGETS             Spheral_DataOutput
        EXPORT              ${PROJECT_NAME}-targets
        DESTINATION         lib
        INCLUDES DESTINATION include
        )

install(FILES       ${DataOutput_headers}
        DESTINATION include/DataOutput
        )
