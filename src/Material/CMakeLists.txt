set(Material_inst
  EquationOfState
  GammaLawGas
  PolytropicEquationOfState
  IsothermalEquationOfState
)

set(Material_sources
  PhysicalConstants.cc
  )
instantiate(Material_inst Material_sources)

set(Material_headers
  )

set(Material_includes
  ${SPHERAL_SOURCE_DIR}
  ${SPHERAL_SOURCE_DIR}/Material
  ${EIGEN_LIB_LOCATION}
  )


blt_add_library(NAME        Spheral_Material
                HEADERS     ${Material_headers}
                SOURCES     ${Material_sources}
                DEPENDS_ON  ${spheral_depends}
                INCLUDES    ${Material_includes}
                )

install(TARGETS             Spheral_Material
        EXPORT              ${PROJECT_NAME}-targets
        DESTINATION         lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Material_headers}
        DESTINATION include/Material
        )

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_Material)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
