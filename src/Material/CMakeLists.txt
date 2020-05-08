include_directories(.)
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
    EquationOfState.hh
    EquationOfStateInline.hh
    PhysicalConstants.hh
    PhysicalConstantsInline.hh
    )

set(Material_python_targets
  MaterialUnits.py
  MaterialEquationsOfState.py
  )

spheral_add_cxx_library(Material)

if(NOT ENABLE_CXXONLY)
  install(FILES       ${Material_python_targets}
          DESTINATION Spheral
          )
endif()

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_Material)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
