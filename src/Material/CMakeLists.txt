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

spheral_install_python_files(
  MaterialUnits.py
  MaterialEquationsOfState.py
  )

spheral_add_cxx_library(Material)

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_Material)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
