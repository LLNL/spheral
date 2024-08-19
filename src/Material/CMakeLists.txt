include_directories(.)

set(Material_sources
    PhysicalConstants.cc
    )

# Add explicit instantiations
list(APPEND Material_sources
     EquationOfStateInst.cc
     GammaLawGasInst.cc
     PolytropicEquationOfStateInst.cc
     IsothermalEquationOfStateInst.cc
     StiffenedGasInst.cc
     )

set(Material_headers
    EquationOfState.hh
    EquationOfStateInline.hh
    GammaLawGas.hh
    StiffenedGas.hh
    HelmholtzEquationOfState.hh
    IsothermalEquationOfState.hh
    IsothermalEquationOfStateInline.hh
    PhysicalConstants.hh
    PhysicalConstantsInline.hh
    PolytropicEquationOfState.hh
    PolytropicEquationOfStateInline.hh
    )

spheral_install_python_files(
  MaterialUnits.py
  MaterialEquationsOfState.py
  computeHugoniot.py
  )

# Helmholtz
if (ENABLE_HELMHOLTZ)
  message("-- Helmholtz equation of state enabled")
  list(APPEND Material_inst      HelmholtzEquationOfState)
  list(APPEND Material_headers   HelmholtzEquationOfState.hh)
  list(APPEND Material_sources   public_helm.f90 invert_helm.f90)
else()
  message("-- Helmholtz equation of state disabled")
endif()

spheral_add_obj_library(Material SPHERAL_OBJ_LIBS)

