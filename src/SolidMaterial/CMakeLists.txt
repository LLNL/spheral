include_directories(.)
set(SolidMaterial_inst
    LinearPolynomialEquationOfState
    GruneisenEquationOfState
    OsborneEquationOfState
    MurnahanEquationOfState
    TillotsonEquationOfState
    StrengthModel
    SteinbergGuinanStrength
    JohnsonCookStrength
    CollinsStrength
    PorousEquationOfState
    PorousStrengthModel
    StrainPorosity
    PhysicsEvolvingMaterialLibrary
   )

set(SolidMaterial_sources )

instantiate(SolidMaterial_inst SolidMaterial_sources)

set(SolidMaterial_headers
    ANEOS.hh
    CollinsStrength.hh
    ConstantStrength.hh
    Geodyn.hh
    GruneisenEquationOfState.hh
    GruneisenEquationOfStateInline.hh
    JohnsonCookStrength.hh
    LinearPolynomialEquationOfState.hh
    LinearPolynomialEquationOfStateInline.hh
    MurnahanEquationOfState.hh
    MurnahanEquationOfStateInline.hh
    NullStrength.hh
    OsborneEquationOfState.hh
    OsborneEquationOfStateInline.hh
    PhysicsEvolvingMaterialLibrary.hh
    PolynomialFit.hh
    PolynomialFitInline.hh
    PorousEquationOfState.hh
    PorousStrengthModel.hh
    SolidEquationOfState.hh
    SolidEquationOfStateInline.hh
    SteinbergGuinanLundStrength.hh
    SteinbergGuinanStrength.hh
    StrainPorosity.hh
    StrainPorosityInline.hh
    StrengthModel.hh
    TillotsonEquationOfState.hh
    TillotsonEquationOfStateInline.hh
    )

set(SolidMaterial_python_targets
    MaterialPropertiesLib.py
    SolidMaterialUnits.py
    SolidMaterialEquationsOfState.py
    ShadowTillotsonEquationOfState.py
    ShadowConstantStrength.py
    )


spheral_add_cxx_library(SolidMaterial)

if(NOT ENABLE_CXXONLY)
  install(FILES       ${SolidMaterial_python_targets}
          DESTINATION Spheral
          )
endif()

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_SolidMaterial)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
