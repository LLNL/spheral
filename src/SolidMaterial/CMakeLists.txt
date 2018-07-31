set(SolidMaterial_sources2d
    CollinsStrengthInst.cc.2.cc
    GruneisenEquationOfStateInst.cc.2.cc
    JohnsonCookStrengthInst.cc.2.cc
    LinearPolynomialEquationOfStateInst.cc.2.cc
    MurnahanEquationOfStateInst.cc.2.cc
    OsborneEquationOfStateInst.cc.2.cc
    PhysicsEvolvingMaterialLibraryInst.cc.2.cc
    PorousEquationOfStateInst.cc.2.cc
    PorousStrengthModelInst.cc.2.cc
    SteinbergGuinanStrengthInst.cc.2.cc
    StrainPorosityInst.cc.2.cc
    TillotsonEquationOfStateInst.cc.2.cc
    )

set(SolidMaterial_sources3d
    CollinsStrengthInst.cc.3.cc
    GruneisenEquationOfStateInst.cc.3.cc
    JohnsonCookStrengthInst.cc.3.cc
    LinearPolynomialEquationOfStateInst.cc.3.cc
    MurnahanEquationOfStateInst.cc.3.cc
    OsborneEquationOfStateInst.cc.3.cc
    PhysicsEvolvingMaterialLibraryInst.cc.3.cc
    PorousEquationOfStateInst.cc.3.cc
    PorousStrengthModelInst.cc.3.cc
    SteinbergGuinanStrengthInst.cc.3.cc
    StrainPorosityInst.cc.3.cc
    TillotsonEquationOfStateInst.cc.3.cc
    )

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

INCLUDE_DIRECTORIES(../)

set(spheral_depends BOOST Eigen qhull silo polytope mpi)

blt_add_library(NAME        Spheral_SolidMaterial
                SOURCES     ${SolidMaterial_sources2d} ${SolidMaterial_sources3d}
                HEADERS     ${SolidMaterial_headers}
                DEPENDS_ON  ${spheral_depends}
                )

install(TARGETS              Spheral_SolidMaterial
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${SolidMaterial_headers}
        DESTINATION include/SolidMaterial
        )
