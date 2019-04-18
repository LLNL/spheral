set(SolidMaterial_inst
    LinearPolynomialEquationOfState
    StrengthModel
   )

set(SolidMaterial_sources )

instantiate(SolidMaterial_inst SolidMaterial_sources)

set(SolidMaterial_headers
    LinearPolynomialEquationOfState.hh
    LinearPolynomialEquationOfStateInline.hh
    NullStrength.hh
    SolidEquationOfState.hh
    SolidEquationOfStateInline.hh
    StrengthModel.hh
    )

blt_add_library(NAME        Spheral_SolidMaterial
                SOURCES     ${SolidMaterial_sources}
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
