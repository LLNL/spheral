set(Material_sources2d
    EquationOfStateInst.cc.2.cc
    GammaLawGasInst.cc.2.cc
    HelmholtzEquationOfStateInst.cc.2.cc
    IsothermalEquationOfStateInst.cc.2.cc
    PolytropicEquationOfStateInst.cc.2.cc
    )

set(Material_sources3d
    EquationOfStateInst.cc.3.cc
    GammaLawGasInst.cc.3.cc
    HelmholtzEquationOfStateInst.cc.3.cc
    IsothermalEquationOfStateInst.cc.3.cc
    PolytropicEquationOfStateInst.cc.3.cc
    )

set(Material_sources
    PhysicalConstants.cc
    )

set(Material_headers
    EquationOfState.hh
    EquationOfStateInline.hh
    GammaLawGas.hh
    HelmholtzEquationOfState.hh
    IsothermalEquationOfState.hh
    IsothermalEquationOfStateInline.hh
    PhysicalConstants.hh
    PhysicalConstantsInline.hh
    PolytropicEquationOfState.hh
    PolytropicEquationOfStateInline.hh
    )

INCLUDE_DIRECTORIES(../)

set(spheral_depends BOOST Eigen qhull silo polytope mpi)

blt_add_library(NAME        Spheral_Material
                SOURCES     ${Material_sources} ${Material_sources2d} ${Material_sources3d}
                HEADERS     ${Material_headers}
                DEPENDS_ON  ${spheral_depends}
                )

install(TARGETS              Spheral_Material
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Material_headers}
        DESTINATION include/Material
        )
