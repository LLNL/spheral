include_directories(.)

set(Porosity_sources
    )

# Add explicit instantiations
list(APPEND Porosity_sources
     PorosityModelInst.cc
     StrainPorosityInst.cc
     PalphaPorosityInst.cc
     PorositySolidMassDensityPolicyInst.cc
     )

set(Porosity_headers
    PorosityModel.hh
    PorosityModelInline.hh
    StrainPorosity.hh
    StrainPorosityInline.hh
    PalphaPorosity.hh
    )

spheral_install_python_files(
    ShadowStrainPorosity.py
    ShadowPalphaPorosity.py
    computeHugoniotWithPorosity.py
    PorousEquationOfState.py
    PorousStrengthModel.py
    )

spheral_add_obj_library(Porosity SPHERAL_OBJ_LIBS)
