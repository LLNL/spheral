include_directories(.)
set(Porosity_inst
    StrainPorosity
    PalphaPorosity
    PorositySolidMassDensityPolicy
    PorousPressurePolicy
    PorousSoundSpeedPolicy
    PorousGammaPolicy
    PorousBulkModulusPolicy
    PorousEntropyPolicy
    PorousShearModulusPolicy
    PorousYieldStrengthPolicy
   )

set(Porosity_sources
    )

instantiate(Porosity_inst Porosity_sources)

set(Porosity_headers
    StrainPorosity.hh
    StrainPorosityInline.hh
    PalphaPorosity.hh
    PorousPressurePolicy.hh
    PorousSoundSpeedPolicy.hh
    PorousGammaPolicy.hh
    PorousBulkModulusPolicy.hh
    PorousEntropyPolicy.hh
    PorousShearModulusPolicy.hh
    PorousYieldStrengthPolicy.hh
    )

spheral_install_python_files(
    ShadowPalphaPorosity.py
     )

spheral_add_obj_library(Porosity SPHERAL_OBJ_LIBS)
