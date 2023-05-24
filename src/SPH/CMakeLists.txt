include_directories(.)
set(SPH_inst
    SPHHydroBase
    PSPHHydroBase
    SolidSPHHydroBase
    computeSPHOmegaGradhCorrection
    computeSPHSumMassDensity
    correctSPHSumMassDensity
    computePSPHCorrections
    computeSumVoronoiCellMassDensity
    )


set(SPH_sources )

if(ENABLE_1D)
  list(APPEND SPH_sources
       SphericalSPHInst.cc
       SphericalSPHHydroBase.cc
       SolidSphericalSPHHydroBase.cc
      )
endif()

if(ENABLE_2D)
  list(APPEND SPH_sources
       SPHHydroBaseRZ.cc
       SolidSPHHydroBaseRZ.cc
       EffectiveRadiusPolicy.cc
      )
endif()

if(ENABLE_3D)
  list(APPEND SPH_sources
       computeHydrostaticEquilibriumPressure.cc
      )
endif()

instantiate(SPH_inst SPH_sources)

set(SPH_headers
    OmegaGradhPolicy.hh
    PSPHHydroBase.hh
    PSPHHydroBaseInline.hh
    SPHHydroBase.hh
    SPHHydroBaseInline.hh
    SPHHydroBaseRZ.hh
    SphericalSPHHydroBase.hh
    SolidSPHHydroBase.hh
    SolidSPHHydroBaseInline.hh
    SolidSPHHydroBaseRZ.hh
    SolidSphericalSPHHydroBase.hh
    computeHydrostaticEquilibriumPressure.hh
    computePSPHCorrections.hh
    computeSPHOmegaGradhCorrection.hh
    computeSPHSumMassDensity.hh
    computeSumVoronoiCellMassDensity.hh
    correctSPHSumMassDensity.hh
    )

spheral_install_python_files(
    SPHHydros.py
    PSPHHydros.py
    SPHUtilities.py
    FacetedSurfaceASPHHydro.py
    )

spheral_add_obj_library(SPH)

