include_directories(.)

set(SPH_sources )

if(ENABLE_1D)
  list(APPEND SPH_sources
       SphericalSPHHydroBase.cc
       SolidSphericalSPHHydroBase.cc
      )
endif()

if(ENABLE_2D)
  list(APPEND SPH_sources
       SPHHydroBaseRZ.cc
       SolidSPHHydroBaseRZ.cc
      )
endif()

if(ENABLE_3D)
  list(APPEND SPH_sources
       computeHydrostaticEquilibriumPressure.cc
      )
endif()

# Add explicit instantiations
list(APPEND SPH_sources
    SPHHydroBaseInst.cc
    PSPHHydroBaseInst.cc
    SolidSPHHydroBaseInst.cc
    computeSPHOmegaGradhCorrectionInst.cc
    computeSPHSumMassDensityInst.cc
    correctSPHSumMassDensityInst.cc
    computePSPHCorrectionsInst.cc
    computeSumVoronoiCellMassDensityInst.cc
    )

if(ENABLE_1D)
  list(APPEND SPH_sources
       SphericalSPHInst.cc
      )
endif()

set(SPH_headers
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

spheral_add_obj_library(SPH SPHERAL_OBJ_LIBS)

