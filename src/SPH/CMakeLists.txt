include_directories(.)
set(SPH_inst
    SPHHydroBase
    PSPHHydroBase
    SolidSPHHydroBase
    SolidFSISPHHydroBase
    FSISPHHydroBase
    computeSPHOmegaGradhCorrection
    computeSPHSumMassDensity
    correctSPHSumMassDensity
    computePSPHCorrections
    computeSumVoronoiCellMassDensity
   )


set(SPH_sources )

if(ENABLE_2D)
  list(APPEND SPH_sources
       SPHHydroBaseGSRZ.cc
       SPHHydroBaseRZ.cc
       SolidSPHHydroBaseRZ.cc
       FSISPHHydroBaseRZ.cc
       SolidFSISPHHydroBaseRZ.cc
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
    SPHHydroBase.hh
    SPHHydroBaseGSRZ.hh
    SPHHydroBaseInline.hh
    SPHHydroBaseRZ.hh
    SolidFSISPHHydroBase.hh
    SolidFSISPHHydroBaseRZ.hh
    FSISPHHydroBase.hh
    FSISPHHydroBaseRZ.hh
    SolidSPHHydroBase.hh
    SolidSPHHydroBaseInline.hh
    SolidSPHHydroBaseRZ.hh
    #SolidSPHHydroBaseRZInline.hh
    computeHydrostaticEquilibriumPressure.hh
    computeSPHOmegaGradhCorrection.hh
    computeSPHSumMassDensity.hh
    computeSumVoronoiCellMassDensity.hh
    correctSPHSumMassDensity.hh
    )

spheral_install_python_files(
  SPHHydros.py
  PSPHHydros.py
  FSISPHHydros.py
  SPHUtilities.py
  FacetedSurfaceASPHHydro.py
  )

spheral_add_cxx_library(SPH)

