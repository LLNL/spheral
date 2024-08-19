include_directories(.)

set(CRKSPH_sources
    centerOfMass.cc
    )

if(ENABLE_2D)
  list(APPEND CRKSPH_sources
       CRKSPHHydroBaseRZ.cc
       SolidCRKSPHHydroBaseRZ.cc
       )
endif()

# Add explicit instantiations
list(APPEND CRKSPH_sources
     computeCRKSPHSumMassDensityInst.cc
     computeSolidCRKSPHSumMassDensityInst.cc
     CRKSPHHydroBaseInst.cc
     detectSurfaceInst.cc
     SolidCRKSPHHydroBaseInst.cc
     editMultimaterialSurfaceTopologyInst.cc
     zerothOrderSurfaceCorrectionsInst.cc
     )

set(CRKSPH_headers
    CRKSPHHydroBase.hh
    CRKSPHHydroBaseInline.hh
    CRKSPHHydroBaseRZ.hh
    CRKSPHVariant.hh
    SolidCRKSPHHydroBase.hh
    SolidCRKSPHHydroBaseInline.hh
    SolidCRKSPHHydroBaseRZ.hh
    centerOfMass.hh
    computeCRKSPHEvaluation.hh
    computeCRKSPHIntegral.hh
    computeCRKSPHSumMassDensity.hh
    computeHullSumMassDensity.hh
    computeSolidCRKSPHSumMassDensity.hh
    computeVoronoiCentroids.hh
    detectSurface.hh
    editMultimaterialSurfaceTopology.hh
    volumeSpacing.hh
    zerothOrderSurfaceCorrections.hh
    )

spheral_add_obj_library(CRKSPH SPHERAL_OBJ_LIBS)

spheral_install_python_files(CRKSPHHydros.py)

