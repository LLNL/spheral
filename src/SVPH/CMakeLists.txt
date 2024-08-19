include_directories(.)

set(SVPH_sources 
  SVPHFieldNames.cc
  )

# Add explicit instantiations
list(APPEND SVPH_sources
     computeSVPHCorrectionsInst.cc
     computeSVPHCorrectionsOnFacesInst.cc
     computeSumVoronoiCellMassDensityFromFacesInst.cc
     sampleFieldListSVPHInst.cc
     gradientFieldListSVPHInst.cc
     SVPHFacetedHydroBaseInst.cc
     SVPHCorrectionsPolicyInst.cc
     SVPHMassDensityPolicyInst.cc
     CompatibleFaceSpecificThermalEnergyPolicyInst.cc
     MeshIdealHPolicyInst.cc
     SpecificThermalEnergyVolumePolicyInst.cc
     CellPressurePolicyInst.cc
     )

set(SVPH_headers
    CellPressurePolicy.hh
    CompatibleFaceSpecificThermalEnergyPolicy.hh
    MeshIdealHPolicy.hh
    SVPHCorrectionsPolicy.hh
    SVPHFacetedHydroBase.hh
    SVPHFacetedHydroBaseInline.hh
    SVPHFieldNames.hh
    SVPHHydroBase.hh
    SVPHHydroBaseInline.hh
    SVPHMassDensityPolicy.hh
    SpecificThermalEnergyVolumePolicy.hh
    computeSVPHCorrections.hh
    computeSVPHCorrectionsOnFaces.hh
    computeSumVoronoiCellMassDensityFromFaces.hh
    gradientFieldListSVPH.hh
    sampleFieldListSVPH.hh
    )

spheral_install_python_files(SVPHHydros.py)

spheral_add_obj_library(SVPH SPHERAL_OBJ_LIBS)

