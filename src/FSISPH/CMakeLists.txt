include_directories(.)

set(FSISPH_sources
    FSIFieldNames.cc )

# Add explicit instantiations
list(APPEND FSISPH_sources
    SolidFSISPHHydroBaseInst.cc
    computeFSISPHSumMassDensityInst.cc
    computeHWeightedFSISPHSumMassDensityInst.cc
    computeInterfacePressureCorrectedSumMassDensityInst.cc
    SlideSurfaceInst.cc
    )

set(FSISPH_headers
    SolidFSISPHHydroBase.hh
    SolidFSISPHHydroBaseInline.hh
    computeFSISPHSumMassDensity.hh
    computeHWeightedFSISPHSumMassDensity.hh
    computeInterfacePressureCorrectedSumMassDensity.hh
    FSIFieldNames.hh
    SlideSurface.hh
    SlideSurfaceInline.hh
    )

spheral_install_python_files(
  FSISPHHydros.py
  SlideSurfaces.py
  )

spheral_add_obj_library(FSISPH SPHERAL_OBJ_LIBS)

