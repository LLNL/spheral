include_directories(.)
set(FSISPH_inst
    SolidFSISPHHydroBase
    computeFSISPHSumMassDensity
    computeHWeightedFSISPHSumMassDensity
    SlideSurface
    )


set(FSISPH_sources
    FSIFieldNames.cc )


instantiate(FSISPH_inst FSISPH_sources)

set(FSISPH_headers
    SolidFSISPHHydroBase.hh
    SolidFSISPHHydroBaseInline.hh
    computeFSISPHSumMassDensity.hh
    computeHWeightedFSISPHSumMassDensity.hh
    FSIFieldNames.hh
    SlideSurface.hh
    SlideSurfaceInline.hh
    )

spheral_install_python_files(
  FSISPHHydros.py
  SlideSurfaces.py
  )

spheral_add_obj_library(FSISPH "SPHERAL_OBJ_LIBS")

