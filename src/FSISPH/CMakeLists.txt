include_directories(.)
set(FSISPH_inst
    SolidFSISPHHydroBase
    computeSurfaceSmoothness
    computeSurfaceNormals
    computeFSISPHSumMassDensity
    SlideSurface
    )


set(FSISPH_sources
    FSIFieldNames.cc )


instantiate(FSISPH_inst FSISPH_sources)

set(FSISPH_headers
    SolidFSISPHHydroBase.hh
    computeSurfaceSmoothness.hh
    computeSurfaceNormals.hh
    computeFSISPHSumMassDensity.hh
    FSIFieldNames.hh
    SlideSurface.hh
    )

spheral_install_python_files(
  FSISPHHydros.py
  SlideSurfaces.py
  )

spheral_add_cxx_library(FSISPH)

