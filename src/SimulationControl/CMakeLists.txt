set(spheralversion "v${CMAKE_PROJECT_VERSION}")
if (EXISTS ${CMAKE_SOURCE_DIR}/.git)
  execute_process(
    COMMAND git --git-dir ${CMAKE_SOURCE_DIR}/.git describe --tags --abbrev=0
    OUTPUT_VARIABLE tagcomm
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  execute_process(
    COMMAND git --git-dir ${CMAKE_SOURCE_DIR}/.git rev-parse --short HEAD
    OUTPUT_VARIABLE githash
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  execute_process(
    COMMAND git --git-dir ${CMAKE_SOURCE_DIR}/.git rev-parse --abbrev-ref HEAD
    OUTPUT_VARIABLE abbref
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  execute_process(
    COMMAND git --git-dir ${CMAKE_SOURCE_DIR}/.git branch --show-current
    OUTPUT_VARIABLE gitbranch
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  set(SPHERAL_VERSION_LIST ${tagcomm} ${githash} ${abbref})
  string(REPLACE ";" " " spheralversion "${SPHERAL_VERSION_LIST}")
endif()
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Spheral_banner.py
  ${CMAKE_CURRENT_BINARY_DIR}/Spheral_banner.py)

# Configure SpheralConfigs.py.in
set(SPHERAL_COMP_CONFIG)
if (SPHERAL_ENABLE_FSISPH)
  list(APPEND SPHERAL_COMP_CONFIG "FSISPH")
endif()
if (SPHERAL_ENABLE_GSPH)
  list(APPEND SPHERAL_COMP_CONFIG "GSPH")
endif()
if (SPHERAL_ENABLE_SVPH)
  list(APPEND SPHERAL_COMP_CONFIG "SVPH")
endif()
set(SPHERAL_HYDRO_IMPORTS "SPHHydros;PSPHHydro;CRKSPHHydros;DEM;SPHUtilities")
if (SPHERAL_ENABLE_FSISPH)
  list(APPEND SPHERAL_HYDRO_IMPORTS "FSISPHHydros;SlideSurfaces")
endif()
if (SPHERAL_ENABLE_GSPH)
  list(APPEND SPHERAL_HYDRO_IMPORTS "GSPHHydros")
endif()
if (SPHERAL_ENABLE_SVPH)
  list(APPEND SPHERAL_HYDRO_IMPORTS "SVPHHydros")
endif()
set(SPHERAL_DIMS 1)
if(ENABLE_2D)
   list(APPEND SPHERAL_DIMS 2)
endif()
if(ENABLE_3D)
   list(APPEND SPHERAL_DIMS 3)
endif()
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/SpheralConfigs.py.in
  ${CMAKE_CURRENT_BINARY_DIR}/SpheralConfigs.py)

spheral_install_python_files(
  Spheral.py
  Spheral1d.py
  Spheral2d.py
  Spheral3d.py
  SphericalSpheral.py
  SpheralRZ.py
  SolidSpheral1d.py
  SolidSpheral2d.py
  SolidSpheral3d.py
  SolidSpheralRZ.py
  spheralDimensions.py
  SpheralTestUtilities.py
  SpheralTimer.py
  SpheralConservation.py
  SpheralController.py
  SpheralOptionParser.py
  SpheralVisitDump.py
  SpheralVoronoiSiloDump.py
  SpheralPolytopeSiloDump.py
  Spheral1dVizDump.py
  SpheralMatplotlib.py
  SpheralTimingParser.py
  findLastRestart.py
  Pnorm.py
  filearraycmp.py
  numericalIntegration.py
  newtonRaphson.py
  bisectFunction.py
  NodeHistory.py
  AverageStrain.py
  EulerianTracerHistory.py
  SpheralHadesDump.py
  siloPointmeshDump.py
  SpheralPointmeshSiloDump.py
  resampleNodeList.py
  overlayNodeList.py
  writeSiloQuadMesh.py
  buildOrderedDict.py
  CaptureStdout.py
  ${CMAKE_CURRENT_BINARY_DIR}/Spheral_banner.py
  ${CMAKE_CURRENT_BINARY_DIR}/SpheralConfigs.py
  )
