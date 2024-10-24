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

  set(SPHERAL_VERSION_LIST ${tagcomm} ${githash} ${abbref})
  string(REPLACE ";" " " spheralversion "${SPHERAL_VERSION_LIST}")
endif()
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Spheral_banner.py
  ${CMAKE_CURRENT_BINARY_DIR}/Spheral_banner.py)

set(_dims 1)
if(ENABLE_2D)
   list(APPEND _dims 2)
endif()
if(ENABLE_3D)
   list(APPEND _dims 3)
endif()

string(REPLACE ";" " " DIMS "${_dims}")
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/spheralDimensions.py.in
  ${CMAKE_CURRENT_BINARY_DIR}/spheralDimensions.py)


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
  ${CMAKE_CURRENT_BINARY_DIR}/spheralDimensions.py
  )
