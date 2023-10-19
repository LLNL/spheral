set (_python_packages
  CXXTypes
  CXXTests
  Geometry
  Silo
  DataOutput
  NodeList
  Field
  FieldList
  Kernel
  Neighbor
  Material
  FileIO
  Utilities
  RK
  DataBase
  Boundary
  Physics
  Hydro
  ExternalForce
  Gravity
  Integrator
  NodeGenerators
  FieldOperations
  SPH
  CRKSPH
  SVPH
  DEM
  GSPH
  FSISPH
  ArtificialViscosity
  Mesh
  Damage
  SolidMaterial
  Strength
  Porosity
  ArtificialConduction
  KernelIntegrator
  polytope
  OpenMP
  )

if(ENABLE_MPI)
  list(APPEND _python_packages Distributed)
endif()
if (ENABLE_ANEOS)
  list(APPEND _python_packages ANEOS)
endif()
if (ENABLE_HELMHOLTZ)
  list(APPEND _python_packages Helmholtz)
endif()

foreach(_python_package ${_python_packages})
  add_subdirectory(${_python_package})
endforeach()

string(REPLACE ";" " " PYTHONPKGS "${_python_packages}")
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/SpheralCompiledPackages.py.in
  ${CMAKE_CURRENT_BINARY_DIR}/SpheralCompiledPackages.py)

set(PYB11_python_targets
  SpheralPickle.py
  PolyhedronFileUtilities.py
  ${CMAKE_CURRENT_BINARY_DIR}/SpheralCompiledPackages.py
  )

install(FILES       ${PYB11_python_targets}
        DESTINATION Spheral
        )

install(FILES       Spheral.pth 
        DESTINATION .
        )
