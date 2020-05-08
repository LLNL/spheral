set (_python_packages
  CXXTypes
  Geometry
  PolyClipper
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
  ArtificialViscosity
  Mesh
  Damage
  SolidMaterial
  Strength
  ArtificialConduction
  polytope
  Distributed
  OpenMP
  )

if(NOT ENABLE_MPI)
  list(REMOVE_ITEM _python_packages
    Distributed
    )
endif()

foreach(_python_package ${_python_packages})
  add_subdirectory(${_python_package})
endforeach()

string(REPLACE ";" " " PYTHONPKGS "${_python_packages}")
configure_file(
  SpheralCompiledPackages.py.in
  ${CMAKE_INSTALL_PREFIX}/Spheral/SpheralCompiledPackages.py)

set(Pybind11Wraps_python_targets
  SpheralPickle.py
  PolyhedronFileUtilities.py
  )

install(FILES       ${Pybind11Wraps_python_targets}
        DESTINATION Spheral
        )

install(FILES       Spheral.pth 
        DESTINATION .
        )
