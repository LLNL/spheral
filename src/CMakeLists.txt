# Add package subdirectories.
set(_packages
  ArtificialConduction
  ArtificialViscosity
  Boundary
  CRKSPH
  #CXXTests
  Damage
  DataBase
  DataOutput
  ExternalForce
  Field
  FieldOperations
  FileIO
  FVCRKH
  Geometry
  Gravity
  Hydro
  Integrator
  Kernel
  Material
  Mesh
  #MHD
  Neighbor
  #NodeGenerators
  NodeList
  Physics
  Utilities
  )

foreach(_package ${_packages})
  add_subdirectory(${_package})
endforeach()
message("${SPHERAL_PYTHON_DEPENDS}")


add_subdirectory(Pybind11Wraps/ArtificialConduction)
add_subdirectory(Pybind11Wraps/CXXTypes)
add_subdirectory(Pybind11Wraps/Geometry)
