# Add package subdirectories.
set(_packages
  Physics
  ArtificialConduction
  Hydro
  Geometry
  FileIO
  NodeList
  Field
  FieldOperations
  CRKSPH
  Neighbor
  DataOutput
  DataBase
  Mesh
  #Utilities
  Boundary
  )

foreach(_package ${_packages})
  add_subdirectory(${_package})
endforeach()
message("${SPHERAL_PYTHON_DEPENDS}")


add_subdirectory(Pybind11Wraps/ArtificialConduction)
