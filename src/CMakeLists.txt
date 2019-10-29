# Add package subdirectories.
set(_packages
  Physics
  ArtificialConduction
  Hydro
  Geometry
  NodeList
  Field
  FieldOperations
  CRKSPH
  #Neighbor
  #Utilities
  )

foreach(_package ${_packages})
  add_subdirectory(${_package})
endforeach()

add_subdirectory(Pybind11Wraps/ArtificialConduction)
