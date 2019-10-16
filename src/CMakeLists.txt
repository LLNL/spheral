# Add package subdirectories.
set(_packages
  ArtificialConduction
  Geometry
  Neighbor
  Utilities
  )

foreach(_package ${_packages})
  add_subdirectory(${_package})
endforeach()
