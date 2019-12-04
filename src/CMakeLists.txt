# Libraries
include_directories(.)
set(SPHERAL_PYTHON_DEPENDS )
#set(spheral_depends BOOST Eigen qhull silo polytope mpi openmp)
set(spheral_depends BOOST Eigen qhull silo polytope )

#add_subdirectory(CInterface)
set(_packages 
  ArtificialViscosity
  Boundary
  CRKSPH
  Damage
  DataBase
  DataOutput
  Distributed
  Field
  FieldOperations
  FileIO
  Geometry
  Hydro
  Integrator
  Kernel
  Material
  Mesh
  Neighbor
  NodeGenerators
  NodeList
  Physics
  SolidMaterial
  SPH
  Strength
  Utilities
  )

foreach(_package ${_packages})
   add_subdirectory(${_package})
endforeach()

add_subdirectory(Pybind11Wraps/polytope)
#add_subdirectory(Pybind11Wraps/Distributed)
add_subdirectory(Pybind11Wraps/OpenMP)
add_subdirectory(Pybind11Wraps/CXXTypes)
add_subdirectory(Pybind11Wraps/Geometry)
add_subdirectory(Pybind11Wraps/PolyClipper)
add_subdirectory(Pybind11Wraps/Silo)
#add_subdirectory(Pybind11Wraps/DataOutput)
add_subdirectory(Pybind11Wraps/NodeList)
add_subdirectory(Pybind11Wraps/Field)
add_subdirectory(Pybind11Wraps/Kernel)
add_subdirectory(Pybind11Wraps/Neighbor)
add_subdirectory(Pybind11Wraps/Material)
#add_subdirectory(Pybind11Wraps/FileIO)
#add_subdirectory(Pybind11Wraps/DataBase)
add_subdirectory(Pybind11Wraps/Boundary)
add_subdirectory(Pybind11Wraps/Physics)
#add_subdirectory(Pybind11Wraps/ArtificialViscosity)
