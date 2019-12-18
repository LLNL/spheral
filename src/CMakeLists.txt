# Libraries
include_directories(.)
set(SPHERAL_PYTHON_DEPENDS )
#set(spheral_depends BOOST Eigen qhull silo polytope mpi openmp)
set(spheral_depends BOOST Eigen qhull silo polytope mpi pybind11)

set(SPHERAL_PYTHON_INSTALL ${PROJECT_BINARY_DIR}/lib)

#add_subdirectory(CInterface)
set(_packages 
  ArtificialViscosity
  ArtificialConduction
  Boundary
  CRKSPH
  Damage
  DataBase
  DataOutput
  Distributed
  ExternalForce
  Field
  FieldOperations
  FileIO
  Geometry
  Gravity
  Hydro
  Integrator
  Kernel
  Material
  Mesh
  Neighbor
  NodeGenerators
  NodeList
  Physics
  #PBGWraps
  Pybind11Wraps
  PythonMPIInterfaces
  SimulationControl
  SolidMaterial
  SPH
  Strength
  SVPH
  Utilities
  )

foreach(_package ${_packages})
   add_subdirectory(${_package})
endforeach()

add_subdirectory(Pybind11Wraps/polytope)
add_subdirectory(Pybind11Wraps/Distributed)
add_subdirectory(Pybind11Wraps/OpenMP)
add_subdirectory(Pybind11Wraps/CXXTypes)
add_subdirectory(Pybind11Wraps/Geometry)
add_subdirectory(Pybind11Wraps/PolyClipper)
add_subdirectory(Pybind11Wraps/Silo)
add_subdirectory(Pybind11Wraps/DataOutput)
add_subdirectory(Pybind11Wraps/NodeList)
add_subdirectory(Pybind11Wraps/Field)
add_subdirectory(Pybind11Wraps/Kernel)
add_subdirectory(Pybind11Wraps/Neighbor)
add_subdirectory(Pybind11Wraps/Material)
add_subdirectory(Pybind11Wraps/FileIO)
add_subdirectory(Pybind11Wraps/DataBase)
add_subdirectory(Pybind11Wraps/Boundary)
add_subdirectory(Pybind11Wraps/Physics)
add_subdirectory(Pybind11Wraps/Hydro)
add_subdirectory(Pybind11Wraps/ExternalForce)
add_subdirectory(Pybind11Wraps/Gravity)
add_subdirectory(Pybind11Wraps/Integrator)
add_subdirectory(Pybind11Wraps/Utilities)
add_subdirectory(Pybind11Wraps/NodeGenerators)
add_subdirectory(Pybind11Wraps/FieldOperations)
add_subdirectory(Pybind11Wraps/SPH)
add_subdirectory(Pybind11Wraps/CRKSPH)
add_subdirectory(Pybind11Wraps/ArtificialViscosity)
add_subdirectory(Pybind11Wraps/SVPH)
add_subdirectory(Pybind11Wraps/Mesh)
add_subdirectory(Pybind11Wraps/Damage)
add_subdirectory(Pybind11Wraps/SolidMaterial)
add_subdirectory(Pybind11Wraps/Strength)
add_subdirectory(Pybind11Wraps/ArtificialConduction)
