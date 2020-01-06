# Libraries
include_directories(.)
set(SPHERAL_PYTHON_DEPENDS )
#set(spheral_depends BOOST Eigen qhull silo polytope mpi openmp)
set(spheral_depends BOOST Eigen qhull silo polytope pybind11)
if(ENABLE_MPI)
  list(APPEND spheral_depends mpi)
endif()

set(SPHERAL_PYTHON_INSTALL ${PROJECT_BINARY_DIR}/lib)

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
  RK
  SolidMaterial
  SPH
  Strength
  SVPH
  Utilities
  )

if(NOT ENABLE_CXXONLY)
  list(APPEND _packages 
    PythonMPIInterfaces
    SimulationControl
    Pybind11Wraps
    )
endif()

if(NOT ENABLE_MPI)
  list(REMOVE_ITEM _packages
    Distributed
    )
endif()

foreach(_package ${_packages})
   add_subdirectory(${_package})
endforeach()


