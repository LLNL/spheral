# Libraries
include_directories(.)
set(spheral_depends BOOST Eigen qhull silo polytope mpi openmp)

add_subdirectory(CInterface)
set(_packages SPH
              CRKSPH
              ArtificialViscosity
              Boundary
              Damage
              DataBase
              DataOutput
              Distributed
              Integrator
              Kernel
              Strength
              SolidMaterial
              Material
              Mesh
              Neighbor
              Hydro
              Physics
              FileIO
              Field
              FieldOperations
              NodeGenerators
              NodeList
              Geometry
              Utilities
     )

foreach(_package ${_packages})
   add_subdirectory(${_package})
endforeach()
