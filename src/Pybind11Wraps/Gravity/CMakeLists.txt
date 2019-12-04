set(PYB11_MODULE_NAME Gravity)

# Include directrory to source file destination
include_directories("${PROJECT_SOURCE_DIR}/src/${PYB11_MODULE_NAME}")
include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${PYBIND11_INCLUDE_DIRS}")
include_directories("${EIGEN_INCLUDE_DIRS}")
include_directories("${POLYTOPE_INCLUDE_DIRS}")

# Include directory to location of pybind include header files
include_directories("${PROJECT_SOURCE_DIR}/src/thirdPartyLibs/pybind11/include")

# Create a static library to compile .cc files
#set(STATIC_MODULE_NAME ${PYB11_MODULE_NAME}static)
#add_library(${STATIC_MODULE_NAME} my_src.cc)
set(MODULE_SHARED_LIB 
  "${PROJECT_BINARY_DIR}/lib/libSpheral_Geometry${CMAKE_SHARED_LIBRARY_SUFFIX}"
  "${PROJECT_BINARY_DIR}/lib/libSpheral_NodeList${CMAKE_SHARED_LIBRARY_SUFFIX}"
  "${PROJECT_BINARY_DIR}/lib/libSpheral_Field${CMAKE_SHARED_LIBRARY_SUFFIX}"
  "${PROJECT_BINARY_DIR}/lib/libSpheral_FieldOperations${CMAKE_SHARED_LIBRARY_SUFFIX}"
  "${PROJECT_BINARY_DIR}/lib/libSpheral_Kernel${CMAKE_SHARED_LIBRARY_SUFFIX}"
  "${PROJECT_BINARY_DIR}/lib/libSpheral_Material${CMAKE_SHARED_LIBRARY_SUFFIX}"
  "${PROJECT_BINARY_DIR}/lib/libSpheral_Neighbor${CMAKE_SHARED_LIBRARY_SUFFIX}"
  "${PROJECT_BINARY_DIR}/lib/libSpheral_DataBase${CMAKE_SHARED_LIBRARY_SUFFIX}"
  "${PROJECT_BINARY_DIR}/lib/libSpheral_Boundary${CMAKE_SHARED_LIBRARY_SUFFIX}"
  "${PROJECT_BINARY_DIR}/lib/libSpheral_Physics${CMAKE_SHARED_LIBRARY_SUFFIX}"
  "${PROJECT_BINARY_DIR}/lib/libSpheral_ArtificialViscosity${CMAKE_SHARED_LIBRARY_SUFFIX}"
  "${PROJECT_BINARY_DIR}/lib/libSpheral_Hydro${CMAKE_SHARED_LIBRARY_SUFFIX}"
  #"${PROJECT_BINARY_DIR}/lib/libSpheral_ExternalForce${CMAKE_SHARED_LIBRARY_SUFFIX}"
  "${PROJECT_BINARY_DIR}/lib/libSpheral_Gravity${CMAKE_SHARED_LIBRARY_SUFFIX}"
  "${PROJECT_BINARY_DIR}/lib/libSpheral_Integrator${CMAKE_SHARED_LIBRARY_SUFFIX}"
  "${PROJECT_BINARY_DIR}/lib/libSpheral_FileIO${CMAKE_SHARED_LIBRARY_SUFFIX}"
  "${PROJECT_BINARY_DIR}/lib/libSpheral_DataOutput${CMAKE_SHARED_LIBRARY_SUFFIX}"
  "${PROJECT_BINARY_DIR}/lib/libSpheral_Utilities${CMAKE_SHARED_LIBRARY_SUFFIX}"
  "${PROJECT_BINARY_DIR}/lib/libSpheral_NodeGenerators${CMAKE_SHARED_LIBRARY_SUFFIX}"
  "${PROJECT_BINARY_DIR}/lib/libSpheral_SPH${CMAKE_SHARED_LIBRARY_SUFFIX}"
  "${PROJECT_BINARY_DIR}/lib/libSpheral_CRKSPH${CMAKE_SHARED_LIBRARY_SUFFIX}"
  #"${PROJECT_BINARY_DIR}/lib/libSpheral_FVCRKH${CMAKE_SHARED_LIBRARY_SUFFIX}"
  #"${PROJECT_BINARY_DIR}/lib/libSpheral_SVPH${CMAKE_SHARED_LIBRARY_SUFFIX}"
  "${PROJECT_BINARY_DIR}/lib/libSpheral_Mesh${CMAKE_SHARED_LIBRARY_SUFFIX}"
  "${PROJECT_BINARY_DIR}/lib/libSpheral_Damage${CMAKE_SHARED_LIBRARY_SUFFIX}"
  "${PROJECT_BINARY_DIR}/lib/libSpheral_SolidMaterial${CMAKE_SHARED_LIBRARY_SUFFIX}"
  "${PROJECT_BINARY_DIR}/lib/libSpheral_Strength${CMAKE_SHARED_LIBRARY_SUFFIX}"
  #"${PROJECT_BINARY_DIR}/lib/libSpheral_ArtificialConduction${CMAKE_SHARED_LIBRARY_SUFFIX}"
  "${PROJECT_BINARY_DIR}/lib/libSpheral_Distributed${CMAKE_SHARED_LIBRARY_SUFFIX}"

  ${QHULL_LIBRARY}
  )

# Generate binding, this creates the <module name>.cc file that is compiled to our final .so
# It looks for <PYB11_MODULE_NAME>MOD.py in this directory.
#   The <>MOD.py file lists the header files in PYB11include list and imports the PYBGEN python files 
#   the user creates for each .cc file.
message(STATUS "Generate ${PYB11_MODULE_NAME} PYB11 Bindings")
include(UsePYB11Generator)
PYB11_GENERATE_BINDINGS()

# Create SHARED library with python bindings, need to set the name and prefix to "" 
# otherwise python will throw <>init module errors on import.
set(PYTHON_MODULE_LIB "Spheral${PYB11_MODULE_NAME}")
message("${PYTHON_MODULE_LIB}")
add_library(${PYTHON_MODULE_LIB} SHARED Spheral${PYB11_GENERATED_SOURCE})
add_dependencies(${PYTHON_MODULE_LIB} ${SPHERAL_PYTHON_DEPENDS})
set_target_properties(${PYTHON_MODULE_LIB} PROPERTIES OUTPUT_NAME Spheral${PYB11_MODULE_NAME})
set_target_properties(${PYTHON_MODULE_LIB} PROPERTIES PREFIX "")
target_link_libraries(${PYTHON_MODULE_LIB} ${MODULE_SHARED_LIB}) #Include our static libs for implemented code.
