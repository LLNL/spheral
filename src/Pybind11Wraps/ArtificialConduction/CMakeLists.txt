set(PYB11_MODULE_NAME ArtificialConduction)

# Include directrory to source file destination
include_directories("${PROJECT_SOURCE_DIR}/src")

# Include directory to location of pybind include header files
include_directories("${PROJECT_SOURCE_DIR}/pybind11/include")

# Create a static library to compile .cc files
set(STATIC_MODULE_NAME ${PYB11_MODULE_NAME}static)
add_library(${STATIC_MODULE_NAME} my_src.cc)
set(MODULE_STATIC_LIB "${PROJECT_BINARY_DIR}/src/lib${STATIC_MODULE_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX}")

# Generate binding, this creates the <module name>.cc file that is compiled to our final .so
# It looks for <PYB11_MODULE_NAME>MOD.py in this directory.
#   The <>MOD.py file lists the header files in PYB11include list and imports the PYBGEN python files 
#   the user creates for each .cc file.
message(STATUS "Generate PYB11 Bindings")
include(UsePYB11Generator)
PYB11_GENERATE_BINDINGS()

# Create SHARED library with python bindings, need to set the name and prefix to "" 
# otherwise python will throw <>init module errors on import.
add_library(pymy_lib SHARED ${PYB11_GENERATED_SOURCE})
set_target_properties(pymy_lib PROPERTIES OUTPUT_NAME ${PYB11_MODULE_NAME})
set_target_properties(pymy_lib PROPERTIES PREFIX "")
target_link_libraries(pymy_lib ${MODULE_STATIC_LIB}) #Include our static libs for implemented code.
