###############################################################################
# Include project's CMake files
###############################################################################
# Basic

# Instantiate and add source files to the project
# _inst_var   : *name* of list variable containing base names to instantiate
#               The file ${_inst_var}Inst.cc.py must exist.
#               If instantiation is disabled, ${_inst_var}.cc will be added
#               to the source files, if it exists.
# _source_var : *name* of list variable to append source files to.
function(instantiate _inst_var _source_var)
  set(_dims 1)
  if(ENABLE_2D)
     list(APPEND _dims 2)
  endif()
  if(ENABLE_3D)
     list(APPEND _dims 3)
  endif()
     
  set(_tmp_source)
  foreach(_inst ${${_inst_var}})
     if(ENABLE_INSTANTIATIONS)
        foreach(_dim ${_dims})
           set(_inst_py ${CMAKE_CURRENT_SOURCE_DIR}/${_inst}Inst.cc.py)
           set(_inst_file ${_inst}Inst${_dim}d.cc)
           # Run the python script to generate the instantiation files
           execute_process(COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/helpers/InstantiationGenerator.py
                                       ${_inst_py} ${_inst_file} ${_dim}
                           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
               
           # Add the instantiation files to the sources.
           list(APPEND _tmp_source ${CMAKE_CURRENT_BINARY_DIR}/${_inst_file})
        endforeach()
     else()
        # If the base source file exists, add it to the source files.
        # These may not exist for header-only implementations.
        if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${_inst}.cc)
           list(APPEND _tmp_source ${_inst}.cc)
        endif()
     endif()
  endforeach()
  set(${_source_var} ${${_source_var}} ${_tmp_source} PARENT_SCOPE)
endfunction()

####################
# Other config defaults
####################

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_definitions("-DDEBUG=1")
else()
    add_definitions("-DDEBUG=0")
endif()

if(ENABLE_CXXONLY)
  add_definitions(-DCXXONLY=1)
endif()
add_definitions(-DUSE_TETGEN=0)
add_definitions(-DUSE_TRIANGLE=0)
add_definitions(-DNOPOLYTOPE=1)
add_definitions(-DSPHERAL1D=1)
if(ENABLE_2D)
  add_definitions(-DSPHERAL2D=1)
endif()
if(ENABLE_3D)
  add_definitions(-DNOR3D=1)
  add_definitions(-DSPHERAL3D=1)
endif()
if(ENABLE_TIMER)
  add_definitions(-DTIMER=1)
endif()

if (ENABLE_MPI)
    add_definitions(-DUSE_MPI=1)
endif()
