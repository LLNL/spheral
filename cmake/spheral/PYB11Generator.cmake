#-----------------------------------------------------------------------------------
# PYB11_GENERATE_BINDINGS
#     - Generates the Python bindings for each module in the list
#     - Generates python stamp files for listing python dependency file to help
#       detecting changes in the pyb11 python files at build time
#
# Variables that must be set before calling PYB11_GENERATE_BINDINGS:
#   PYB11_MODULE_NAME
#     - Pyb11 module to be generated
#   PYTHON_EXE
#     - Python executable
#   <PYB11_MODULE_NAME>_DEPENDS
#     - Any target dependencies that must be built before generating the module
#
# To get the names of the generated source
# use: ${PYB11_GENERATED_SOURCE}
#-----------------------------------------------------------------------------------

macro(PYB11_GENERATE_BINDINGS)
  set(PYB11_SOURCE "${PYB11_MODULE_NAME}MOD.py")
  set(PYB11_GENERATED_SOURCE "${PYB11_MODULE_NAME}.cc")

  # List directories in which spheral .py files can be found.
  set(PYTHON_ENV 
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/${PYB11_MODULE_NAME}:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/polytope:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/Distributed:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/OpenMP:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/CXXTypes:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/Geometry:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/PolyClipper:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/Silo:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/DataOutput:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/NodeList:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/Field:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/FieldList:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/Kernel:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/Neighbor:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/Material:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/FileIO:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/DataBase:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/Boundary:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/Physics:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/Hydro:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/ExternalForce:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/Gravity:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/Integrator:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/Utilities:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/NodeGenerators:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/FieldOperations:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/SPH:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/RK:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/CRKSPH:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/ArtificialViscosity:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/SVPH:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/Mesh:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/Damage:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/SolidMaterial:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/Strength:"
      "${PROJECT_SOURCE_DIR}/src/Pybind11Wraps/ArtificialConduction:"
      "${PROJECT_SOURCE_DIR}/src/SimulationControl")

  # Format list into a one line shell friendly format
  STRING(REPLACE ";" "<->" PYTHON_ENV_STR ${PYTHON_ENV})


  # Generating python stamp files to detect changes in PYB11_SOURCE and
  # its included modules
  if(EXISTS ${PYTHON_EXE})
    # Python must exist to generate at config time
    if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/${PYB11_MODULE_NAME}_stamp.cmake")
      # Generate stamp files at config time
      execute_process(COMMAND env PYTHONPATH=\"${PYTHON_ENV_STR}\"
                      ${PYTHON_EXE} ${PROJECT_SOURCE_DIR}/src/helpers/moduleCheck.py 
                      ${PYB11_MODULE_NAME}
                      ${CMAKE_CURRENT_SOURCE_DIR}/${PYB11_SOURCE}
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                      )
    endif()

    # Include list of dependent python files
    include(${CMAKE_CURRENT_BINARY_DIR}/${PYB11_MODULE_NAME}_stamp.cmake)
  endif()

  # Always regenerate the stamp files at build time. Any change in the stamp file
  # will trigger a rebuild of the target pyb11 module
  add_custom_target(${PYB11_MODULE_NAME}_stamp ALL
                    COMMAND env PYTHONPATH=\"${PYTHON_ENV_STR}\"
                    ${PYTHON_EXE} ${PROJECT_SOURCE_DIR}/src/helpers/moduleCheck.py
                    ${PYB11_MODULE_NAME}
                    ${CMAKE_CURRENT_SOURCE_DIR}/${PYB11_SOURCE}
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                    DEPENDS python-install
                    )

  # Generate the actual pyb11 module cpp source file
  add_custom_command(OUTPUT Spheral${PYB11_GENERATED_SOURCE}
                     COMMAND env PYTHONPATH=\"${PYTHON_ENV_STR}\"
                     ${PYTHON_EXE} -c
                     'from PYB11Generator import * \; 
                     import ${PYB11_MODULE_NAME}MOD \;
                     PYB11generateModule(${PYB11_MODULE_NAME}MOD, \"Spheral${PYB11_MODULE_NAME}\") '
                     DEPENDS ${PYB11_MODULE_NAME}_stamp ${${PYB11_MODULE_NAME}_DEPENDS} ${PYB11_SOURCE}
                     )

endmacro()
