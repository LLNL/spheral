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

macro(PYB11_GENERATE_BINDINGS PYB11_MODULE_NAME)
  set(PYB11_SOURCE "${PYB11_MODULE_NAME}MOD.py")
  set(PYB11_GENERATED_SOURCE "Spheral${PYB11_MODULE_NAME}.cc")

  # List directories in which spheral .py files can be found.
  set(PYTHON_ENV 
      ${EXTRA_PYB11_SPHERAL_ENV_VARS}
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/${PYB11_MODULE_NAME}:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/polytope:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/Distributed:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/OpenMP:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/CXXTypes:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/Geometry:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/PolyClipper:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/Silo:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/DataOutput:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/NodeList:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/Field:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/FieldList:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/Kernel:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/Neighbor:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/Material:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/FileIO:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/DataBase:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/Boundary:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/Physics:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/Hydro:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/ExternalForce:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/Gravity:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/Integrator:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/Utilities:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/NodeGenerators:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/FieldOperations:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/SPH:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/RK:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/CRKSPH:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/ArtificialViscosity:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/SVPH:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/Mesh:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/Damage:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/SolidMaterial:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/Strength:"
      "${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/ArtificialConduction:"
      "${CMAKE_BINARY_DIR}/src/SimulationControl")

  # Format list into a one line shell friendly format
  STRING(REPLACE ";" "<->" PYTHON_ENV_STR ${PYTHON_ENV})


  # Generating python stamp files to detect changes in PYB11_SOURCE and
  # its included modules
  if(EXISTS ${PYTHON_EXE})
    # Python must exist to generate at config time
    if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/${PYB11_MODULE_NAME}_stamp.cmake")
      # Generate stamp files at config time
      execute_process(COMMAND env PYTHONPATH=\"${PYTHON_ENV_STR}\"
                      ${PYTHON_EXE} ${SPHERAL_ROOT_DIR}/src/helpers/moduleCheck.py 
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
                    ${PYTHON_EXE} ${SPHERAL_ROOT_DIR}/src/helpers/moduleCheck.py
                    ${PYB11_MODULE_NAME}
                    ${CMAKE_CURRENT_SOURCE_DIR}/${PYB11_SOURCE}
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                    DEPENDS python-install
                    )

  # Generate the actual pyb11 module cpp source file
  add_custom_command(OUTPUT ${PYB11_GENERATED_SOURCE}
                     COMMAND env PYTHONPATH=\"${PYTHON_ENV_STR}\"
                     ${PYTHON_EXE} -c
                     'from PYB11Generator import * \; 
                     import ${PYB11_MODULE_NAME}MOD \;
                     PYB11generateModule(${PYB11_MODULE_NAME}MOD, \"Spheral${PYB11_MODULE_NAME}\") '
                     DEPENDS ${PYB11_MODULE_NAME}_stamp ${${PYB11_MODULE_NAME}_DEPENDS} ${PYB11_SOURCE}
                     )

endmacro()
