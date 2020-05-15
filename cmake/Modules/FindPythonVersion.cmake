include(${CMAKE_CURRENT_LIST_DIR}/FindPackageHandleStandardArgs.cmake)
function(find_python_version)
  if(PYTHON_EXECUTABLE)
    set(CMDSTRING "import sys; print sys.version.split(' ')[0].split('rc')[0]")
    execute_process(COMMAND ${PYTHON_EXECUTABLE} "-c" "${CMDSTRING}"
                    RESULT_VARIABLE _python_version_status
		    OUTPUT_VARIABLE _python_version_string
		    ERROR_QUIET
		    OUTPUT_STRIP_TRAILING_WHITESPACE)
    if (NOT _python_version_status)
      # Break full version into major.minor.patch based on how cmake does version checking
      STRING(REGEX REPLACE "^([0-9]+)\\.[0-9]+\\.[0-9]+" "\\1" _major "${_python_version_string}")
      STRING(REGEX REPLACE "^[0-9]+\\.([0-9])+\\.[0-9]+" "\\1" _minor "${_python_version_string}")
      STRING(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9])+" "\\1" _patch "${_python_version_string}")

      # Store the major.minor version for later use
      set(PYTHON_VERSION "${_major}.${_minor}" 
	  CACHE STRING "Python major.minor version number" FORCE)

      # For backward compatibility to cmake's find_package(PythonInterp)
      set(PYTHON_VERSION_STRING ${_python_version_string} 
	  CACHE STRING "Python major.minor.patch version number" FORCE)

      message(STATUS "Python major/minor version number is ${PYTHON_VERSION}")
    endif()
  endif()
endfunction(find_python_version)