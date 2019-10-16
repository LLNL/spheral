include(${CMAKE_CURRENT_LIST_DIR}/FindPackageHandleStandardArgs.cmake)
function(find_python_site_packages)
  if(PYTHON_EXECUTABLE)
    set (CMDSTRING "from distutils.sysconfig import get_python_lib; print get_python_lib()")
    execute_process(COMMAND ${PYTHON_EXECUTABLE} "-c" "${CMDSTRING}"
                    RESULT_VARIABLE _site_packages_status 
		    OUTPUT_VARIABLE _site_packages_location
		    ERROR_QUIET 
		    OUTPUT_STRIP_TRAILING_WHITESPACE)
    if (NOT _site_packages_status)
      set(PYTHON_SITE_PACKAGES ${_site_packages_location}
          CACHE STRING "Python site-packages directory lcoation" FORCE)
      message(STATUS "Python site-packages located at ${PYTHON_SITE_PACKAGES}")
    endif()
  endif()
endfunction(find_python_site_packages)
