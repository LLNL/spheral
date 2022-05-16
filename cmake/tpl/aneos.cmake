set(ANEOS_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${lib_name})
set(ANEOS_DIST "M-ANEOS-v1.0.tar.gz")
set(ANEOS_CACHE "${CACHE_DIR}/${ANEOS_DIST}")
set(ANEOS_URL "https://github.com/isale-code/M-ANEOS/releases/download/v1.0beta/${ANEOS_DIST}")
set(ANEOS_SRC_DIR "${ANEOS_PREFIX}/src/aneos/src")
set(ANEOS_DEST_DIR "${${lib_name}_DIR}/lib")
set(ANEOS_INPUT_SRC_DIR "${ANEOS_PREFIX}/src/aneos/input")
set(ANEOS_INPUT_DEST_DIR "${${lib_name}_DIR}/input")
set(ANEOS_NO_INCLUDES On)   # ANEOS does not produce any include header files

#set(${lib_name}_INCLUDES aneos.h)
set(${lib_name}_libs libaneos.a)

if (EXISTS ${ANEOS_CACHE})
  set(ANEOS_URL ${ANEOS_CACHE})
endif()


set(ANEOS_INPUT_DEST_DIR ${ANEOS_INPUT_DEST_DIR} PARENT_SCOPE)
