include(FindPackageHandleStandardArgs)

find_path(Silo_INCLUDE_DIR
  NAMES silo.h
  PATHS ${Silo_ROOT}/include
  PATH_SUFFIXES Silo
  )
find_library(Silo_LIBRARY
  NAMES libsiloh5.a libsiloh5.so
  PATHS ${Silo_ROOT}/lib
  )

find_package_handle_standard_args(Silo
  FOUND_VAR Silo_FOUND
  REQUIRED_VARS
  Silo_LIBRARY
  Silo_INCLUDE_DIR
  )



