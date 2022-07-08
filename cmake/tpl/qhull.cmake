# Setting this to just the release library until we support TPL debug builds on LC
set(${lib_name}_libs libqhullstatic.a)
if(ENABLE_STATIC_TPL)
  set(${lib_name}_libs libqhullstatic.a)
endif()
set(QHULL_BUILD_TYPE Release)
# if (CMAKE_BUILD_TYPE STREQUAL "Debug")
#   set(${lib_name}_libs libqhullstatic_d.a)
# else()
#   set(${lib_name}_libs libqhullstatic.a)
# endif()
# set(QHULL_BUILD_TYPE ${CMAKE_BUILD_TYPE})

