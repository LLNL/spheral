set(BOOST_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${lib_name})
set(BOOST_SRC_DIR ${BOOST_PREFIX}/src/boost)
set(BOOST_BUILD_DIR ${BOOST_PREFIX}/src/boost-build)
set(BOOST_URL "https://sourceforge.net/projects/boost/files/boost/1.74.0/boost_1_74_0.tar.bz2/download")
set(BOOST_CACHE "${CACHE_DIR}/boost_1_74_0.tar.bz2")

if (NOT BOOST_HEADER_ONLY)
  set(${lib_name}_libs libboost_filesystem.so libboost_system.so)
endif()

set(BOOST_INSTALL_DIR ${${lib_name}_DIR} PARENT_SCOPE)
