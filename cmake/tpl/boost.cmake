if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  set(TOOLSET "gcc")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
  set(TOOLSET "intel-linux")
else()
  set(TOOLSET ${CMAKE_CXX_COMPILER_ID})
endif()

set(BOOST_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${lib_name})
set(BOOST_SRC_DIR ${BOOST_PREFIX}/src/boost)
set(BOOST_BUILD_DIR ${BOOST_PREFIX}/src/boost-build)
set(BOOST_URL "https://sourceforge.net/projects/boost/files/boost/1.69.0/boost_1_69_0.tar.bz2")
set(BOOST_CACHE "${CACHE_DIR}/boost_1_69_0.tar.bz2")

set(${lib_name}_libs libboost_filesystem.so libboost_system.so)

if(${lib_name}_BUILD)

  if (EXISTS ${BOOST_CACHE})
    set(BOOST_URL ${BOOST_CACHE})
  endif()

  ExternalProject_add(${lib_name}
    PREFIX ${BOOST_PREFIX}
    URL ${BOOST_URL} 
    DOWNLOAD_DIR ${CACHE_DIR}
    CONFIGURE_COMMAND env CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} CFLAGS=-fPIC CXXFLAGS=-fPIC FCFLAGS=-fPIC ${BOOST_SRC_DIR}/bootstrap.sh
    --with-toolset=${TOOLSET}
    --without-libraries=atomic,container,coroutine,log,chrono,context,date_time,exception,fiber,graph,graph_parallel,iostreams,locale,math,mpi,program_options,python,random,regex,serialization,test,thread,timer,type_erasure,wave
    --prefix=${${lib_name}_DIR}
    BUILD_IN_SOURCE 1
    BUILD_COMMAND ${BOOST_SRC_DIR}/b2 install
    INSTALL_COMMAND sleep 1

    LOG_DOWNLOAD ${OUT_PROTOCOL_EP}
    LOG_CONFIGURE ${OUT_PROTOCOL_EP}
    LOG_BUILD ${OUT_PROTOCOL_EP}
    LOG_INSTALL ${OUT_PROTOCOL_EP}
  )

endif()
set(BOOST_INSTALL_DIR ${${lib_name}_DIR} PARENT_SCOPE)
