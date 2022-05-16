set(HDF5_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${lib_name})
set(HDF5_CACHE "${CACHE_DIR}/hdf5-1.10.4.tar.bz2")
set(HDF5_URL "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.4/src/hdf5-1.10.4.tar.bz2")

set(${lib_name}_libs libhdf5.a libhdf5_hl.a)

set(HDF5_INSTALL_DIR ${${lib_name}_DIR} PARENT_SCOPE)
