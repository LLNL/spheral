set(${lib_name}_libs libpolytope.a)

include(${SPHERAL_ROOT_DIR}/cmake/spheral/SpheralHandleExt.cmake)
Spheral_Handle_Ext(${lib_name} libpolytope APPLE)
