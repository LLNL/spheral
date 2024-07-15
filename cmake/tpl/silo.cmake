set(${lib_name}_libs libsiloh5.a)

include(${SPHERAL_ROOT_DIR}/cmake/spheral/SpheralHandleExt.cmake)
Spheral_Handle_Ext(${lib_name} libsiloh5 APPLE)
