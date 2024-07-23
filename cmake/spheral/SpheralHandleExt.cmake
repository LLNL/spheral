
#----------------------------------------------------------------------------------------
#                                   Spheral_Handle_Ext
#----------------------------------------------------------------------------------------

# ----------------------
# INPUT-OUTPUT VARIABLES
# ----------------------
# <lib_name>      : REQUIRED : name of target TPL
# <libname>       : REQUIRED : library name to change extension
# APPLE           : REQUIRED : flag for Mac OSX

# -----------------------
# OUTPUT VARIABLES TO USE
# -----------------------
# <lib_name>_libs : list of library names with modified extension
#----------------------------------------------------------------------------------------

function(Spheral_Handle_Ext lib_name libname APPLE)

  if(APPLE)
    set(SHARED_EXT "dylib")
    set(${lib_name}_libs ${libname}.${SHARED_EXT})
  endif()

endfunction()
