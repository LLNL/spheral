
#----------------------------------------------------------------------------------------
#                                   Spheral_Handle_Ext
#----------------------------------------------------------------------------------------

# -------------------------------------------
# VARIABLES THAT NEED TO BE PREVIOUSLY DEFINED
# -------------------------------------------
# <lib_name>_DIR        : REQUIRED : The installation location of the TPL
# <lib_name>_INCLUDES   : OPTIONAL : Specific includes for the TPL

# ----------------------
# INPUT-OUTPUT VARIABLES
# ----------------------
# <lib_name>     : REQUIRED : The name of the target TPL
# TPL_CMAKE_DIR  : REQUIRED : Directory containing files for each TPL
#                             listing their library names

# -----------------------
# OUTPUT VARIABLES TO USE - Made available implicitly after function call
# -----------------------
# <lib_name> : Exportable target for the TPL
#----------------------------------------------------------------------------------------

function(Spheral_Handle_Ext lib_name libname APPLE)

  if(APPLE)
    set(SHARED_EXT "dylib")
    set(${lib_name}_libs ${libname}.${SHARED_EXT})
  endif()

endfunction()
