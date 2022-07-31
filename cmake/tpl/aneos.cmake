set(${lib_name}_NO_INCLUDES On)
set(${lib_name}_libs libaneos.a)

set(ANEOS_INPUT_DEST_DIR "${${lib_name}_DIR}/input")
set(ANEOS_INPUT_DEST_DIR ${ANEOS_INPUT_DEST_DIR} PARENT_SCOPE)
