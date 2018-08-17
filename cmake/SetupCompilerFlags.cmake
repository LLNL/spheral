#
# Should not have to do this but cmake doesnt support enabling 
# intel compilers for c++11
#
if (COMPILER_FAMILY_IS_INTEL AND BLT_CXX_STD STREQUAL "c++11")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()


#######################################################################
#                              C++ Flags
#######################################################################

# enable the 'restrict' keyword for disambiguating pointers
blt_append_custom_compiler_flag(
    FLAGS_VAR restrict_flag
    DEFAULT  "-restrict"
    CLANG    " "
    GNU      "-restrict"
    INTEL    "-restrict"
    )

# Set optimization level
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    blt_append_custom_compiler_flag(
        FLAGS_VAR CMAKE_CXX_FLAGS
        DEFAULT  "-O0"
        CLANG    "-O0"
        GNU      "-O0"
        INTEL    "-O0"
        )
else()
    blt_append_custom_compiler_flag(
        FLAGS_VAR CMAKE_CXX_FLAGS
        DEFAULT  "-O3"
        CLANG    "-O3"
        GNU      "-O3"
        INTEL    "-O3"
        )
endif()

# enable debug symbols
# NOTE: We do this on all builds and then strip out symbols for the performance gain
#  because it forces the compiler to byte-align everything
blt_append_custom_compiler_flag(
    FLAGS_VAR CMAKE_CXX_FLAGS
    DEFAULT  "-g"
    )

# May generate Intel(R) SSE3, SSE2, and SSE instructions
blt_append_custom_compiler_flag(
    FLAGS_VAR CMAKE_CXX_FLAGS
    DEFAULT  "-msse3"
    CLANG    " "
    GNU      "-msse3"
    INTEL    "-msse3"
    )

# Floating point flags
blt_append_custom_compiler_flag(
    FLAGS_VAR CMAKE_CXX_FLAGS
    DEFAULT  ""
    CLANG    "-ffp-contract=off"
    GNU      "-ffp-contract=off -ffloat-store"
    INTEL    "-fp-model strict -fp-model source -prec-div -prec-sqrt -no-ftz"
    )

# Warning flags
# ToDo: add -Werror-all back
blt_append_custom_compiler_flag(
    FLAGS_VAR CMAKE_CXX_FLAGS
    DEFAULT  " "
    CLANG    "-Wno-unused-variable -Wno-unused-argument -Wno-unused-value -Wno-unused-parameter -Wno-unused-local-typedef -Wno-sign-compare"
    GNU      " "
    INTEL    " "
    )

#######################################################################
#                           CUDA Flags
#######################################################################

# Set optimization level
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    list (APPEND CUDA_NVCC_FLAGS -G -O0 -g -Xcompiler '-g' -Xcompiler '-O0')
else()
    list (APPEND CUDA_NVCC_FLAGS -O3 -Xcompiler '-O3' )
endif()

