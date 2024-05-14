include_directories(.)
set(Kernel_inst
    BSplineKernel
    W4SplineKernel
    GaussianKernel
    SuperGaussianKernel
    PiGaussianKernel
    HatKernel
    SincKernel
    NSincPolynomialKernel
    NBSplineKernel
    QuarticSplineKernel
    QuinticSplineKernel
    VolumeIntegrationFunctions
    WendlandC2Kernel
    WendlandC4Kernel
    WendlandC6Kernel
    ExpInvKernel
    )


set(Kernel_sources
    SphericalKernel.cc
    SphericalBiCubicSplineKernel.cc)


instantiate(Kernel_inst Kernel_sources)


set(ENABLE_1D_OLD ${ENABLE_1D})
if (NOT ${ENABLE_1D_OLD})
   set(ENABLE_1D TRUE CACHE BOOL "" FORCE)
endif()

set(ENABLE_2D_OLD ${ENABLE_2D})
if (NOT ${ENABLE_2D_OLD})
   set(ENABLE_2D TRUE CACHE BOOL "" FORCE)
endif()

#higher dimensional things can depend on the lower dimension TableKernel instantiations
set(Kernel_inst TableKernel)
instantiate(Kernel_inst Kernel_sources)

set(ENABLE_1D ${ENABLE_1D_OLD} CACHE BOOL "" FORCE)
set(ENABLE_2D ${ENABLE_2D_OLD} CACHE BOOL "" FORCE)

set(Kernel_headers
    BSplineKernel.hh
    BSplineKernelInline.hh
    ExpInvKernel.hh
    ExpInvKernelInline.hh
    GaussianKernel.hh
    GaussianKernelInline.hh
    HatKernel.hh
    HatKernelInline.hh
    Kernel.hh
    KernelInline.hh
    KernelTraits.hh
    NBSplineKernel.hh
    NBSplineKernelInline.hh
    NSincPolynomialKernel.hh
    NSincPolynomialKernelInline.hh
    PiGaussianKernel.hh
    PiGaussianKernelInline.hh
    QuarticSplineKernel.hh
    QuarticSplineKernelInline.hh
    QuinticSplineKernel.hh
    SincKernel.hh
    SincKernelInline.hh
    SuperGaussianKernel.hh
    SuperGaussianKernelInline.hh
    TableKernel.hh
    TableKernelInline.hh
    SphericalKernel.hh
    SphericalKernelInline.hh
    SphericalBiCubicSplineKernel.hh
    SphericalBiCubicSplineKernelInline.hh
    VolumeIntegrationFunctions.hh
    W4SplineKernel.hh
    W4SplineKernelInline.hh
    WendlandC2Kernel.hh
    WendlandC2KernelInline.hh
    WendlandC4Kernel.hh
    WendlandC4KernelInline.hh
    WendlandC6Kernel.hh
    WendlandC6KernelInline.hh
    )

spheral_add_obj_library(Kernel SPHERAL_OBJ_LIBS)

