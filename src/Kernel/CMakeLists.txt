include_directories(.)

set(Kernel_sources
    SphericalKernel.cc
    SphericalBiCubicSplineKernel.cc)

# Add explicit instantiations
list(APPEND Kernel_sources
     BSplineKernelInst.cc
     W4SplineKernelInst.cc
     GaussianKernelInst.cc
     SuperGaussianKernelInst.cc
     PiGaussianKernelInst.cc
     HatKernelInst.cc
     SincKernelInst.cc
     NSincPolynomialKernelInst.cc
     NBSplineKernelInst.cc
     QuarticSplineKernelInst.cc
     QuinticSplineKernelInst.cc
     TableKernelInst.cc
     VolumeIntegrationFunctionsInst.cc
     WendlandC2KernelInst.cc
     WendlandC4KernelInst.cc
     WendlandC6KernelInst.cc
     ExpInvKernelInst.cc
     )

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

