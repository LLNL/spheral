set(Kernel_sources2d
    BSplineKernelInst.cc.2.cc
    ExpInvKernelInst.cc.2.cc
    GaussianKernelInst.cc.2.cc
    HatKernelInst.cc.2.cc
    NBSplineKernelInst.cc.2.cc
    NSincPolynomialKernelInst.cc.2.cc
    PiGaussianKernelInst.cc.2.cc
    QuarticSplineKernelInst.cc.2.cc
    QuinticSplineKernelInst.cc.2.cc
    SincKernelInst.cc.2.cc
    SuperGaussianKernelInst.cc.2.cc
    TableKernelInst.cc.2.cc
    VolumeIntegrationFunctionsInst.cc.2.cc
    W4SplineKernelInst.cc.2.cc
    WendlandC2KernelInst.cc.2.cc
    WendlandC4KernelInst.cc.2.cc
    WendlandC6KernelInst.cc.2.cc
    )

set(Kernel_sources3d
    BSplineKernelInst.cc.3.cc
    ExpInvKernelInst.cc.3.cc
    GaussianKernelInst.cc.3.cc
    HatKernelInst.cc.3.cc
    NBSplineKernelInst.cc.3.cc
    NSincPolynomialKernelInst.cc.3.cc
    PiGaussianKernelInst.cc.3.cc
    QuarticSplineKernelInst.cc.3.cc
    QuinticSplineKernelInst.cc.3.cc
    SincKernelInst.cc.3.cc
    SuperGaussianKernelInst.cc.3.cc
    TableKernelInst.cc.3.cc
    VolumeIntegrationFunctionsInst.cc.3.cc
    W4SplineKernelInst.cc.3.cc
    WendlandC2KernelInst.cc.3.cc
    WendlandC4KernelInst.cc.3.cc
    WendlandC6KernelInst.cc.3.cc
    )

set(Kernel_headers
    BSplineKernel.hh
    BSplineKernelInline.hh
    ExpInvKernel.hh
    ExpInvKernelInline.hh
    GaussianKernel.hh
    GaussianKernelInline.hh
    HKernel.hh
    HKernelInline.hh
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

INCLUDE_DIRECTORIES(../)

set(spheral_depends BOOST Eigen qhull silo polytope mpi)

blt_add_library(NAME        Spheral_Kernel
                SOURCES     ${Kernel_sources2d} ${Kernel_sources3d}
                HEADERS     ${Kernel_headers}
                DEPENDS_ON  ${spheral_depends}
                )

install(TARGETS              Spheral_Kernel
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Kernel_headers}
        DESTINATION include/Kernel
        )
