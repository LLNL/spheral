set(Kernel_inst
    BSplineKernel
    ExpInvKernel
    GaussianKernel
    HatKernel
    NBSplineKernel
    NSincPolynomialKernel
    PiGaussianKernel
    QuarticSplineKernel
    QuinticSplineKernel
    SincKernel
    SuperGaussianKernel
    TableKernel
    VolumeIntegrationFunctions
    W4SplineKernel
    WendlandC2Kernel
    WendlandC4Kernel
    WendlandC6Kernel
   )


set(Kernel_sources )

instantiate(Kernel_inst Kernel_sources)

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

INCLUDE_DIRECTORIES(. ../)

set(spheral_depends BOOST Eigen qhull silo polytope mpi)

blt_add_library(NAME        Spheral_Kernel
                SOURCES     ${Kernel_sources}
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
