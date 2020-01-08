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
    TableKernel
    VolumeIntegrationFunctions
    WendlandC2Kernel
    WendlandC4Kernel
    WendlandC6Kernel
    ExpInvKernel
    )


set(Kernel_sources )

instantiate(Kernel_inst Kernel_sources)

set(Kernel_headers
    GaussianKernel.hh
    GaussianKernelInline.hh
    Kernel.hh
    KernelInline.hh
    KernelTraits.hh
    NBSplineKernel.hh
    NBSplineKernelInline.hh
    PiGaussianKernel.hh
    PiGaussianKernelInline.hh
    TableKernel.hh
    TableKernelInline.hh
    VolumeIntegrationFunctions.hh
    )

spheral_add_cxx_library(Kernel)

install(TARGETS              Spheral_Kernel
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Kernel_headers}
        DESTINATION include/Kernel
        )

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_Kernel)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
