set(Kernel_inst
    GaussianKernel
    NBSplineKernel
    PiGaussianKernel
    TableKernelAle
    VolumeIntegrationFunctionsAle
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
