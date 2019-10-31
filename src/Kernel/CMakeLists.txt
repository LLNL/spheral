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
  BSplineKernel.hh
  BSplineKernelInline.hh
  ExpInvKernel.hh
  ExpInvKernelInline.hh
  GaussianKernel.hh
  GaussianKernelInline.hh
  HatKernel.hh
  HatKernelInline.hh
  HKernel.hh
  HKernelInline.hh
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

set(Kernel_includes
  ${SPHERAL_SOURCE_DIR}
  ${SPHERAL_SOURCE_DIR}/Kernel
  ${EIGEN_LIB_LOCATION}
  )


blt_add_library(NAME        Spheral_Kernel
                HEADERS     ${Kernel_headers}
                SOURCES     ${Kernel_sources}
                DEPENDS_ON  ${spheral_depends}
                INCLUDES    ${Kernel_includes}
                )

install(TARGETS             Spheral_Kernel
        EXPORT              ${PROJECT_NAME}-targets
        DESTINATION         lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Kernel_headers}
        DESTINATION include/Kernel
        )

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_Kernel)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
