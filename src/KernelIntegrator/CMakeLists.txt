include_directories(.)

set(KernelIntegrator_sources
  GaussLegendreValues.cc
  SymmetricTetrahedralValues.cc
  SymmetricTriangularValues.cc
  )

# Add explicit instantiations
list(APPEND KernelIntegrator_sources
     BilinearIndexInst.cc
     FlatConnectivityInst.cc
     KernelIntegralInst.cc
     KernelIntegratorInst.cc
     ManufacturedSolutionInst.cc
     RKIntegrationKernelInst.cc
     SPHIntegrationKernelInst.cc
     )

set(KernelIntegrator_headers
  BilinearIndex.hh
  BilinearIndexInline.hh
  DataTypeFunctions.hh
  FlatConnectivity.hh
  FlatConnectivityInline.hh
  GaussLegendreValues.hh
  IntegrationCoefficient.hh
  IntegrationCoefficientInline.hh
  IntegrationKernel.hh
  KernelIntegral.hh
  KernelIntegralInline.hh
  KernelIntegrationData.hh
  KernelIntegrator.hh
  KernelIntegratorInline.hh
  ManufacturedSolution.hh
  RKIntegrationKernel.hh
  RKIntegrationKernelInline.hh
  SPHIntegrationKernel.hh
  SymmetricTetrahedralValues.hh
  SymmetricTriangularValues.hh
  )

spheral_add_obj_library(KernelIntegrator SPHERAL_OBJ_LIBS)
