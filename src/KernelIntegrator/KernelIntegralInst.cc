//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "KernelIntegral.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class LinearKernel<Dim<1>>;
  template class LinearGrad<Dim<1>>;
  template class LinearKernelVector<Dim<1>>;
  template class LinearKernelStdVector<Dim<1>>;
  template class LinearGradStdVector<Dim<1>>;
  template class BilinearKernelKernel<Dim<1>>;
  template class BilinearGradKernel<Dim<1>>;
  template class BilinearKernelGrad<Dim<1>>;
  template class BilinearGradDotGrad<Dim<1>>;
  template class BilinearGradProdGrad<Dim<1>>;
  template class BilinearSurfaceNormalKernelKernelFromGrad<Dim<1>>;
  template class LinearSurfaceKernel<Dim<1>>;
  template class LinearSurfaceNormalKernel<Dim<1>>;
  template class LinearSurfaceNormalKernelStdVector<Dim<1>>;
  template class BilinearSurfaceKernelKernel<Dim<1>>;
  template class BilinearSurfaceNormalKernelKernel<Dim<1>>;
  template class BilinearSurfaceNormalKernelDotGrad<Dim<1>>;
  template class CellCoefficient<Dim<1>>;
  template class SurfaceNormalCoefficient<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class LinearKernel<Dim<2>>;
  template class LinearGrad<Dim<2>>;
  template class LinearKernelVector<Dim<2>>;
  template class LinearKernelStdVector<Dim<2>>;
  template class LinearGradStdVector<Dim<2>>;
  template class BilinearKernelKernel<Dim<2>>;
  template class BilinearGradKernel<Dim<2>>;
  template class BilinearKernelGrad<Dim<2>>;
  template class BilinearGradDotGrad<Dim<2>>;
  template class BilinearGradProdGrad<Dim<2>>;
  template class BilinearSurfaceNormalKernelKernelFromGrad<Dim<2>>;
  template class LinearSurfaceKernel<Dim<2>>;
  template class LinearSurfaceNormalKernel<Dim<2>>;
  template class LinearSurfaceNormalKernelStdVector<Dim<2>>;
  template class BilinearSurfaceKernelKernel<Dim<2>>;
  template class BilinearSurfaceNormalKernelKernel<Dim<2>>;
  template class BilinearSurfaceNormalKernelDotGrad<Dim<2>>;
  template class CellCoefficient<Dim<2>>;
  template class SurfaceNormalCoefficient<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class LinearKernel<Dim<3>>;
  template class LinearGrad<Dim<3>>;
  template class LinearKernelVector<Dim<3>>;
  template class LinearKernelStdVector<Dim<3>>;
  template class LinearGradStdVector<Dim<3>>;
  template class BilinearKernelKernel<Dim<3>>;
  template class BilinearGradKernel<Dim<3>>;
  template class BilinearKernelGrad<Dim<3>>;
  template class BilinearGradDotGrad<Dim<3>>;
  template class BilinearGradProdGrad<Dim<3>>;
  template class BilinearSurfaceNormalKernelKernelFromGrad<Dim<3>>;
  template class LinearSurfaceKernel<Dim<3>>;
  template class LinearSurfaceNormalKernel<Dim<3>>;
  template class LinearSurfaceNormalKernelStdVector<Dim<3>>;
  template class BilinearSurfaceKernelKernel<Dim<3>>;
  template class BilinearSurfaceNormalKernelKernel<Dim<3>>;
  template class BilinearSurfaceNormalKernelDotGrad<Dim<3>>;
  template class CellCoefficient<Dim<3>>;
  template class SurfaceNormalCoefficient<Dim<3>>;
#endif
}