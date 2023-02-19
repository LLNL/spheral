text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "KernelIntegral.cc"

namespace Spheral {
  template class LinearKernel<Dim<%(ndim)s>>;
  template class LinearGrad<Dim<%(ndim)s>>;
  template class LinearKernelVector<Dim<%(ndim)s>>;
  template class LinearKernelStdVector<Dim<%(ndim)s>>;
  template class LinearGradStdVector<Dim<%(ndim)s>>;
  template class BilinearKernelKernel<Dim<%(ndim)s>>;
  template class BilinearGradKernel<Dim<%(ndim)s>>;
  template class BilinearKernelGrad<Dim<%(ndim)s>>;
  template class BilinearGradDotGrad<Dim<%(ndim)s>>;
  template class BilinearGradProdGrad<Dim<%(ndim)s>>;
  template class BilinearSurfaceNormalKernelKernelFromGrad<Dim<%(ndim)s>>;
  template class LinearSurfaceKernel<Dim<%(ndim)s>>;
  template class LinearSurfaceNormalKernel<Dim<%(ndim)s>>;
  template class LinearSurfaceNormalKernelStdVector<Dim<%(ndim)s>>;
  template class BilinearSurfaceKernelKernel<Dim<%(ndim)s>>;
  template class BilinearSurfaceNormalKernelKernel<Dim<%(ndim)s>>;
  template class BilinearSurfaceNormalKernelDotGrad<Dim<%(ndim)s>>;
  template class CellCoefficient<Dim<%(ndim)s>>;
  template class SurfaceNormalCoefficient<Dim<%(ndim)s>>;
}
"""
