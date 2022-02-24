text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "SPH/SPHHydroBase.cc"
#include "SPH/SPHEvaluateDerivativesImpl.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SPHHydroBase<%(Dim)s>;
  template void SPHHydroBase<%(Dim)s>::evaluateDerivativesImpl(const %(Scalar)s,
                                                               const %(Scalar)s,
                                                               const DataBase<%(Dim)s>&,
                                                               const State<%(Dim)s>&,
                                                               StateDerivatives<%(Dim)s>&,
                                                               const TableKernel<%(Dim)s>&,
                                                               const TableKernel<%(Dim)s>&,
                                                               const TableKernel<%(Dim)s>&) const;
}
"""
