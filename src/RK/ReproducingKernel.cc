#include "RK/ReproducingKernel.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
ReproducingKernel<Dimension>::
ReproducingKernel(const TableKernel<Dimension>& W,
                  const RKOrder order):
  mW(W),
  mOrder(order) {
  switch(order) {
  case RKOrder::ZerothOrder:
    mEvaluateBaseKernel = &RKUtilities<Dimension, RKOrder::ZerothOrder>::evaluateBaseKernel;
    mEvaluateBaseGradient(&RKUtilities<Dimension, RKOrder::ZerothOrder>::evaluateBaseGradient);
    mEvaluateBaseHessian(&RKUtilities<Dimension, RKOrder::ZerothOrder>::evaluateBaseHessian);
    mEvaluateKernelAndGradient(&RKUtilities<Dimension, RKOrder::ZerothOrder>::evaluateKernelAndGradient);
    mEvaluateKernel(&RKUtilities<Dimension, RKOrder::ZerothOrder>::evaluateKernel);
    mEvaluateGradient(&RKUtilities<Dimension, RKOrder::ZerothOrder>::evaluateGradient);
    mEvaluateHessian(&RKUtilities<Dimension, RKOrder::ZerothOrder>::evaluateHessian);
    mEvaluateKernelAndGradient(&RKUtilities<Dimension, RKOrder::ZerothOrder>::evaluateKernelAndGradient);
    mEvaluateKernelAndGradients(&RKUtilities<Dimension, RKOrder::ZerothOrder>::evaluateKernelAndGradients);
    mComputeCorrections(&RKUtilities<Dimension, RKOrder::ZerothOrder>::computeCorrections);
    mComputeNormal(&RKUtilities<Dimension, RKOrder::ZerothOrder>::computeNormal);
    break;
  case RKOrder::LinearOrder:
  case RKOrder::QuadraticOrder:
  case RKOrder::CubicOrder:
  case RKOrder::QuarticOrder:
  case RKOrder::QuinticOrder:
  case RKOrder::SexticOrder:
  case RKOrder::SepticOrder:
  default:
    VERIFY2("Unknown order passed to ReproducingKernel", false);
  }
}

}
