//---------------------------------Spheral++----------------------------------//
// SPHSmoothingScale
//
// Implements the standard SPH scalar smoothing scale algorithm.
//
// Created by JMO, Wed Sep 14 13:50:49 PDT 2005
//----------------------------------------------------------------------------//
#include "SPHSmoothingScale.hh"
#include "Kernel/TableKernel.hh"

namespace Spheral {
namespace NodeSpace {

using KernelSpace::TableKernel;
using std::min;
using std::max;
using std::abs;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHSmoothingScale<Dimension>::
SPHSmoothingScale():
  SmoothingScaleBase<Dimension>() {
}

//------------------------------------------------------------------------------
// Copy constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHSmoothingScale<Dimension>::
SPHSmoothingScale(const SPHSmoothingScale<Dimension>& rhs):
  SmoothingScaleBase<Dimension>(rhs) {
}

//------------------------------------------------------------------------------
// Assignment.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHSmoothingScale<Dimension>&
SPHSmoothingScale<Dimension>::
operator=(const SPHSmoothingScale& rhs) {
  SmoothingScaleBase<Dimension>::operator=(rhs);
  return *this;
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHSmoothingScale<Dimension>::
~SPHSmoothingScale<Dimension>() {
}

//------------------------------------------------------------------------------
// Time derivative of the smoothing scale.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::SymTensor
SPHSmoothingScale<Dimension>::
smoothingScaleDerivative(const SymTensor& H,
                         const Tensor& DvDx,
                         const Scalar hmin,
                         const Scalar hmax,
                         const Scalar hminratio,
                         const Scalar nPerh,
                         const int maxNumNeighbors) const {
  return -H/(Dimension::nDim)*DvDx.Trace();
}

//------------------------------------------------------------------------------
// Compute an idealized new H based on the given moments.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::SymTensor
SPHSmoothingScale<Dimension>::
idealSmoothingScale(const SymTensor& H,
                    const Scalar zerothMoment,
                    const SymTensor& secondMoment,
                    const int numNeighbors,
                    const TableKernel<Dimension>& W,
                    const Scalar hmin,
                    const Scalar hmax,
                    const Scalar hminratio,
                    const Scalar nPerh,
                    const int maxNumNeighbors) const {

  // Pre-conditions.
  REQUIRE2(fuzzyEqual(H.Trace(), Dimension::nDim*H.xx(), 1.0e-5), H << " : " << H.Trace() << " " << Dimension::nDim*H.xx());
  REQUIRE2(zerothMoment >= 0.0, zerothMoment);
  REQUIRE2(numNeighbors >= 0, numNeighbors);

  // Determine the current effective number of nodes per smoothing scale.
  Scalar currentNodesPerSmoothingScale;
  if (fuzzyEqual(zerothMoment, 0.0)) {

    // This node appears to be in isolation.  It's not clear what to do here --
    // for now we'll punt and say you should double the current smoothing scale.
    currentNodesPerSmoothingScale = 0.5*nPerh;

  } else {

    // Query from the kernel the equivalent nodes per smoothing scale
    // for the observed sum.
    currentNodesPerSmoothingScale = W.equivalentNodesPerSmoothingScale(zerothMoment);
  }
  CHECK(currentNodesPerSmoothingScale > 0.0);

  // Determine if we should limit the new h by the total number of neighbors.
  // number of neighbors.
  const Scalar maxNeighborLimit = Dimension::rootnu(double(maxNumNeighbors)/double(max(1, numNeighbors)));

  // The ratio of the desired to current nodes per smoothing scale.
  const Scalar s = min(maxNeighborLimit, min(4.0, max(0.25, nPerh/(currentNodesPerSmoothingScale + 1.0e-30))));
  CHECK2(s > 0.0, "Bad scaling:  " << s << " " << maxNumNeighbors << " " << numNeighbors << " " << maxNeighborLimit << " " << nPerh << " " << currentNodesPerSmoothingScale);

  // Now determine how to scale the current H to the desired value.
  Scalar a;
  if (s < 1.0) {
    a = 0.4*(1.0 + s*s);
  } else {
    a = 0.4*(1.0 + 1.0/(s*s*s));
  }
  CHECK(1.0 - a + a*s > 0.0);
  const Scalar hi0 = 1.0/H.xx();
  const Scalar hi = min(hmax, max(hmin, hi0*(1.0 - a + a*s)));
  CHECK(hi > 0.0);

  return 1.0/hi * SymTensor::one;
}

//------------------------------------------------------------------------------
// Determine a new smoothing scale as a replacement for the old, using assorted
// limiting on the ideal H measurement.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::SymTensor
SPHSmoothingScale<Dimension>::
newSmoothingScale(const SymTensor& H,
                  const Scalar zerothMoment,
                  const SymTensor& secondMoment,
                  const int numNeighbors,
                  const TableKernel<Dimension>& W,
                  const Scalar hmin,
                  const Scalar hmax,
                  const Scalar hminratio,
                  const Scalar nPerh,
                  const int maxNumNeighbors) const {
  return idealSmoothingScale(H, 
                             zerothMoment,
                             secondMoment,
                             numNeighbors,
                             W,
                             hmin,
                             hmax,
                             hminratio,
                             nPerh,
                             maxNumNeighbors);
}

}
}

