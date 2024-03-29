//---------------------------------Spheral++----------------------------------//
// SPHSmoothingScale
//
// Implements the standard SPH scalar smoothing scale algorithm.
//
// Created by JMO, Wed Sep 14 13:50:49 PDT 2005
//----------------------------------------------------------------------------//
#include "SPHSmoothingScale.hh"
#include "Geometry/Dimension.hh"
#include "Kernel/TableKernel.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Mesh/Mesh.hh"

#include <cmath>
#include <vector>

namespace Spheral {

using std::min;
using std::max;
using std::abs;
using std::vector;

namespace {

//------------------------------------------------------------------------------
// Convert a given number of neighbors to the equivalent 1D "radius" in nodes.
//------------------------------------------------------------------------------
template<typename Dimension> inline double equivalentRadius(const double n);

// 1D
template<>
inline double
equivalentRadius<Dim<1> >(const double n) {
  return 0.5*n;
}

// 2D
template<>
inline double
equivalentRadius<Dim<2> >(const double n) {
  return std::sqrt(n/M_PI);
}

// 3D
template<>
inline double
equivalentRadius<Dim<3> >(const double n) {
  return Dim<3>::rootnu(3.0*n/(4.0*M_PI));
}

}

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
                         const Vector& /*pos*/,
                         const Tensor& DvDx,
                         const Scalar /*hmin*/,
                         const Scalar /*hmax*/,
                         const Scalar /*hminratio*/,
                         const Scalar /*nPerh*/) const {
  return -H/(Dimension::nDim)*DvDx.Trace();
}

//------------------------------------------------------------------------------
// Compute an idealized new H based on the given moments.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::SymTensor
SPHSmoothingScale<Dimension>::
idealSmoothingScale(const SymTensor& H,
                    const Vector& /*pos*/,
                    const Scalar zerothMoment,
                    const Vector& firstMoment,
                    const SymTensor& /*secondMomentEta*/,
                    const SymTensor& /*secondMomentLab*/,
                    const TableKernel<Dimension>& W,
                    const Scalar hmin,
                    const Scalar hmax,
                    const Scalar /*hminratio*/,
                    const Scalar nPerh,
                    const ConnectivityMap<Dimension>& /*connectivityMap*/,
                    const unsigned /*nodeListi*/,
                    const unsigned /*i*/) const {

  // Pre-conditions.
  // REQUIRE2(fuzzyEqual(H.Trace(), Dimension::nDim*H.xx(), 1.0e-5), H << " : " << H.Trace() << " " << Dimension::nDim*H.xx());
  REQUIRE2(zerothMoment >= 0.0, zerothMoment);

  // // Count how many neighbors we currently sample by gather.
  // unsigned n0 = 0;
  // const double kernelExtent = W.kernelExtent();
  // const vector<const NodeList<Dimension>*> nodeLists = connectivityMap.nodeLists();
  // const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
  // const unsigned numNodeLists = nodeLists.size();
  // for (unsigned nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
  //   const Field<Dimension, Vector>& posj = nodeLists[nodeListj]->positions();
  //   for (vector<int>::const_iterator jItr = fullConnectivity[nodeListj].begin();
  //        jItr != fullConnectivity[nodeListj].end();
  //        ++jItr) {
  //     const unsigned j = *jItr;
  //     const double etai = (H*(pos - posj[j])).magnitude();
  //     if (etai <= kernelExtent) ++n0;
  //   }
  // }

  // // We compute an upper-bound for h depending on if we're getting too many neighbors.
  // const double targetRadius = kernelExtent*nPerh;
  // double currentActualRadius = equivalentRadius<Dimension>(double(n0));  // This is radius in number of nodes.
  // const double maxNeighborLimit = 1.25*targetRadius/(currentActualRadius + 1.0e-30);

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
  CHECK2(currentNodesPerSmoothingScale > 0.0, "Bad estimate for nPerh effective from kernel: " << currentNodesPerSmoothingScale);

  // The ratio of the desired to current nodes per smoothing scale.
  const Scalar s = std::min(4.0, std::max(0.25, nPerh/(currentNodesPerSmoothingScale + 1.0e-30)));
  // const Scalar s = min(4.0, max(0.25, min(maxNeighborLimit, nPerh/(currentNodesPerSmoothingScale + 1.0e-30))));
  CHECK(s > 0.0);

  // Now determine how to scale the current H to the desired value.
  Scalar a;
  if (s < 1.0) {
    a = 0.4*(1.0 + s*s);
  } else {
    a = 0.4*(1.0 + 1.0/(s*s*s));
  }
  CHECK(1.0 - a + a*s > 0.0);
  const double hi0 = 1.0/H.xx();
  const double hi1 = std::min(hmax, std::max(hmin, hi0*(1.0 - a + a*s)));

  // Turn the new vote into the SPH tensor and we're done.
  CHECK(hi1 > 0.0);
  return 1.0/hi1 * SymTensor::one;
}

//------------------------------------------------------------------------------
// Determine a new smoothing scale as a replacement for the old, using assorted
// limiting on the ideal H measurement.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::SymTensor
SPHSmoothingScale<Dimension>::
newSmoothingScale(const SymTensor& H,
                  const Vector& pos,
                  const Scalar zerothMoment,
                  const Vector& firstMoment,
                  const SymTensor& secondMomentEta,
                  const SymTensor& secondMomentLab,
                  const TableKernel<Dimension>& W,
                  const Scalar hmin,
                  const Scalar hmax,
                  const Scalar hminratio,
                  const Scalar nPerh,
                  const ConnectivityMap<Dimension>& connectivityMap,
                  const unsigned nodeListi,
                  const unsigned i) const {
  return idealSmoothingScale(H, 
                             pos,
                             zerothMoment,
                             firstMoment,
                             secondMomentEta,
                             secondMomentLab,
                             W,
                             hmin,
                             hmax,
                             hminratio,
                             nPerh,
                             connectivityMap,
                             nodeListi,
                             i);
}

//------------------------------------------------------------------------------
// Use the volumes of tessellation to set the new Hs.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::SymTensor
SPHSmoothingScale<Dimension>::
idealSmoothingScale(const SymTensor& /*H*/,
                    const Mesh<Dimension>& /*mesh*/,
                    const typename Mesh<Dimension>::Zone& zone,
                    const Scalar hmin,
                    const Scalar hmax,
                    const Scalar /*hminratio*/,
                    const Scalar nPerh) const {
  const Scalar vol = zone.volume();
  CHECK(vol > 0.0);
  const Scalar hi = std::max(hmin, std::min(hmax, nPerh * Dimension::rootnu(vol)));
  return 1.0/hi * SymTensor::one;
}

}
