//------------------------------------------------------------------------------
// Compute the SPH mass density summation.
//------------------------------------------------------------------------------

#include "DEM/computeParticleRadius.hh"
#include "Field/FieldList.hh"
#include "NodeList/NodeList.hh"

namespace Spheral {

template<typename Dimension>
void
computeParticleRadius(const FieldList<Dimension, typename Dimension::SymTensor>& H,
                            FieldList<Dimension, typename Dimension::Scalar>& particleRadius) {

  // Pre-conditions.
  const auto numNodeLists = particleRadius.size();
  REQUIRE(H.size() == numNodeLists);

  // First the self contribution.
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto n = particleRadius[nodeListi]->numInternalElements();
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      particleRadius(nodeListi, i) =  Dimension::rootnu(std::max(0.0,H(nodeListi, i).Determinant()));
    }
  }

}

}
