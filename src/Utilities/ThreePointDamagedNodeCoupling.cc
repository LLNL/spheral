//---------------------------------Spheral++----------------------------------//
// ThreePointDamagedNodeCoupling
//
// A functor class encapsulating how we couple solid nodes in the presence of
// multiple materials and damage.
// 
// This for uses the "three point" formalism, which allows damaged points to
// cut communication between pairs that talk across them.
//
// Unlike other NodeCouping objects, ThreePointDamagedNodeCoupling directly
// sets the pairwise f_couple in the NodePairList, so all work is done in the
// constructor.
//
// Created by JMO, Fri Jan  1 15:18:38 PST 2021
//----------------------------------------------------------------------------//
#include "Utilities/ThreePointDamagedNodeCoupling.hh"
#include "Utilities/pointDistances.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
ThreePointDamagedNodeCoupling<Dimension>::
ThreePointDamagedNodeCoupling(const FieldList<Dimension, Vector>& position,
                              const FieldList<Dimension, SymTensor>& H,
                              const FieldList<Dimension, SymTensor>& damage,
                              const ConnectivityMap<Dimension>& connectivity,
                              const TableKernel<Dimension>& W,
                              NodePairList& pairs) {

  const auto W0 = W.kernelValue(0.0, 1.0);

  // For each interacting pair we need to compute the effective damage shielding, expressed
  // as the f_couple parameter in the NodePairIdxType.
  const auto np = pairs.size();
  const auto numNodeLists = position.size();
  Vector b;
#pragma omp parallel for private(b)
  for (auto kk = 0u; kk < np; ++kk) {
    auto& pair = pairs[kk];
    const auto& xi = position(pair.i_list, pair.i_node);
    const auto& xj = position(pair.j_list, pair.j_node);
    const auto  xji = xj - xi;
    const auto  xhatji = xji.unitVector();

    // Initialize the coupling coefficient
    auto& fij = pair.f_couple;
    fij = 1.0;

    // Find the common neighbors for this pair.
    const auto intersection_list = connectivity.connectivityIntersectionForNodes(pair.i_list, pair.i_node,
                                                                                 pair.j_list, pair.j_node);
    for (auto nodeListk = 0u; nodeListk < numNodeLists; ++nodeListk) {
      for (const auto k: intersection_list[nodeListk]) {

        // State for node k
        const auto& xk = position(nodeListk, k);
        const auto& Hk = H(nodeListk, k);
        const auto& Dk = damage(nodeListk, k);

        // We only proceed if the closest point to k on (i,j) is bounded by (i,j)
        if (closestPointOnSegment(xk, xi, xj, b)) {
          fij *= std::max(0.0, std::min(1.0, 1.0 - (Dk*xhatji).magnitude() * W.kernelValue((Hk*(b - xk)).magnitude(), 1.0)/W0));
        }
      }
    }
    CHECK(fij >= 0.0 and fij <= 1.0);
  }
}

}
