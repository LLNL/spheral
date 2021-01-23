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
#include "Utilities/Timer.hh"

// Declare timers
extern Timer TIME_Damage;
extern Timer TIME_ThreePointCoupling;

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
ThreePointDamagedNodeCoupling<Dimension>::
ThreePointDamagedNodeCoupling(const FieldList<Dimension, Vector>& position,
                              const FieldList<Dimension, SymTensor>& H,
                              const FieldList<Dimension, SymTensor>& damage,
                              const TableKernel<Dimension>& W,
                              const ConnectivityMap<Dimension>& connectivity,
                              const bool useIntersectConnectivity,
                              NodePairList& pairs) {

  TIME_Damage.start();
  TIME_ThreePointCoupling.start();
  const auto W0 = W.kernelValue(0.0, 1.0);
  const auto etamax2 = W.kernelExtent() * W.kernelExtent();
  const auto npairs = pairs.size();
  const auto numNodeLists = position.numFields();

  // For each interacting pair we need to compute the effective damage shielding, expressed
  // as the f_couple parameter in the NodePairIdxType.
  // Everyone starts out fully coupled.
#pragma omp parallel for
  for (auto i = 0u; i < npairs; ++i) {
    pairs[i].f_couple = 1.0;
  }

  // This branch computes the same thing in different ways for efficiency.
  // If a small fraction of the points are damaged it's better not to make
  // ConnectivityMap compute the neighbor intersetions for each pair, and
  // instead do the Nneigh^2 work here for the points that are damaged.
  if (not useIntersectConnectivity) {

    // Walk all the points with non-negligible damage.
    Vector b;
    for (auto kl = 0u; kl < numNodeLists; ++kl) {
      const auto nk = position[kl]->numElements();
#pragma omp parallel for private(b)
      for (auto k = 0u; k < nk; ++k) {
        const auto& Dk = damage(kl, k);
        if (Dk.Trace() > 1e-3) {

          // This point has damage, so damage the pair interactions for all pairs that k is in
          // the intersection set of.
          const auto& xk = position(kl, k);
          const auto& Hk = H(kl, k);
          const auto& connectivity_k = connectivity.connectivityForNode(kl, k);
          CHECK(fullConnectivity.size() == numNodeLists);
          for (auto il = 0u; il  < numNodeLists; ++il) {
            const auto ni = connectivity_k[il].size();
            for (auto ii = 0u; ii < ni; ++ii) {
              const auto  i = connectivity_k[il][ii];
              const auto& xi = position(il, i);
              const auto& Hi = H(il, i);
              for (auto jl = il; jl < numNodeLists; ++jl) {
                const auto nj = connectivity_k[jl].size();
                for (auto jj = (jl == il ? ii + 1u : 0u); jj < nj; ++jj) {
                  const auto  j = connectivity_k[jl][jj];
                  const auto& xj = position(jl, j);
                  const auto& Hj = H(jl, j);
                  const auto  xji = xj - xi;
                  if (std::min((Hi*xji).magnitude2(), (Hj*xji).magnitude2()) < etamax2) {
                    // k is in the intersection of (i,j) connectivity, but is it geometrically between (i,j)?
                    if (closestPointOnSegment(xk, xi, xj, b)) {
                      // Yep!
                      NodePairIdxType pair(i, il, j, jl);
                      auto itr = std::find(pairs.begin(), pairs.end(), pair);
                      if (itr != pairs.end()) {                  // Because we don't recompute connectivity during a time step
                        const auto xhatji = xji.unitVector();
#pragma omp atomic
                        itr->f_couple *= std::max(0.0, std::min(1.0, 1.0 - (Dk*xhatji).magnitude() * W.kernelValue((Hk*(b - xk)).magnitude(), 1.0)/W0));
                        // if (k == 48) std::cerr << "      (" << i << " " << j << " " << k << ") " << Dk << " " << itr->f_couple << std::endl;
                      }
                    }
                  }
                  // Damage (k,j)
                  auto pair = std::min(NodePairIdxType(k, kl, j, jl), NodePairIdxType(j, jl, k, kl));
                  auto itr = std::find(pairs.begin(), pairs.end(), pair);
                  CHECK(itr != pairs.end());
                  const auto xhatjk = (xk - xj).unitVector();
#pragma omp atomic
                  itr->f_couple *= std::max(0.0, std::min(1.0, 1.0 - (Dk*xhatjk).magnitude()));
                }
              }
              // Damage (k,i)
              auto pair = std::min(NodePairIdxType(k, kl, i, il), NodePairIdxType(i, il, k, kl));
              auto itr = std::find(pairs.begin(), pairs.end(), pair);
              CHECK(itr != pairs.end());
              const auto xhatik = (xk - xi).unitVector();
#pragma omp atomic
              itr->f_couple *= std::max(0.0, std::min(1.0, 1.0 - (Dk*xhatik).magnitude()));
            }
          }
        }
      }
    }

  } else {

    // This branch is much faster for computing the shielding when a non-negligible amount of material is damaged.
    // The (big) downside is we have to have the ConnectivityMap compute the neighbor set intersections for use
    // here.
    Vector b;
#pragma omp parallel for private(b)
    for (auto kk = 0u; kk < npairs; ++kk) {
      auto& pair = pairs[kk];
      auto& fij = pair.f_couple;
      const auto& xi = position(pair.i_list, pair.i_node);
      const auto& xj = position(pair.j_list, pair.j_node);
      const auto  xji = xj - xi;
      const auto  xhatji = xji.unitVector();

      // Find the common neighbors for this pair.
      const auto intersection_list = connectivity.connectivityIntersectionForNodes(pair.i_list, pair.i_node,
                                                                                   pair.j_list, pair.j_node);
      // const auto& intersection_list = connectivity.intersectionConnectivity(pair);
      for (auto nodeListk = 0u; nodeListk < numNodeLists; ++nodeListk) {
        // if (pair.i_node == 46) {
        //   std::cerr << "intersection list " << pair << " : ";
        //   std::cerr << "   [";
        //   std::copy(intersection_list[nodeListk].begin(), intersection_list[nodeListk].end(), std::ostream_iterator<int>(std::cerr, " "));
        //   std::cerr << "]\n";
        // }
        for (const auto k: intersection_list[nodeListk]) {

          // State for node k
          const auto& xk = position(nodeListk, k);
          const auto& Hk = H(nodeListk, k);
          const auto& Dk = damage(nodeListk, k);

          // if (k == 48 and pair.i_node == 46) {
          //   closestPointOnSegment(xk, xi, xj, b);
          //   std::cerr << " --> " << pair << " " << xi << " " << xk << " " << xj << " " << Dk << " " << closestPointOnSegment(xk, xi, xj, b) << " " << b << std::endl;
          // }

          // We only proceed if the closest point to k on (i,j) is bounded by (i,j)
          if (Dk.Trace() > 1e-3 and closestPointOnSegment(xk, xi, xj, b)) {
            fij *= std::max(0.0, std::min(1.0, 1.0 - (Dk*xhatji).magnitude() * W.kernelValue((Hk*(b - xk)).magnitude(), 1.0)/W0));
          }
        }
      }
      CHECK(fij >= 0.0 and fij <= 1.0);
    }
  }
  
  // std::cerr << "At the end:\n"
  //           << "    (46 47) : " << damage(0,47) << " " << std::find(pairs.begin(), pairs.end(), NodePairIdxType(46, 0, 47, 0))->f_couple << "\n"
  //           << "    (46 48) : " << damage(0,48) << " " << std::find(pairs.begin(), pairs.end(), NodePairIdxType(46, 0, 48, 0))->f_couple << "\n"
  //           << "    (46 49) : " << damage(0,49) << " " << std::find(pairs.begin(), pairs.end(), NodePairIdxType(46, 0, 49, 0))->f_couple << "\n"
  //           << "    (46 50) : " << damage(0,50) << " " << std::find(pairs.begin(), pairs.end(), NodePairIdxType(46, 0, 50, 0))->f_couple << "\n";

  TIME_ThreePointCoupling.stop();
  TIME_Damage.stop();
}

}
