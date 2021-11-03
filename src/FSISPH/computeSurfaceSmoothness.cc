#include "FSISPH/computeSurfaceSmoothness.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"
#include "Field/FieldList.hh"

namespace Spheral{

template<typename Dimension>
void
computeSurfaceSmoothness(const ConnectivityMap<Dimension>& connectivityMap,
                            const TableKernel<Dimension>& W,
                            const FieldList<Dimension, typename Dimension::Vector>& position,
                            const FieldList<Dimension, typename Dimension::Scalar>& mass,
                            const FieldList<Dimension, typename Dimension::Scalar>& massDensity,
                            const FieldList<Dimension, typename Dimension::SymTensor>& H,
                            FieldList<Dimension, typename Dimension::Vector>& surfaceNormals,
                            FieldList<Dimension, typename Dimension::Scalar>& surfaceFraction,
                            FieldList<Dimension, typename Dimension::Scalar>& surfaceNeighborFraction,
                            FieldList<Dimension, typename Dimension::Scalar>& surfaceSmoothness) {

  // Pre-conditions.
  const auto numNodeLists = massDensity.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(mass.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);

  // The set of interacting node pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

  // Now the pair contributions.
#pragma omp parallel
  {
    int i, j, nodeListi, nodeListj;
    auto surfaceSmoothness_thread = surfaceSmoothness.threadCopy();
    auto surfaceFraction_thread = surfaceFraction.threadCopy();
    auto surfaceNeighborFraction_thread = surfaceNeighborFraction.threadCopy();
#pragma omp for
    for (auto k = 0u; k < npairs; ++k) {
      i = pairs[k].i_node;
      j = pairs[k].j_node;
      nodeListi = pairs[k].i_list;
      nodeListj = pairs[k].j_list;

      if(nodeListi!=nodeListj){
        
        // State for node i
        const auto& ri = position(nodeListi, i);
        const auto  mi = mass(nodeListi, i);
        const auto  rhoi = massDensity(nodeListi, i);
        const auto& Hi = H(nodeListi, i);
        const auto& ni = surfaceNormals(nodeListi,i);
        const auto  voli = mi/rhoi;

        // State for node j
        const auto& rj = position(nodeListj, j);
        const auto  mj = mass(nodeListj, j);
        const auto  rhoj = massDensity(nodeListj, j);
        const auto& Hj = H(nodeListj, j);
        const auto& nj = surfaceNormals(nodeListj,j);
        const auto  volj = mj/rhoj;

        // Kernel weighting and gradient.
        const auto rij = ri - rj;
        const auto etai = Hi*rij;
        const auto etaj = Hj*rij;

        const auto Wi = W.kernelValue(etai.magnitude(), Hi.Determinant());
        const auto Wj = W.kernelValue(etaj.magnitude(), Hj.Determinant());
        const auto Wij = 0.5*(Wi+Wj);

        surfaceNeighborFraction_thread(nodeListi, i) += 1.0;
        surfaceNeighborFraction_thread(nodeListj, j) += 1.0;

        surfaceFraction_thread(nodeListi, i) += volj*Wij;
        surfaceFraction_thread(nodeListj, j) += voli*Wij;

        const auto nirij =  ni.dot(rij);
        const auto njrij = -nj.dot(rij);

        const auto isSubmerged = (nirij > 0.0 or njrij > 0.0);
        if (!isSubmerged){
          const auto ninja = std::max(ni.dot(-nj),0.0);
          surfaceSmoothness_thread(nodeListi, i) += volj * ninja * Wij;
          surfaceSmoothness_thread(nodeListj, j) += voli * ninja * Wij;
        }
      }
    }

#pragma omp critical
    {
      surfaceFraction_thread.threadReduce();
      surfaceNeighborFraction_thread.threadReduce();
      surfaceSmoothness_thread.threadReduce();
    }
  }

  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto n = surfaceFraction[nodeListi]->numInternalElements();
 #pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      const auto  numNeighborsi = connectivityMap.numNeighborsForNode(nodeListi, i);
      surfaceSmoothness(nodeListi,i) /= max(surfaceFraction(nodeListi,i),1.0e-5);
      surfaceNeighborFraction(nodeListi,i) /= numNeighborsi;
     }
    
   }

} // function

} //spheral namespace
