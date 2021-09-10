#include "FSISPH/smoothSurfaceNormals.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"
#include "Field/FieldList.hh"

namespace Spheral{

template<typename Dimension>
void
smoothSurfaceNormals(const ConnectivityMap<Dimension>& connectivityMap,
                            const TableKernel<Dimension>& W,
                            const FieldList<Dimension, typename Dimension::Vector>& position,
                            const FieldList<Dimension, typename Dimension::Scalar>& mass,
                            const FieldList<Dimension, typename Dimension::Scalar>& massDensity,
                            const FieldList<Dimension, typename Dimension::SymTensor>& H,
                            FieldList<Dimension, typename Dimension::Vector>& interfaceNormals) {

  typedef typename Dimension::Vector Vector;

  // Pre-conditions.
  const auto numNodeLists = massDensity.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(mass.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);

  // The set of interacting node pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

   FieldList<Dimension, Vector> n0(FieldStorageType::CopyFields);
   for (auto nodeListi = 0u; nodeListi != numNodeLists; ++nodeListi) {
     n0.appendNewField("normalsTemp", massDensity[nodeListi]->nodeList(), Vector::zero);
  }

  // Now the pair contributions.
#pragma omp parallel
  {
    int i, j, nodeListi, nodeListj;
    auto interfaceNormals_thread = n0.threadCopy();

#pragma omp for
    for (auto k = 0u; k < npairs; ++k) {
      i = pairs[k].i_node;
      j = pairs[k].j_node;
      nodeListi = pairs[k].i_list;
      nodeListj = pairs[k].j_list;

      if(nodeListi==nodeListj){
        // State for node i
        const auto& ri = position(nodeListi, i);
        const auto  mi = mass(nodeListi, i);
        const auto  rhoi = massDensity(nodeListi, i);
        const auto& Hi = H(nodeListi, i);
        const auto& ni = interfaceNormals(nodeListi,i);

        // State for node j
        const auto& rj = position(nodeListj, j);
        const auto  mj = mass(nodeListj, j);
        const auto  rhoj = massDensity(nodeListj, j);
        const auto& Hj = H(nodeListj, j);
        const auto& nj = interfaceNormals(nodeListj,j);

        // Kernel weighting and gradient.
        const auto rij = ri - rj;
        const auto etai = Hi*rij;
        const auto etaj = Hj*rij;

        const auto Wi = W.kernelValue(etai.magnitude(), Hi.Determinant());
        const auto Wj = W.kernelValue(etaj.magnitude(), Hj.Determinant());

        interfaceNormals_thread(nodeListi, i) += mj/rhoj * nj * Wi;
        interfaceNormals_thread(nodeListj, j) += mi/rhoi * ni * Wj;
      }
    }

#pragma omp critical
    {
      interfaceNormals_thread.threadReduce();
    }
  }

  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
     const auto n = interfaceNormals[nodeListi]->numInternalElements();
 #pragma omp parallel for
     for (auto i = 0u; i < n; ++i) {
       interfaceNormals(nodeListi,i) = n0(nodeListi,i)/max(n0(nodeListi,i).magnitude(),1e-30);
     }
    
   }

} // function

} //spheral namespace
