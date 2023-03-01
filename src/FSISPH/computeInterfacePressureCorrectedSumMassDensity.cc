#include "FSISPH/computeInterfacePressureCorrectedSumMassDensity.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"
#include "Utilities/safeInv.hh"

#include <limits.h>

namespace Spheral{

template<typename Dimension>
void
computeInterfacePressureCorrectedSumMassDensity(const ConnectivityMap<Dimension>& connectivityMap,
                                     const TableKernel<Dimension>& W,
                                     const std::vector<int>& sumDensityNodeLists,
                                     const FieldList<Dimension, typename Dimension::Vector>& position,
                                     const FieldList<Dimension, typename Dimension::Scalar>& mass,
                                     const FieldList<Dimension, typename Dimension::SymTensor>& H,
                                     const FieldList<Dimension, typename Dimension::Scalar>& volume,
                                     const FieldList<Dimension, typename Dimension::Scalar>& pressure,
                                     const FieldList<Dimension, typename Dimension::Scalar>& soundSpeed,
                                           FieldList<Dimension, typename Dimension::Scalar>& massDensity) {

  // Pre-conditions.
  const auto numNodeLists = massDensity.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(volume.size() == numNodeLists);
  REQUIRE(pressure.size() == numNodeLists);
  REQUIRE(soundSpeed.size() == numNodeLists);
  REQUIRE(mass.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);

  // Some useful variables.
  const auto W0 = W.kernelValue(0.0, 1.0);
  const auto linEtaMax = 1.0;
  const auto linEtaMin = 0.5;
  const auto tiny = std::numeric_limits<double>::epsilon();

  // The set of interacting node pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto n = massDensity[nodeListi]->numInternalElements();
    if (sumDensityNodeLists[nodeListi]==1){
#pragma omp parallel for
      for (auto i = 0u; i < n; ++i) {
        const auto  mi = mass(nodeListi, i);
        const auto& Hi = H(nodeListi, i);
        const auto  Hdeti = Hi.Determinant();
        massDensity(nodeListi,i) = mi*Hdeti*W0;
      } // for node
    }   // if 
  }     // for nodelist

  // Now the pair contributions.
#pragma omp parallel
  {
    int i, j, nodeListi, nodeListj;
    auto massDensity_thread = massDensity.threadCopy();
#pragma omp for
    for (auto k = 0u; k < npairs; ++k) {
      i = pairs[k].i_node;
      j = pairs[k].j_node;
      nodeListi = pairs[k].i_list;
      nodeListj = pairs[k].j_list;

      const auto  mi = mass(nodeListi, i);
      const auto  mj = mass(nodeListj, j);
      const auto& ri = position(nodeListi, i);
      const auto& rj = position(nodeListj, j);
      const auto& Hi = H(nodeListi, i);
      const auto& Hj = H(nodeListj, j);

      const auto  rij = ri - rj;

      const auto  Hdeti = Hi.Determinant();
      const auto  etai = (Hi*rij).magnitude();
      const auto  Wi = W.kernelValue(etai, Hdeti);

      const auto  Hdetj = Hj.Determinant();
      const auto  etaj = (Hj*rij).magnitude();
      const auto  Wj = W.kernelValue(etaj, Hdetj);
      
      auto meffi = mi;
      auto meffj = mj;

      if(nodeListi!=nodeListj){
        const auto   Vi = volume(nodeListi, i);
        const auto   Vj = volume(nodeListj, j);
        const auto   Pi = pressure(nodeListi, i);
        const auto   Pj = pressure(nodeListj, j);
        const auto   ci = soundSpeed(nodeListi, i);
        const auto   cj = soundSpeed(nodeListj, j);
        const auto rhoi = mi*safeInv(Vi,tiny);
        const auto rhoj = mj*safeInv(Vj,tiny);

        const auto drhodPi = std::min(std::max((Pj-Pi)*safeInv(ci*ci,tiny),rhoi*linEtaMax),-rhoi*linEtaMin);
        const auto drhodPj = std::min(std::max((Pi-Pj)*safeInv(cj*cj,tiny),rhoj*linEtaMax),-rhoi*linEtaMin);
        meffi = Vi*(rhoj+drhodPj);
        meffj = Vj*(rhoi+drhodPi);
      }
      if(sumDensityNodeLists[nodeListi]==1){
        massDensity_thread(nodeListi, i) += meffj*Wi;
      } 
      
      if(sumDensityNodeLists[nodeListj]==1){
        massDensity_thread(nodeListj, j) += meffi*Wj;
      }
    }

#pragma omp critical
    {
      massDensity_thread.threadReduce();
    }
  }


} // function

} //spheral namespace