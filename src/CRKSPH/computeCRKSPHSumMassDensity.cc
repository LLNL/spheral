//------------------------------------------------------------------------------
// Compute the CRKSPH mass density summation.
//------------------------------------------------------------------------------
#include "computeCRKSPHSumMassDensity.hh"
#include "interpolateCRKSPH.hh"
#include "computeCRKSPHMoments.hh"
#include "computeCRKSPHCorrections.hh"
#include "CRKSPHCorrectionParams.hh"
#include "NodeList/NodeList.hh"
#include "Hydro/HydroFieldNames.hh"

namespace Spheral {

using std::min;
using std::max;
using std::abs;

template<typename Dimension>
void
computeCRKSPHSumMassDensity(const ConnectivityMap<Dimension>& connectivityMap,
                            const TableKernel<Dimension>& W,
                            const FieldList<Dimension, typename Dimension::Vector>& position,
                            const FieldList<Dimension, typename Dimension::Scalar>& mass,
                            const FieldList<Dimension, typename Dimension::Scalar>& vol,
                            const FieldList<Dimension, typename Dimension::SymTensor>& H,
                            const FieldList<Dimension, int>& voidPoint,
                            FieldList<Dimension, typename Dimension::Scalar>& massDensity) {

  // Pre-conditions.
  const size_t numNodeLists = massDensity.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(mass.size() == numNodeLists);
  REQUIRE(vol.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef typename Dimension::FourthRankTensor FourthRankTensor;
  typedef typename Dimension::FifthRankTensor FifthRankTensor;

  // Some scratch variables.
  Scalar Wi, Wj;
  Vector rij, etai, etaj;
  Vector Bi = Vector::zero, Bj = Vector::zero;
  Tensor Ci = Tensor::zero, Cj = Tensor::zero;

  FieldList<Dimension, Scalar> wsum(FieldStorageType::CopyFields), vol1(FieldStorageType::CopyFields);
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    wsum.appendNewField("weight sum", position[nodeListi]->nodeList(), 0.0);
    vol1.appendNewField("sampled volume", position[nodeListi]->nodeList(), 0.0);
  }

  massDensity = 0.0;
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const FluidNodeList<Dimension>& nodeList = dynamic_cast<const FluidNodeList<Dimension>&>(massDensity[nodeListi]->nodeList());
    const Scalar rhoMin = nodeList.rhoMin();
    const Scalar rhoMax = nodeList.rhoMax();

    // Iterate over the nodes in this node list.
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      // Get the state for node i.
      const Vector& ri = position(nodeListi, i);
      const Scalar mi = mass(nodeListi, i);
      const Scalar Vi = vol(nodeListi, i);
      const Scalar rhoi = massDensity(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();

      // Loop over the neighbors in just point i's NodeList.
      const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
      const vector<int>& connectivity = fullConnectivity[nodeListi];
      const int firstGhostNodej = massDensity[nodeListi]->nodeList().firstGhostNode();
      for (vector<int>::const_iterator jItr = connectivity.begin();
           jItr != connectivity.end();
           ++jItr) {
        const int j = *jItr;

        // Check if this node pair has already been calculated.
        if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                     nodeListi, j,
                                                     firstGhostNodej) and
            voidPoint(nodeListi, j) == 0) {
          const Vector& rj = position(nodeListi, j);
          const Scalar mj = mass(nodeListi, j);
          const Scalar Vj = vol(nodeListi, j);
          const Scalar rhoj = massDensity(nodeListi, j);
          const SymTensor& Hj = H(nodeListi, j);
          const Scalar Hdetj = Hj.Determinant();

          // Kernel weighting and gradient.
          rij = ri - rj;
          etai = Hi*rij;
          etaj = Hj*rij;
          Wi = W.kernelValue(etai.magnitude(), Hdeti);
          Wj = W.kernelValue(etaj.magnitude(), Hdetj);

          // Sum the pair-wise contributions.
          wsum(nodeListi, i) += Vj*Wj;
          wsum(nodeListi, j) += Vi*Wi;
          massDensity(nodeListi, i) += mj * Vj*Wj;
          massDensity(nodeListi, j) += mi * Vi*Wi;
          vol1(nodeListi, i) += Vj * Vj*Wj;
          vol1(nodeListi, j) += Vi * Vi*Wi;
        }
      }
  
      // Finalize the density and volume for node i.
      const Scalar W0 = W.kernelValue(0.0, Hdeti);
      wsum(nodeListi, i) += Vi*W0;
      CHECK(wsum(nodeListi, i) > 0.0);
      vol1(nodeListi, i) = (vol1(nodeListi, i) + Vi*Vi*W0)/wsum(nodeListi, i);
      massDensity(nodeListi, i) = max(max(rhoMin, 0.1*mi*H(nodeListi, i).Determinant()),
                                      min(rhoMax,
                                          (massDensity(nodeListi, i) + mi*Vi*W0)/
                                          (wsum(nodeListi, i)*vol1(nodeListi, i))));
      ENSURE(vol1(nodeListi, i) > 0.0);
      ENSURE(massDensity(nodeListi, i) > 0.0);
    }
  }
}

}
