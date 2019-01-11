//------------------------------------------------------------------------------
// Compute the Voronoi cell mass density summation.
//------------------------------------------------------------------------------
#include "computeSumVoronoiCellMassDensity.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"
#include "NodeList/FluidNodeList.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Utilities/safeInv.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

template<typename Dimension>
void
computeSumVoronoiCellMassDensity(const ConnectivityMap<Dimension>& connectivityMap,
                                 const TableKernel<Dimension>& W,
                                 const FieldList<Dimension, typename Dimension::Vector>& position,
                                 const FieldList<Dimension, typename Dimension::Scalar>& mass,
                                 const FieldList<Dimension, typename Dimension::Scalar>& volume,
                                 const FieldList<Dimension, typename Dimension::SymTensor>& H,
                                 FieldList<Dimension, typename Dimension::Scalar>& massDensity) {

  // Pre-conditions.
  const size_t numNodeLists = massDensity.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(mass.size() == numNodeLists);
  REQUIRE(volume.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Zero out the result, and prepare a FieldList to hold the effective volume.
  massDensity = 0.0;
  FieldList<Dimension, Scalar> Veff(FieldStorageType::CopyFields);
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const NodeList<Dimension>& nodeList = massDensity[nodeListi]->nodeList();
    Veff.appendNewField("effective volume", nodeList, 0.0);
  }

  // Walk the FluidNodeLists.
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const FluidNodeList<Dimension>& nodeList = dynamic_cast<const FluidNodeList<Dimension>&>(massDensity[nodeListi]->nodeList());
    const int firstGhostNodei = nodeList.firstGhostNode();
    const Scalar rhoMin = nodeList.rhoMin();
    const Scalar rhoMax = nodeList.rhoMax();

    // Iterate over the nodes in this node list.
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      // Get the state for node i.
      const Vector& ri = position(nodeListi, i);
      const Scalar& mi = mass(nodeListi, i);
      const Scalar& Vi = volume(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();

      // Get the neighbors for this node (in this NodeList).  We use the approximation here
      // that nodes from other NodeLists do not contribute to the density of this one.
      const vector<int>& connectivity = connectivityMap.connectivityForNode(nodeListi, i)[nodeListi];
      for (vector<int>::const_iterator jItr = connectivity.begin();
           jItr != connectivity.end();
           ++jItr) {
        const int j = *jItr;

        // Check if this node pair has already been calculated.
        if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                     nodeListi, j,
                                                     firstGhostNodei)) {
          const Vector& rj = position(nodeListi, j);
          const Scalar& mj = mass(nodeListi, j);
          const Scalar& Vj = volume(nodeListi, j);
          const SymTensor& Hj = H(nodeListi, j);
          const Scalar Hdetj = Hj.Determinant();

          // Kernel weighting and gradient.
          const Vector rij = ri - rj;
          const Scalar etai = (Hi*rij).magnitude();
          const Scalar etaj = (Hj*rij).magnitude();
          const Scalar Wi = W.kernelValue(etai, Hdeti);
          const Scalar Wj = W.kernelValue(etaj, Hdetj);

          // Sum the pair-wise contributions.
          Veff(nodeListi, i) += Vj*Wi;
          massDensity(nodeListi, i) += mj*Wi;

          Veff(nodeListi, j) += Vi*Wj;
          massDensity(nodeListi, j) += mi*Wj;
        }
      }
      // Finalize the density for node i.
      const Scalar W0 = W.kernelValue(0.0, Hdeti);
      massDensity(nodeListi, i) = max(rhoMin, 
                                      min(rhoMax,
                                          (massDensity(nodeListi, i) + mi*W0) * 
                                          safeInv(Veff(nodeListi, i) + Vi*W0)));

      CHECK(massDensity(nodeListi, i) > 0.0);
    }
  }
}

}
