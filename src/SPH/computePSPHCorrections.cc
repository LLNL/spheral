//------------------------------------------------------------------------------
// Compute the PSPH grad h corrections due to Hopkins 2013
//------------------------------------------------------------------------------
#include "computePSPHCorrections.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"
#include "Hydro/HydroFieldNames.hh"

namespace Spheral {
namespace SPHSpace {

using namespace std;
using std::min;
using std::max;
using std::abs;

using FieldSpace::Field;
using FieldSpace::FieldList;
using NeighborSpace::ConnectivityMap;
using KernelSpace::TableKernel;
using NodeSpace::NodeList;

template<typename Dimension>
void
computePSPHCorrections(const ConnectivityMap<Dimension>& connectivityMap,
                               const TableKernel<Dimension>& W,
                               const FieldList<Dimension, typename Dimension::Scalar>& mass,
                               const FieldList<Dimension, typename Dimension::Vector>& position,
                               const FieldList<Dimension, typename Dimension::Scalar>& specificThermalEnergy,
                               const FieldList<Dimension, typename Dimension::SymTensor>& H,
                               FieldList<Dimension, typename Dimension::Scalar>& PSPHpbar,
                               FieldList<Dimension, typename Dimension::Scalar>& PSPHcorrection) {

  // Pre-conditions.
  const size_t numNodeLists = PSPHpbar.size();
  REQUIRE(mass.size() == numNodeLists);
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(specificThermalEnergy.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);
  REQUIRE(PSPHcorrection.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Zero out the result.
  PSPHpbar = 0.0;
  PSPHcorrection = 0.0;
  //const double gamma = 5.0/3.0;//NEEDS TO COME FROM THE INTERFACE!
  const double gamma = 1.5;//NEEDS TO COME FROM THE INTERFACE!

  // Walk the FluidNodeLists.
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const NodeList<Dimension>& nodeList = PSPHpbar[nodeListi]->nodeList();

    // Stuff we're going to accumulate.
    Field<Dimension, Scalar> gradsum("sum of the gradient", nodeList);

    // Iterate over the nodes in this node list.
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      // Get the state for node i.
      const Vector& ri = position(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();
      const Scalar invhi = (Hi.Trace()/Dimension::nDim);

      // Neighbors!
      const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
      CHECK(fullConnectivity.size() == numNodeLists);

      
      Scalar gradPbari=0.0;//DpbarDh
      Scalar Nbari=0.0;//Averaged Number of particles
      Scalar gradNbari=0.0;//DnbarDh

      // Iterate over the neighbor NodeLists.
      for (int nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {

        // Iterate over the neighbors for in this NodeList.
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        for (vector<int>::const_iterator jItr = connectivity.begin();
             jItr != connectivity.end();
             ++jItr) {
          const int j = *jItr;

          const Vector& rj = position(nodeListj, j);
          const SymTensor& Hj = H(nodeListj, j);
          const Scalar Hdetj = Hj.Determinant();
          const Scalar& mj = mass(nodeListj, j);
          const Scalar& epsj = specificThermalEnergy(nodeListj, j);

          // Kernel weighting and gradient.
          const Vector rij = ri - rj;
          const Scalar etai = (Hi*rij).magnitude();
          const Scalar etaj = (Hj*rij).magnitude();
          const std::pair<double, double> WWi = W.kernelAndGradValue(etai, Hdeti);
          const Scalar& Wi = WWi.first;
          const Scalar& gWi = WWi.second;
          const std::pair<double, double> WWj = W.kernelAndGradValue(etaj, Hdetj);
          const Scalar& Wj = WWj.first;
          const Scalar& gWj = WWj.second;
          const Scalar xj=(gamma-1)*mj*epsj;
          const Scalar gradh=invhi*(Dimension::nDim*Wi+etai*gWi);
          PSPHpbar(nodeListi, i) += xj*Wi;
          Nbari += Wi;
          gradPbari -= xj*gradh;
          gradNbari -= gradh;

        }
      }
      const Scalar fi=1.0+gradNbari*safeInv(Dimension::nDim*Nbari*invhi);
      PSPHcorrection(nodeListi, i)=gradPbari*safeInv(Dimension::nDim*(gamma-1)*Nbari*invhi*fi);
    }
  }
}

}
}

