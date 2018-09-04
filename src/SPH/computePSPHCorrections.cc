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
computePSPHCorrections(const ConnectivityMap<Dimension>& connectivityMap,
                       const TableKernel<Dimension>& W,
                       const FieldList<Dimension, typename Dimension::Scalar>& mass,
                       const FieldList<Dimension, typename Dimension::Vector>& position,
                       const FieldList<Dimension, typename Dimension::Scalar>& specificThermalEnergy,
                       const FieldList<Dimension, typename Dimension::Scalar>& gamma,
                       const FieldList<Dimension, typename Dimension::SymTensor>& H,
                       const bool computeMassDensity,
                       FieldList<Dimension, typename Dimension::Scalar>& PSPHmassDensity,
                       FieldList<Dimension, typename Dimension::Scalar>& PSPHpbar,
                       FieldList<Dimension, typename Dimension::Scalar>& PSPHsoundSpeed,
                       FieldList<Dimension, typename Dimension::Scalar>& PSPHcorrection) {

  // Pre-conditions.
  const size_t numNodeLists = PSPHpbar.size();
  REQUIRE(mass.size() == numNodeLists);
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(specificThermalEnergy.size() == numNodeLists);
  REQUIRE(gamma.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);
  REQUIRE(PSPHmassDensity.size() == numNodeLists);
  REQUIRE(PSPHsoundSpeed.size() == numNodeLists);
  REQUIRE(PSPHcorrection.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Zero out the result.
  if (computeMassDensity) PSPHmassDensity = 0.0;
  PSPHpbar = 0.0;
  PSPHcorrection = 0.0;
  const double tiny = 1.0e-30;

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
      const Scalar mi = mass(nodeListi, i);
      const Scalar invhi = (Hi.Trace()/Dimension::nDim);
      const Scalar epsi = specificThermalEnergy(nodeListi, i);
      const Scalar gammai = gamma(nodeListi, i);
     
      // Self-contribution!!
      Scalar W0  = W(0.0, Hdeti);
      const Scalar xi = (gammai-1.0)*mi*epsi;
      Scalar gradh0 = invhi*(Dimension::nDim*W0);
      PSPHpbar(nodeListi, i) = xi*W0;
      Scalar Nbari = W0;//Averaged Number of particles 
      Scalar gradPbari = -xi*gradh0;//DpbarDh
      Scalar gradNbari = -gradh0;//DnbarDh
      if (computeMassDensity) PSPHmassDensity(nodeListi, i) += mi*W0;

      // Neighbors!
      const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
      CHECK(fullConnectivity.size() == numNodeLists);

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
          const Scalar mj = mass(nodeListj, j);
          const Scalar epsj = specificThermalEnergy(nodeListj, j);
          const Scalar gammaj = gamma(nodeListj, j);

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
          const Scalar xj=(gammaj-1.0)*mj*epsj;
          const Scalar gradh=invhi*(Dimension::nDim*Wi+etai*gWi);
          if (computeMassDensity) PSPHmassDensity(nodeListi, i) += (nodeListi == nodeListj ? mj : mi)*Wj;
          PSPHpbar(nodeListi, i) += xj*Wi;
          Nbari += Wi;
          gradPbari -= xj*gradh;
          gradNbari -= gradh;

        }
      }
      //const Scalar fi=1.0+gradNbari*safeInv(Dimension::nDim*Nbari*invhi);
      //PSPHcorrection(nodeListi, i)=gradPbari*safeInv(Dimension::nDim*(gammai-1.0)*Nbari*invhi*fi);
      const Scalar fi=1.0+gradNbari/max(Dimension::nDim*Nbari*invhi,tiny);
      PSPHcorrection(nodeListi, i)=gradPbari/max(Dimension::nDim*(gammai-1.0)*Nbari*invhi*fi,tiny);
      CHECK2((gammai-1.0)*epsi >= 0.0, i << " " << gammai << " " << epsi);
      PSPHsoundSpeed(nodeListi, i) = sqrt(std::max(0.0, gammai*(gammai - 1.0)*epsi));
    }
  }
}

}
