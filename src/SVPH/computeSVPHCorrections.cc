//------------------------------------------------------------------------------
// Modified ideas from CSPH using the SVPH normalized kernel estimates.
//------------------------------------------------------------------------------
#include "computeSVPHCorrections.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Geometry/outerProduct.hh"
#include "Geometry/innerProduct.hh"

namespace Spheral {
namespace SVPHSpace {

using namespace std;
using std::min;
using std::max;
using std::abs;

using FieldSpace::Field;
using FieldSpace::FieldList;
using NeighborSpace::ConnectivityMap;
using KernelSpace::TableKernel;
using NodeSpace::NodeList;
using Geometry::outerProduct;
using Geometry::innerProduct;

template<typename Dimension>
void
computeSVPHCorrections(const ConnectivityMap<Dimension>& connectivityMap,
                       const TableKernel<Dimension>& W,
                       const FieldList<Dimension, typename Dimension::Scalar>& volume,
                       const FieldList<Dimension, typename Dimension::Vector>& position,
                       const FieldList<Dimension, typename Dimension::SymTensor>& H,
                       FieldList<Dimension, typename Dimension::Vector>& B,
                       FieldList<Dimension, typename Dimension::Tensor>& gradB) {

  // Pre-conditions.
  const size_t numNodeLists = B.size();
  REQUIRE(volume.size() == numNodeLists);
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);
  REQUIRE(gradB.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  // Zero out the result.
  B = Vector::zero;
  gradB = Tensor::zero;

  // Walk the FluidNodeLists.
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const NodeList<Dimension>& nodeList = B[nodeListi]->nodeList();
    const int firstGhostNodei = nodeList.firstGhostNode();

    // We can derive everything in terms of the first and second moments 
    // of the local positions.
    Field<Dimension, Vector> m1("first moment", nodeList);
    Field<Dimension, SymTensor> m2("second moment", nodeList);
    Field<Dimension, Tensor> gradm1("gradient of the first moment", nodeList);
    Field<Dimension, ThirdRankTensor> gradm2("gradient of the second moment", nodeList);

    // Iterate over the nodes in this node list.
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      // Get the state for node i.
      const Scalar Vi = volume(nodeListi, i);
      const Vector& ri = position(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);

      // Self contribution.
      gradm1(i) += Vi*W(0.0, 1.0);

      // Neighbors!
      const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
      CHECK(fullConnectivity.size() == numNodeLists);
      const vector<int>& connectivity = fullConnectivity[nodeListi];

      // Iterate over the neighbors for in this NodeList.
      for (vector<int>::const_iterator jItr = connectivity.begin();
           jItr != connectivity.end();
           ++jItr) {
        const int j = *jItr;

        // // Check if this node pair has already been calculated.
        // if (connectivityMap.calculatePairInteraction(nodeListi, i, 
        //                                              nodeListi, j,
        //                                              firstGhostNodei)) {

          // State of node j.
          const Scalar Vj = volume(nodeListi, j);
          const Vector& rj = position(nodeListi, j);
          const SymTensor& Hj = H(nodeListi, j);

          // Kernel weighting and gradient.
          const Vector rij = ri - rj;
          const Vector etai = Hi*rij;
          const Vector etaj = Hj*rij;
          const std::pair<double, double> WWi = W.kernelAndGradValue(etai.magnitude(), 1.0);
          const Scalar& Wi = WWi.first;
          const Vector gradWi = -(Hi*etai.unitVector())*WWi.second;
          const std::pair<double, double> WWj = W.kernelAndGradValue(etaj.magnitude(), 1.0);
          const Scalar& Wj = WWj.first;
          const Vector gradWj = (Hj*etaj.unitVector())*WWj.second;

          // First moment. 
          m1(i) += Vj*Wj * rij;
          gradm1(i) += Vj*( rij*gradWj + Tensor::one*Wj);

          // Second moment.
          const SymTensor thpt = rij.selfdyad();
          m2(i) += Vj*Wj * thpt;
          // gradm2(i) += Vj*(outerProduct<Dimension>(gradWj, thpt) +
          //                  Wj*(outerProduct<Dimension>(rij, Tensor::one) + outerProduct<Dimension>(Tensor::one, rij)));
          gradm2(i) += Vj*(outerProduct<Dimension>(thpt, gradWj) +
                           Wj*(outerProduct<Dimension>(rij, Tensor::one) + outerProduct<Dimension>(Tensor::one, rij)));

          // // First moment. 
          // m1(i) += Vj*Wj * rij;
          // m1(j) -= Vi*Wi * rij;
          // gradm1(i) += Vj*( rij*gradWj + Tensor::one*Wj);
          // gradm1(j) += Vi*(-rij*gradWi + Tensor::one*Wi);

          // // Second moment.
          // const SymTensor thpt = rij.selfdyad();
          // m2(i) += Vj*Wj * thpt;
          // m2(j) += Vi*Wi * thpt;
          // gradm2(i) += Vj*outerProduct<Dimension>(thpt, gradWj);
          // gradm2(j) += Vi*outerProduct<Dimension>(thpt, gradWi);
          // for (size_t ii = 0; ii != Dimension::nDim; ++ii) {
          //   for (size_t jj = 0; jj != Dimension::nDim; ++jj) {
          //     gradm2(i)(ii, jj, jj) += Vj*Wj * rij(ii);
          //     gradm2(j)(ii, jj, jj) -= Vi*Wi * rij(ii);

          //     gradm2(i)(ii, jj, ii) += Vj*Wj * rij(jj);
          //     gradm2(j)(ii, jj, ii) -= Vi*Wi * rij(jj);
          //   }
          // }

        // }
      }

      // Based on the moments we can calculate the SVPH corrections terms and their gradients.
      if (i < firstGhostNodei) {
        CHECK2(abs(m2(i).Determinant()) > 1.0e-30, i << " " << m2(i).Determinant());
        const SymTensor m2inv = m2(i).Inverse();
        B(nodeListi, i) = -(m2inv*m1(i));

        // gradB(nodeListi, i) = -innerProduct<Dimension>(m2inv, gradm1(i)) + 
        //   innerProduct<Dimension>(innerProduct<Dimension>(innerProduct<Dimension>(m2inv, gradm2(i)), m2inv), m1(i));

        // gradB(nodeListi, i) = -innerProduct<Dimension>(m2inv, gradm1(i)) + 
        //   innerProduct<Dimension>(innerProduct<Dimension>(innerProduct<Dimension>(m2inv, gradm2(i)), m2inv), m1(i));

        gradB(nodeListi, i) = -innerProduct<Dimension>(gradm1(i), m2inv) + 
          innerProduct<Dimension>(innerProduct<Dimension>(innerProduct<Dimension>(m1(i), m2inv), gradm2(i)), m2inv);
      }
    }
  }
}

}
}

