//---------------------------------Spheral++------------------------------------
// Compute the CSPH corrections.
// This version combines the reproducing kernel interpolation of Liu et al.
// and the linear gradient correctios of Randles & Libersy.
//------------------------------------------------------------------------------
#include "computeCSPHCorrections2.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Geometry/outerProduct.hh"
#include "Geometry/innerProduct.hh"

namespace Spheral {
namespace CSPHSpace {

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
computeCSPHCorrections2(const ConnectivityMap<Dimension>& connectivityMap,
                        const TableKernel<Dimension>& W,
                        const FieldList<Dimension, typename Dimension::Scalar>& weight,
                        const FieldList<Dimension, typename Dimension::Vector>& position,
                        const FieldList<Dimension, typename Dimension::SymTensor>& H,
                        const bool coupleNodeLists,
                        FieldList<Dimension, typename Dimension::Scalar>& A0,
                        FieldList<Dimension, typename Dimension::Scalar>& A,
                        FieldList<Dimension, typename Dimension::Vector>& B,
                        FieldList<Dimension, typename Dimension::Tensor>& M) {

  // Pre-conditions.
  const size_t numNodeLists = A.size();
  REQUIRE(weight.size() == numNodeLists);
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);
  REQUIRE(A0.size() == numNodeLists);
  REQUIRE(B.size() == numNodeLists);
  REQUIRE(M.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  // Zero out the result.
  A0 = 0.0;
  A = 0.0;
  B = Vector::zero;
  M = Tensor::zero;

  // We can derive everything in terms of the zeroth, first, and second moments 
  // of the local positions.
  FieldList<Dimension, Scalar> m0(FieldSpace::Copy);
  FieldList<Dimension, Vector> m1(FieldSpace::Copy);
  FieldList<Dimension, SymTensor> m2(FieldSpace::Copy);
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const NodeList<Dimension>& nodeList = A[nodeListi]->nodeList();
    m0.appendNewField("zeroth moment", nodeList, 0.0);
    m1.appendNewField("first moment", nodeList, Vector::zero);
    m2.appendNewField("second moment", nodeList, SymTensor::zero);
  }

  // Walk the FluidNodeLists.
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const NodeList<Dimension>& nodeList = A[nodeListi]->nodeList();
    const int firstGhostNodei = nodeList.firstGhostNode();

    // Iterate over the nodes in this node list.
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      // Get the state for node i.
      const Scalar wi = weight(nodeListi, i);
      const Vector& ri = position(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();

      // Self contribution.
      const Scalar wwi = wi*W(0.0, Hdeti);
      m0(nodeListi, i) += wwi;

      // Neighbors!
      const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
      CHECK(fullConnectivity.size() == numNodeLists);

      for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
        if (coupleNodeLists or nodeListi == nodeListj) {
          const vector<int>& connectivity = fullConnectivity[nodeListj];
          const int firstGhostNodej = A[nodeListj]->nodeList().firstGhostNode();

          // Iterate over the neighbors for in this NodeList.
          for (vector<int>::const_iterator jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            const int j = *jItr;

            // Check if this node pair has already been calculated.
            if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                         nodeListj, j,
                                                         firstGhostNodej)) {

              // State of node j.
              const Scalar wj = weight(nodeListj, j);
              const Vector& rj = position(nodeListj, j);
              const SymTensor& Hj = H(nodeListj, j);
              const Scalar Hdetj = Hj.Determinant();

              // Kernel weighting and gradient.
              const Vector rij = ri - rj;
              const Vector etai = Hi*rij;
              const Vector etaj = Hj*rij;
              const std::pair<double, double> WWi = W.kernelAndGradValue(etai.magnitude(), Hdeti);
              const Scalar& Wi = WWi.first;
              const Vector gradWi = -(Hi*etai.unitVector())*WWi.second;
              const std::pair<double, double> WWj = W.kernelAndGradValue(etaj.magnitude(), Hdetj);
              const Scalar& Wj = WWj.first;
              const Vector gradWj = (Hj*etaj.unitVector())*WWj.second;

              // Zeroth moment. 
              const Scalar wwi = wi*Wi;
              const Scalar wwj = wj*Wj;
              m0(nodeListi, i) += wwj;
              m0(nodeListj, j) += wwi;

              // First moment. 
              m1(nodeListi, i) += wwj * rij;
              m1(nodeListj, j) -= wwi * rij;

              // Second moment.
              const SymTensor thpt = rij.selfdyad();
              m2(nodeListi, i) += wwj*thpt;
              m2(nodeListj, j) += wwi*thpt;

              // The derivative correction M.
              M(nodeListi, i) -= wj*rij.dyad(gradWj);
              M(nodeListj, j) += wi*rij.dyad(gradWi);
            }
          }
        }
      }

      // Based on the moments we can calculate the CSPH corrections terms.
      if (i < firstGhostNodei) {
        CHECK2(abs(m2(nodeListi, i).Determinant()) > 1.0e-30, i << " " << m0(nodeListi, i) << " " << m2(nodeListi, i).Determinant());
        const SymTensor m2inv = m2(nodeListi, i).Inverse();
        const Vector m2invm1 = m2inv*m1(nodeListi, i);
        const Scalar Ainv = m0(nodeListi, i) - m2invm1.dot(m1(nodeListi, i));
        CHECK(Ainv != 0.0);
        A0(nodeListi, i) = 1.0/m0(nodeListi, i);
        A(nodeListi, i) = 1.0/Ainv;
        B(nodeListi, i) = -m2invm1;
        M(nodeListi, i) = M(nodeListi, i).Inverse();

        // // BLAGO!
        // // Force only zeroth corrections.
        // A(nodeListi, i) = A0(nodeListi, i);
        // B(nodeListi, i) = Vector::zero;
        // gradA(nodeListi, i) = -FastMath::square(A(nodeListi, i))*gradm0(nodeListi, i);
        // gradB(nodeListi, i) = Tensor::zero;
        // // BLAGO!

      }
    }
  }

}

}
}

