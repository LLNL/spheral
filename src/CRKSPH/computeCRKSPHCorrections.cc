//------------------------------------------------------------------------------
// Compute the CRKSPH corrections.
//------------------------------------------------------------------------------
#include <stdio.h>
#include "computeCRKSPHCorrections.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Geometry/outerProduct.hh"
#include "Geometry/innerProduct.hh"

namespace Spheral {
namespace CRKSPHSpace {

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
computeCRKSPHCorrections(const ConnectivityMap<Dimension>& connectivityMap,
                         const TableKernel<Dimension>& W,
                         const FieldList<Dimension, typename Dimension::Scalar>& weight,
                         const FieldList<Dimension, typename Dimension::Vector>& position,
                         const FieldList<Dimension, typename Dimension::SymTensor>& H,
                         const FieldList<Dimension, typename Dimension::Scalar>& damage,
                         const bool coupleNodeLists,
                         FieldList<Dimension, typename Dimension::Scalar>& m0,
                         FieldList<Dimension, typename Dimension::Vector>& m1,
                         FieldList<Dimension, typename Dimension::SymTensor>& m2,
                         FieldList<Dimension, typename Dimension::Scalar>& A0,
                         FieldList<Dimension, typename Dimension::Scalar>& A,
                         FieldList<Dimension, typename Dimension::Vector>& B,
                         FieldList<Dimension, typename Dimension::Vector>& gradA0,
                         FieldList<Dimension, typename Dimension::Vector>& gradA,
                         FieldList<Dimension, typename Dimension::Tensor>& gradB) {

  // Pre-conditions.
  const size_t numNodeLists = A.size();
  REQUIRE(weight.size() == numNodeLists);
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);
  REQUIRE(damage.size() == numNodeLists or damage.size() == 0);
  REQUIRE(m0.size() == numNodeLists);
  REQUIRE(m1.size() == numNodeLists);
  REQUIRE(m2.size() == numNodeLists);
  REQUIRE(A0.size() == numNodeLists);
  REQUIRE(B.size() == numNodeLists);
  REQUIRE(gradA.size() == numNodeLists);
  REQUIRE(gradB.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  // Zero out the result.
  m0 = 0.0;
  m1 = Vector::zero;
  m2 = SymTensor::zero;
  A0 = 0.0;
  A = 0.0;
  B = Vector::zero;
  gradA = Vector::zero;
  gradB = Tensor::zero;

  // We can derive everything in terms of the zeroth, first, and second moments 
  // of the local positions.
  // Field<Dimension, Scalar> m0("zeroth moment", nodeList);
  // Field<Dimension, Vector> m1("first moment", nodeList);
  // Field<Dimension, SymTensor> m2("second moment", nodeList);
  FieldList<Dimension, Vector> gradm0(FieldSpace::Copy);
  FieldList<Dimension, Tensor> gradm1(FieldSpace::Copy);
  FieldList<Dimension, ThirdRankTensor> gradm2(FieldSpace::Copy);
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const NodeList<Dimension>& nodeList = A[nodeListi]->nodeList();
    gradm0.appendNewField("grad zeroth moment", nodeList, Vector::zero);
    gradm1.appendNewField("grad first moment", nodeList, Tensor::zero);
    gradm2.appendNewField("grad second moment", nodeList, ThirdRankTensor::zero);
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
      gradm1(nodeListi, i) += Tensor::one*wwi;

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
              gradm0(nodeListi, i) += wj*gradWj;
              gradm0(nodeListj, j) += wi*gradWi;

              // First moment. 
              m1(nodeListi, i) += wwj * rij;
              m1(nodeListj, j) -= wwi * rij;
              gradm1(nodeListi, i) += wj*(outerProduct<Dimension>( rij, gradWj) + Tensor::one*Wj);
              gradm1(nodeListj, j) += wi*(outerProduct<Dimension>(-rij, gradWi) + Tensor::one*Wi);

              // Second moment.
              const SymTensor thpt = rij.selfdyad();
              m2(nodeListi, i) += wwj*thpt;
              m2(nodeListj, j) += wwi*thpt;
              gradm2(nodeListi, i) += wj*outerProduct<Dimension>(thpt, gradWj);
              gradm2(nodeListj, j) += wi*outerProduct<Dimension>(thpt, gradWi);

              
              for (size_t ii = 0; ii != Dimension::nDim; ++ii) {
                for (size_t jj = 0; jj != Dimension::nDim; ++jj) {
                  gradm2(nodeListi, i)(ii, jj, jj) += wwj*rij(ii);
                  gradm2(nodeListi, i)(jj, ii, jj) += wwj*rij(ii);

                  gradm2(nodeListj, j)(ii, jj, jj) -= wwi*rij(ii);
                  gradm2(nodeListj, j)(jj, ii, jj) -= wwi*rij(ii);
                }
              }
            }
          }
        }
      }

      // Based on the moments we can calculate the CRKSPH corrections terms and their gradients.
      if (i < firstGhostNodei) {
        // CHECK2(abs(m2(nodeListi, i).Determinant()) > 1.0e-30, i << " " << m0(nodeListi, i) << " " << m2(nodeListi, i) << " " << m2(nodeListi, i).Determinant());
        const SymTensor m2inv = abs(m2(nodeListi, i).Determinant()) > 0.0 ? m2(nodeListi, i).Inverse() : SymTensor::zero;
        const Vector m2invm1 = m2inv*m1(nodeListi, i);
        const Scalar Ainv = m0(nodeListi, i) - m2invm1.dot(m1(nodeListi, i));
        CHECK(Ainv != 0.0);
        A0(nodeListi, i) = 1.0/m0(nodeListi, i);
        A(nodeListi, i) = 1.0/Ainv;
        B(nodeListi, i) = -m2invm1;
        gradA0(nodeListi, i) = -FastMath::square(A0(nodeListi, i))*gradm0(nodeListi, i);
        gradA(nodeListi, i) = -A(nodeListi, i)*A(nodeListi, i)*gradm0(nodeListi, i);
        for (size_t ii = 0; ii != Dimension::nDim; ++ii) {
          for (size_t jj = 0; jj != Dimension::nDim; ++jj) {
            for (size_t kk = 0; kk != Dimension::nDim; ++kk) {
              gradA(nodeListi, i)(ii) += A(nodeListi, i)*A(nodeListi, i)*m2inv(jj,kk)*m1(nodeListi, i)(kk)*gradm1(nodeListi, i)(jj,ii);
              gradA(nodeListi, i)(ii) += A(nodeListi, i)*A(nodeListi, i)*m2inv(jj,kk)*gradm1(nodeListi, i)(kk,ii)*m1(nodeListi, i)(jj);
              gradB(nodeListi, i)(ii,jj) -= m2inv(ii,kk)*gradm1(nodeListi, i)(kk,jj);
              for (size_t ll = 0; ll != Dimension::nDim; ++ll) {
                for (size_t mm = 0; mm != Dimension::nDim; ++mm) {
                  gradA(nodeListi, i)(ii) -= A(nodeListi, i)*A(nodeListi, i)*m2inv(jj,kk)*gradm2(nodeListi, i)(kk,ll,ii)*m2inv(ll,mm)*m1(nodeListi, i)(mm)*m1(nodeListi, i)(jj);
                  gradB(nodeListi, i)(ii,jj) += m2inv(ii,kk)*gradm2(nodeListi, i)(kk,ll,jj)*m2inv(ll,mm)*m1(nodeListi, i)(mm);
                }
              }
            }
          }
        }

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

  // Note -- we are suspending the next order corrections for now!

}

}
}

