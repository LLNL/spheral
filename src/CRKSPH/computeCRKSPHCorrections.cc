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
#include "Geometry/invertRankNTensor.hh"

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

//------------------------------------------------------------------------------
// CRKSPH Zeroth Order Corrections Method 
//------------------------------------------------------------------------------
template<typename Dimension>
void
computeZerothCRKSPHCorrections(const ConnectivityMap<Dimension>& connectivityMap,
                         const TableKernel<Dimension>& W,
                         const FieldList<Dimension, typename Dimension::Scalar>& weight,
                         const FieldList<Dimension, typename Dimension::Vector>& position,
                         const FieldList<Dimension, typename Dimension::SymTensor>& H,
                         FieldList<Dimension, typename Dimension::Scalar>& A,
                         FieldList<Dimension, typename Dimension::Vector>& gradA) {

  // Pre-conditions.
  const size_t numNodeLists = A.size();
  REQUIRE(weight.size() == numNodeLists);
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);
  REQUIRE(gradA.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  // We can derive everything in terms of the zeroth, first, and second moments 
  // of the local positions.
  FieldList<Dimension, Scalar> m0(FieldSpace::Copy);
  FieldList<Dimension, Vector> gradm0(FieldSpace::Copy);
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const NodeList<Dimension>& nodeList = A[nodeListi]->nodeList();
    m0.appendNewField("zeroth moment", nodeList, 0.0);
    gradm0.appendNewField("grad zeroth moment", nodeList, Vector::zero);
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
      const Vector& ri = position(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();

      // Self contribution.
      const Scalar wwi =  weight(nodeListi, i)*W(0.0, Hdeti);
      m0(nodeListi, i) += wwi;

      // Neighbors!
      const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
      CHECK(fullConnectivity.size() == numNodeLists);

      for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
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
            const Vector& rj = position(nodeListj, j);
            const SymTensor& Hj = H(nodeListj, j);
            const Scalar Hdetj = Hj.Determinant();

            // Node weighting with pair-wise coupling.
            const Scalar wi = weight(nodeListi, i);
            const Scalar wj = weight(nodeListj, j);

            // Kernel weighting and gradient.
            const Vector rij = ri - rj;
            const Vector etai = Hi*rij;
            const Vector etaj = Hj*rij;
            const std::pair<double, double> WWi = W.kernelAndGradValue(etai.magnitude(), Hdeti);
            const Scalar Wi = WWi.first;
            const Vector gradWi = -(Hi*etai.unitVector())*WWi.second;
            const std::pair<double, double> WWj = W.kernelAndGradValue(etaj.magnitude(), Hdetj);
            const Scalar Wj = WWj.first;
            const Vector gradWj = (Hj*etaj.unitVector())*WWj.second;

            // Zeroth moment. 
            const Scalar wwi = wi*Wi;
            const Scalar wwj = wj*Wj;
            m0(nodeListi, i) += wwj;
            m0(nodeListj, j) += wwi;
            gradm0(nodeListi, i) += wj*gradWj;
            gradm0(nodeListj, j) += wi*gradWi;

          }
        }
      }

      // Based on the moments we can calculate the CRKSPH corrections terms and their gradients.
      if (i < firstGhostNodei) {
        CHECK(m0(nodeListi, i) != 0.0);
        A(nodeListi, i) = 1.0/m0(nodeListi, i);
        gradA(nodeListi, i) = -FastMath::square(A(nodeListi, i))*gradm0(nodeListi, i);
      }
    }
  }
}

//------------------------------------------------------------------------------
// CRKSPH Linear Corrections Method 
//------------------------------------------------------------------------------
template<typename Dimension>
void
computeLinearCRKSPHCorrections(const ConnectivityMap<Dimension>& connectivityMap,
                         const TableKernel<Dimension>& W,
                         const FieldList<Dimension, typename Dimension::Scalar>& weight,
                         const FieldList<Dimension, typename Dimension::Vector>& position,
                         const FieldList<Dimension, typename Dimension::SymTensor>& H,
                         FieldList<Dimension, typename Dimension::Scalar>& A,
                         FieldList<Dimension, typename Dimension::Vector>& B,
                         FieldList<Dimension, typename Dimension::Vector>& gradA,
                         FieldList<Dimension, typename Dimension::Tensor>& gradB) {

  // Pre-conditions.
  const size_t numNodeLists = A.size();
  REQUIRE(weight.size() == numNodeLists);
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);
  REQUIRE(B.size() == numNodeLists);
  REQUIRE(gradA.size() == numNodeLists);
  REQUIRE(gradB.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  // We can derive everything in terms of the zeroth, first, and second moments 
  // of the local positions.
  FieldList<Dimension, Scalar> m0(FieldSpace::Copy);
  FieldList<Dimension, Vector> m1(FieldSpace::Copy);
  FieldList<Dimension, SymTensor> m2(FieldSpace::Copy);
  FieldList<Dimension, Vector> gradm0(FieldSpace::Copy);
  FieldList<Dimension, Tensor> gradm1(FieldSpace::Copy);
  FieldList<Dimension, ThirdRankTensor> gradm2(FieldSpace::Copy);
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const NodeList<Dimension>& nodeList = A[nodeListi]->nodeList();
    m0.appendNewField("zeroth moment", nodeList, 0.0);
    m1.appendNewField("first moment", nodeList, Vector::zero);
    m2.appendNewField("second moment", nodeList, SymTensor::zero);
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
      const Vector& ri = position(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();

      // Self contribution.
      const Scalar wwi =  weight(nodeListi, i)*W(0.0, Hdeti);
      m0(nodeListi, i) += wwi;
      gradm1(nodeListi, i) += Tensor::one*wwi;

      // Neighbors!
      const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
      CHECK(fullConnectivity.size() == numNodeLists);

      for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
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
            const Vector& rj = position(nodeListj, j);
            const SymTensor& Hj = H(nodeListj, j);
            const Scalar Hdetj = Hj.Determinant();

            // Node weighting with pair-wise coupling.
            const Scalar wi = weight(nodeListi, i);
            const Scalar wj = weight(nodeListj, j);

            // Kernel weighting and gradient.
            const Vector rij = ri - rj;
            const Vector etai = Hi*rij;
            const Vector etaj = Hj*rij;
            const std::pair<double, double> WWi = W.kernelAndGradValue(etai.magnitude(), Hdeti);
            const Scalar Wi = WWi.first;
            const Vector gradWi = -(Hi*etai.unitVector())*WWi.second;
            const std::pair<double, double> WWj = W.kernelAndGradValue(etaj.magnitude(), Hdetj);
            const Scalar Wj = WWj.first;
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

      // Based on the moments we can calculate the CRKSPH corrections terms and their gradients.
      if (i < firstGhostNodei) {
        const SymTensor m2inv = abs(m2(nodeListi, i).Determinant()) > 1.0e-10 ? m2(nodeListi, i).Inverse() : SymTensor::zero;
        const Vector m2invm1 = m2inv*m1(nodeListi, i);
        const Scalar Ainv = m0(nodeListi, i) - m2invm1.dot(m1(nodeListi, i));
        CHECK(Ainv != 0.0);
        A(nodeListi, i) = 1.0/Ainv;
        B(nodeListi, i) = -m2invm1;
        gradA(nodeListi, i) = -A(nodeListi, i)*A(nodeListi, i)*gradm0(nodeListi, i);
        gradB(nodeListi, i) = Tensor::zero;
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
      }
    }
  }
}

//------------------------------------------------------------------------------
// CRKSPH Quadratic Corrections Method 
//------------------------------------------------------------------------------
template<typename Dimension>
void
computeQuadraticCRKSPHCorrections(const ConnectivityMap<Dimension>& connectivityMap,
                                  const TableKernel<Dimension>& W,
                                  const FieldList<Dimension, typename Dimension::Scalar>& weight,
                                  const FieldList<Dimension, typename Dimension::Vector>& position,
                                  const FieldList<Dimension, typename Dimension::SymTensor>& H,
                                  FieldList<Dimension, typename Dimension::Scalar>& A,
                                  FieldList<Dimension, typename Dimension::Vector>& B,
                                  FieldList<Dimension, typename Dimension::Tensor>& C,
                                  FieldList<Dimension, typename Dimension::Vector>& gradA,
                                  FieldList<Dimension, typename Dimension::Tensor>& gradB,
                                  FieldList<Dimension, typename Dimension::ThirdRankTensor>& gradC) {

  // Pre-conditions.
  const size_t numNodeLists = A.size();
  REQUIRE(weight.size() == numNodeLists);
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);
  REQUIRE(B.size() == numNodeLists);
  REQUIRE(C.size() == numNodeLists);
  REQUIRE(gradA.size() == numNodeLists);
  REQUIRE(gradB.size() == numNodeLists);
  REQUIRE(gradC.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef typename Dimension::FourthRankTensor FourthRankTensor;

  const double tiny = 1.0e-30;

  // We can derive everything in terms of the zeroth, first, and second moments 
  // of the local positions.
  FieldList<Dimension, Scalar> m0(FieldSpace::Copy);
  FieldList<Dimension, Vector> m1(FieldSpace::Copy);
  FieldList<Dimension, SymTensor> m2(FieldSpace::Copy);
  FieldList<Dimension, ThirdRankTensor> m3(FieldSpace::Copy);
  FieldList<Dimension, FourthRankTensor> m4(FieldSpace::Copy);
  FieldList<Dimension, Vector> gradm0(FieldSpace::Copy);
  FieldList<Dimension, Tensor> gradm1(FieldSpace::Copy);
  FieldList<Dimension, ThirdRankTensor> gradm2(FieldSpace::Copy);
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const NodeList<Dimension>& nodeList = A[nodeListi]->nodeList();
    m0.appendNewField("zeroth moment", nodeList, 0.0);
    m1.appendNewField("first moment", nodeList, Vector::zero);
    m2.appendNewField("second moment", nodeList, SymTensor::zero);
    m3.appendNewField("third moment", nodeList, ThirdRankTensor::zero);
    m4.appendNewField("fourth moment", nodeList, FourthRankTensor::zero);
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
      const Vector& ri = position(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();

      // Self contribution.
      const Scalar wwi =  weight(nodeListi, i)*W(0.0, Hdeti);
      m0(nodeListi, i) += wwi;
      gradm1(nodeListi, i) += Tensor::one*wwi;

      // Neighbors!
      const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
      CHECK(fullConnectivity.size() == numNodeLists);

      for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
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
            const Vector& rj = position(nodeListj, j);
            const SymTensor& Hj = H(nodeListj, j);
            const Scalar Hdetj = Hj.Determinant();

            // Node weighting with pair-wise coupling.
            const Scalar wi = weight(nodeListi, i);
            const Scalar wj = weight(nodeListj, j);

            // Kernel weighting and gradient.
            const Vector rij = ri - rj;
            const Vector etai = Hi*rij;
            const Vector etaj = Hj*rij;
            const std::pair<double, double> WWi = W.kernelAndGradValue(etai.magnitude(), Hdeti);
            const Scalar Wi = WWi.first;
            const Vector gradWi = -(Hi*etai.unitVector())*WWi.second;
            const std::pair<double, double> WWj = W.kernelAndGradValue(etaj.magnitude(), Hdetj);
            const Scalar Wj = WWj.first;
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

            // Third and fourth moment.
            ThirdRankTensor thpt3;
            FourthRankTensor thpt4;
            for (size_t ii = 0; ii != Dimension::nDim; ++ii) {
              for (size_t jj = 0; jj != Dimension::nDim; ++jj) {
                for (size_t kk = 0; kk != Dimension::nDim; ++kk) {
                  thpt3(ii,jj,kk) = rij(ii)*rij(jj)*rij(kk);
                  for (size_t mm = 0; mm != Dimension::nDim; ++mm) {
                    thpt4(ii,jj,kk,mm) = rij(ii)*rij(jj)*rij(kk)*rij(mm);
                  }
                }
              }
            }
            m3(nodeListi, i) += wwj*thpt3;
            m3(nodeListj, j) += wwi*thpt3;
            m4(nodeListi, i) += wwj*thpt4;
            m4(nodeListj, j) += wwi*thpt4;

          }
        }
      }

      // Based on the moments we can calculate the CRKSPH corrections terms and their gradients.
      // if (i < firstGhostNodei) {
      //   const SymTensor m2inv = abs(m2(nodeListi, i).Determinant()) > 1.0e-10 ? m2(nodeListi, i).Inverse() : SymTensor::zero;
      //   const FourthRankTensor m4inv = invertRankNTensor(m4(nodeListi, i));
      //   const Scalar Bden = 1.0 - innerProduct<Dimension>(m3(nodeListi, i), innerProduct<Dimension>(m4inv, innerProduct<Dimension>(m3(nodeListi, i), m2inv)));
      //   B(nodeListi, i) = (-innerProduct<Dimension>(m1(nodeListi, i), m2inv) + 
      //                      innerProduct<Dimension>(m2(nodeListi, i), innerProduct<Dimension>(m4inv, innerProduct<Dimension>(m3(nodeListi, i), m2inv))))/max(tiny, abs(Bden))*sgn(Bden);
      //   C(nodeListi, i) = -innerProduct<Dimension>(m2(nodeListi, i) + innerProduct<Dimension>(B(nodeListi, i), m3(nodeListi, i)), m4inv);
      //   const Scalar Ainv = m0(nodeListi, i) - innerProduct<Dimension>(B(nodeListi, i), m1(nodeListi, i)) + innerProduct<Dimension>(C(nodeListi, i), m2(nodeListi, i));
      //   CHECK(Ainv != 0.0);
      //   A(nodeListi, i) = 1.0/Ainv;
      // }
    }
  }
}
//------------------------------------------------------------------------------
// A modified form where we can specify a pair-wise node coupling functor.
// This method returns both the corrections assuming full coupling (A,B,gradA,gradB)
// as well as corrections employing the specified coupling (Ac,Bc,gradAc,gradBc).
//------------------------------------------------------------------------------
template<typename Dimension>
void
computeCRKSPHCorrections(const ConnectivityMap<Dimension>& connectivityMap,
                         const TableKernel<Dimension>& W,
                         const FieldList<Dimension, typename Dimension::Scalar>& weight,
                         const FieldList<Dimension, typename Dimension::Vector>& position,
                         const FieldList<Dimension, typename Dimension::SymTensor>& H,
                         const NodeCoupling& nodeCoupling,
                         FieldList<Dimension, typename Dimension::Scalar>& A,
                         FieldList<Dimension, typename Dimension::Vector>& B,
                         FieldList<Dimension, typename Dimension::Vector>& gradA,
                         FieldList<Dimension, typename Dimension::Tensor>& gradB,
                         FieldList<Dimension, typename Dimension::Scalar>& Ac,
                         FieldList<Dimension, typename Dimension::Vector>& Bc,
                         FieldList<Dimension, typename Dimension::Vector>& gradAc,
                         FieldList<Dimension, typename Dimension::Tensor>& gradBc) {

  // Pre-conditions.
  const size_t numNodeLists = A.size();
  REQUIRE(weight.size() == numNodeLists);
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);
  REQUIRE(B.size() == numNodeLists);
  REQUIRE(gradA.size() == numNodeLists);
  REQUIRE(gradB.size() == numNodeLists);
  REQUIRE(Ac.size() == numNodeLists);
  REQUIRE(Bc.size() == numNodeLists);
  REQUIRE(gradAc.size() == numNodeLists);
  REQUIRE(gradBc.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  // We can derive everything in terms of the zeroth, first, and second moments 
  // of the local positions.
  FieldList<Dimension, Scalar> m0(FieldSpace::Copy), m0c(FieldSpace::Copy);
  FieldList<Dimension, Vector> m1(FieldSpace::Copy), m1c(FieldSpace::Copy);
  FieldList<Dimension, SymTensor> m2(FieldSpace::Copy), m2c(FieldSpace::Copy);
  FieldList<Dimension, Vector> gradm0(FieldSpace::Copy), gradm0c(FieldSpace::Copy);
  FieldList<Dimension, Tensor> gradm1(FieldSpace::Copy), gradm1c(FieldSpace::Copy);
  FieldList<Dimension, ThirdRankTensor> gradm2(FieldSpace::Copy), gradm2c(FieldSpace::Copy);
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const NodeList<Dimension>& nodeList = A[nodeListi]->nodeList();
    m0.appendNewField("zeroth moment", nodeList, 0.0);
    m1.appendNewField("first moment", nodeList, Vector::zero);
    m2.appendNewField("second moment", nodeList, SymTensor::zero);
    gradm0.appendNewField("grad zeroth moment", nodeList, Vector::zero);
    gradm1.appendNewField("grad first moment", nodeList, Tensor::zero);
    gradm2.appendNewField("grad second moment", nodeList, ThirdRankTensor::zero);
    m0c.appendNewField("zeroth moment coupling", nodeList, 0.0);
    m1c.appendNewField("first moment coupling", nodeList, Vector::zero);
    m2c.appendNewField("second moment coupling", nodeList, SymTensor::zero);
    gradm0c.appendNewField("grad zeroth moment coupling", nodeList, Vector::zero);
    gradm1c.appendNewField("grad first moment coupling", nodeList, Tensor::zero);
    gradm2c.appendNewField("grad second moment coupling", nodeList, ThirdRankTensor::zero);
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
      const Vector& ri = position(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();

      // Self contribution.
      const Scalar wwi = weight(nodeListi, i)*W(0.0, Hdeti);
      m0(nodeListi, i) += wwi;
      gradm1(nodeListi, i) += Tensor::one*wwi;
      m0c(nodeListi, i) += wwi;
      gradm1c(nodeListi, i) += Tensor::one*wwi;

      // Neighbors!
      const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
      CHECK(fullConnectivity.size() == numNodeLists);

      for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
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
            const Vector& rj = position(nodeListj, j);
            const SymTensor& Hj = H(nodeListj, j);
            const Scalar Hdetj = Hj.Determinant();

            // Find the pair weighting scaling.
            const double fij = nodeCoupling(nodeListi, i, nodeListj, j);
            CHECK(fij >= 0.0 and fij <= 1.0);

            // Node weighting with pair-wise coupling.
            const Scalar wi = weight(nodeListi, i);
            const Scalar wj = weight(nodeListj, j);

            // Kernel weighting and gradient.
            const Vector rij = ri - rj;
            const Vector etai = Hi*rij;
            const Vector etaj = Hj*rij;
            const std::pair<double, double> WWi = W.kernelAndGradValue(etai.magnitude(), Hdeti);
            const Scalar Wi = WWi.first;
            const Vector gradWi = -(Hi*etai.unitVector())*WWi.second;
            const std::pair<double, double> WWj = W.kernelAndGradValue(etaj.magnitude(), Hdetj);
            const Scalar Wj = WWj.first;
            const Vector gradWj = (Hj*etaj.unitVector())*WWj.second;

            // Zeroth moment. 
            const Scalar wwi = wi*Wi;
            const Scalar wwj = wj*Wj;
            m0(nodeListi, i) += wwj;
            m0(nodeListj, j) += wwi;
            gradm0(nodeListi, i) += wj*gradWj;
            gradm0(nodeListj, j) += wi*gradWi;

            m0c(nodeListi, i) += fij*wwj;
            m0c(nodeListj, j) += fij*wwi;
            gradm0c(nodeListi, i) += fij*wj*gradWj;
            gradm0c(nodeListj, j) += fij*wi*gradWi;

            // First moment. 
            m1(nodeListi, i) += wwj * rij;
            m1(nodeListj, j) -= wwi * rij;
            gradm1(nodeListi, i) += wj*(outerProduct<Dimension>( rij, gradWj) + Tensor::one*Wj);
            gradm1(nodeListj, j) += wi*(outerProduct<Dimension>(-rij, gradWi) + Tensor::one*Wi);

            m1c(nodeListi, i) += fij*wwj * rij;
            m1c(nodeListj, j) -= fij*wwi * rij;
            gradm1c(nodeListi, i) += fij*wj*(outerProduct<Dimension>( rij, gradWj) + Tensor::one*Wj);
            gradm1c(nodeListj, j) += fij*wi*(outerProduct<Dimension>(-rij, gradWi) + Tensor::one*Wi);

            // Second moment.
            const SymTensor thpt = rij.selfdyad();
            m2(nodeListi, i) += wwj*thpt;
            m2(nodeListj, j) += wwi*thpt;
            gradm2(nodeListi, i) += wj*outerProduct<Dimension>(thpt, gradWj);
            gradm2(nodeListj, j) += wi*outerProduct<Dimension>(thpt, gradWi);

            m2c(nodeListi, i) += fij*wwj*thpt;
            m2c(nodeListj, j) += fij*wwi*thpt;
            gradm2c(nodeListi, i) += fij*wj*outerProduct<Dimension>(thpt, gradWj);
            gradm2c(nodeListj, j) += fij*wi*outerProduct<Dimension>(thpt, gradWi);

            for (size_t ii = 0; ii != Dimension::nDim; ++ii) {
              for (size_t jj = 0; jj != Dimension::nDim; ++jj) {
                gradm2(nodeListi, i)(ii, jj, jj) += wwj*rij(ii);
                gradm2(nodeListi, i)(jj, ii, jj) += wwj*rij(ii);

                gradm2(nodeListj, j)(ii, jj, jj) -= wwi*rij(ii);
                gradm2(nodeListj, j)(jj, ii, jj) -= wwi*rij(ii);

                gradm2c(nodeListi, i)(ii, jj, jj) += fij*wwj*rij(ii);
                gradm2c(nodeListi, i)(jj, ii, jj) += fij*wwj*rij(ii);

                gradm2c(nodeListj, j)(ii, jj, jj) -= fij*wwi*rij(ii);
                gradm2c(nodeListj, j)(jj, ii, jj) -= fij*wwi*rij(ii);
              }
            }
          }
        }
      }

      // Based on the moments we can calculate the CRKSPH corrections terms and their gradients.
      if (i < firstGhostNodei) {
        const Scalar m2det = abs(m2(nodeListi, i).Determinant());
        const SymTensor m2inv = m2det > 1.0e-10 ? m2(nodeListi, i).Inverse() * m2det/max(1e-8, m2det) : SymTensor::zero;
        const Vector m2invm1 = m2inv*m1(nodeListi, i);
        const Scalar Ainv = m0(nodeListi, i) - m2invm1.dot(m1(nodeListi, i));
        CHECK(Ainv != 0.0);
        A(nodeListi, i) = 1.0/Ainv;
        B(nodeListi, i) = -m2invm1;
        gradA(nodeListi, i) = -A(nodeListi, i)*A(nodeListi, i)*gradm0(nodeListi, i);
        gradB(nodeListi, i) = Tensor::zero;
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

        const Scalar m2cdet = abs(m2c(nodeListi, i).Determinant());
        const SymTensor m2cinv = m2cdet > 1.0e-10 ? m2c(nodeListi, i).Inverse() * m2cdet/max(1e-8, m2cdet) : SymTensor::zero;
        const Vector m2cinvm1c = m2cinv*m1c(nodeListi, i);
        const Scalar Acinv = m0c(nodeListi, i) - m2cinvm1c.dot(m1c(nodeListi, i));
        CHECK(Acinv != 0.0);
        Ac(nodeListi, i) = 1.0/Acinv;
        Bc(nodeListi, i) = -m2cinvm1c;
        gradAc(nodeListi, i) = -Ac(nodeListi, i)*Ac(nodeListi, i)*gradm0c(nodeListi, i);
        gradBc(nodeListi, i) = Tensor::zero;
        for (size_t ii = 0; ii != Dimension::nDim; ++ii) {
          for (size_t jj = 0; jj != Dimension::nDim; ++jj) {
            for (size_t kk = 0; kk != Dimension::nDim; ++kk) {
              gradAc(nodeListi, i)(ii) += Ac(nodeListi, i)*Ac(nodeListi, i)*m2cinv(jj,kk)*m1c(nodeListi, i)(kk)*gradm1c(nodeListi, i)(jj,ii);
              gradAc(nodeListi, i)(ii) += Ac(nodeListi, i)*Ac(nodeListi, i)*m2cinv(jj,kk)*gradm1c(nodeListi, i)(kk,ii)*m1c(nodeListi, i)(jj);
              gradBc(nodeListi, i)(ii,jj) -= m2cinv(ii,kk)*gradm1c(nodeListi, i)(kk,jj);
              for (size_t ll = 0; ll != Dimension::nDim; ++ll) {
                for (size_t mm = 0; mm != Dimension::nDim; ++mm) {
                  gradAc(nodeListi, i)(ii) -= Ac(nodeListi, i)*Ac(nodeListi, i)*m2cinv(jj,kk)*gradm2c(nodeListi, i)(kk,ll,ii)*m2cinv(ll,mm)*m1c(nodeListi, i)(mm)*m1c(nodeListi, i)(jj);
                  gradBc(nodeListi, i)(ii,jj) += m2cinv(ii,kk)*gradm2c(nodeListi, i)(kk,ll,jj)*m2cinv(ll,mm)*m1c(nodeListi, i)(mm);
                }
              }
            }
          }
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
// The basic method, assuming all nodes are fully coupled.
//------------------------------------------------------------------------------
template<typename Dimension>
void
computeCRKSPHCorrections(const ConnectivityMap<Dimension>& connectivityMap,
                         const TableKernel<Dimension>& W,
                         const FieldList<Dimension, typename Dimension::Scalar>& weight,
                         const FieldList<Dimension, typename Dimension::Vector>& position,
                         const FieldList<Dimension, typename Dimension::SymTensor>& H,
                         const int correctionOrder,
                         FieldList<Dimension, typename Dimension::Scalar>& A,
                         FieldList<Dimension, typename Dimension::Vector>& B,
                         FieldList<Dimension, typename Dimension::Tensor>& C,
                         FieldList<Dimension, typename Dimension::Vector>& gradA,
                         FieldList<Dimension, typename Dimension::Tensor>& gradB,
                         FieldList<Dimension, typename Dimension::ThirdRankTensor>& gradC) {
                         if(correctionOrder == 0){
                           computeZerothCRKSPHCorrections(connectivityMap,W,weight,position,H,A,gradA);
                         }else if(correctionOrder == 1){
                           computeLinearCRKSPHCorrections(connectivityMap,W,weight,position,H,A,B,gradA,gradB);
                         }else if(correctionOrder == 2){
                           computeQuadraticCRKSPHCorrections(connectivityMap,W,weight,position,H,A,B,C,gradA,gradB,gradC);
                         }
      
}


}//End Namespace
}

