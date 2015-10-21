//------------------------------------------------------------------------------
// Compute the CRKSPH corrections.
//------------------------------------------------------------------------------
#include <stdio.h>
#include "Eigen/Dense"

#include "computeCRKSPHCorrections.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Geometry/outerProduct.hh"
#include "Geometry/innerProduct.hh"
#include "Geometry/innerDoubleProduct.hh"
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
using Geometry::innerDoubleProduct;

namespace {

//------------------------------------------------------------------------------
// This is a specialized solver for C used in the quadratic RK "Mike" method.
// 1D
void solveC(const Dim<1>::FourthRankTensor& L,
            const Dim<1>::Tensor& Q,
            Dim<1>::Tensor& C,
            const double& tiny = 1.0e-15) {
  C(0,0) = std::abs(L(0,0,0,0)) > 1.0e-15 ? Q(0,0)/L(0,0,0,0) : 0.0;
}

void solveGradC(const Dim<1>::FourthRankTensor& L,
            const Dim<1>::ThirdRankTensor& gQ,
            Dim<1>::ThirdRankTensor& gC,
            const double& tiny = 1.0e-15) {
  gC(0,0,0) = std::abs(L(0,0,0,0)) > 1.0e-15 ? gQ(0,0,0)/L(0,0,0,0) : 0.0;
}


// 2D
void solveC(const Dim<2>::FourthRankTensor& L,
            const Dim<2>::Tensor& Q,
            Dim<2>::Tensor& C) {
  typedef Eigen::Matrix<double, 3, 3, Eigen::RowMajor> MatrixType;
  typedef Eigen::Matrix<double, 3, 1> VectorType;
  MatrixType A;
  VectorType b;
  A << L(0,0,0,0), L(1,0,0,0) + L(0,1,0,0), L(1,1,0,0),
       L(0,0,1,0), L(1,0,1,0) + L(0,1,1,0), L(1,1,1,0),
       L(0,0,1,1), L(1,0,1,1) + L(0,1,1,1), L(1,1,1,1);
  b << Q(0,0), Q(1,0),  Q(1,1);
  VectorType x = A.colPivHouseholderQr().solve(b);
  C(0,0) = x(0);
  C(0,1) = x(1);
  C(1,0) = x(1);
  C(1,1) = x(2);
}
void solveGradC(const Dim<2>::FourthRankTensor& L,
            const Dim<2>::ThirdRankTensor& gQ,
            Dim<2>::ThirdRankTensor& gC) {
  typedef Eigen::Matrix<double, 3, 3, Eigen::RowMajor> MatrixType;
  typedef Eigen::Matrix<double, 3, 1> VectorType;
  MatrixType A;
  VectorType b1;
  VectorType b2;
  A << L(0,0,0,0), L(1,0,0,0) + L(0,1,0,0), L(1,1,0,0),
       L(0,0,1,0), L(1,0,1,0) + L(0,1,1,0), L(1,1,1,0),
       L(0,0,1,1), L(1,0,1,1) + L(0,1,1,1), L(1,1,1,1);
  b1 << gQ(0,0,0), gQ(1,0,0),  gQ(1,1,0);
  b2 << gQ(0,0,1), gQ(1,0,1),  gQ(1,1,1);
  VectorType x1 = A.colPivHouseholderQr().solve(b1);
  VectorType x2 = A.colPivHouseholderQr().solve(b2);
  gC(0,0,0) = x1(0);
  gC(0,1,0) = x1(1);
  gC(1,0,0) = x1(1);
  gC(1,1,0) = x1(2);
  gC(0,0,1) = x2(0);
  gC(0,1,1) = x2(1);
  gC(1,0,1) = x2(1);
  gC(1,1,1) = x2(2);
 
}


// 3D
void solveC(const Dim<3>::FourthRankTensor& L,
            const Dim<3>::Tensor& Q,
            Dim<3>::Tensor& C) {
  typedef Eigen::Matrix<double, 6, 6, Eigen::RowMajor> MatrixType;
  typedef Eigen::Matrix<double, 6, 1> VectorType;
  MatrixType A;
  VectorType b;
  A << L(0,0,0,0), L(1,0,0,0) + L(0,1,0,0), L(2,0,0,0) + L(0,2,0,0), L(1,1,0,0), L(2,1,0,0) + L(1,2,0,0), L(2,2,0,0), 
       L(0,0,1,0), L(1,0,1,0) + L(0,1,1,0), L(2,0,1,0) + L(0,2,1,0), L(1,1,1,0), L(2,1,1,0) + L(1,2,1,0), L(2,2,1,0), 
       L(0,0,2,0), L(1,0,2,0) + L(0,1,2,0), L(2,0,2,0) + L(0,2,2,0), L(1,1,2,0), L(2,1,2,0) + L(1,2,2,0), L(2,2,2,0), 
       L(0,0,1,1), L(1,0,1,1) + L(0,1,1,1), L(2,0,1,1) + L(0,2,1,1), L(1,1,1,1), L(2,1,1,1) + L(1,2,1,1), L(2,2,1,1),
       L(0,0,2,1), L(1,0,2,1) + L(0,1,2,1), L(2,0,2,1) + L(0,2,2,1), L(1,1,2,1), L(2,1,2,1) + L(1,2,2,1), L(2,2,2,1),
       L(0,0,2,2), L(1,0,2,2) + L(0,1,2,2), L(2,0,2,2) + L(0,2,2,2), L(1,1,2,2), L(2,1,2,2) + L(1,2,2,2), L(2,2,2,2);
  b << Q(0,0), Q(1,0), Q(2,0), Q(1,1), Q(2,1), Q(2,2);
  VectorType x = A.colPivHouseholderQr().solve(b);
  C(0,0) = x(0);
  C(0,1) = x(1);
  C(0,2) = x(2);
  C(1,0) = x(1);
  C(1,1) = x(3);
  C(1,2) = x(4);
  C(2,0) = x(2);
  C(2,1) = x(4);
  C(2,2) = x(5);
}
void solveGradC(const Dim<3>::FourthRankTensor& L,
            const Dim<3>::ThirdRankTensor& gQ,
            Dim<3>::ThirdRankTensor& gC) {
  typedef Eigen::Matrix<double, 6, 6, Eigen::RowMajor> MatrixType;
  typedef Eigen::Matrix<double, 6, 1> VectorType;
  MatrixType A;
  VectorType b1;
  VectorType b2;
  A << L(0,0,0,0), L(1,0,0,0) + L(0,1,0,0), L(2,0,0,0) + L(0,2,0,0), L(1,1,0,0), L(2,1,0,0) + L(1,2,0,0), L(2,2,0,0),
       L(0,0,1,0), L(1,0,1,0) + L(0,1,1,0), L(2,0,1,0) + L(0,2,1,0), L(1,1,1,0), L(2,1,1,0) + L(1,2,1,0), L(2,2,1,0),
       L(0,0,2,0), L(1,0,2,0) + L(0,1,2,0), L(2,0,2,0) + L(0,2,2,0), L(1,1,2,0), L(2,1,2,0) + L(1,2,2,0), L(2,2,2,0),
       L(0,0,1,1), L(1,0,1,1) + L(0,1,1,1), L(2,0,1,1) + L(0,2,1,1), L(1,1,1,1), L(2,1,1,1) + L(1,2,1,1), L(2,2,1,1),
       L(0,0,2,1), L(1,0,2,1) + L(0,1,2,1), L(2,0,2,1) + L(0,2,2,1), L(1,1,2,1), L(2,1,2,1) + L(1,2,2,1), L(2,2,2,1),
       L(0,0,2,2), L(1,0,2,2) + L(0,1,2,2), L(2,0,2,2) + L(0,2,2,2), L(1,1,2,2), L(2,1,2,2) + L(1,2,2,2), L(2,2,2,2);
  b1 << gQ(0,0,0), gQ(1,0,0), gQ(2,0,0), gQ(1,1,0), gQ(2,1,0), gQ(2,2,0);
  b2 << gQ(0,0,1), gQ(1,0,1), gQ(2,0,1), gQ(1,1,1), gQ(2,1,1), gQ(2,2,1);
  VectorType x1 = A.colPivHouseholderQr().solve(b1);
  VectorType x2 = A.colPivHouseholderQr().solve(b2);
  gC(0,0,0) = x1(0);
  gC(0,1,0) = x1(1);
  gC(0,2,0) = x1(2);
  gC(1,0,0) = x1(1);
  gC(1,1,0) = x1(3);
  gC(1,2,0) = x1(4);
  gC(2,0,0) = x1(2);
  gC(2,1,0) = x1(4);
  gC(2,2,0) = x1(5);
  gC(0,0,1) = x2(0);
  gC(0,1,1) = x2(1);
  gC(0,2,1) = x2(2);
  gC(1,0,1) = x2(1);
  gC(1,1,1) = x2(3);
  gC(1,2,1) = x2(4);
  gC(2,0,1) = x2(2);
  gC(2,1,1) = x2(4);
  gC(2,2,1) = x2(5);
}



//------------------------------------------------------------------------------

}

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
// CRKSPH Quadratic Corrections Method (tensor based Mike version).
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
  typedef typename Dimension::FifthRankTensor FifthRankTensor;

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
  FieldList<Dimension, FourthRankTensor> gradm3(FieldSpace::Copy);
  FieldList<Dimension, FifthRankTensor> gradm4(FieldSpace::Copy);
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
    gradm3.appendNewField("grad third moment", nodeList, FourthRankTensor::zero);
    gradm4.appendNewField("grad fourth moment", nodeList, FifthRankTensor::zero);
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

            // Third Moment
            ThirdRankTensor thpt3 = outerProduct<Dimension>(thpt,rij);
            m3(nodeListi, i) += wwj*thpt3;
            m3(nodeListj, j) -= wwi*thpt3;
            gradm3(nodeListi, i) += wj*outerProduct<Dimension>(thpt3, gradWj);
            gradm3(nodeListj, j) -= wi*outerProduct<Dimension>(thpt3, gradWi);

            for (size_t ii = 0; ii != Dimension::nDim; ++ii) {
              for (size_t jj = 0; jj != Dimension::nDim; ++jj) {
               for (size_t kk = 0; kk != Dimension::nDim; ++kk) {
                gradm3(nodeListi, i)(ii, jj, kk, kk) += wwj*rij(ii)*rij(jj);
                gradm3(nodeListi, i)(ii, jj, kk, jj) += wwj*rij(ii)*rij(kk);
                gradm3(nodeListi, i)(ii, jj, kk, ii) += wwj*rij(jj)*rij(kk);

                gradm3(nodeListj, j)(ii, jj, kk, kk) += wwi*rij(ii)*rij(jj);
                gradm3(nodeListj, j)(ii, jj, kk, jj) += wwi*rij(ii)*rij(kk);
                gradm3(nodeListj, j)(ii, jj, kk, ii) += wwi*rij(jj)*rij(kk);

               }
              }
            }

            // Fourth Moment
            FourthRankTensor thpt4 = outerProduct<Dimension>(thpt3,rij);
            m4(nodeListi, i) += wwj*thpt4;
            m4(nodeListj, j) += wwi*thpt4;
            gradm4(nodeListi, i) += wj*outerProduct<Dimension>(thpt4, gradWj);
            gradm4(nodeListj, j) += wi*outerProduct<Dimension>(thpt4, gradWi);

            for (size_t ii = 0; ii != Dimension::nDim; ++ii) {
             for (size_t jj = 0; jj != Dimension::nDim; ++jj) {
              for (size_t kk = 0; kk != Dimension::nDim; ++kk) {
               for (size_t ll = 0; ll != Dimension::nDim; ++ll) {
                  gradm4(nodeListi, i)(ii, jj, kk, ll, ll) += wwj*rij(ii)*rij(jj)*rij(kk);
                  gradm4(nodeListi, i)(ii, jj, kk, ll, kk) += wwj*rij(ii)*rij(jj)*rij(ll);
                  gradm4(nodeListi, i)(ii, jj, kk, ll, jj) += wwj*rij(ii)*rij(kk)*rij(ll);
                  gradm4(nodeListi, i)(ii, jj, kk, ll, ii) += wwj*rij(jj)*rij(kk)*rij(ll);

                  gradm4(nodeListj, j)(ii, jj, kk, ll, ll) -= wwi*rij(ii)*rij(jj)*rij(kk);
                  gradm4(nodeListj, j)(ii, jj, kk, ll, kk) -= wwi*rij(ii)*rij(jj)*rij(ll);
                  gradm4(nodeListj, j)(ii, jj, kk, ll, jj) -= wwi*rij(ii)*rij(kk)*rij(ll);
                  gradm4(nodeListj, j)(ii, jj, kk, ll, ii) -= wwi*rij(jj)*rij(kk)*rij(ll);
               }
              }
             }
            }
          }
        }
      }

      // Based on the moments we can calculate the CRKSPH corrections terms and their gradients.
      if (i < firstGhostNodei) {
        // const Scalar m0i = m0(nodeListi, i);
        // const Scalar m1i = m1(nodeListi, i).x();
        // const Scalar m2i = m2(nodeListi, i).xx();
        // const Scalar m3i = m3(nodeListi, i)(0,0,0);
        // const Scalar m4i = m4(nodeListi, i)(0,0,0,0);
        // const Scalar m2inv = abs(m2i) > 1.0e-30 ? 1.0/m2i : 0.0;
        // const Scalar L = m3i*m2inv*m3i - m4i;
        // const Scalar Linv = abs(L) > 1.0e-30 ? 1.0/L : 0.0;
        // C(nodeListi, i).xx((m2i - m1i*m2inv*m3i)*Linv);
        // B(nodeListi, i).x(-(m1i + C(nodeListi, i).xx()*m3i)*m2inv);
        // const Scalar Ainv = m0i + B(nodeListi, i).x()*m1i + C(nodeListi, i).xx()*m2i;
        // CHECK(Ainv != 0.0);
        // A(nodeListi, i) = 1.0/Ainv;

        // cerr << " --> " << i << " " 
        //      << A(nodeListi,i)*(m0i + B(nodeListi, i).x()*m1i + C(nodeListi, i).xx()*m2i) << " "
        //      << A(nodeListi,i)*(m1i + B(nodeListi, i).x()*m2i + C(nodeListi, i).xx()*m3i) << " "
        //      << A(nodeListi,i)*(m2i + B(nodeListi, i).x()*m3i + C(nodeListi, i).xx()*m4i) << endl
        //      << "     " << m0i << " " << m1i << " " << m2i << " " << m3i << " " << m4i << endl
        //      << "     " << A(nodeListi,i) << " " << B(nodeListi,i) << " " << C(nodeListi, i) << endl;

        const Scalar  m0i = m0(nodeListi, i);
        const Vector& m1i = m1(nodeListi, i);
        const SymTensor& m2i = m2(nodeListi, i);
        const ThirdRankTensor& m3i = m3(nodeListi, i);
        const FourthRankTensor& m4i = m4(nodeListi, i);
        const Vector&  gm0i = gradm0(nodeListi, i);
        const Tensor& gm1i = gradm1(nodeListi, i);
        const ThirdRankTensor& gm2i = gradm2(nodeListi, i);
        const FourthRankTensor& gm3i = gradm3(nodeListi, i);
        const FifthRankTensor& gm4i = gradm4(nodeListi, i);
        const SymTensor m2inv = abs(m2i.Determinant()) > 1.0e-15 ? m2i.Inverse() : SymTensor::zero;
        const FourthRankTensor L = innerProduct<Dimension>(m3i, innerProduct<Dimension>(m2inv, m3i)) - m4i;

        // const FourthRankTensor Linv = invertRankNTensor(L);
        // C(nodeListi, i) = innerDoubleProduct<Dimension>(m2i - innerProduct<Dimension>(m1i, innerProduct<Dimension>(m2inv, m3i)), Linv);

        // const Tensor L2 = innerDoubleProduct<Dimension>(Tensor::one, L);
        // const Tensor Q = m2i - innerProduct<Dimension>(m1i, innerProduct<Dimension>(m2inv, m3i));
        // const Tensor Qinv = abs(Q.Determinant()) > 1.0e-15 ? Q.Inverse() : Tensor::zero;
        // const Tensor Cinv = innerProduct<Dimension>(L2, Qinv);
        // C(nodeListi, i) = abs(Cinv.Determinant()) > 1.0e-15 ? Cinv.Inverse() : Tensor::zero;

        const Tensor Q = m2i - innerProduct<Dimension>(m1i, innerProduct<Dimension>(m2inv, m3i));
        solveC(L, Q, C(nodeListi, i));
        const Vector coef = (m1i + innerDoubleProduct<Dimension>(C(nodeListi, i), m3i));

        B(nodeListi, i) = -innerProduct<Dimension>(coef, m2inv);
        const Scalar Ainv = m0i + innerProduct<Dimension>(B(nodeListi, i), m1i) + innerDoubleProduct<Dimension>(C(nodeListi, i), m2i);
        CHECK(Ainv != 0.0);
        A(nodeListi, i) = 1.0/Ainv;
        // Gradients (unfortunately for some gradient terms need to explicitly do for loops (or would need inner product that specified which indices).

        ThirdRankTensor gradQ = gradm2(nodeListi, i) - innerProduct<Dimension>(m1i, innerProduct<Dimension>(m2inv, gm3i));
        FifthRankTensor gradL = innerProduct<Dimension>(m3i, innerProduct<Dimension>(m2inv, gm3i)) - gm4i;
        for (size_t ii = 0; ii != Dimension::nDim; ++ii) {
         for (size_t jj = 0; jj != Dimension::nDim; ++jj) {
          for (size_t kk = 0; kk != Dimension::nDim; ++kk) {
           for (size_t ll = 0; ll != Dimension::nDim; ++ll) {
            for (size_t mm = 0; mm != Dimension::nDim; ++mm) {
             gradQ(ii,jj,kk) -= gm1i(ll,kk)*m2inv(ll,mm)*m3i(mm,ii,jj);
             for (size_t nn = 0; nn != Dimension::nDim; ++nn) {
              for (size_t oo = 0; oo != Dimension::nDim; ++oo) {
               gradQ(ii,jj,kk) += m1i(ll)*m2inv(ll,mm)*gm2i(mm,nn,kk)*m2inv(nn,oo)*m3i(oo,ii,jj);
               gradL(ii,jj,kk,ll,mm) += gm3i(ii,jj,nn,mm)*m2inv(nn,oo)*m3i(oo,kk,ll);
               for (size_t pp = 0; pp != Dimension::nDim; ++pp) {
                for (size_t qq = 0; qq != Dimension::nDim; ++qq) {
                 gradL(ii,jj,kk,ll,mm) -= m3i(ii,jj,nn)*m2inv(nn,oo)*gm2i(oo,pp,mm)*m2inv(pp,qq)*m3i(qq,kk,ll);
                }
               }
              }
             }
            }
           }
          }
         }
        }
        gradQ -= innerDoubleProduct<Dimension>(C(nodeListi, i),gradL);//Not really gradQ but gradQ - CgradL 
        solveGradC(L,gradQ,gradC(nodeListi, i));
        gradB(nodeListi, i) = -innerProduct<Dimension>(m2inv,(gm1i+innerDoubleProduct<Dimension>(m3i,gradC(nodeListi, i))+innerDoubleProduct<Dimension>(C(nodeListi, i),gm3i)));
        for (size_t ii = 0; ii != Dimension::nDim; ++ii) {
         for (size_t jj = 0; jj != Dimension::nDim; ++jj) {
          for (size_t kk = 0; kk != Dimension::nDim; ++kk) {
           for (size_t ll = 0; ll != Dimension::nDim; ++ll) {
            for (size_t mm = 0; mm != Dimension::nDim; ++mm) {
              gradB(nodeListi, i)(ii,jj) += coef(kk)*m2inv(kk,ll)*gm2i(ll,mm,jj)*m2inv(mm,ii);
            }
           }
          }
         }
        }
        gradA(nodeListi, i) = -A(nodeListi, i)*A(nodeListi, i)*(gm0i + innerProduct<Dimension>(m1i,gradB(nodeListi, i)) + innerProduct<Dimension>(B(nodeListi, i),gm1i)+innerDoubleProduct<Dimension>(C(nodeListi, i),gm2i)+innerDoubleProduct<Dimension>(m2i,gradC(nodeListi, i)));

      

        // cerr << " --> " << i << " " 
        //      << A(nodeListi,i)*(m0i + B(nodeListi, i).dot(m1i) + innerDoubleProduct<Dimension>(C(nodeListi, i), m2i)) << " "
        //      << A(nodeListi,i)*(m1i + innerProduct<Dimension>(B(nodeListi, i), m2i) + innerDoubleProduct<Dimension>(C(nodeListi, i), m3i)) << " "
        //      << A(nodeListi,i)*(m2i + innerProduct<Dimension>(B(nodeListi, i), m3i) + innerDoubleProduct<Dimension>(C(nodeListi, i), m4i)) << endl
        //      << "     " << m0i << " " << m1i << " " << m2i << " " << m3i << " " << m4i << endl
        //      << "     " << A(nodeListi,i) << " " << B(nodeListi,i) << " " << C(nodeListi, i) << endl;

      }

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
                         const CRKOrder correctionOrder,
                         FieldList<Dimension, typename Dimension::Scalar>& A,
                         FieldList<Dimension, typename Dimension::Vector>& B,
                         FieldList<Dimension, typename Dimension::Tensor>& C,
                         FieldList<Dimension, typename Dimension::Vector>& gradA,
                         FieldList<Dimension, typename Dimension::Tensor>& gradB,
                         FieldList<Dimension, typename Dimension::ThirdRankTensor>& gradC) {
                         if(correctionOrder == ZerothOrder){
                           computeZerothCRKSPHCorrections(connectivityMap,W,weight,position,H,A,gradA);
                         }else if(correctionOrder == LinearOrder){
                           computeLinearCRKSPHCorrections(connectivityMap,W,weight,position,H,A,B,gradA,gradB);
                         }else if(correctionOrder == QuadraticOrder){
                           computeQuadraticCRKSPHCorrections(connectivityMap,W,weight,position,H,A,B,C,gradA,gradB,gradC);
                         }
      
}


}//End Namespace
}

