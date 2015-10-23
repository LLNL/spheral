//------------------------------------------------------------------------------
// Compute the moments necessary for CRKSPH corrections.
//------------------------------------------------------------------------------
#include <stdio.h>

#include "computeCRKSPHMoments.hh"
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

//------------------------------------------------------------------------------
// Compute the moments.
//------------------------------------------------------------------------------
template<typename Dimension>
void
computeCRKSPHMoments(const ConnectivityMap<Dimension>& connectivityMap,
                     const TableKernel<Dimension>& W,
                     const FieldList<Dimension, typename Dimension::Scalar>& weight,
                     const FieldList<Dimension, typename Dimension::Vector>& position,
                     const FieldList<Dimension, typename Dimension::SymTensor>& H,
                     const CRKOrder correctionOrder,
                     FieldList<Dimension, typename Dimension::Scalar>& m0,
                     FieldList<Dimension, typename Dimension::Vector>& m1,
                     FieldList<Dimension, typename Dimension::SymTensor>& m2,
                     FieldList<Dimension, typename Dimension::ThirdRankTensor>& m3,
                     FieldList<Dimension, typename Dimension::FourthRankTensor>& m4,
                     FieldList<Dimension, typename Dimension::Vector>& gradm0,
                     FieldList<Dimension, typename Dimension::Tensor>& gradm1,
                     FieldList<Dimension, typename Dimension::ThirdRankTensor>& gradm2,
                     FieldList<Dimension, typename Dimension::FourthRankTensor>& gradm3,
                     FieldList<Dimension, typename Dimension::FifthRankTensor>& gradm4) {

  // Pre-conditions.
  const size_t numNodeLists = weight.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);
  REQUIRE(m0.size() == numNodeLists);
  REQUIRE(m1.size() == numNodeLists);
  REQUIRE(m2.size() == numNodeLists);
  REQUIRE(gradm0.size() == numNodeLists);
  REQUIRE(gradm1.size() == numNodeLists);
  REQUIRE(gradm2.size() == numNodeLists);
  if (correctionOrder == QuadraticOrder) {
    REQUIRE(m3.size() == numNodeLists);
    REQUIRE(m4.size() == numNodeLists);
    REQUIRE(gradm3.size() == numNodeLists);
    REQUIRE(gradm4.size() == numNodeLists);
  }

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef typename Dimension::FourthRankTensor FourthRankTensor;
  typedef typename Dimension::FifthRankTensor FifthRankTensor;

  // Zero things out.
  m0 = 0.0;
  m1 = Vector::zero;
  m2 = SymTensor::zero;
  gradm0 = Vector::zero;
  gradm1 = Tensor::zero;
  gradm2 = ThirdRankTensor::zero;
  if (correctionOrder == QuadraticOrder) {
    m3 = ThirdRankTensor::zero;
    m4 = FourthRankTensor::zero;
    gradm3 = FourthRankTensor::zero;
    gradm4 = FifthRankTensor::zero;
  }

  // Walk the FluidNodeLists.
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const NodeList<Dimension>& nodeList = weight[nodeListi]->nodeList();
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

      Scalar thpt;
      for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        const int firstGhostNodej = weight[nodeListj]->nodeList().firstGhostNode();

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
            for (size_t ii = 0; ii != Dimension::nDim; ++ii) {
              gradm0(nodeListi, i)(ii) += wj*gradWj(ii);
              gradm0(nodeListj, j)(ii) += wi*gradWi(ii);
            }

            // First moment. 
            for (size_t ii = 0; ii != Dimension::nDim; ++ii) {
              m1(nodeListi, i)(ii) += wwj * rij(ii);
              m1(nodeListj, j)(ii) -= wwi * rij(ii);
              for (size_t jj = 0; jj != Dimension::nDim; ++jj) {
                gradm1(nodeListi, i)(ii,jj) += wj*rij(ii)*gradWj(jj);
                gradm1(nodeListj, j)(ii,jj) -= wi*rij(ii)*gradWi(jj);
              }
              gradm1(nodeListi, i)(ii,ii) += wj*Wj;
              gradm1(nodeListj, j)(ii,ii) += wi*Wi;
            }

            // Second moment.
            for (size_t ii = 0; ii != Dimension::nDim; ++ii) {
              for (size_t jj = ii; jj != Dimension::nDim; ++jj) {   // 'cause m2 is a symmetric tensor.
                thpt = rij(ii)*rij(jj);
                m2(nodeListi,i)(ii,jj) += wwj*thpt;
                m2(nodeListj,j)(ii,jj) += wwi*thpt;
              }
              for (size_t jj = 0; jj != Dimension::nDim; ++jj) {
                thpt = rij(ii)*rij(jj);
                for (size_t kk = 0; kk != Dimension::nDim; ++kk) {
                  gradm2(nodeListi, i)(ii,jj,kk) += wj*thpt*gradWj(kk);
                  gradm2(nodeListj, j)(ii,jj,kk) += wi*thpt*gradWi(kk);
                }
                gradm2(nodeListi, i)(ii, jj, jj) += wwj*rij(ii);
                gradm2(nodeListi, i)(jj, ii, jj) += wwj*rij(ii);
                gradm2(nodeListj, j)(ii, jj, jj) -= wwi*rij(ii);
                gradm2(nodeListj, j)(jj, ii, jj) -= wwi*rij(ii);
              }
            }

            // We only need the next moments if doing quadratic CRK.  We avoid it otherwise
            // since this is a lot of memory and expense.
            if (correctionOrder == QuadraticOrder) {

              // Third Moment
              for (size_t ii = 0; ii != Dimension::nDim; ++ii) {
                for (size_t jj = 0; jj != Dimension::nDim; ++jj) {
                  for (size_t kk = 0; kk != Dimension::nDim; ++kk) {
                    thpt = rij(ii)*rij(jj)*rij(kk);
                    m3(nodeListi,i)(ii,jj,kk) += wwj*thpt;
                    m3(nodeListj,j)(ii,jj,kk) -= wwi*thpt;
                    for (size_t mm = 0; mm != Dimension::nDim; ++mm) {
                      gradm3(nodeListi,i)(ii,jj,kk,mm) += wj*thpt*gradWj(mm);
                      gradm3(nodeListj,j)(ii,jj,kk,mm) -= wi*thpt*gradWi(mm);
                    }
                    gradm3(nodeListi,i)(ii, jj, kk, kk) += wwj*rij(ii)*rij(jj);
                    gradm3(nodeListi,i)(ii, jj, kk, jj) += wwj*rij(ii)*rij(kk);
                    gradm3(nodeListi,i)(ii, jj, kk, ii) += wwj*rij(jj)*rij(kk);
                    gradm3(nodeListj,j)(ii, jj, kk, kk) += wwi*rij(ii)*rij(jj);
                    gradm3(nodeListj,j)(ii, jj, kk, jj) += wwi*rij(ii)*rij(kk);
                    gradm3(nodeListj,j)(ii, jj, kk, ii) += wwi*rij(jj)*rij(kk);
                  }
                }
              }

              // Fourth Moment
              for (size_t ii = 0; ii != Dimension::nDim; ++ii) {
                for (size_t jj = 0; jj != Dimension::nDim; ++jj) {
                  for (size_t kk = 0; kk != Dimension::nDim; ++kk) {
                    for (size_t mm = 0; mm != Dimension::nDim; ++mm) {
                      thpt = rij(ii)*rij(jj)*rij(kk)*rij(mm);
                      m4(nodeListi,i)(ii,jj,kk,mm) += wwj*thpt;
                      m4(nodeListj,j)(ii,jj,kk,mm) += wwi*thpt;
                      for (size_t nn = 0; nn != Dimension::nDim; ++nn) {
                        gradm4(nodeListi,i)(ii,jj,kk,mm,nn) += wj*thpt*gradWj(nn);
                        gradm4(nodeListj,j)(ii,jj,kk,mm,nn) += wi*thpt*gradWi(nn);
                      }
                      gradm4(nodeListi,i)(ii, jj, kk, mm, mm) += wwj*rij(ii)*rij(jj)*rij(kk);
                      gradm4(nodeListi,i)(ii, jj, kk, mm, kk) += wwj*rij(ii)*rij(jj)*rij(mm);
                      gradm4(nodeListi,i)(ii, jj, kk, mm, jj) += wwj*rij(ii)*rij(kk)*rij(mm);
                      gradm4(nodeListi,i)(ii, jj, kk, mm, ii) += wwj*rij(jj)*rij(kk)*rij(mm);

                      gradm4(nodeListj,j)(ii, jj, kk, mm, mm) -= wwi*rij(ii)*rij(jj)*rij(kk);
                      gradm4(nodeListj,j)(ii, jj, kk, mm, kk) -= wwi*rij(ii)*rij(jj)*rij(mm);
                      gradm4(nodeListj,j)(ii, jj, kk, mm, jj) -= wwi*rij(ii)*rij(kk)*rij(mm);
                      gradm4(nodeListj,j)(ii, jj, kk, mm, ii) -= wwi*rij(jj)*rij(kk)*rij(mm);
                    }
                  }
                }
              }
            }

          }
        }
      }
    }
  }
}

}//End Namespace
}

