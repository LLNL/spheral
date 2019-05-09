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

using std::min;
using std::max;
using std::abs;

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
                     const NodeCoupling& nodeCoupling,
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
  if (correctionOrder == CRKOrder::QuadraticOrder) {
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
  if (correctionOrder == CRKOrder::QuadraticOrder) {
    m3 = ThirdRankTensor::zero;
    m4 = FourthRankTensor::zero;
    gradm3 = FourthRankTensor::zero;
    gradm4 = FifthRankTensor::zero;
  }

  // Walk the FluidNodeLists.
  for (auto nodeListi = 0; nodeListi < numNodeLists; ++nodeListi) {
    const auto ni = connectivityMap.numNodes(nodeListi);

    // Iterate over the nodes in this node list.
#pragma omp parallel for
    for (auto iItr = 0; iItr < ni; ++iItr) {
      const auto i = connectivityMap.ithNode(nodeListi, iItr);

      // Get the state for node i.
      const auto  weighti = weight(nodeListi, i);
      const auto& ri = position(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();

      // Self contribution.
      const auto wwi = weighti*W(0.0, Hdeti);
      m0(nodeListi, i) += wwi;
      gradm1(nodeListi, i) += Tensor::one*wwi;

      // Neighbors!
      const auto& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
      CHECK(fullConnectivity.size() == numNodeLists);

      Scalar thpt;
      for (auto nodeListj = 0; nodeListj < numNodeLists; ++nodeListj) {
        const auto& connectivity = fullConnectivity[nodeListj];

        // Iterate over the neighbors for in this NodeList.
        for (auto jItr = connectivity.begin(); jItr != connectivity.end(); ++jItr) {
          const auto j = *jItr;

          // State of node j.
          const auto  weightj = weight(nodeListj, j);
          const auto& rj = position(nodeListj, j);
          const auto& Hj = H(nodeListj, j);
          const auto  Hdetj = Hj.Determinant();

          // Find the effective weights of i->j and j->i.
          // const Scalar wi = 2.0*weighti*weightj/(weighti + weightj);
          const auto wi = 0.5*(weighti + weightj);
          const auto wj = wi;
          // const auto wi = weighti;
          // const auto wj = weightj;

          // Find the pair weighting scaling.
          const auto fij = nodeCoupling(nodeListi, i, nodeListj, j);
          CHECK(fij >= 0.0 and fij <= 1.0);

          // Kernel weighting and gradient.
          const auto rij = ri - rj;
          auto etai = Hi*rij;
          auto etaj = Hj*rij;

          const auto WWi = W.kernelAndGradValue(etai.magnitude(), Hdeti);
          const auto WWj = W.kernelAndGradValue(etaj.magnitude(), Hdetj);

          // // j
          // const Scalar Wi = WWi.first;
          // const Scalar Wj = WWj.first;
          // const Vector gradWj =  (Hj*etaj.unitVector())*WWj.second;
          // const Vector gradWi = -(Hi*etai.unitVector())*WWi.second;

          // // i
          // // const Scalar Wi = WWj.first;
          // // const Scalar Wj = WWi.first;
          // // const Vector gradWj = -(Hi*etai.unitVector())*WWi.second;
          // // const Vector gradWi =  (Hj*etaj.unitVector())*WWj.second;

          // ij
          const auto Wi = 0.5*(WWi.first + WWj.first);
          const auto Wj = Wi;
          const auto gradWj = 0.5*((Hj*etaj.unitVector())*WWj.second +
                                   (Hi*etai.unitVector())*WWi.second);
          const auto gradWi = -gradWj;

          // Zeroth moment. 
          const auto wwi = wi*Wi;
          const auto wwj = wj*Wj;
          m0(nodeListi, i) += fij*wwj;
          for (auto ii = 0; ii != Dimension::nDim; ++ii) {
            gradm0(nodeListi, i)(ii) += fij*wj*gradWj(ii);
          }

          // First moment. 
          for (auto ii = 0; ii != Dimension::nDim; ++ii) {
            m1(nodeListi, i)(ii) += fij*wwj * rij(ii);
            for (auto jj = 0; jj != Dimension::nDim; ++jj) {
              gradm1(nodeListi, i)(ii,jj) += fij*wj*rij(ii)*gradWj(jj);
            }
            gradm1(nodeListi, i)(ii,ii) += fij*wj*Wj;
          }

          // Second moment.
          for (auto ii = 0; ii != Dimension::nDim; ++ii) {
            for (auto jj = ii; jj != Dimension::nDim; ++jj) {   // 'cause m2 is a symmetric tensor.
              thpt = rij(ii)*rij(jj);
              m2(nodeListi,i)(ii,jj) += fij*wwj*thpt;
            }
            for (auto jj = 0; jj != Dimension::nDim; ++jj) {
              thpt = rij(ii)*rij(jj);
              for (auto kk = 0; kk != Dimension::nDim; ++kk) {
                gradm2(nodeListi, i)(ii,jj,kk) += fij*wj*thpt*gradWj(kk);
              }
              gradm2(nodeListi, i)(ii, jj, jj) += fij*wwj*rij(ii);
              gradm2(nodeListi, i)(jj, ii, jj) += fij*wwj*rij(ii);
            }
          }

          // We only need the next moments if doing quadratic CRK.  We avoid it otherwise
          // since this is a lot of memory and expense.
          if (correctionOrder == CRKOrder::QuadraticOrder) {

            // Third Moment
            for (auto ii = 0; ii != Dimension::nDim; ++ii) {
              for (auto jj = 0; jj != Dimension::nDim; ++jj) {
                for (auto kk = 0; kk != Dimension::nDim; ++kk) {
                  thpt = rij(ii)*rij(jj)*rij(kk);
                  m3(nodeListi,i)(ii,jj,kk) += fij*wwj*thpt;
                  for (auto mm = 0; mm != Dimension::nDim; ++mm) {
                    gradm3(nodeListi,i)(ii,jj,kk,mm) += fij*wj*thpt*gradWj(mm);
                  }
                  gradm3(nodeListi,i)(ii, jj, kk, kk) += fij*wwj*rij(ii)*rij(jj);
                  gradm3(nodeListi,i)(ii, jj, kk, jj) += fij*wwj*rij(ii)*rij(kk);
                  gradm3(nodeListi,i)(ii, jj, kk, ii) += fij*wwj*rij(jj)*rij(kk);
                }
              }
            }

            // Fourth Moment
            for (auto ii = 0; ii != Dimension::nDim; ++ii) {
              for (auto jj = 0; jj != Dimension::nDim; ++jj) {
                for (auto kk = 0; kk != Dimension::nDim; ++kk) {
                  for (auto mm = 0; mm != Dimension::nDim; ++mm) {
                    thpt = rij(ii)*rij(jj)*rij(kk)*rij(mm);
                    m4(nodeListi,i)(ii,jj,kk,mm) += fij*wwj*thpt;
                    for (auto nn = 0; nn != Dimension::nDim; ++nn) {
                      gradm4(nodeListi,i)(ii,jj,kk,mm,nn) += fij*wj*thpt*gradWj(nn);
                    }
                    gradm4(nodeListi,i)(ii, jj, kk, mm, mm) += fij*wwj*rij(ii)*rij(jj)*rij(kk);
                    gradm4(nodeListi,i)(ii, jj, kk, mm, kk) += fij*wwj*rij(ii)*rij(jj)*rij(mm);
                    gradm4(nodeListi,i)(ii, jj, kk, mm, jj) += fij*wwj*rij(ii)*rij(kk)*rij(mm);
                    gradm4(nodeListi,i)(ii, jj, kk, mm, ii) += fij*wwj*rij(jj)*rij(kk)*rij(mm);
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
