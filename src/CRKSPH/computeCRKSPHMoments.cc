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

using std::vector;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

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

  // Connectivity
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

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

  const auto W0 = W(0.0, 1.0);

  // Self contributions.
  for (auto nodeListi = 0; nodeListi < numNodeLists; ++nodeListi) {
    const auto ni = m0[nodeListi]->numInternalElements();
#pragma omp parallel for
    for (auto i = 0; i < ni; ++i) {
      const auto  weighti = weight(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      const auto  wwi = weighti*Hdeti*W0;
      m0(nodeListi, i) += wwi;
      gradm1(nodeListi, i) += Tensor::one*wwi;
    }
  }

  // Walk all the interacting pairs.
#pragma omp parallel
  {
    // Thread private scratch variables
    int i, j, nodeListi, nodeListj;
    Scalar Wi, gWi, Wj, gWj, thpt;
    Vector gradWi, gradWj;

    typename SpheralThreads<Dimension>::FieldListStack threadStack, quadThreadStack;
    auto m0_thread = m0.threadCopy(threadStack);
    auto m1_thread = m1.threadCopy(threadStack);
    auto m2_thread = m2.threadCopy(threadStack);
    auto m3_thread = m3.threadCopy(quadThreadStack);
    auto m4_thread = m4.threadCopy(quadThreadStack);
    auto gradm0_thread = gradm0.threadCopy(threadStack);
    auto gradm1_thread = gradm1.threadCopy(threadStack);
    auto gradm2_thread = gradm2.threadCopy(threadStack);
    auto gradm3_thread = gradm3.threadCopy(quadThreadStack);
    auto gradm4_thread = gradm4.threadCopy(quadThreadStack);

#pragma omp for
    for (auto kk = 0; kk < npairs; ++kk) {
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      // Get the state for node i.
      const auto  weighti = weight(nodeListi, i);
      const auto& ri = position(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();

      // State of node j.
      const auto  weightj = weight(nodeListj, j);
      const auto& rj = position(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto  Hdetj = Hj.Determinant();

      // Find the pair weighting scaling.
      const double fij = nodeCoupling(nodeListi, i, nodeListj, j);
      CHECK(fij >= 0.0 and fij <= 1.0);

      // Find the effective weights of i->j and j->i.
      // const Scalar wi = 2.0*weighti*weightj/(weighti + weightj);
      // const Scalar wi = 0.5*(weighti + weightj);
      // const Scalar wj = wi;
      const auto wi = fij*weighti;
      const auto wj = fij*weightj;

      // Kernel weighting and gradient.
      const auto rij = ri - rj;
      const auto etai = Hi*rij;
      const auto etaj = Hj*rij;

      std::tie(Wi, gWi) = W.kernelAndGradValue(etai.magnitude(), Hdeti);
      std::tie(Wj, gWj) = W.kernelAndGradValue(etaj.magnitude(), Hdetj);

      // j
      gradWj =  (Hj*etaj.unitVector())*gWj;
      gradWi = -(Hi*etai.unitVector())*gWi;
            
      // // i
      // std::swap(Wi, Wj);
      // gradWj = -(Hi*etai.unitVector())*gWi;
      // gradWi =  (Hj*etaj.unitVector())*gWj

      // // ij
      // gradWj = 0.5*((Hj*etaj.unitVector())*gWj +
      //               (Hi*etai.unitVector())*gWi);
      // gradWi = -gradWj;

      // Zeroth moment. 
      const auto wwi = wi*Wi;
      const auto wwj = wj*Wj;
      m0_thread(nodeListi, i) += wwj;
      m0_thread(nodeListj, j) += wwi;
      for (auto ii = 0; ii != Dimension::nDim; ++ii) {
        gradm0_thread(nodeListi, i)(ii) += wj*gradWj(ii);
        gradm0_thread(nodeListj, j)(ii) += wi*gradWi(ii);
      }

      // First moment. 
      for (auto ii = 0; ii != Dimension::nDim; ++ii) {
        m1_thread(nodeListi, i)(ii) += wwj * rij(ii);
        m1_thread(nodeListj, j)(ii) -= wwi * rij(ii);
        for (auto jj = 0; jj != Dimension::nDim; ++jj) {
          gradm1_thread(nodeListi, i)(ii,jj) += wj*rij(ii)*gradWj(jj);
          gradm1_thread(nodeListj, j)(ii,jj) -= wi*rij(ii)*gradWi(jj);
        }
        gradm1_thread(nodeListi, i)(ii,ii) += wj*Wj;
        gradm1_thread(nodeListj, j)(ii,ii) += wi*Wi;
      }

      // Second moment.
      for (auto ii = 0; ii != Dimension::nDim; ++ii) {
        for (auto jj = ii; jj != Dimension::nDim; ++jj) {   // 'cause m2 is a symmetric tensor.
          thpt = rij(ii)*rij(jj);
          m2_thread(nodeListi,i)(ii,jj) += wwj*thpt;
          m2_thread(nodeListj,j)(ii,jj) += wwi*thpt;
        }
        for (auto jj = 0; jj != Dimension::nDim; ++jj) {
          thpt = rij(ii)*rij(jj);
          for (auto kk = 0; kk != Dimension::nDim; ++kk) {
            gradm2_thread(nodeListi, i)(ii,jj,kk) += wj*thpt*gradWj(kk);
            gradm2_thread(nodeListj, j)(ii,jj,kk) += wi*thpt*gradWi(kk);
          }
          gradm2_thread(nodeListi, i)(ii, jj, jj) += wwj*rij(ii);
          gradm2_thread(nodeListi, i)(jj, ii, jj) += wwj*rij(ii);
          gradm2_thread(nodeListj, j)(ii, jj, jj) -= wwi*rij(ii);
          gradm2_thread(nodeListj, j)(jj, ii, jj) -= wwi*rij(ii);
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
              m3_thread(nodeListi,i)(ii,jj,kk) += wwj*thpt;
              m3_thread(nodeListj,j)(ii,jj,kk) -= wwi*thpt;
              for (auto mm = 0; mm != Dimension::nDim; ++mm) {
                gradm3_thread(nodeListi,i)(ii,jj,kk,mm) += wj*thpt*gradWj(mm);
                gradm3_thread(nodeListj,j)(ii,jj,kk,mm) -= wi*thpt*gradWi(mm);
              }
              gradm3_thread(nodeListi,i)(ii, jj, kk, kk) += wwj*rij(ii)*rij(jj);
              gradm3_thread(nodeListi,i)(ii, jj, kk, jj) += wwj*rij(ii)*rij(kk);
              gradm3_thread(nodeListi,i)(ii, jj, kk, ii) += wwj*rij(jj)*rij(kk);
              gradm3_thread(nodeListj,j)(ii, jj, kk, kk) += wwi*rij(ii)*rij(jj);
              gradm3_thread(nodeListj,j)(ii, jj, kk, jj) += wwi*rij(ii)*rij(kk);
              gradm3_thread(nodeListj,j)(ii, jj, kk, ii) += wwi*rij(jj)*rij(kk);
            }
          }
        }

        // Fourth Moment
        for (auto ii = 0; ii != Dimension::nDim; ++ii) {
          for (auto jj = 0; jj != Dimension::nDim; ++jj) {
            for (auto kk = 0; kk != Dimension::nDim; ++kk) {
              for (auto mm = 0; mm != Dimension::nDim; ++mm) {
                thpt = rij(ii)*rij(jj)*rij(kk)*rij(mm);
                m4_thread(nodeListi,i)(ii,jj,kk,mm) += wwj*thpt;
                m4_thread(nodeListj,j)(ii,jj,kk,mm) += wwi*thpt;
                for (auto nn = 0; nn != Dimension::nDim; ++nn) {
                  gradm4_thread(nodeListi,i)(ii,jj,kk,mm,nn) += wj*thpt*gradWj(nn);
                  gradm4_thread(nodeListj,j)(ii,jj,kk,mm,nn) += wi*thpt*gradWi(nn);
                }
                gradm4_thread(nodeListi,i)(ii, jj, kk, mm, mm) += wwj*rij(ii)*rij(jj)*rij(kk);
                gradm4_thread(nodeListi,i)(ii, jj, kk, mm, kk) += wwj*rij(ii)*rij(jj)*rij(mm);
                gradm4_thread(nodeListi,i)(ii, jj, kk, mm, jj) += wwj*rij(ii)*rij(kk)*rij(mm);
                gradm4_thread(nodeListi,i)(ii, jj, kk, mm, ii) += wwj*rij(jj)*rij(kk)*rij(mm);

                gradm4_thread(nodeListj,j)(ii, jj, kk, mm, mm) -= wwi*rij(ii)*rij(jj)*rij(kk);
                gradm4_thread(nodeListj,j)(ii, jj, kk, mm, kk) -= wwi*rij(ii)*rij(jj)*rij(mm);
                gradm4_thread(nodeListj,j)(ii, jj, kk, mm, jj) -= wwi*rij(ii)*rij(kk)*rij(mm);
                gradm4_thread(nodeListj,j)(ii, jj, kk, mm, ii) -= wwi*rij(jj)*rij(kk)*rij(mm);
              }
            }
          }
        }
      }
    }

    // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);
    threadReduceFieldLists<Dimension>(quadThreadStack);

  }   // OMP parallel
}

}//End Namespace
