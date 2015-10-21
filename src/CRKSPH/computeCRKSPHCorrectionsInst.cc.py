text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "computeCRKSPHCorrections.cc"

namespace Spheral {
namespace CRKSPHSpace {

//------------------------------------------------------------------------------
// CRKSPH Quadratic Corrections Method (tensor based Mike version).
//------------------------------------------------------------------------------
void
computeQuadraticCRKSPHCorrectionsMike(const ConnectivityMap<Dim< %(ndim)s > >& connectivityMap,
                                      const TableKernel<Dim< %(ndim)s > >& W,
                                      const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                      const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                      const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& H,
                                      FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& A,
                                      FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& B,
                                      FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>& C,
                                      FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& gradA,
                                      FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>& gradB,
                                      FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::ThirdRankTensor>& gradC) {

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

  typedef Dim<%(ndim)s> Dimension;
  typedef Dimension::Scalar Scalar;
  typedef Dimension::Vector Vector;
  typedef Dimension::Tensor Tensor;
  typedef Dimension::SymTensor SymTensor;
  typedef Dimension::ThirdRankTensor ThirdRankTensor;
  typedef Dimension::FourthRankTensor FourthRankTensor;
  typedef Dimension::FifthRankTensor FifthRankTensor;

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
    for (ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
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

//            gradm1(nodeListi, i) += wj*(outerProduct<Dimension>( rij, gradWj) + Tensor::one*Wj);
//            gradm1(nodeListj, j) += wi*(outerProduct<Dimension>(-rij, gradWi) + Tensor::one*Wi);
"""
for ii in xrange(idim):
    for jj in xrange(idim):
        if ii == jj:
            text += """
            gradm1(nodeListi, i)(%(ii)i,%(jj)i) += wj*( rij(%(ii)i)*gradWj(%(jj)i) + Wj);
            gradm1(nodeListj, j)(%(ii)i,%(jj)i) += wi*(-rij(%(ii)i)*gradWi(%(jj)i) + Wi);
""" % {"ii" : ii, "jj" : jj}
        else:
            text += """
            gradm1(nodeListi, i)(%(ii)i,%(jj)i) +=  wj*rij(%(ii)i)*gradWj(%(jj)i);
            gradm1(nodeListj, j)(%(ii)i,%(jj)i) += -wi*rij(%(ii)i)*gradWi(%(jj)i);
""" % {"ii" : ii, "jj" : jj}

text += """
            // Second moment.
            const SymTensor thpt = rij.selfdyad();
            m2(nodeListi, i) += wwj*thpt;
            m2(nodeListj, j) += wwi*thpt;
"""

for ii in xrange(idim):
    for jj in xrange(idim):
        for kk in xrange(idim):
            text += """
            gradm2(nodeListi, i)(%(ii)i,%(jj)i,%(kk)i) += wj*thpt(%(ii)i,%(jj)i)*gradWj(%(kk)i);
            gradm2(nodeListj, j)(%(ii)i,%(jj)i,%(kk)i) += wi*thpt(%(ii)i,%(jj)i)*gradWi(%(kk)i);
""" % {"ii": ii, "jj": jj, "kk": kk}
        text += """
            gradm2(nodeListi, i)(%(ii)i, %(jj)i, %(jj)i) += wwj*rij(%(ii)i);
            gradm2(nodeListi, i)(%(jj)i, %(ii)i, %(jj)i) += wwj*rij(%(ii)i);

            gradm2(nodeListj, j)(%(ii)i, %(jj)i, %(jj)i) -= wwi*rij(%(ii)i);
            gradm2(nodeListj, j)(%(jj)i, %(ii)i, %(jj)i) -= wwi*rij(%(ii)i);
""" % {"ii": ii, "jj": jj}

text += """
            // Third Moment
            ThirdRankTensor thpt3 = outerProduct<Dimension>(thpt,rij);
            m3(nodeListi, i) += wwj*thpt3;
            m3(nodeListj, j) -= wwi*thpt3;
"""

for ii in xrange(idim):
    for jj in xrange(idim):
        for kk in xrange(idim):
            for mm in xrange(idim):
                text += """
            gradm3(nodeListi, i)(%(ii)i,%(jj)i,%(kk)i,%(mm)i) += wj*thpt3(%(ii)i,%(jj)i,%(kk)i)*gradWj(%(mm)i);
            gradm3(nodeListj, j)(%(ii)i,%(jj)i,%(kk)i,%(mm)i) -= wi*thpt3(%(ii)i,%(jj)i,%(kk)i)*gradWi(%(mm)i);
""" % {"ii": ii, "jj": jj, "kk": kk, "mm": mm}

            text += """
            gradm3(nodeListi, i)(%(ii)i, %(jj)i, %(kk)i, %(kk)i) += wwj*rij(%(ii)i)*rij(%(jj)i);
            gradm3(nodeListi, i)(%(ii)i, %(jj)i, %(kk)i, %(jj)i) += wwj*rij(%(ii)i)*rij(%(kk)i);
            gradm3(nodeListi, i)(%(ii)i, %(jj)i, %(kk)i, %(ii)i) += wwj*rij(%(jj)i)*rij(%(kk)i);

            gradm3(nodeListj, j)(%(ii)i, %(jj)i, %(kk)i, %(kk)i) += wwi*rij(%(ii)i)*rij(%(jj)i);
            gradm3(nodeListj, j)(%(ii)i, %(jj)i, %(kk)i, %(jj)i) += wwi*rij(%(ii)i)*rij(%(kk)i);
            gradm3(nodeListj, j)(%(ii)i, %(jj)i, %(kk)i, %(ii)i) += wwi*rij(%(jj)i)*rij(%(kk)i);
""" % {"ii": ii, "jj": jj, "kk": kk}

text += """
            // Fourth Moment
            FourthRankTensor thpt4 = outerProduct<Dimension>(thpt3,rij);
            m4(nodeListi, i) += wwj*thpt4;
            m4(nodeListj, j) += wwi*thpt4;
"""

for ii in xrange(idim):
    for jj in xrange(idim):
        for kk in xrange(idim):
            for mm in xrange(idim):
                for qq in xrange(idim):
                    text += """
            gradm4(nodeListi, i)(%(ii)i,%(jj)i,%(kk)i,%(mm)i,%(qq)i) += wj*thpt4(%(ii)i,%(jj)i,%(kk)i,%(mm)i)*gradWj(%(qq)i);
            gradm4(nodeListj, j)(%(ii)i,%(jj)i,%(kk)i,%(mm)i,%(qq)i) += wi*thpt4(%(ii)i,%(jj)i,%(kk)i,%(mm)i)*gradWi(%(qq)i);
""" % {"ii": ii, "jj": jj, "kk": kk, "mm": mm, "qq": qq}

                text += """
            gradm4(nodeListi, i)(%(ii)i, %(jj)i, %(kk)i, %(mm)i, %(mm)i) += wwj*rij(%(ii)i)*rij(%(jj)i)*rij(%(kk)i);
            gradm4(nodeListi, i)(%(ii)i, %(jj)i, %(kk)i, %(mm)i, %(kk)i) += wwj*rij(%(ii)i)*rij(%(jj)i)*rij(%(mm)i);
            gradm4(nodeListi, i)(%(ii)i, %(jj)i, %(kk)i, %(mm)i, %(jj)i) += wwj*rij(%(ii)i)*rij(%(kk)i)*rij(%(mm)i);
            gradm4(nodeListi, i)(%(ii)i, %(jj)i, %(kk)i, %(mm)i, %(ii)i) += wwj*rij(%(jj)i)*rij(%(kk)i)*rij(%(mm)i);

            gradm4(nodeListj, j)(%(ii)i, %(jj)i, %(kk)i, %(mm)i, %(mm)i) -= wwi*rij(%(ii)i)*rij(%(jj)i)*rij(%(kk)i);
            gradm4(nodeListj, j)(%(ii)i, %(jj)i, %(kk)i, %(mm)i, %(kk)i) -= wwi*rij(%(ii)i)*rij(%(jj)i)*rij(%(mm)i);
            gradm4(nodeListj, j)(%(ii)i, %(jj)i, %(kk)i, %(mm)i, %(jj)i) -= wwi*rij(%(ii)i)*rij(%(kk)i)*rij(%(mm)i);
            gradm4(nodeListj, j)(%(ii)i, %(jj)i, %(kk)i, %(mm)i, %(ii)i) -= wwi*rij(%(jj)i)*rij(%(kk)i)*rij(%(mm)i);
""" % {"ii": ii, "jj": jj, "kk": kk, "mm": mm}

text += """
          }
        }
      }

      // Based on the moments we can calculate the CRKSPH corrections terms and their gradients.
      if (i < firstGhostNodei) {
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

}
}


namespace Spheral {
  namespace CRKSPHSpace {
    template void computeCRKSPHCorrections(const ConnectivityMap<Dim< %(ndim)s > >&, 
                                           const TableKernel<Dim< %(ndim)s > >&, 
                                           const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                           const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                           const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>&,
                                           const CRKOrder correctionOrder,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::ThirdRankTensor>&);

    template void computeZerothCRKSPHCorrections(const ConnectivityMap<Dim< %(ndim)s > >&, 
                                           const TableKernel<Dim< %(ndim)s > >&, 
                                           const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                           const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                           const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&);

    template void computeLinearCRKSPHCorrections(const ConnectivityMap<Dim< %(ndim)s > >&, 
                                           const TableKernel<Dim< %(ndim)s > >&, 
                                           const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                           const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                           const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>&);

    template void computeCRKSPHCorrections(const ConnectivityMap<Dim< %(ndim)s > >&, 
                                           const TableKernel<Dim< %(ndim)s > >&, 
                                           const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                           const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                           const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>&,
                                           const NodeCoupling& nodeCoupling,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                                           FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>&);
  }
}

"""
