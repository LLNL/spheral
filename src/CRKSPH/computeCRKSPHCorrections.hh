//---------------------------------Spheral++------------------------------------
// Compute the CRKSPH corrections.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeCRKSPHCorrections__
#define __Spheral__computeCRKSPHCorrections__

#include "CRKSPHCorrectionParams.hh"
namespace Spheral {

  // Forward declarations.
  namespace NeighborSpace {
    template<typename Dimension> class ConnectivityMap;
  }
  namespace KernelSpace {
    template<typename Dimension> class TableKernel;
  }
  namespace FieldSpace {
    template<typename Dimension, typename DataType> class FieldList;
  }

  namespace CRKSPHSpace {

    // Function to compute CRK corrections based on the moments.
    template<typename Dimension>
    void
    computeCRKSPHCorrections(const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& m0,
                             const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& m1,
                             const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& m2,
                             const FieldSpace::FieldList<Dimension, typename Dimension::ThirdRankTensor>& m3,
                             const FieldSpace::FieldList<Dimension, typename Dimension::FourthRankTensor>& m4,
                             const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& gradm0,
                             const FieldSpace::FieldList<Dimension, typename Dimension::Tensor>& gradm1,
                             const FieldSpace::FieldList<Dimension, typename Dimension::ThirdRankTensor>& gradm2,
                             const FieldSpace::FieldList<Dimension, typename Dimension::FourthRankTensor>& gradm3,
                             const FieldSpace::FieldList<Dimension, typename Dimension::FifthRankTensor>& gradm4,
                             const CRKOrder correctionOrder,
                             FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& A,
                             FieldSpace::FieldList<Dimension, typename Dimension::Vector>& B,
                             FieldSpace::FieldList<Dimension, typename Dimension::Tensor>& C,
                             FieldSpace::FieldList<Dimension, typename Dimension::Vector>& gradA,
                             FieldSpace::FieldList<Dimension, typename Dimension::Tensor>& gradB,
                             FieldSpace::FieldList<Dimension, typename Dimension::ThirdRankTensor>& gradC);

    // Zeroth Correction 
    template<typename Dimension>
    void
    computeZerothCRKSPHCorrections(const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& m0,
                                   const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& gradm0,
                                   FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& A,
                                   FieldSpace::FieldList<Dimension, typename Dimension::Vector>& gradA);

    // Linear Correction 
    template<typename Dimension>
    void
    computeLinearCRKSPHCorrections(const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& m0,
                                   const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& m1,
                                   const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& m2,
                                   const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& gradm0,
                                   const FieldSpace::FieldList<Dimension, typename Dimension::Tensor>& gradm1,
                                   const FieldSpace::FieldList<Dimension, typename Dimension::ThirdRankTensor>& gradm2,
                                   FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& A,
                                   FieldSpace::FieldList<Dimension, typename Dimension::Vector>& B,
                                   FieldSpace::FieldList<Dimension, typename Dimension::Vector>& gradA,
                                   FieldSpace::FieldList<Dimension, typename Dimension::Tensor>& gradB);

    // Quadratic Correction 
    template<typename Dimension>
    void
    computeQuadraticCRKSPHCorrections(const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& m0,
                                      const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& m1,
                                      const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& m2,
                                      const FieldSpace::FieldList<Dimension, typename Dimension::ThirdRankTensor>& m3,
                                      const FieldSpace::FieldList<Dimension, typename Dimension::FourthRankTensor>& m4,
                                      const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& gradm0,
                                      const FieldSpace::FieldList<Dimension, typename Dimension::Tensor>& gradm1,
                                      const FieldSpace::FieldList<Dimension, typename Dimension::ThirdRankTensor>& gradm2,
                                      const FieldSpace::FieldList<Dimension, typename Dimension::FourthRankTensor>& gradm3,
                                      const FieldSpace::FieldList<Dimension, typename Dimension::FifthRankTensor>& gradm4,
                                      FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& A,
                                      FieldSpace::FieldList<Dimension, typename Dimension::Vector>& B,
                                      FieldSpace::FieldList<Dimension, typename Dimension::Tensor>& C,
                                      FieldSpace::FieldList<Dimension, typename Dimension::Vector>& gradA,
                                      FieldSpace::FieldList<Dimension, typename Dimension::Tensor>& gradB,
                                      FieldSpace::FieldList<Dimension, typename Dimension::ThirdRankTensor>& gradC);
  }
}

#endif
