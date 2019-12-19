//---------------------------------Spheral++------------------------------------
// Compute the CRKSPH corrections.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeCRKSPHCorrections__
#define __Spheral__computeCRKSPHCorrections__

#include "CRKSPHCorrectionParams.hh"
namespace Spheral {

// Forward declarations.
template<typename Dimension> class ConnectivityMap;
template<typename Dimension> class TableKernel;
template<typename Dimension, typename DataType> class FieldList;

// Function to compute CRK corrections based on the moments.
template<typename Dimension>
void
computeCRKSPHCorrections(const FieldList<Dimension, typename Dimension::Scalar>& m0,
                         const FieldList<Dimension, typename Dimension::Vector>& m1,
                         const FieldList<Dimension, typename Dimension::SymTensor>& m2,
                         const FieldList<Dimension, typename Dimension::ThirdRankTensor>& m3,
                         const FieldList<Dimension, typename Dimension::FourthRankTensor>& m4,
                         const FieldList<Dimension, typename Dimension::Vector>& gradm0,
                         const FieldList<Dimension, typename Dimension::Tensor>& gradm1,
                         const FieldList<Dimension, typename Dimension::ThirdRankTensor>& gradm2,
                         const FieldList<Dimension, typename Dimension::FourthRankTensor>& gradm3,
                         const FieldList<Dimension, typename Dimension::FifthRankTensor>& gradm4,
                         const FieldList<Dimension, typename Dimension::SymTensor>& H,
                         const FieldList<Dimension, int>& surfacePoint,
                         const CRKOrder correctionOrder,
                         FieldList<Dimension, typename Dimension::Scalar>& A,
                         FieldList<Dimension, typename Dimension::Vector>& B,
                         FieldList<Dimension, typename Dimension::Tensor>& C,
                         FieldList<Dimension, typename Dimension::Vector>& gradA,
                         FieldList<Dimension, typename Dimension::Tensor>& gradB,
                         FieldList<Dimension, typename Dimension::ThirdRankTensor>& gradC);

// Zeroth Correction 
template<typename Dimension>
void
computeZerothCRKSPHCorrections(const FieldList<Dimension, typename Dimension::Scalar>& m0,
                               const FieldList<Dimension, typename Dimension::Vector>& gradm0,
                               FieldList<Dimension, typename Dimension::Scalar>& A,
                               FieldList<Dimension, typename Dimension::Vector>& gradA);

// Linear Correction 
template<typename Dimension>
void
computeLinearCRKSPHCorrections(const FieldList<Dimension, typename Dimension::Scalar>& m0,
                               const FieldList<Dimension, typename Dimension::Vector>& m1,
                               const FieldList<Dimension, typename Dimension::SymTensor>& m2,
                               const FieldList<Dimension, typename Dimension::Vector>& gradm0,
                               const FieldList<Dimension, typename Dimension::Tensor>& gradm1,
                               const FieldList<Dimension, typename Dimension::ThirdRankTensor>& gradm2,
                               const FieldList<Dimension, typename Dimension::SymTensor>& H,
                               const FieldList<Dimension, int>& surfacePoint,
                               FieldList<Dimension, typename Dimension::Scalar>& A,
                               FieldList<Dimension, typename Dimension::Vector>& B,
                               FieldList<Dimension, typename Dimension::Vector>& gradA,
                               FieldList<Dimension, typename Dimension::Tensor>& gradB);

// Quadratic Correction 
template<typename Dimension>
void
computeQuadraticCRKSPHCorrections(const FieldList<Dimension, typename Dimension::Scalar>& m0,
                                  const FieldList<Dimension, typename Dimension::Vector>& m1,
                                  const FieldList<Dimension, typename Dimension::SymTensor>& m2,
                                  const FieldList<Dimension, typename Dimension::ThirdRankTensor>& m3,
                                  const FieldList<Dimension, typename Dimension::FourthRankTensor>& m4,
                                  const FieldList<Dimension, typename Dimension::Vector>& gradm0,
                                  const FieldList<Dimension, typename Dimension::Tensor>& gradm1,
                                  const FieldList<Dimension, typename Dimension::ThirdRankTensor>& gradm2,
                                  const FieldList<Dimension, typename Dimension::FourthRankTensor>& gradm3,
                                  const FieldList<Dimension, typename Dimension::FifthRankTensor>& gradm4,
                                  const FieldList<Dimension, typename Dimension::SymTensor>& H,
                                  const FieldList<Dimension, int>& surfacePoint,
                                  FieldList<Dimension, typename Dimension::Scalar>& A,
                                  FieldList<Dimension, typename Dimension::Vector>& B,
                                  FieldList<Dimension, typename Dimension::Tensor>& C,
                                  FieldList<Dimension, typename Dimension::Vector>& gradA,
                                  FieldList<Dimension, typename Dimension::Tensor>& gradB,
                                  FieldList<Dimension, typename Dimension::ThirdRankTensor>& gradC);
}

#endif
