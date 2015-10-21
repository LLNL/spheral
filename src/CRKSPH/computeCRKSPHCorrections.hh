//---------------------------------Spheral++------------------------------------
// Compute the CRKSPH corrections.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeCRKSPHCorrections__
#define __Spheral__computeCRKSPHCorrections__

#include "SolidSPH/NodeCoupling.hh"
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
    // Function to compute CRK corrections assuming full pair-wise node coupling
    template<typename Dimension>
    void
    computeCRKSPHCorrections(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                             const KernelSpace::TableKernel<Dimension>& W,
                             const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& weight,
                             const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                             const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
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
    computeZerothCRKSPHCorrections(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                             const KernelSpace::TableKernel<Dimension>& W,
                             const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& weight,
                             const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                             const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                             FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& A,
                             FieldSpace::FieldList<Dimension, typename Dimension::Tensor>& gradA);


    // Linear Correction 
    template<typename Dimension>
    void
    computeLinearCRKSPHCorrections(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                             const KernelSpace::TableKernel<Dimension>& W,
                             const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& weight,
                             const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                             const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                             FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& A,
                             FieldSpace::FieldList<Dimension, typename Dimension::Vector>& B,
                             FieldSpace::FieldList<Dimension, typename Dimension::Vector>& gradA,
                             FieldSpace::FieldList<Dimension, typename Dimension::Tensor>& gradB);
    // Quadratic Correction 
    template<typename Dimension>
    void
    computeQuadraticCRKSPHCorrections(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                             const KernelSpace::TableKernel<Dimension>& W,
                             const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& weight,
                             const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                             const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                             FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& A,
                             FieldSpace::FieldList<Dimension, typename Dimension::Vector>& B,
                             FieldSpace::FieldList<Dimension, typename Dimension::Tensor>& C,
                             FieldSpace::FieldList<Dimension, typename Dimension::Vector>& gradA,
                             FieldSpace::FieldList<Dimension, typename Dimension::Tensor>& gradB,
                             FieldSpace::FieldList<Dimension, typename Dimension::ThirdRankTensor>& gradC);

    // Version allowing arbitrary function for pair-wise node coupling.
    template<typename Dimension>
    void
    computeCRKSPHCorrections(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                             const KernelSpace::TableKernel<Dimension>& W,
                             const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& weight,
                             const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                             const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                             const NodeCoupling& nodeCoupling,
                             FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& A,
                             FieldSpace::FieldList<Dimension, typename Dimension::Vector>& B,
                             FieldSpace::FieldList<Dimension, typename Dimension::Vector>& gradA,
                             FieldSpace::FieldList<Dimension, typename Dimension::Tensor>& gradB,
                             FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& Ac,
                             FieldSpace::FieldList<Dimension, typename Dimension::Vector>& Bc,
                             FieldSpace::FieldList<Dimension, typename Dimension::Vector>& gradAc,
                             FieldSpace::FieldList<Dimension, typename Dimension::Tensor>& gradBc);

  }
}

#endif
