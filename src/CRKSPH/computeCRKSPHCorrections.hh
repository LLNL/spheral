//---------------------------------Spheral++------------------------------------
// Compute the CRKSPH corrections.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeCRKSPHCorrections__
#define __Spheral__computeCRKSPHCorrections__

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

    // Full version with arbitrary function for pair-wise node coupling.
    template<typename Dimension>
    void
    computeCRKSPHCorrections(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                             const KernelSpace::TableKernel<Dimension>& W,
                             const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& weight,
                             const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                             const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                             double (*pairWeightFunctionPtr)(const unsigned, const unsigned, const unsigned, const unsigned),
                             FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& m0,
                             FieldSpace::FieldList<Dimension, typename Dimension::Vector>& m1,
                             FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& m2,
                             FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& A0,
                             FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& A,
                             FieldSpace::FieldList<Dimension, typename Dimension::Vector>& B,
                             FieldSpace::FieldList<Dimension, typename Dimension::Vector>& gradA0,
                             FieldSpace::FieldList<Dimension, typename Dimension::Vector>& gradA,
                             FieldSpace::FieldList<Dimension, typename Dimension::Tensor>& gradB);

    // Version assuming full pair-wise node coupling.
    template<typename Dimension>
    void
    computeCRKSPHCorrections(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                             const KernelSpace::TableKernel<Dimension>& W,
                             const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& weight,
                             const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                             const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                             FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& m0,
                             FieldSpace::FieldList<Dimension, typename Dimension::Vector>& m1,
                             FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& m2,
                             FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& A0,
                             FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& A,
                             FieldSpace::FieldList<Dimension, typename Dimension::Vector>& B,
                             FieldSpace::FieldList<Dimension, typename Dimension::Vector>& gradA0,
                             FieldSpace::FieldList<Dimension, typename Dimension::Vector>& gradA,
                             FieldSpace::FieldList<Dimension, typename Dimension::Tensor>& gradB);

  }
}

#endif
