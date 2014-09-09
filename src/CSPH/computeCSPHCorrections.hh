//---------------------------------Spheral++------------------------------------
// Compute the CSPH corrections.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeCSPHCorrections__
#define __Spheral__computeCSPHCorrections__

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

  namespace CSPHSpace {

    template<typename Dimension>
    void
    computeCSPHCorrections(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                           const KernelSpace::TableKernel<Dimension>& W,
                           const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& weight,
                           const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                           const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                           FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& A0,
                           FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& A,
                           FieldSpace::FieldList<Dimension, typename Dimension::Vector>& B,
                           FieldSpace::FieldList<Dimension, typename Dimension::Vector>& C,
                           FieldSpace::FieldList<Dimension, typename Dimension::Tensor>& D,
                           FieldSpace::FieldList<Dimension, typename Dimension::Vector>& gradA,
                           FieldSpace::FieldList<Dimension, typename Dimension::Tensor>& gradB);
  }
}

#endif
