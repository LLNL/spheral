//---------------------------------Spheral++------------------------------------
// Compute the TaylorSPH corrections.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeTaylorSPHCorrections__
#define __Spheral__computeTaylorSPHCorrections__

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

  namespace TaylorSPHSpace {

    template<typename Dimension>
    void
    computeTaylorSPHCorrections(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                           const KernelSpace::TableKernel<Dimension>& W,
                           const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& weight,
                           const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                           const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                           FieldSpace::FieldList<Dimension, typename Dimension::Tensor>& D);
  }
}

#endif
