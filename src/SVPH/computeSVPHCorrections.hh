//---------------------------------Spheral++------------------------------------
// Compute the SVPH corrections.
// Based on the reproducing kernel ideas, similar to CSPH.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeSVPHCorrections__
#define __Spheral__computeSVPHCorrections__

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

  namespace SVPHSpace {

    template<typename Dimension>
    void
    computeSVPHCorrections(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                           const KernelSpace::TableKernel<Dimension>& W,
                           const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& volume,
                           const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                           const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                           FieldSpace::FieldList<Dimension, typename Dimension::Vector>& B,
                           FieldSpace::FieldList<Dimension, typename Dimension::Tensor>& gradB);
  }
}

#endif
