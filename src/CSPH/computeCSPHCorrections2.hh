//---------------------------------Spheral++------------------------------------
// Compute the CSPH corrections.
// This version combines the reproducing kernel interpolation of Liu et al.
// and the linear gradient correctios of Randles & Libersy.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeCSPHCorrections2__
#define __Spheral__computeCSPHCorrections2__

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
    computeCSPHCorrections2(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                            const KernelSpace::TableKernel<Dimension>& W,
                            const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& weight,
                            const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                            const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                            const bool coupleNodeLists,
                            FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& A0,
                            FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& A,
                            FieldSpace::FieldList<Dimension, typename Dimension::Vector>& B,
                            FieldSpace::FieldList<Dimension, typename Dimension::Tensor>& M);

  }
}

#endif
