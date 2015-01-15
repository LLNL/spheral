//---------------------------------Spheral++------------------------------------
// Compute the CSPH corrections.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeCSPHEvaluation__
#define __Spheral__computeCSPHEvaluation__

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
    computeCSPHEvaluation(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                           const KernelSpace::TableKernel<Dimension>& W,
                           const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& weight,
                           const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                           const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                           size_t nodeListi, const int i, typename Dimension::Vector reval,
                           const bool coupleNodeLists, typename Dimension::Scalar& WCSPH, typename Dimension::Vector& gradWCSPH);
  }
}

#endif
