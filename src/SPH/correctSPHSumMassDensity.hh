//---------------------------------Spheral++------------------------------------
// Apply a single corrective pass to the SPH mass density summation.
//------------------------------------------------------------------------------
#ifndef __Spheral__correctSPHSumMassDensity__
#define __Spheral__correctSPHSumMassDensity__

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

  namespace SPHSpace {

    template<typename Dimension>
    void
    correctSPHSumMassDensity(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                             const KernelSpace::TableKernel<Dimension>& W,
                             const bool sumOverAllNodeLists,
                             const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                             const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& mass,
                             const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                             FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& massDensity);
  }
}

#endif
