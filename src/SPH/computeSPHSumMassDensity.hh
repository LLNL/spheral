//---------------------------------Spheral++------------------------------------
// Compute the SPH mass density summation.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeSPHSumMassDensity__
#define __Spheral__computeSPHSumMassDensity__

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
    computeSPHSumMassDensity(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                             const KernelSpace::TableKernel<Dimension>& W,
                             const bool sumOverAllNodeLists,
                             const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                             const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& mass,
                             const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                             FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& massDensity);
  }
}

#endif
