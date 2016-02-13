//---------------------------------Spheral++------------------------------------
// Compute the hull mass density summation.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeHullSumMassDensity__
#define __Spheral__computeHullSumMassDensity__

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
  namespace BoundarySpace {
    template<typename Dimension> class Boundary;
  }
  class NodeCoupling;

  namespace CRKSPHSpace {
    template<typename Dimension>
    void
    computeHullSumMassDensity(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                              const KernelSpace::TableKernel<Dimension>& W,
                              const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                              const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& mass,
                              const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                              const NodeCoupling& nodeCoupling,
                              FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& massDensity);
  }
}

#endif
