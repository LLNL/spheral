//---------------------------------Spheral++------------------------------------
// Compute the CRKSPH mass density summation.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeSolidCRKSPHSumMassDensity__
#define __Spheral__computeSolidCRKSPHSumMassDensity__

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
    computeSolidCRKSPHSumMassDensity(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                                     const KernelSpace::TableKernel<Dimension>& W,
                                     const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                                     const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& mass,
                                     const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                                     const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& massDensity0,
                                     const NodeCoupling& nodeCoupling,
                                     const bool correctSum,
                                     FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& massDensity);
  }
}

#endif
