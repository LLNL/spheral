//---------------------------------Spheral++------------------------------------
// Compute the CRKSPH mass density summation.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeCRKSPHSumMassDensity__
#define __Spheral__computeCRKSPHSumMassDensity__

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

  namespace CRKSPHSpace {

    template<typename Dimension>
    void
    computeCRKSPHSumMassDensity(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                              const KernelSpace::TableKernel<Dimension>& W,
                              const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                              const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& mass,
                              const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                              const typename std::vector<BoundarySpace::Boundary<Dimension>*>::const_iterator& boundaryBegin,
                              const typename std::vector<BoundarySpace::Boundary<Dimension>*>::const_iterator& boundaryEnd,
                              FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& massDensity);
  }
}

#endif
