//---------------------------------Spheral++------------------------------------
// Compute the CSPH mass density summation.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeCSPHSumMassDensity__
#define __Spheral__computeCSPHSumMassDensity__

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

  namespace CSPHSpace {

    template<typename Dimension>
    void
    computeCSPHSumMassDensity(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
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
