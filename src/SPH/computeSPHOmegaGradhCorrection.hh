//---------------------------------Spheral++------------------------------------
// Compute the SPH grad h correction due to Springel et al.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeSPHOmegaGradhCorrection__
#define __Spheral__computeSPHOmegaGradhCorrection__

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
    computeSPHOmegaGradhCorrection(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                                   const KernelSpace::TableKernel<Dimension>& W,
                                   const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                                   const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                                   FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& omegaGradh);
  }
}

#endif
