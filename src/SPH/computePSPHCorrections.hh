//---------------------------------Spheral++------------------------------------
// Compute the PSPH correction due to Hopkins 2013
//------------------------------------------------------------------------------
#ifndef __Spheral__computePSPHCorrections__
#define __Spheral__computePSPHCorrections__

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
    computePSPHCorrections(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                           const KernelSpace::TableKernel<Dimension>& W,
                           const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& mass,
                           const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                           const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& specificThermalEnergy,
                           const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& gamma,
                           const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                           FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& PSPHpbar,
                           FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& PSPHsoundSpeed,
                           FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& PSPHcorrection);
  }
}

#endif
