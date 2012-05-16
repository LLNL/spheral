//---------------------------------Spheral++------------------------------------
// Compute the SPH pressure per point that will produce the input accelerations
// given our standard SPH momentum equation.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeSPHHydrostaticEquilibriumPressure__
#define __Spheral__computeSPHHydrostaticEquilibriumPressure__

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
    computeSPHHydrostaticEquilibriumPressure(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                                             const KernelSpace::TableKernel<Dimension>& W,
                                             const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                                             const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& mass,
                                             const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                                             const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& acceleration,
                                             const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& massDensity,
                                             FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& pressure);

  }
}

#endif
