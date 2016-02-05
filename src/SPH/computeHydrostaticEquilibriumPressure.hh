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
  namespace DataBaseSpace {
    template<typename Dimension> class DataBase;
  }
  namespace KernelSpace {
    template<typename Dimension> class TableKernel;
  }
  namespace FieldSpace {
    template<typename Dimension, typename DataType> class FieldList;
  }

  namespace SPHSpace {

    void
    computeSPHHydrostaticEquilibriumPressure(const DataBaseSpace::DataBase<Dim<3> >& db,
                                             const KernelSpace::TableKernel<Dim<3> >& W,
                                             const FieldSpace::FieldList<Dim<3>, Dim<3>::Vector>& acceleration,
                                             FieldSpace::FieldList<Dim<3>, Dim<3>::Scalar>& pressure);

  }
}

#endif
