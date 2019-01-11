//---------------------------------Spheral++------------------------------------
// Compute the SPH pressure per point that will produce the input accelerations
// given our standard SPH momentum equation.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeSPHHydrostaticEquilibriumPressure__
#define __Spheral__computeSPHHydrostaticEquilibriumPressure__

namespace Spheral {

  // Forward declarations.
  template<typename Dimension> class ConnectivityMap;
  template<typename Dimension> class DataBase;
  template<typename Dimension> class TableKernel;
  template<typename Dimension, typename DataType> class FieldList;

  void computeSPHHydrostaticEquilibriumPressure(const DataBase<Dim<3> >& db,
                                                const TableKernel<Dim<3> >& W,
                                                const FieldList<Dim<3>, Dim<3>::Vector>& acceleration,
                                                FieldList<Dim<3>, Dim<3>::Scalar>& pressure);

}

#endif
