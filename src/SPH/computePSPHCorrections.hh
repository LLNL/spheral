//---------------------------------Spheral++------------------------------------
// Compute the PSPH correction due to Hopkins 2013
//------------------------------------------------------------------------------
#ifndef __Spheral__computePSPHCorrections__
#define __Spheral__computePSPHCorrections__

namespace Spheral {

  // Forward declarations.
  template<typename Dimension> class ConnectivityMap;
  template<typename Dimension> class TableKernel;
  template<typename Dimension, typename DataType> class FieldList;

  template<typename Dimension>
  void
  computePSPHCorrections(const ConnectivityMap<Dimension>& connectivityMap,
                         const TableKernel<Dimension>& W,
                         const FieldList<Dimension, typename Dimension::Scalar>& mass,
                         const FieldList<Dimension, typename Dimension::Vector>& position,
                         const FieldList<Dimension, typename Dimension::Scalar>& specificThermalEnergy,
                         const FieldList<Dimension, typename Dimension::Scalar>& gamma,
                         const FieldList<Dimension, typename Dimension::SymTensor>& H,
                         const bool computeMassDensity,
                         FieldList<Dimension, typename Dimension::Scalar>& PSPHmassDensity,
                         FieldList<Dimension, typename Dimension::Scalar>& PSPHpbar,
                         FieldList<Dimension, typename Dimension::Scalar>& PSPHsoundSpeed,
                         FieldList<Dimension, typename Dimension::Scalar>& PSPHcorrection);

}

#endif
