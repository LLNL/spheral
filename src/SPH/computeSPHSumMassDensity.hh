//---------------------------------Spheral++------------------------------------
// Compute the SPH mass density summation.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeSPHSumMassDensity__
#define __Spheral__computeSPHSumMassDensity__

namespace Spheral {

  // Forward declarations.
  template<typename Dimension> class ConnectivityMap;
  template<typename Dimension, typename DataType> class FieldList;

  template<typename Dimension, typename KernelType>
  void
  computeSPHSumMassDensity(const ConnectivityMap<Dimension>& connectivityMap,
                           const KernelType& W,
                           const bool sumOverAllNodeLists,
                           const FieldList<Dimension, typename Dimension::Vector>& position,
                           const FieldList<Dimension, typename Dimension::Scalar>& mass,
                           const FieldList<Dimension, typename Dimension::SymTensor>& H,
                           FieldList<Dimension, typename Dimension::Scalar>& massDensity);

}

#endif
