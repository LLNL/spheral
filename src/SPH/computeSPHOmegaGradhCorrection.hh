//---------------------------------Spheral++------------------------------------
// Compute the SPH grad h correction due to Springel et al.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeSPHOmegaGradhCorrection__
#define __Spheral__computeSPHOmegaGradhCorrection__

namespace Spheral {

  // Forward declarations.
  template<typename Dimension> class ConnectivityMap;
  template<typename Dimension> class TableKernel;
  template<typename Dimension, typename DataType> class FieldList;

  template<typename Dimension>
  void
  computeSPHOmegaGradhCorrection(const ConnectivityMap<Dimension>& connectivityMap,
                                 const TableKernel<Dimension>& W,
                                 const FieldList<Dimension, typename Dimension::Vector>& position,
                                 const FieldList<Dimension, typename Dimension::SymTensor>& H,
                                 FieldList<Dimension, typename Dimension::Scalar>& omegaGradh);

}

#endif
