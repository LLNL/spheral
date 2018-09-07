//---------------------------------Spheral++------------------------------------
// Compute the SVPH corrections.
// Based on the reproducing kernel ideas, similar to CRKSPH.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeSVPHCorrections__
#define __Spheral__computeSVPHCorrections__

#include "Geometry/Dimension.hh"

namespace Spheral {

  // Forward declarations.
  template<typename Dimension> class ConnectivityMap;
  template<typename Dimension> class TableKernel;
  template<typename Dimension, typename DataType> class Field;
  template<typename Dimension, typename DataType> class FieldList;

  // Single NodeList version.
  template<typename Dimension>
  void
  computeSVPHCorrections(const ConnectivityMap<Dimension>& connectivityMap,
                         const TableKernel<Dimension>& W,
                         const FieldList<Dimension, typename Dimension::Scalar>& volume,
                         const FieldList<Dimension, typename Dimension::Vector>& position,
                         const FieldList<Dimension, typename Dimension::SymTensor>& H,
                         Field<Dimension, typename Dimension::Scalar>& A,
                         Field<Dimension, typename Dimension::Vector>& B,
                         Field<Dimension, typename Dimension::Tensor>& gradB);

  // Full FieldList version.
  template<typename Dimension>
  void
  computeSVPHCorrections(const ConnectivityMap<Dimension>& connectivityMap,
                         const TableKernel<Dimension>& W,
                         const FieldList<Dimension, typename Dimension::Scalar>& volume,
                         const FieldList<Dimension, typename Dimension::Vector>& position,
                         const FieldList<Dimension, typename Dimension::SymTensor>& H,
                         FieldList<Dimension, typename Dimension::Scalar>& A,
                         FieldList<Dimension, typename Dimension::Vector>& B,
                         FieldList<Dimension, typename Dimension::Tensor>& gradB);

  // // Specializations.
  // template<>
  // void
  // computeSVPHCorrections<Dim<2> >(const ConnectivityMap<Dim<2> >& connectivityMap,
  //                                 const TableKernel<Dim<2> >& W,
  //                                 const FieldList<Dim<2> ,  Dim<2>::Scalar>& volume,
  //                                 const FieldList<Dim<2> ,  Dim<2>::Vector>& position,
  //                                 const FieldList<Dim<2> ,  Dim<2>::SymTensor>& H,
  //                                 FieldList<Dim<2> ,  Dim<2>::Scalar>& A,
  //                                 FieldList<Dim<2> ,  Dim<2>::Vector>& B,
  //                                 FieldList<Dim<2> ,  Dim<2>::Tensor>& gradB);
}

#endif
