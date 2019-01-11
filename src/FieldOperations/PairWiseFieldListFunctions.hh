//---------------------------------Spheral++----------------------------------//
// PairWiseFieldListFunctions -- A set of global functions which can be applied
// to FieldLists using low order pairwise operations.
//
// Created by JMO, Sat Apr 12 15:51:01 PDT 2003
//----------------------------------------------------------------------------//
#include "Geometry/MathTraits.hh"

namespace Spheral {

template<typename Dimension> class TableKernel;

template<typename Dimension, typename DataType> class FieldList;

// Calculate the gradient of the given FieldList.
template<typename Dimension, typename DataType>
FieldList<Dimension, typename MathTraits<Dimension, DataType>::GradientType>
gradientPairWise
(const FieldList<Dimension, DataType>& fieldList,
 const FieldList<Dimension, typename Dimension::Vector>& position,
 const FieldList<Dimension, typename Dimension::Scalar>& weight,
 const FieldList<Dimension, typename Dimension::Scalar>& mass,
 const FieldList<Dimension, typename Dimension::Scalar>& rho,
 const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
 const TableKernel<Dimension>& kernel);

// Calculate the divergence of the given FieldList.
template<typename Dimension, typename DataType>
FieldList<Dimension, typename MathTraits<Dimension, DataType>::DivergenceType>
divergencePairWise
(const FieldList<Dimension, DataType>& fieldList,
 const FieldList<Dimension, typename Dimension::Vector>& position,
 const FieldList<Dimension, typename Dimension::Scalar>& weight,
 const FieldList<Dimension, typename Dimension::Scalar>& mass,
 const FieldList<Dimension, typename Dimension::Scalar>& rho,
 const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
 const TableKernel<Dimension>& kernel);

}
