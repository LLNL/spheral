//---------------------------------Spheral++----------------------------------//
// computeEigenValues
//
// A convenience compiled method to compute the eigen values & eigen vectors
// for a while field of symmetric tensors in one pass.
//----------------------------------------------------------------------------//

#include "computeEigenValues.hh"
#include "Field/Field.hh"

namespace Spheral {

template<typename Dimension>
void
computeEigenValues(const FieldView<Dimension, typename Dimension::SymTensor>& field,
                   FieldView<Dimension, typename Dimension::Vector>& eigenValues,
                   FieldView<Dimension, typename Dimension::Tensor>& eigenVectors) {

  // Pre-conditions.
  VERIFY(eigenValues->nodeListPtr() == field->nodeListPtr());
  VERIFY(eigenVectors->nodeListPtr() == field->nodeListPtr());

  // Do them all one by one.
  typename Dimension::SymTensor::EigenStructType eigen;
  for (unsigned i = 0; i != field->numElements(); ++i) {
    eigen = field(i).eigenVectors();
    eigenValues(i) = eigen.eigenValues;
    eigenVectors(i) = eigen.eigenVectors;
  }
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template void computeEigenValues(const FieldView<Dim<1>, Dim<1>::SymTensor>&, FieldView<Dim<1>, Dim<1>::Vector>&, FieldView<Dim<1>, Dim<1>::Tensor>&);
template void computeEigenValues(const FieldView<Dim<2>, Dim<2>::SymTensor>&, FieldView<Dim<2>, Dim<2>::Vector>&, FieldView<Dim<2>, Dim<2>::Tensor>&);
template void computeEigenValues(const FieldView<Dim<3>, Dim<3>::SymTensor>&, FieldView<Dim<3>, Dim<3>::Vector>&, FieldView<Dim<3>, Dim<3>::Tensor>&);

}
