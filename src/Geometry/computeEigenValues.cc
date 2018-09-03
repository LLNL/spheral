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
computeEigenValues(const Field<Dimension, typename Dimension::SymTensor>& field,
                   Field<Dimension, typename Dimension::Vector>& eigenValues,
                   Field<Dimension, typename Dimension::Tensor>& eigenVectors) {

  // Pre-conditions.
  VERIFY(eigenValues.nodeListPtr() == field.nodeListPtr());
  VERIFY(eigenVectors.nodeListPtr() == field.nodeListPtr());

  // Do them all one by one.
  typename Dimension::SymTensor::EigenStructType eigen;
  for (unsigned i = 0; i != field.numElements(); ++i) {
    eigen = field(i).eigenVectors();
    eigenValues(i) = eigen.eigenValues;
    eigenVectors(i) = eigen.eigenVectors;
  }
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template void computeEigenValues(const Field<Dim<1>, Dim<1>::SymTensor>&, Field<Dim<1>, Dim<1>::Vector>&, Field<Dim<1>, Dim<1>::Tensor>&);
template void computeEigenValues(const Field<Dim<2>, Dim<2>::SymTensor>&, Field<Dim<2>, Dim<2>::Vector>&, Field<Dim<2>, Dim<2>::Tensor>&);
template void computeEigenValues(const Field<Dim<3>, Dim<3>::SymTensor>&, Field<Dim<3>, Dim<3>::Vector>&, Field<Dim<3>, Dim<3>::Tensor>&);

}
