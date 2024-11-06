//---------------------------------Spheral++----------------------------------//
// computeEigenValues
//
// A convenience compiled method to compute the eigen values & eigen vectors
// for a while field of symmetric tensors in one pass.
//----------------------------------------------------------------------------//
#ifndef __Spheral_Geometry_computeEigenValues__
#define __Spheral_Geometry_computeEigenValues__

namespace Spheral {

template<typename Dimension, typename Value> class Field;
template<typename Dimension, typename Value> class FieldView;

template<typename Dimension>
void
computeEigenValues(const Field<Dimension, typename Dimension::SymTensor>& field,
                   Field<Dimension, typename Dimension::Vector>& eigenValues,
                   Field<Dimension, typename Dimension::Tensor>& eigenVectors);

}

#endif
