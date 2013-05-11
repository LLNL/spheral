//---------------------------------Spheral++----------------------------------//
// computeEigenValues
//
// A convenience compiled method to compute the eigen values & eigen vectors
// for a while field of symmetric tensors in one pass.
//----------------------------------------------------------------------------//
#ifndef __Spheral_Geometry_computeEigenValues__
#define __Spheral_Geometry_computeEigenValues__

namespace Spheral {

namespace FieldSpace {
  template<typename Dimension, typename Value> class Field;
}

template<typename Dimension>
void
computeEigenValues(const FieldSpace::Field<Dimension, typename Dimension::SymTensor>& field,
                   FieldSpace::Field<Dimension, typename Dimension::Vector>& eigenValues,
                   FieldSpace::Field<Dimension, typename Dimension::Tensor>& eigenVectors);

}

#endif
