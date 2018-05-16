//------------------------------------------------------------------------------
// Use geometric clipping to remap a set of conserved fields.
// Currently only works single NodeList -> single NodeList, no boundaries.
//------------------------------------------------------------------------------
#ifndef __Spheral_overlayRemapFields__
#define __Spheral_overlayRemapFields__

#include <vector>
#include "Boundary/Boundary.hh"
#include "Field/Field.hh"

namespace Spheral {

template<typename Dimension>
void
overlayRemapFields(const std::vector<BoundarySpace::Boundary<Dimension>*>& boundaries,
                   const std::vector<FieldSpace::Field<Dimension, typename Dimension::Scalar>*>& scalarDonorFields,
                   const std::vector<FieldSpace::Field<Dimension, typename Dimension::Vector>*>& vectorDonorFields,
                   const std::vector<FieldSpace::Field<Dimension, typename Dimension::Tensor>*>& tensorDonorFields,
                   const std::vector<FieldSpace::Field<Dimension, typename Dimension::SymTensor>*>& symTensorDonorFields,
                   std::vector<FieldSpace::Field<Dimension, typename Dimension::Scalar>*>& scalarAcceptorFields,
                   std::vector<FieldSpace::Field<Dimension, typename Dimension::Vector>*>& vectorAcceptorFields,
                   std::vector<FieldSpace::Field<Dimension, typename Dimension::Tensor>*>& tensorAcceptorFields,
                   std::vector<FieldSpace::Field<Dimension, typename Dimension::SymTensor>*>& symTensorAcceptorFields);

}

#endif
