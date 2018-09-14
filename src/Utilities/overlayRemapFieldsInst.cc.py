text = ""
if ndim != "1":
    text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Utilities/overlayRemapFields.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
template void 
overlayRemapFields(const std::vector<BoundarySpace::Boundary<Dim<%(ndim)s>>*>& boundaries,
                   const std::vector<FieldSpace::Field<Dim<%(ndim)s>, Dim<%(ndim)s>::Scalar>*>& scalarDonorFields,
                   const std::vector<FieldSpace::Field<Dim<%(ndim)s>, Dim<%(ndim)s>::Vector>*>& vectorDonorFields,
                   const std::vector<FieldSpace::Field<Dim<%(ndim)s>, Dim<%(ndim)s>::Tensor>*>& tensorDonorFields,
                   const std::vector<FieldSpace::Field<Dim<%(ndim)s>, Dim<%(ndim)s>::SymTensor>*>& symTensorDonorFields,
                   std::vector<FieldSpace::Field<Dim<%(ndim)s>, Dim<%(ndim)s>::Scalar>*>& scalarAcceptorFields,
                   std::vector<FieldSpace::Field<Dim<%(ndim)s>, Dim<%(ndim)s>::Vector>*>& vectorAcceptorFields,
                   std::vector<FieldSpace::Field<Dim<%(ndim)s>, Dim<%(ndim)s>::Tensor>*>& tensorAcceptorFields,
                   std::vector<FieldSpace::Field<Dim<%(ndim)s>, Dim<%(ndim)s>::SymTensor>*>& symTensorAcceptorFields);
}

"""
