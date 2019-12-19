text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "CRKSPH/interpolateCRKSPH.cc"
#include "Geometry/Dimension.hh"
#include "SPH/NodeCoupling.hh"

namespace Spheral {

template 
std::vector<boost::variant<FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Scalar>,
                           FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Vector>,
                           FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Tensor>,
                           FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::SymTensor>,
                           FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::ThirdRankTensor>>>
interpolateCRKSPH<Dim<%(ndim)s>>(const std::vector<boost::variant<FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Scalar>,
                                                                  FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Vector>,
                                                                  FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Tensor>,
                                                                  FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::SymTensor>,
                                                                  FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::ThirdRankTensor>>>& fieldLists,
                                          const FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Vector>& position,
                                          const FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Scalar>& weight,
                                          const FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::SymTensor>& Hfield,
                                          const FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Scalar>& A,
                                          const FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Vector>& B,
                                          const FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Tensor>& C,
                                          const ConnectivityMap<Dim<%(ndim)s> >& connectivityMap,
                                          const CRKOrder correctionOrder,
                                          const TableKernel< Dim<%(ndim)s> >& kernel,
                                          const NodeCoupling& nodeCoupling);

}

"""
