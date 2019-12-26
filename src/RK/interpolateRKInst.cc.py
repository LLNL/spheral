text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "RK/interpolateRK.cc"
#include "Geometry/Dimension.hh"
#include "SPH/NodeCoupling.hh"

namespace Spheral {

template 
std::vector<boost::variant<FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Scalar>,
                           FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Vector>,
                           FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Tensor>,
                           FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::SymTensor>,
                           FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::ThirdRankTensor>>>
interpolateRK<Dim<%(ndim)s>>(const std::vector<boost::variant<FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Scalar>,
                                                              FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Vector>,
                                                              FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Tensor>,
                                                              FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::SymTensor>,
                                                              FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::ThirdRankTensor>>>& fieldLists,
                                          const FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Vector>& position,
                                          const FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Scalar>& weight,
                                          const FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::SymTensor>& H,
                                          const ConnectivityMap<Dim<%(ndim)s> >& connectivityMap,
                                          const TableKernel< Dim<%(ndim)s> >& kernel,
                                          const RKOrder correctionOrder,
                                          const FieldList<Dim<%(ndim)s>, std::vector<double>>& corrections,
                                          const NodeCoupling& nodeCoupling);

}

"""
