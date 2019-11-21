text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"

#include "computeRKCorrections.cc"

namespace Spheral {

template void computeRKCorrections<Dim<%(ndim)s>>(const ConnectivityMap<Dim<%(ndim)s>>& connectivityMap,
                                                  const TableKernel<Dim<%(ndim)s>>& W,
                                                  const FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Scalar>& volume,
                                                  const FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Vector>& position,
                                                  const FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::SymTensor>& H,
                                                  const CRKOrder correctionOrder,
                                                  FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Scalar>& A,
                                                  FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Vector>& B,
                                                  FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Tensor>& C,
                                                  FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::ThirdRankTensor>& D,
                                                  FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Vector>& gradA,
                                                  FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Tensor>& gradB,
                                                  FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::ThirdRankTensor>& gradC,
                                                  FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::FourthRankTensor>& gradD,
                                                  FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Tensor>& hessA,
                                                  FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::ThirdRankTensor>& hessB,
                                                  FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::FourthRankTensor>& hessC,
                                                  FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::FifthRankTensor>& hessD);

}

"""
