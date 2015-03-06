//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "interpolateCRKSPH.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace CRKSPHSpace {

    using KernelSpace::TableKernel;
    using NeighborSpace::ConnectivityMap;

    //--------------------------------------------------------------------------
    // 1D
    //--------------------------------------------------------------------------
    template 
    FieldList<Dim<1>, Dim<1>::Scalar> 
    interpolateCRKSPH<Dim<1>, Dim<1>::Scalar>(const FieldList<Dim<1>, Dim<1>::Scalar>& fieldList,
                                            const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                            const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                            const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                            const bool coupleNodeLists,
                                            const FieldList<Dim<1>, Dim<1>::Scalar>& A,
                                            const FieldList<Dim<1>, Dim<1>::Vector>& B,
                                            const ConnectivityMap<Dim<1> >& connectivityMap,
                                            const TableKernel< Dim<1> >& kernel);
    template 
    FieldList<Dim<1>, Dim<1>::Vector> 
    interpolateCRKSPH<Dim<1>, Dim<1>::Vector>(const FieldList<Dim<1>, Dim<1>::Vector>& fieldList,
                                            const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                            const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                            const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                            const bool coupleNodeLists,
                                            const FieldList<Dim<1>, Dim<1>::Scalar>& A,
                                            const FieldList<Dim<1>, Dim<1>::Vector>& B,
                                            const ConnectivityMap<Dim<1> >& connectivityMap,
                                            const TableKernel< Dim<1> >& kernel);

    template 
    FieldList<Dim<1>, Dim<1>::Tensor> 
    interpolateCRKSPH<Dim<1>, Dim<1>::Tensor>(const FieldList<Dim<1>, Dim<1>::Tensor>& fieldList,
                                            const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                            const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                            const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                            const bool coupleNodeLists,
                                            const FieldList<Dim<1>, Dim<1>::Scalar>& A,
                                            const FieldList<Dim<1>, Dim<1>::Vector>& B,
                                            const ConnectivityMap<Dim<1> >& connectivityMap,
                                            const TableKernel< Dim<1> >& kernel);

    template 
    FieldList<Dim<1>, Dim<1>::SymTensor> 
    interpolateCRKSPH<Dim<1>, Dim<1>::SymTensor>(const FieldList<Dim<1>, Dim<1>::SymTensor>& fieldList,
                                               const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                               const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                               const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                               const bool coupleNodeLists,
                                               const FieldList<Dim<1>, Dim<1>::Scalar>& A,
                                               const FieldList<Dim<1>, Dim<1>::Vector>& B,
                                               const ConnectivityMap<Dim<1> >& connectivityMap,
                                               const TableKernel< Dim<1> >& kernel);

    template 
    FieldList<Dim<1>, Dim<1>::ThirdRankTensor> 
    interpolateCRKSPH<Dim<1>, Dim<1>::ThirdRankTensor>(const FieldList<Dim<1>, Dim<1>::ThirdRankTensor>& fieldList,
                                                     const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                                     const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                                     const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                                     const bool coupleNodeLists,
                                                     const FieldList<Dim<1>, Dim<1>::Scalar>& A,
                                                     const FieldList<Dim<1>, Dim<1>::Vector>& B,
                                                     const ConnectivityMap<Dim<1> >& connectivityMap,
                                                     const TableKernel< Dim<1> >& kernel);

    //--------------------------------------------------------------------------
    // 2D
    //--------------------------------------------------------------------------
    template 
    FieldList<Dim<2>, Dim<2>::Scalar> 
    interpolateCRKSPH<Dim<2>, Dim<2>::Scalar>(const FieldList<Dim<2>, Dim<2>::Scalar>& fieldList,
                                            const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                            const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                            const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                            const bool coupleNodeLists,
                                            const FieldList<Dim<2>, Dim<2>::Scalar>& A,
                                            const FieldList<Dim<2>, Dim<2>::Vector>& B,
                                            const ConnectivityMap<Dim<2> >& connectivityMap,
                                            const TableKernel< Dim<2> >& kernel);
    template 
    FieldList<Dim<2>, Dim<2>::Vector> 
    interpolateCRKSPH<Dim<2>, Dim<2>::Vector>(const FieldList<Dim<2>, Dim<2>::Vector>& fieldList,
                                            const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                            const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                            const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                            const bool coupleNodeLists,
                                            const FieldList<Dim<2>, Dim<2>::Scalar>& A,
                                            const FieldList<Dim<2>, Dim<2>::Vector>& B,
                                            const ConnectivityMap<Dim<2> >& connectivityMap,
                                            const TableKernel< Dim<2> >& kernel);

    template 
    FieldList<Dim<2>, Dim<2>::Tensor> 
    interpolateCRKSPH<Dim<2>, Dim<2>::Tensor>(const FieldList<Dim<2>, Dim<2>::Tensor>& fieldList,
                                            const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                            const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                            const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                            const bool coupleNodeLists,
                                            const FieldList<Dim<2>, Dim<2>::Scalar>& A,
                                            const FieldList<Dim<2>, Dim<2>::Vector>& B,
                                            const ConnectivityMap<Dim<2> >& connectivityMap,
                                            const TableKernel< Dim<2> >& kernel);

    template 
    FieldList<Dim<2>, Dim<2>::SymTensor> 
    interpolateCRKSPH<Dim<2>, Dim<2>::SymTensor>(const FieldList<Dim<2>, Dim<2>::SymTensor>& fieldList,
                                               const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                               const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                               const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                               const bool coupleNodeLists,
                                               const FieldList<Dim<2>, Dim<2>::Scalar>& A,
                                               const FieldList<Dim<2>, Dim<2>::Vector>& B,
                                               const ConnectivityMap<Dim<2> >& connectivityMap,
                                               const TableKernel< Dim<2> >& kernel);

    template 
    FieldList<Dim<2>, Dim<2>::ThirdRankTensor> 
    interpolateCRKSPH<Dim<2>, Dim<2>::ThirdRankTensor>(const FieldList<Dim<2>, Dim<2>::ThirdRankTensor>& fieldList,
                                                     const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                                     const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                                     const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                                     const bool coupleNodeLists,
                                                     const FieldList<Dim<2>, Dim<2>::Scalar>& A,
                                                     const FieldList<Dim<2>, Dim<2>::Vector>& B,
                                                     const ConnectivityMap<Dim<2> >& connectivityMap,
                                                     const TableKernel< Dim<2> >& kernel);

    //--------------------------------------------------------------------------
    // 3D
    //--------------------------------------------------------------------------
    template 
    FieldList<Dim<3>, Dim<3>::Scalar> 
    interpolateCRKSPH<Dim<3>, Dim<3>::Scalar>(const FieldList<Dim<3>, Dim<3>::Scalar>& fieldList,
                                            const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                            const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                            const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                            const bool coupleNodeLists,
                                            const FieldList<Dim<3>, Dim<3>::Scalar>& A,
                                            const FieldList<Dim<3>, Dim<3>::Vector>& B,
                                            const ConnectivityMap<Dim<3> >& connectivityMap,
                                            const TableKernel< Dim<3> >& kernel);
    template 
    FieldList<Dim<3>, Dim<3>::Vector> 
    interpolateCRKSPH<Dim<3>, Dim<3>::Vector>(const FieldList<Dim<3>, Dim<3>::Vector>& fieldList,
                                            const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                            const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                            const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                            const bool coupleNodeLists,
                                            const FieldList<Dim<3>, Dim<3>::Scalar>& A,
                                            const FieldList<Dim<3>, Dim<3>::Vector>& B,
                                            const ConnectivityMap<Dim<3> >& connectivityMap,
                                            const TableKernel< Dim<3> >& kernel);

    template 
    FieldList<Dim<3>, Dim<3>::Tensor> 
    interpolateCRKSPH<Dim<3>, Dim<3>::Tensor>(const FieldList<Dim<3>, Dim<3>::Tensor>& fieldList,
                                            const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                            const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                            const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                            const bool coupleNodeLists,
                                            const FieldList<Dim<3>, Dim<3>::Scalar>& A,
                                            const FieldList<Dim<3>, Dim<3>::Vector>& B,
                                            const ConnectivityMap<Dim<3> >& connectivityMap,
                                            const TableKernel< Dim<3> >& kernel);

    template 
    FieldList<Dim<3>, Dim<3>::SymTensor> 
    interpolateCRKSPH<Dim<3>, Dim<3>::SymTensor>(const FieldList<Dim<3>, Dim<3>::SymTensor>& fieldList,
                                               const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                               const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                               const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                               const bool coupleNodeLists,
                                               const FieldList<Dim<3>, Dim<3>::Scalar>& A,
                                               const FieldList<Dim<3>, Dim<3>::Vector>& B,
                                               const ConnectivityMap<Dim<3> >& connectivityMap,
                                               const TableKernel< Dim<3> >& kernel);

    template 
    FieldList<Dim<3>, Dim<3>::ThirdRankTensor> 
    interpolateCRKSPH<Dim<3>, Dim<3>::ThirdRankTensor>(const FieldList<Dim<3>, Dim<3>::ThirdRankTensor>& fieldList,
                                                     const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                                     const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                                     const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                                     const bool coupleNodeLists,
                                                     const FieldList<Dim<3>, Dim<3>::Scalar>& A,
                                                     const FieldList<Dim<3>, Dim<3>::Vector>& B,
                                                     const ConnectivityMap<Dim<3> >& connectivityMap,
                                                     const TableKernel< Dim<3> >& kernel);

  }
}
