//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "gradientCRKSPH.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace CRKSPHSpace {

    using KernelSpace::TableKernel;
    using NeighborSpace::ConnectivityMap;

    // 1D
    template 
    FieldList<Dim<1>, MathTraits<Dim<1>, Dim<1>::Scalar>::GradientType> 
    gradientCRKSPH<Dim<1>, Dim<1>::Scalar>(const FieldList<Dim<1>, Dim<1>::Scalar>& fieldList,
					 const FieldList<Dim<1>, Dim<1>::Vector>& position,
					 const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
					 const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
					 const FieldList<Dim<1>, Dim<1>::Scalar>& A,
					 const FieldList<Dim<1>, Dim<1>::Vector>& B,
					 const FieldList<Dim<1>, Dim<1>::Vector>& C,
					 const FieldList<Dim<1>, Dim<1>::Tensor>& D,
					 const FieldList<Dim<1>, Dim<1>::Vector>& gradA,
					 const FieldList<Dim<1>, Dim<1>::Tensor>& gradB,
					 const ConnectivityMap<Dim<1> >& connectivityMap,
					 const TableKernel< Dim<1> >& kernel);
    template 
    FieldList<Dim<1>, MathTraits<Dim<1>, Dim<1>::Vector>::GradientType> 
    gradientCRKSPH<Dim<1>, Dim<1>::Vector>(const FieldList<Dim<1>, Dim<1>::Vector>& fieldList,
					 const FieldList<Dim<1>, Dim<1>::Vector>& position,
					 const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
					 const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
					 const FieldList<Dim<1>, Dim<1>::Scalar>& A,
					 const FieldList<Dim<1>, Dim<1>::Vector>& B,
					 const FieldList<Dim<1>, Dim<1>::Vector>& C,
					 const FieldList<Dim<1>, Dim<1>::Tensor>& D,
					 const FieldList<Dim<1>, Dim<1>::Vector>& gradA,
					 const FieldList<Dim<1>, Dim<1>::Tensor>& gradB,
					 const ConnectivityMap<Dim<1> >& connectivityMap,
					 const TableKernel< Dim<1> >& kernel);

    // 2D
    template 
    FieldList<Dim<2>, MathTraits<Dim<2>, Dim<2>::Scalar>::GradientType> 
    gradientCRKSPH<Dim<2>, Dim<2>::Scalar>(const FieldList<Dim<2>, Dim<2>::Scalar>& fieldList,
					 const FieldList<Dim<2>, Dim<2>::Vector>& position,
					 const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
					 const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
					 const FieldList<Dim<2>, Dim<2>::Scalar>& A,
					 const FieldList<Dim<2>, Dim<2>::Vector>& B,
					 const FieldList<Dim<2>, Dim<2>::Vector>& C,
					 const FieldList<Dim<2>, Dim<2>::Tensor>& D,
					 const FieldList<Dim<2>, Dim<2>::Vector>& gradA,
					 const FieldList<Dim<2>, Dim<2>::Tensor>& gradB,
					 const ConnectivityMap<Dim<2> >& connectivityMap,
					 const TableKernel< Dim<2> >& kernel);
    template 
    FieldList<Dim<2>, MathTraits<Dim<2>, Dim<2>::Vector>::GradientType> 
    gradientCRKSPH<Dim<2>, Dim<2>::Vector>(const FieldList<Dim<2>, Dim<2>::Vector>& fieldList,
					 const FieldList<Dim<2>, Dim<2>::Vector>& position,
					 const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
					 const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
					 const FieldList<Dim<2>, Dim<2>::Scalar>& A,
					 const FieldList<Dim<2>, Dim<2>::Vector>& B,
					 const FieldList<Dim<2>, Dim<2>::Vector>& C,
					 const FieldList<Dim<2>, Dim<2>::Tensor>& D,
					 const FieldList<Dim<2>, Dim<2>::Vector>& gradA,
					 const FieldList<Dim<2>, Dim<2>::Tensor>& gradB,
					 const ConnectivityMap<Dim<2> >& connectivityMap,
					 const TableKernel< Dim<2> >& kernel);

    // 3D
    template 
    FieldList<Dim<3>, MathTraits<Dim<3>, Dim<3>::Scalar>::GradientType> 
    gradientCRKSPH<Dim<3>, Dim<3>::Scalar>(const FieldList<Dim<3>, Dim<3>::Scalar>& fieldList,
					 const FieldList<Dim<3>, Dim<3>::Vector>& position,
					 const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
					 const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
					 const FieldList<Dim<3>, Dim<3>::Scalar>& A,
					 const FieldList<Dim<3>, Dim<3>::Vector>& B,
					 const FieldList<Dim<3>, Dim<3>::Vector>& C,
					 const FieldList<Dim<3>, Dim<3>::Tensor>& D,
					 const FieldList<Dim<3>, Dim<3>::Vector>& gradA,
					 const FieldList<Dim<3>, Dim<3>::Tensor>& gradB,
					 const ConnectivityMap<Dim<3> >& connectivityMap,
					 const TableKernel< Dim<3> >& kernel);
    template 
    FieldList<Dim<3>, MathTraits<Dim<3>, Dim<3>::Vector>::GradientType> 
    gradientCRKSPH<Dim<3>, Dim<3>::Vector>(const FieldList<Dim<3>, Dim<3>::Vector>& fieldList,
					 const FieldList<Dim<3>, Dim<3>::Vector>& position,
					 const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
					 const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
					 const FieldList<Dim<3>, Dim<3>::Scalar>& A,
					 const FieldList<Dim<3>, Dim<3>::Vector>& B,
					 const FieldList<Dim<3>, Dim<3>::Vector>& C,
					 const FieldList<Dim<3>, Dim<3>::Tensor>& D,
					 const FieldList<Dim<3>, Dim<3>::Vector>& gradA,
					 const FieldList<Dim<3>, Dim<3>::Tensor>& gradB,
					 const ConnectivityMap<Dim<3> >& connectivityMap,
					 const TableKernel< Dim<3> >& kernel);

  }
}
