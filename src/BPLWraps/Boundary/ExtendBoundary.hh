#include "Boundary/Boundary.hh"
#include "Geometry/Dimension.hh"

#ifndef __GCCXML__
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "NodeList/NodeList.hh"
#include "DataBase/DataBase.hh"
#endif

using namespace Spheral::BoundarySpace;
using namespace Spheral::FieldSpace;

//------------------------------------------------------------------------------
// Provide instances for each of the templated FieldList Boundary methods.
//------------------------------------------------------------------------------
void
applyScalarFieldListGhostBoundary1d(Boundary<Spheral::Dim<1> >* self,
                                    FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Scalar>& fieldList);

void
applyVectorFieldListGhostBoundary1d(Boundary<Spheral::Dim<1> >* self,
                                    FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Vector>& fieldList);

void
applyTensorFieldListGhostBoundary1d(Boundary<Spheral::Dim<1> >* self,
                                    FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Tensor>& fieldList);

void
applySymTensorFieldListGhostBoundary1d(Boundary<Spheral::Dim<1> >* self,
                                       FieldList<Spheral::Dim<1>, Spheral::Dim<1>::SymTensor>& fieldList);

//------------------------------------------------------------------------------
void
enforceScalarFieldListBoundary1d(Boundary<Spheral::Dim<1> >* self,
                                 FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Scalar>& fieldList);

void
enforceVectorFieldListBoundary1d(Boundary<Spheral::Dim<1> >* self,
                                 FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Vector>& fieldList);

void
enforceTensorFieldListBoundary1d(Boundary<Spheral::Dim<1> >* self,
                                 FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Tensor>& fieldList);

void
enforceSymTensorFieldListBoundary1d(Boundary<Spheral::Dim<1> >* self,
                                    FieldList<Spheral::Dim<1>, Spheral::Dim<1>::SymTensor>& fieldList);

//------------------------------------------------------------------------------
void
applyScalarFieldListGhostBoundary2d(Boundary<Spheral::Dim<2> >* self,
                                    FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Scalar>& fieldList);

void
applyVectorFieldListGhostBoundary2d(Boundary<Spheral::Dim<2> >* self,
                                    FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Vector>& fieldList);

void
applyTensorFieldListGhostBoundary2d(Boundary<Spheral::Dim<2> >* self,
                                    FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Tensor>& fieldList);

void
applySymTensorFieldListGhostBoundary2d(Boundary<Spheral::Dim<2> >* self,
                                       FieldList<Spheral::Dim<2>, Spheral::Dim<2>::SymTensor>& fieldList);

//------------------------------------------------------------------------------
void
enforceScalarFieldListBoundary2d(Boundary<Spheral::Dim<2> >* self,
                                 FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Scalar>& fieldList);

void
enforceVectorFieldListBoundary2d(Boundary<Spheral::Dim<2> >* self,
                                 FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Vector>& fieldList);

void
enforceTensorFieldListBoundary2d(Boundary<Spheral::Dim<2> >* self,
                                 FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Tensor>& fieldList);

void
enforceSymTensorFieldListBoundary2d(Boundary<Spheral::Dim<2> >* self,
                                    FieldList<Spheral::Dim<2>, Spheral::Dim<2>::SymTensor>& fieldList);

//------------------------------------------------------------------------------
void
applyScalarFieldListGhostBoundary3d(Boundary<Spheral::Dim<3> >* self,
                                    FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Scalar>& fieldList);

void
applyVectorFieldListGhostBoundary3d(Boundary<Spheral::Dim<3> >* self,
                                    FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Vector>& fieldList);

void
applyTensorFieldListGhostBoundary3d(Boundary<Spheral::Dim<3> >* self,
                                    FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Tensor>& fieldList);

void
applySymTensorFieldListGhostBoundary3d(Boundary<Spheral::Dim<3> >* self,
                                       FieldList<Spheral::Dim<3>, Spheral::Dim<3>::SymTensor>& fieldList);

//------------------------------------------------------------------------------
void
enforceScalarFieldListBoundary3d(Boundary<Spheral::Dim<3> >* self,
                                 FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Scalar>& fieldList);

void
enforceVectorFieldListBoundary3d(Boundary<Spheral::Dim<3> >* self,
                                 FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Vector>& fieldList);

void
enforceTensorFieldListBoundary3d(Boundary<Spheral::Dim<3> >* self,
                                 FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Tensor>& fieldList);

void
enforceSymTensorFieldListBoundary3d(Boundary<Spheral::Dim<3> >* self,
                                    FieldList<Spheral::Dim<3>, Spheral::Dim<3>::SymTensor>& fieldList);
