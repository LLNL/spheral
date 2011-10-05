#include "ExtendBoundary.hh"

//------------------------------------------------------------------------------
// Provide instances for each of the templated FieldList Boundary methods.
//------------------------------------------------------------------------------
void
applyScalarFieldListGhostBoundary1d(Boundary<Spheral::Dim<1> >* self,
                                    FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Scalar>& fieldList) {
  self->applyFieldListGhostBoundary(fieldList);
}

void
applyVectorFieldListGhostBoundary1d(Boundary<Spheral::Dim<1> >* self,
                                    FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Vector>& fieldList) {
  self->applyFieldListGhostBoundary(fieldList);
}

void
applyTensorFieldListGhostBoundary1d(Boundary<Spheral::Dim<1> >* self,
                                    FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Tensor>& fieldList) {
  self->applyFieldListGhostBoundary(fieldList);
}

void
applySymTensorFieldListGhostBoundary1d(Boundary<Spheral::Dim<1> >* self,
                                    FieldList<Spheral::Dim<1>, Spheral::Dim<1>::SymTensor>& fieldList) {
  self->applyFieldListGhostBoundary(fieldList);
}

//------------------------------------------------------------------------------
void
enforceScalarFieldListBoundary1d(Boundary<Spheral::Dim<1> >* self,
                                      FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Scalar>& fieldList) {
  self->enforceFieldListBoundary(fieldList);
}

void
enforceVectorFieldListBoundary1d(Boundary<Spheral::Dim<1> >* self,
                                      FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Vector>& fieldList) {
  self->enforceFieldListBoundary(fieldList);
}

void
enforceTensorFieldListBoundary1d(Boundary<Spheral::Dim<1> >* self,
                                      FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Tensor>& fieldList) {
  self->enforceFieldListBoundary(fieldList);
}

void
enforceSymTensorFieldListBoundary1d(Boundary<Spheral::Dim<1> >* self,
                                      FieldList<Spheral::Dim<1>, Spheral::Dim<1>::SymTensor>& fieldList) {
  self->enforceFieldListBoundary(fieldList);
}

//------------------------------------------------------------------------------
void
applyScalarFieldListGhostBoundary2d(Boundary<Spheral::Dim<2> >* self,
                                    FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Scalar>& fieldList) {
  self->applyFieldListGhostBoundary(fieldList);
}

void
applyVectorFieldListGhostBoundary2d(Boundary<Spheral::Dim<2> >* self,
                                    FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Vector>& fieldList) {
  self->applyFieldListGhostBoundary(fieldList);
}

void
applyTensorFieldListGhostBoundary2d(Boundary<Spheral::Dim<2> >* self,
                                    FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Tensor>& fieldList) {
  self->applyFieldListGhostBoundary(fieldList);
}

void
applySymTensorFieldListGhostBoundary2d(Boundary<Spheral::Dim<2> >* self,
                                    FieldList<Spheral::Dim<2>, Spheral::Dim<2>::SymTensor>& fieldList) {
  self->applyFieldListGhostBoundary(fieldList);
}

//------------------------------------------------------------------------------
void
enforceScalarFieldListBoundary2d(Boundary<Spheral::Dim<2> >* self,
                                      FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Scalar>& fieldList) {
  self->enforceFieldListBoundary(fieldList);
}

void
enforceVectorFieldListBoundary2d(Boundary<Spheral::Dim<2> >* self,
                                      FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Vector>& fieldList) {
  self->enforceFieldListBoundary(fieldList);
}

void
enforceTensorFieldListBoundary2d(Boundary<Spheral::Dim<2> >* self,
                                      FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Tensor>& fieldList) {
  self->enforceFieldListBoundary(fieldList);
}

void
enforceSymTensorFieldListBoundary2d(Boundary<Spheral::Dim<2> >* self,
                                      FieldList<Spheral::Dim<2>, Spheral::Dim<2>::SymTensor>& fieldList) {
  self->enforceFieldListBoundary(fieldList);
}

//------------------------------------------------------------------------------
void
applyScalarFieldListGhostBoundary3d(Boundary<Spheral::Dim<3> >* self,
                                    FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Scalar>& fieldList) {
  self->applyFieldListGhostBoundary(fieldList);
}

void
applyVectorFieldListGhostBoundary3d(Boundary<Spheral::Dim<3> >* self,
                                    FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Vector>& fieldList) {
  self->applyFieldListGhostBoundary(fieldList);
}

void
applyTensorFieldListGhostBoundary3d(Boundary<Spheral::Dim<3> >* self,
                                    FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Tensor>& fieldList) {
  self->applyFieldListGhostBoundary(fieldList);
}

void
applySymTensorFieldListGhostBoundary3d(Boundary<Spheral::Dim<3> >* self,
                                    FieldList<Spheral::Dim<3>, Spheral::Dim<3>::SymTensor>& fieldList) {
  self->applyFieldListGhostBoundary(fieldList);
}

//------------------------------------------------------------------------------
void
enforceScalarFieldListBoundary3d(Boundary<Spheral::Dim<3> >* self,
                                      FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Scalar>& fieldList) {
  self->enforceFieldListBoundary(fieldList);
}

void
enforceVectorFieldListBoundary3d(Boundary<Spheral::Dim<3> >* self,
                                      FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Vector>& fieldList) {
  self->enforceFieldListBoundary(fieldList);
}

void
enforceTensorFieldListBoundary3d(Boundary<Spheral::Dim<3> >* self,
                                      FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Tensor>& fieldList) {
  self->enforceFieldListBoundary(fieldList);
}

void
enforceSymTensorFieldListBoundary3d(Boundary<Spheral::Dim<3> >* self,
                                      FieldList<Spheral::Dim<3>, Spheral::Dim<3>::SymTensor>& fieldList) {
  self->enforceFieldListBoundary(fieldList);
}
