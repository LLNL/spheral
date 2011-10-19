#include "ExtendDataBase.hh"

//------------------------------------------------------------------------------
// Provide a list interface for the NodeLists in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
boost::python::list
dataBaseNodeLists(Spheral::DataBaseSpace::DataBase<Dimension>* self) {
  return Spheral::iteratorsAsListByRef(self->nodeListBegin(), 
                                       self->nodeListEnd());
}

boost::python::list
dataBaseNodeLists1d(Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self) {
  return dataBaseNodeLists(self);
}

boost::python::list
dataBaseNodeLists2d(Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self) {
  return dataBaseNodeLists(self);
}

boost::python::list
dataBaseNodeLists3d(Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self) {
  return dataBaseNodeLists(self);
}

//------------------------------------------------------------------------------
// Provide a list interface for the FluidNodeLists in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
boost::python::list
dataBaseFluidNodeLists(Spheral::DataBaseSpace::DataBase<Dimension>* self) {
  return Spheral::iteratorsAsListByRef(self->fluidNodeListBegin(), 
                                       self->fluidNodeListEnd());
}

boost::python::list
dataBaseFluidNodeLists1d(Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self) {
  return dataBaseFluidNodeLists(self);
}

boost::python::list
dataBaseFluidNodeLists2d(Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self) {
  return dataBaseFluidNodeLists(self);
}

boost::python::list
dataBaseFluidNodeLists3d(Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self) {
  return dataBaseFluidNodeLists(self);
}

// //------------------------------------------------------------------------------
// // Provide a list interface for the MashNodeLists in the DataBase.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// boost::python::list
// dataBaseMashNodeLists(Spheral::DataBaseSpace::DataBase<Dimension>* self) {
//   return Spheral::iteratorsAsListByRef(self->mashNodeListBegin(), 
//                                        self->mashNodeListEnd());
// }

// inline
// boost::python::list
// dataBaseMashNodeLists1d(Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self) {
//   return dataBaseMashNodeLists(self);
// }

// inline
// boost::python::list
// dataBaseMashNodeLists2d(Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self) {
//   return dataBaseMashNodeLists(self);
// }

// inline
// boost::python::list
// dataBaseMashNodeLists3d(Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self) {
//   return dataBaseMashNodeLists(self);
// }

//------------------------------------------------------------------------------
// Expose the template methods for creating new global FieldLists.
//------------------------------------------------------------------------------
// 1-D
FieldList<Spheral::Dim<1>, int>
newGlobalIntFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                        const int val) {
  return self->newGlobalFieldList(val);
}

FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Scalar>
newGlobalScalarFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                           const Spheral::Dim<1>::Scalar val) {
  return self->newGlobalFieldList(val);
}

FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Vector>
newGlobalVectorFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                           const Spheral::Dim<1>::Vector val) {
  return self->newGlobalFieldList(val);
}

FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Tensor>
newGlobalTensorFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                           const Spheral::Dim<1>::Tensor val) {
  return self->newGlobalFieldList(val);
}

FieldList<Spheral::Dim<1>, Spheral::Dim<1>::SymTensor>
newGlobalSymTensorFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                              const Spheral::Dim<1>::SymTensor val) {
  return self->newGlobalFieldList(val);
}

// 2-D
FieldList<Spheral::Dim<2>, int>
newGlobalIntFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                        const int val) {
  return self->newGlobalFieldList(val);
}

FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Scalar>
newGlobalScalarFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                           const Spheral::Dim<2>::Scalar val) {
  return self->newGlobalFieldList(val);
}

FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Vector>
newGlobalVectorFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                           const Spheral::Dim<2>::Vector val) {
  return self->newGlobalFieldList(val);
}

FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Tensor>
newGlobalTensorFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                           const Spheral::Dim<2>::Tensor val) {
  return self->newGlobalFieldList(val);
}

FieldList<Spheral::Dim<2>, Spheral::Dim<2>::SymTensor>
newGlobalSymTensorFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                              const Spheral::Dim<2>::SymTensor val) {
  return self->newGlobalFieldList(val);
}

// 3-D
FieldList<Spheral::Dim<3>, int>
newGlobalIntFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                        const int val) {
  return self->newGlobalFieldList(val);
}

FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Scalar>
newGlobalScalarFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                           const Spheral::Dim<3>::Scalar val) {
  return self->newGlobalFieldList(val);
}

FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Vector>
newGlobalVectorFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                           const Spheral::Dim<3>::Vector val) {
  return self->newGlobalFieldList(val);
}

FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Tensor>
newGlobalTensorFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                           const Spheral::Dim<3>::Tensor val) {
  return self->newGlobalFieldList(val);
}

FieldList<Spheral::Dim<3>, Spheral::Dim<3>::SymTensor>
newGlobalSymTensorFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                              const Spheral::Dim<3>::SymTensor val) {
  return self->newGlobalFieldList(val);
}

//------------------------------------------------------------------------------
// Expose the template methods for creating new fluid FieldLists.
//------------------------------------------------------------------------------
// 1-D
FieldList<Spheral::Dim<1>, int>
newFluidIntFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                       const int val) {
  return self->newGlobalFieldList(val);
}

FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Scalar>
newFluidScalarFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                          const Spheral::Dim<1>::Scalar val) {
  return self->newFluidFieldList(val);
}

FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Vector>
newFluidVectorFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                          const Spheral::Dim<1>::Vector val) {
  return self->newFluidFieldList(val);
}

FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Tensor>
newFluidTensorFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                          const Spheral::Dim<1>::Tensor val) {
  return self->newFluidFieldList(val);
}

FieldList<Spheral::Dim<1>, Spheral::Dim<1>::SymTensor>
newFluidSymTensorFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                             const Spheral::Dim<1>::SymTensor val) {
  return self->newFluidFieldList(val);
}

// 2-D
FieldList<Spheral::Dim<2>, int>
newFluidIntFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                       const int val) {
  return self->newGlobalFieldList(val);
}

FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Scalar>
newFluidScalarFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                          const Spheral::Dim<2>::Scalar val) {
  return self->newFluidFieldList(val);
}

FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Vector>
newFluidVectorFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                          const Spheral::Dim<2>::Vector val) {
  return self->newFluidFieldList(val);
}

FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Tensor>
newFluidTensorFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                          const Spheral::Dim<2>::Tensor val) {
  return self->newFluidFieldList(val);
}

FieldList<Spheral::Dim<2>, Spheral::Dim<2>::SymTensor>
newFluidSymTensorFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                             const Spheral::Dim<2>::SymTensor val) {
  return self->newFluidFieldList(val);
}

// 3-D
FieldList<Spheral::Dim<3>, int>
newFluidIntFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                       const int val) {
  return self->newGlobalFieldList(val);
}

FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Scalar>
newFluidScalarFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                          const Spheral::Dim<3>::Scalar val) {
  return self->newFluidFieldList(val);
}

FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Vector>
newFluidVectorFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                          const Spheral::Dim<3>::Vector val) {
  return self->newFluidFieldList(val);
}

FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Tensor>
newFluidTensorFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                          const Spheral::Dim<3>::Tensor val) {
  return self->newFluidFieldList(val);
}

FieldList<Spheral::Dim<3>, Spheral::Dim<3>::SymTensor>
newFluidSymTensorFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                             const Spheral::Dim<3>::SymTensor val) {
  return self->newFluidFieldList(val);
}

// //------------------------------------------------------------------------------
// // Expose the template methods for creating new mash FieldLists.
// //------------------------------------------------------------------------------
// // 1-D
// inline
// FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Scalar>
// newMashScalarFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self) {
//   return self->newMashFieldList(Spheral::Dim<1>::Scalar());
// }

// inline
// FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Vector>
// newMashVectorFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self) {
//   return self->newMashFieldList(Spheral::Dim<1>::Vector());
// }

// inline
// FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Tensor>
// newMashTensorFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self) {
//   return self->newMashFieldList(Spheral::Dim<1>::Tensor());
// }

// inline
// FieldList<Spheral::Dim<1>, Spheral::Dim<1>::SymTensor>
// newMashSymTensorFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self) {
//   return self->newMashFieldList(Spheral::Dim<1>::SymTensor());
// }

// // 2-D
// inline
// FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Scalar>
// newMashScalarFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self) {
//   return self->newMashFieldList(Spheral::Dim<2>::Scalar());
// }

// inline
// FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Vector>
// newMashVectorFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self) {
//   return self->newMashFieldList(Spheral::Dim<2>::Vector());
// }

// inline
// FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Tensor>
// newMashTensorFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self) {
//   return self->newMashFieldList(Spheral::Dim<2>::Tensor());
// }

// inline
// FieldList<Spheral::Dim<2>, Spheral::Dim<2>::SymTensor>
// newMashSymTensorFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self) {
//   return self->newMashFieldList(Spheral::Dim<2>::SymTensor());
// }

// // 3-D
// inline
// FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Scalar>
// newMashScalarFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self) {
//   return self->newMashFieldList(Spheral::Dim<3>::Scalar());
// }

// inline
// FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Vector>
// newMashVectorFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self) {
//   return self->newMashFieldList(Spheral::Dim<3>::Vector());
// }

// inline
// FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Tensor>
// newMashTensorFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self) {
//   return self->newMashFieldList(Spheral::Dim<3>::Tensor());
// }

// inline
// FieldList<Spheral::Dim<3>, Spheral::Dim<3>::SymTensor>
// newMashSymTensorFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self) {
//   return self->newMashFieldList(Spheral::Dim<3>::SymTensor());
// }

//------------------------------------------------------------------------------
// Expose the template methods for resizing global FieldLists.
//------------------------------------------------------------------------------
// 1-D
void
resizeGlobalScalarFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                              FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Scalar>& fieldList) {
  self->resizeGlobalFieldList(fieldList);
}

void
resizeGlobalVectorFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                              FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Vector>& fieldList) {
  self->resizeGlobalFieldList(fieldList);
}

void
resizeGlobalTensorFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                              FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Tensor>& fieldList) {
  self->resizeGlobalFieldList(fieldList);
}

void
resizeGlobalSymTensorFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                                 FieldList<Spheral::Dim<1>, Spheral::Dim<1>::SymTensor>& fieldList) {
  self->resizeGlobalFieldList(fieldList);
}

// 2-D
void
resizeGlobalScalarFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                              FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Scalar>& fieldList) {
  self->resizeGlobalFieldList(fieldList);
}

void
resizeGlobalVectorFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                              FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Vector>& fieldList) {
  self->resizeGlobalFieldList(fieldList);
}

void
resizeGlobalTensorFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                              FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Tensor>& fieldList) {
  self->resizeGlobalFieldList(fieldList);
}

void
resizeGlobalSymTensorFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                                 FieldList<Spheral::Dim<2>, Spheral::Dim<2>::SymTensor>& fieldList) {
  self->resizeGlobalFieldList(fieldList);
}

// 3-D
void
resizeGlobalScalarFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                              FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Scalar>& fieldList) {
  self->resizeGlobalFieldList(fieldList);
}

void
resizeGlobalVectorFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                              FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Vector>& fieldList) {
  self->resizeGlobalFieldList(fieldList);
}

void
resizeGlobalTensorFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                              FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Tensor>& fieldList) {
  self->resizeGlobalFieldList(fieldList);
}

void
resizeGlobalSymTensorFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                                 FieldList<Spheral::Dim<3>, Spheral::Dim<3>::SymTensor>& fieldList) {
  self->resizeGlobalFieldList(fieldList);
}

//------------------------------------------------------------------------------
// Expose the template methods for resizing fluid FieldLists.
//------------------------------------------------------------------------------
// 1-D
void
resizeFluidScalarFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                              FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Scalar>& fieldList) {
  self->resizeFluidFieldList(fieldList);
}

void
resizeFluidVectorFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                              FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Vector>& fieldList) {
  self->resizeFluidFieldList(fieldList);
}

void
resizeFluidTensorFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                              FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Tensor>& fieldList) {
  self->resizeFluidFieldList(fieldList);
}

void
resizeFluidSymTensorFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                                 FieldList<Spheral::Dim<1>, Spheral::Dim<1>::SymTensor>& fieldList) {
  self->resizeFluidFieldList(fieldList);
}

// 2-D
void
resizeFluidScalarFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                              FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Scalar>& fieldList) {
  self->resizeFluidFieldList(fieldList);
}

void
resizeFluidVectorFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                              FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Vector>& fieldList) {
  self->resizeFluidFieldList(fieldList);
}

void
resizeFluidTensorFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                              FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Tensor>& fieldList) {
  self->resizeFluidFieldList(fieldList);
}

void
resizeFluidSymTensorFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                                 FieldList<Spheral::Dim<2>, Spheral::Dim<2>::SymTensor>& fieldList) {
  self->resizeFluidFieldList(fieldList);
}

// 3-D
void
resizeFluidScalarFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                              FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Scalar>& fieldList) {
  self->resizeFluidFieldList(fieldList);
}

void
resizeFluidVectorFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                              FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Vector>& fieldList) {
  self->resizeFluidFieldList(fieldList);
}

void
resizeFluidTensorFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                              FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Tensor>& fieldList) {
  self->resizeFluidFieldList(fieldList);
}

void
resizeFluidSymTensorFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                                 FieldList<Spheral::Dim<3>, Spheral::Dim<3>::SymTensor>& fieldList) {
  self->resizeFluidFieldList(fieldList);
}

// //------------------------------------------------------------------------------
// // Expose the template methods for resizing mash FieldLists.
// //------------------------------------------------------------------------------
// // 1-D
// inline
// void
// resizeMashScalarFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
//                               FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Scalar>& fieldList) {
//   self->resizeMashFieldList(fieldList);
// }

// inline
// void
// resizeMashVectorFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
//                               FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Vector>& fieldList) {
//   self->resizeMashFieldList(fieldList);
// }

// inline
// void
// resizeMashTensorFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
//                               FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Tensor>& fieldList) {
//   self->resizeMashFieldList(fieldList);
// }

// inline
// void
// resizeMashSymTensorFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
//                                  FieldList<Spheral::Dim<1>, Spheral::Dim<1>::SymTensor>& fieldList) {
//   self->resizeMashFieldList(fieldList);
// }

// // 2-D
// inline
// void
// resizeMashScalarFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
//                               FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Scalar>& fieldList) {
//   self->resizeMashFieldList(fieldList);
// }

// inline
// void
// resizeMashVectorFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
//                               FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Vector>& fieldList) {
//   self->resizeMashFieldList(fieldList);
// }

// inline
// void
// resizeMashTensorFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
//                               FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Tensor>& fieldList) {
//   self->resizeMashFieldList(fieldList);
// }

// inline
// void
// resizeMashSymTensorFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
//                                  FieldList<Spheral::Dim<2>, Spheral::Dim<2>::SymTensor>& fieldList) {
//   self->resizeMashFieldList(fieldList);
// }

// // 3-D
// inline
// void
// resizeMashScalarFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
//                               FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Scalar>& fieldList) {
//   self->resizeMashFieldList(fieldList);
// }

// inline
// void
// resizeMashVectorFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
//                               FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Vector>& fieldList) {
//   self->resizeMashFieldList(fieldList);
// }

// inline
// void
// resizeMashTensorFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
//                               FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Tensor>& fieldList) {
//   self->resizeMashFieldList(fieldList);
// }

// inline
// void
// resizeMashSymTensorFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
//                                  FieldList<Spheral::Dim<3>, Spheral::Dim<3>::SymTensor>& fieldList) {
//   self->resizeMashFieldList(fieldList);
// }
