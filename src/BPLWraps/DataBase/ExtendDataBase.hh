#include "DataBase/DataBase.hh"
#include "Geometry/Dimension.hh"

#ifndef __GCCXML__
#include "Kernel/TableKernel.hh"
#include "Field/FieldList.hh"
#endif

#include "BPLWraps/iteratorsAsList.hh"

using Spheral::FieldSpace::FieldList;

//------------------------------------------------------------------------------
// Provide a list interface for the NodeLists in the DataBase.
//------------------------------------------------------------------------------
boost::python::list
dataBaseNodeLists1d(Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self);

boost::python::list
dataBaseNodeLists2d(Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self);

boost::python::list
dataBaseNodeLists3d(Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self);

//------------------------------------------------------------------------------
// Provide a list interface for the FluidNodeLists in the DataBase.
//------------------------------------------------------------------------------
boost::python::list
dataBaseFluidNodeLists1d(Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self);

boost::python::list
dataBaseFluidNodeLists2d(Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self);

boost::python::list
dataBaseFluidNodeLists3d(Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self);

//------------------------------------------------------------------------------
// Expose the template methods for creating new global FieldLists.
//------------------------------------------------------------------------------
// 1-D
FieldList<Spheral::Dim<1>, int>
newGlobalIntFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                        const int val = 0);

FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Scalar>
newGlobalScalarFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                           const Spheral::Dim<1>::Scalar val = 0.0);

FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Vector>
newGlobalVectorFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                           const Spheral::Dim<1>::Vector val = Spheral::Dim<1>::Vector());

FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Tensor>
newGlobalTensorFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                           const Spheral::Dim<1>::Tensor val = Spheral::Dim<1>::Tensor());

FieldList<Spheral::Dim<1>, Spheral::Dim<1>::SymTensor>
newGlobalSymTensorFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                              const Spheral::Dim<1>::SymTensor val = Spheral::Dim<1>::SymTensor());

// 2-D
FieldList<Spheral::Dim<2>, int>
newGlobalIntFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                        const int val = 0);

FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Scalar>
newGlobalScalarFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                           const Spheral::Dim<2>::Scalar val = Spheral::Dim<2>::Scalar());

FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Vector>
newGlobalVectorFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                           const Spheral::Dim<2>::Vector val = Spheral::Dim<2>::Vector());

FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Tensor>
newGlobalTensorFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                           const Spheral::Dim<2>::Tensor val = Spheral::Dim<2>::Tensor());

FieldList<Spheral::Dim<2>, Spheral::Dim<2>::SymTensor>
newGlobalSymTensorFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                              const Spheral::Dim<2>::SymTensor val = Spheral::Dim<2>::SymTensor());

// 3-D
FieldList<Spheral::Dim<3>, int>
newGlobalIntFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                        const int val = 0);

FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Scalar>
newGlobalScalarFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                           const Spheral::Dim<3>::Scalar val = Spheral::Dim<3>::Scalar());

FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Vector>
newGlobalVectorFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                           const Spheral::Dim<3>::Vector val = Spheral::Dim<3>::Vector());

FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Tensor>
newGlobalTensorFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                           const Spheral::Dim<3>::Tensor val = Spheral::Dim<3>::Tensor());

FieldList<Spheral::Dim<3>, Spheral::Dim<3>::SymTensor>
newGlobalSymTensorFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                              const Spheral::Dim<3>::SymTensor val = Spheral::Dim<3>::SymTensor());

//------------------------------------------------------------------------------
// Expose the template methods for creating new fluid FieldLists.
//------------------------------------------------------------------------------
// 1-D
FieldList<Spheral::Dim<1>, int>
newFluidIntFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                       const int val = 0);

FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Scalar>
newFluidScalarFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                          const Spheral::Dim<1>::Scalar val = Spheral::Dim<1>::Scalar());

FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Vector>
newFluidVectorFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                          const Spheral::Dim<1>::Vector val = Spheral::Dim<1>::Vector());

FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Tensor>
newFluidTensorFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                          const Spheral::Dim<1>::Tensor val = Spheral::Dim<1>::Tensor());

FieldList<Spheral::Dim<1>, Spheral::Dim<1>::SymTensor>
newFluidSymTensorFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                             const Spheral::Dim<1>::SymTensor val = Spheral::Dim<1>::SymTensor());

// 2-D
FieldList<Spheral::Dim<2>, int>
newFluidIntFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                       const int val = 0);

FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Scalar>
newFluidScalarFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                          const Spheral::Dim<2>::Scalar val = Spheral::Dim<2>::Scalar());

FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Vector>
newFluidVectorFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                          const Spheral::Dim<2>::Vector val = Spheral::Dim<2>::Vector());

FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Tensor>
newFluidTensorFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                          const Spheral::Dim<2>::Tensor val = Spheral::Dim<2>::Tensor());

FieldList<Spheral::Dim<2>, Spheral::Dim<2>::SymTensor>
newFluidSymTensorFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                             const Spheral::Dim<2>::SymTensor val = Spheral::Dim<2>::SymTensor());

// 3-D
FieldList<Spheral::Dim<3>, int>
newFluidIntFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                       const int val = 0);

FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Scalar>
newFluidScalarFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                          const Spheral::Dim<3>::Scalar val = Spheral::Dim<3>::Scalar());

FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Vector>
newFluidVectorFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                          const Spheral::Dim<3>::Vector val = Spheral::Dim<3>::Vector());

FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Tensor>
newFluidTensorFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                          const Spheral::Dim<3>::Tensor val = Spheral::Dim<3>::Tensor());

FieldList<Spheral::Dim<3>, Spheral::Dim<3>::SymTensor>
newFluidSymTensorFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                             const Spheral::Dim<3>::SymTensor val = Spheral::Dim<3>::SymTensor());

//------------------------------------------------------------------------------
// Expose the template methods for resizing global FieldLists.
//------------------------------------------------------------------------------
// 1-D
void
resizeGlobalScalarFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                              FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Scalar>& fieldList);

void
resizeGlobalVectorFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                              FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Vector>& fieldList);

void
resizeGlobalTensorFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                              FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Tensor>& fieldList);

void
resizeGlobalSymTensorFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                                 FieldList<Spheral::Dim<1>, Spheral::Dim<1>::SymTensor>& fieldList);

// 2-D
void
resizeGlobalScalarFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                              FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Scalar>& fieldList);

void
resizeGlobalVectorFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                              FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Vector>& fieldList);

void
resizeGlobalTensorFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                              FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Tensor>& fieldList);

void
resizeGlobalSymTensorFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                                 FieldList<Spheral::Dim<2>, Spheral::Dim<2>::SymTensor>& fieldList);

// 3-D
void
resizeGlobalScalarFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                              FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Scalar>& fieldList);

void
resizeGlobalVectorFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                              FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Vector>& fieldList);

void
resizeGlobalTensorFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                              FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Tensor>& fieldList);

void
resizeGlobalSymTensorFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                                 FieldList<Spheral::Dim<3>, Spheral::Dim<3>::SymTensor>& fieldList);

//------------------------------------------------------------------------------
// Expose the template methods for resizing fluid FieldLists.
//------------------------------------------------------------------------------
// 1-D
void
resizeFluidScalarFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                             FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Scalar>& fieldList);

void
resizeFluidVectorFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                             FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Vector>& fieldList);

void
resizeFluidTensorFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                             FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Tensor>& fieldList);

void
resizeFluidSymTensorFieldList1d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<1> >* self,
                                FieldList<Spheral::Dim<1>, Spheral::Dim<1>::SymTensor>& fieldList);

// 2-D
void
resizeFluidScalarFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                             FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Scalar>& fieldList);

void
resizeFluidVectorFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                             FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Vector>& fieldList);

void
resizeFluidTensorFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                             FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Tensor>& fieldList);

void
resizeFluidSymTensorFieldList2d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<2> >* self,
                                FieldList<Spheral::Dim<2>, Spheral::Dim<2>::SymTensor>& fieldList);

// 3-D
void
resizeFluidScalarFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                             FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Scalar>& fieldList);

void
resizeFluidVectorFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                             FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Vector>& fieldList);

void
resizeFluidTensorFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                             FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Tensor>& fieldList);

void
resizeFluidSymTensorFieldList3d(const Spheral::DataBaseSpace::DataBase<Spheral::Dim<3> >* self,
                                FieldList<Spheral::Dim<3>, Spheral::Dim<3>::SymTensor>& fieldList);
