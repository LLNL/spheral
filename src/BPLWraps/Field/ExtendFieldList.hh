#ifndef __ExtendFieldList_hh__
#define __ExtendFieldList_hh__

#include "Geometry/Dimension.hh"
#include "Field/FieldList.hh"
// #include "NodeList/NodeList.hh"
#include "Kernel/TableKernel.hh"

#include "BPLWraps/getSetItem.hh"
#include "BPLWraps/iteratorsAsList.hh"

using namespace Spheral::FieldSpace;

typedef unsigned long long ULL;

typedef Spheral::FieldSpace::Field<Spheral::Dim<1>, int> IntField1d;
typedef Spheral::FieldSpace::Field<Spheral::Dim<1>, ULL> ULLField1d;
typedef Spheral::FieldSpace::Field<Spheral::Dim<1>, Spheral::Dim<1>::Scalar> ScalarField1d;
typedef Spheral::FieldSpace::Field<Spheral::Dim<1>, Spheral::Dim<1>::Vector> VectorField1d;
typedef Spheral::FieldSpace::Field<Spheral::Dim<1>, Spheral::Dim<1>::Vector3d> Vector3dField1d;
typedef Spheral::FieldSpace::Field<Spheral::Dim<1>, Spheral::Dim<1>::Tensor> TensorField1d;
typedef Spheral::FieldSpace::Field<Spheral::Dim<1>, Spheral::Dim<1>::SymTensor> SymTensorField1d;
typedef Spheral::FieldSpace::Field<Spheral::Dim<1>, Spheral::Dim<1>::ThirdRankTensor> ThirdRankTensorField1d;

typedef Spheral::FieldSpace::Field<Spheral::Dim<2>, int> IntField2d;
typedef Spheral::FieldSpace::Field<Spheral::Dim<2>, ULL> ULLField2d;
typedef Spheral::FieldSpace::Field<Spheral::Dim<2>, Spheral::Dim<2>::Scalar> ScalarField2d;
typedef Spheral::FieldSpace::Field<Spheral::Dim<2>, Spheral::Dim<2>::Vector> VectorField2d;
typedef Spheral::FieldSpace::Field<Spheral::Dim<2>, Spheral::Dim<2>::Vector3d> Vector3dField2d;
typedef Spheral::FieldSpace::Field<Spheral::Dim<2>, Spheral::Dim<2>::Tensor> TensorField2d;
typedef Spheral::FieldSpace::Field<Spheral::Dim<2>, Spheral::Dim<2>::SymTensor> SymTensorField2d;
typedef Spheral::FieldSpace::Field<Spheral::Dim<2>, Spheral::Dim<2>::ThirdRankTensor> ThirdRankTensorField2d;

typedef Spheral::FieldSpace::Field<Spheral::Dim<3>, int> IntField3d;
typedef Spheral::FieldSpace::Field<Spheral::Dim<3>, ULL> ULLField3d;
typedef Spheral::FieldSpace::Field<Spheral::Dim<3>, Spheral::Dim<3>::Scalar> ScalarField3d;
typedef Spheral::FieldSpace::Field<Spheral::Dim<3>, Spheral::Dim<3>::Vector> VectorField3d;
typedef Spheral::FieldSpace::Field<Spheral::Dim<3>, Spheral::Dim<3>::Vector3d> Vector3dField3d;
typedef Spheral::FieldSpace::Field<Spheral::Dim<3>, Spheral::Dim<3>::Tensor> TensorField3d;
typedef Spheral::FieldSpace::Field<Spheral::Dim<3>, Spheral::Dim<3>::SymTensor> SymTensorField3d;
typedef Spheral::FieldSpace::Field<Spheral::Dim<3>, Spheral::Dim<3>::ThirdRankTensor> ThirdRankTensorField3d;

typedef Spheral::FieldSpace::FieldList<Spheral::Dim<1>, int> IntFieldList1d;
typedef Spheral::FieldSpace::FieldList<Spheral::Dim<1>, ULL> ULLFieldList1d;
typedef Spheral::FieldSpace::FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Scalar> ScalarFieldList1d;
typedef Spheral::FieldSpace::FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Vector> VectorFieldList1d;
typedef Spheral::FieldSpace::FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Vector3d> Vector3dFieldList1d;
typedef Spheral::FieldSpace::FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Tensor> TensorFieldList1d;
typedef Spheral::FieldSpace::FieldList<Spheral::Dim<1>, Spheral::Dim<1>::SymTensor> SymTensorFieldList1d;
typedef Spheral::FieldSpace::FieldList<Spheral::Dim<1>, Spheral::Dim<1>::ThirdRankTensor> ThirdRankTensorFieldList1d;

typedef Spheral::FieldSpace::FieldList<Spheral::Dim<2>, int> IntFieldList2d;
typedef Spheral::FieldSpace::FieldList<Spheral::Dim<2>, ULL> ULLFieldList2d;
typedef Spheral::FieldSpace::FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Scalar> ScalarFieldList2d;
typedef Spheral::FieldSpace::FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Vector> VectorFieldList2d;
typedef Spheral::FieldSpace::FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Vector3d> Vector3dFieldList2d;
typedef Spheral::FieldSpace::FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Tensor> TensorFieldList2d;
typedef Spheral::FieldSpace::FieldList<Spheral::Dim<2>, Spheral::Dim<2>::SymTensor> SymTensorFieldList2d;
typedef Spheral::FieldSpace::FieldList<Spheral::Dim<2>, Spheral::Dim<2>::ThirdRankTensor> ThirdRankTensorFieldList2d;

typedef Spheral::FieldSpace::FieldList<Spheral::Dim<3>, int> IntFieldList3d;
typedef Spheral::FieldSpace::FieldList<Spheral::Dim<3>, ULL> ULLFieldList3d;
typedef Spheral::FieldSpace::FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Scalar> ScalarFieldList3d;
typedef Spheral::FieldSpace::FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Vector> VectorFieldList3d;
typedef Spheral::FieldSpace::FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Vector3d> Vector3dFieldList3d;
typedef Spheral::FieldSpace::FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Tensor> TensorFieldList3d;
typedef Spheral::FieldSpace::FieldList<Spheral::Dim<3>, Spheral::Dim<3>::SymTensor> SymTensorFieldList3d;
typedef Spheral::FieldSpace::FieldList<Spheral::Dim<3>, Spheral::Dim<3>::ThirdRankTensor> ThirdRankTensorFieldList3d;

namespace Spheral {
namespace FieldSpace {

// Declare our methods.
IntField1d* getItemIntFieldList1d(IntFieldList1d& self, const unsigned int index);
ULLField1d* getItemULLFieldList1d(ULLFieldList1d& self, const unsigned int index);
ScalarField1d* getItemScalarFieldList1d(ScalarFieldList1d& self, const unsigned int index);
Field<Dim<1>, Dim<1>::Vector>* getItemVectorFieldList1d(FieldList<Dim<1>, Dim<1>::Vector>& self, const unsigned int index);
Field<Dim<1>, Dim<1>::Vector3d>* getItemVector3dFieldList1d(FieldList<Dim<1>, Dim<1>::Vector3d>& self, const unsigned int index);
Field<Dim<1>, Dim<1>::Tensor>* getItemTensorFieldList1d(FieldList<Dim<1>, Dim<1>::Tensor>& self, const unsigned int index);
Field<Dim<1>, Dim<1>::SymTensor>* getItemSymTensorFieldList1d(FieldList<Dim<1>, Dim<1>::SymTensor>& self, const unsigned int index);
Field<Dim<1>, Dim<1>::ThirdRankTensor>* getItemThirdRankTensorFieldList1d(FieldList<Dim<1>, Dim<1>::ThirdRankTensor>& self, const unsigned int index);
IntField2d* getItemIntFieldList2d(IntFieldList2d& self, const unsigned int index);
ULLField2d* getItemULLFieldList2d(ULLFieldList2d& self, const unsigned int index);
ScalarField2d* getItemScalarFieldList2d(ScalarFieldList2d& self, const unsigned int index);
Field<Dim<2>, Dim<2>::Vector>* getItemVectorFieldList2d(FieldList<Dim<2>, Dim<2>::Vector>& self, const unsigned int index);
Field<Dim<2>, Dim<2>::Vector3d>* getItemVector3dFieldList2d(FieldList<Dim<2>, Dim<2>::Vector3d>& self, const unsigned int index);
Field<Dim<2>, Dim<2>::Tensor>* getItemTensorFieldList2d(FieldList<Dim<2>, Dim<2>::Tensor>& self, const unsigned int index);
Field<Dim<2>, Dim<2>::SymTensor>* getItemSymTensorFieldList2d(FieldList<Dim<2>, Dim<2>::SymTensor>& self, const unsigned int index);
Field<Dim<2>, Dim<2>::ThirdRankTensor>* getItemThirdRankTensorFieldList2d(FieldList<Dim<2>, Dim<2>::ThirdRankTensor>& self, const unsigned int index);
IntField3d* getItemIntFieldList3d(IntFieldList3d& self, const unsigned int index);
ULLField3d* getItemULLFieldList3d(ULLFieldList3d& self, const unsigned int index);
ScalarField3d* getItemScalarFieldList3d(ScalarFieldList3d& self, const unsigned int index);
Field<Dim<3>, Dim<3>::Vector>* getItemVectorFieldList3d(FieldList<Dim<3>, Dim<3>::Vector>& self, const unsigned int index);
Field<Dim<3>, Dim<3>::Vector3d>* getItemVector3dFieldList3d(FieldList<Dim<3>, Dim<3>::Vector3d>& self, const unsigned int index);
Field<Dim<3>, Dim<3>::Tensor>* getItemTensorFieldList3d(FieldList<Dim<3>, Dim<3>::Tensor>& self, const unsigned int index);
Field<Dim<3>, Dim<3>::SymTensor>* getItemSymTensorFieldList3d(FieldList<Dim<3>, Dim<3>::SymTensor>& self, const unsigned int index);
Field<Dim<3>, Dim<3>::ThirdRankTensor>* getItemThirdRankTensorFieldList3d(FieldList<Dim<3>, Dim<3>::ThirdRankTensor>& self, const unsigned int index);
boost::python::list intFields1d(FieldList<Dim<1>, int>* self);
boost::python::list ULLFields1d(FieldList<Dim<1>, ULL>* self);
boost::python::list scalarFields1d(FieldList<Dim<1>, Dim<1>::Scalar>* self);
boost::python::list vectorFields1d(FieldList<Dim<1>, Dim<1>::Vector>* self);
boost::python::list vector3dFields1d(FieldList<Dim<1>, Dim<1>::Vector3d>* self);
boost::python::list tensorFields1d(FieldList<Dim<1>, Dim<1>::Tensor>* self);
boost::python::list symTensorFields1d(FieldList<Dim<1>, Dim<1>::SymTensor>* self);
boost::python::list thirdRankTensorFields1d(FieldList<Dim<1>, Dim<1>::ThirdRankTensor>* self);
boost::python::list intFields2d(FieldList<Dim<2>, int>* self);
boost::python::list ULLFields2d(FieldList<Dim<2>, ULL>* self);
boost::python::list scalarFields2d(FieldList<Dim<2>, Dim<2>::Scalar>* self);
boost::python::list vectorFields2d(FieldList<Dim<2>, Dim<2>::Vector>* self);
boost::python::list vector3dFields2d(FieldList<Dim<2>, Dim<2>::Vector3d>* self);
boost::python::list tensorFields2d(FieldList<Dim<2>, Dim<2>::Tensor>* self);
boost::python::list symTensorFields2d(FieldList<Dim<2>, Dim<2>::SymTensor>* self);
boost::python::list thirdRankTensorFields2d(FieldList<Dim<2>, Dim<2>::ThirdRankTensor>* self);
boost::python::list intFields3d(FieldList<Dim<3>, int>* self);
boost::python::list ULLFields3d(FieldList<Dim<3>, ULL>* self);
boost::python::list scalarFields3d(FieldList<Dim<3>, Dim<3>::Scalar>* self);
boost::python::list vectorFields3d(FieldList<Dim<3>, Dim<3>::Vector>* self);
boost::python::list vector3dFields3d(FieldList<Dim<3>, Dim<3>::Vector3d>* self);
boost::python::list tensorFields3d(FieldList<Dim<3>, Dim<3>::Tensor>* self);
boost::python::list symTensorFields3d(FieldList<Dim<3>, Dim<3>::SymTensor>* self);
boost::python::list thirdRankTensorFields3d(FieldList<Dim<3>, Dim<3>::ThirdRankTensor>* self);
Field<Dim<1>, int>& extractFieldForNodeListInt1d(const FieldList<Dim<1>, int>* self,
                                                 const NodeSpace::NodeList<Dim<1> >& nodeList);
Field<Dim<1>, ULL>& extractFieldForNodeListULL1d(const FieldList<Dim<1>, ULL>* self,
                                                                const NodeSpace::NodeList<Dim<1> >& nodeList);
Field<Dim<1>, Dim<1>::Scalar>& extractFieldForNodeListScalar1d(const FieldList<Dim<1>, Dim<1>::Scalar>* self,
                                                                           const NodeSpace::NodeList<Dim<1> >& nodeList);
Field<Dim<1>, Dim<1>::Vector>& extractFieldForNodeListVector1d(const FieldList<Dim<1>, Dim<1>::Vector>* self,
                                                                           const NodeSpace::NodeList<Dim<1> >& nodeList);
Field<Dim<1>, Dim<1>::Vector3d>& extractFieldForNodeListVector3d1d(const FieldList<Dim<1>, Dim<1>::Vector3d>* self,
                                                                               const NodeSpace::NodeList<Dim<1> >& nodeList);
Field<Dim<1>, Dim<1>::Tensor>& extractFieldForNodeListTensor1d(const FieldList<Dim<1>, Dim<1>::Tensor>* self,
                                                                           const NodeSpace::NodeList<Dim<1> >& nodeList);
Field<Dim<1>, Dim<1>::SymTensor>& extractFieldForNodeListSymTensor1d(const FieldList<Dim<1>, Dim<1>::SymTensor>* self,
                                                                                 const NodeSpace::NodeList<Dim<1> >& nodeList);
Field<Dim<1>, Dim<1>::ThirdRankTensor>& extractFieldForNodeListThirdRankTensor1d(const FieldList<Dim<1>, Dim<1>::ThirdRankTensor>* self,
                                                                                 const NodeSpace::NodeList<Dim<1> >& nodeList);
Field<Dim<2>, int>& extractFieldForNodeListInt2d(const FieldList<Dim<2>, int>* self,
                                                 const NodeSpace::NodeList<Dim<2> >& nodeList);
Field<Dim<2>, ULL>& extractFieldForNodeListULL2d(const FieldList<Dim<2>, ULL>* self,
                                                                const NodeSpace::NodeList<Dim<2> >& nodeList);
Field<Dim<2>, Dim<2>::Scalar>& extractFieldForNodeListScalar2d(const FieldList<Dim<2>, Dim<2>::Scalar>* self,
                                                                           const NodeSpace::NodeList<Dim<2> >& nodeList);
Field<Dim<2>, Dim<2>::Vector>& extractFieldForNodeListVector2d(const FieldList<Dim<2>, Dim<2>::Vector>* self,
                                                                           const NodeSpace::NodeList<Dim<2> >& nodeList);
Field<Dim<2>, Dim<2>::Vector3d>& extractFieldForNodeListVector3d2d(const FieldList<Dim<2>, Dim<2>::Vector3d>* self,
                                                                               const NodeSpace::NodeList<Dim<2> >& nodeList);
Field<Dim<2>, Dim<2>::Tensor>& extractFieldForNodeListTensor2d(const FieldList<Dim<2>, Dim<2>::Tensor>* self,
                                                                           const NodeSpace::NodeList<Dim<2> >& nodeList);
Field<Dim<2>, Dim<2>::SymTensor>& extractFieldForNodeListSymTensor2d(const FieldList<Dim<2>, Dim<2>::SymTensor>* self,
                                                                                 const NodeSpace::NodeList<Dim<2> >& nodeList);
Field<Dim<2>, Dim<2>::ThirdRankTensor>& extractFieldForNodeListThirdRankTensor2d(const FieldList<Dim<2>, Dim<2>::ThirdRankTensor>* self,
                                                                                 const NodeSpace::NodeList<Dim<2> >& nodeList);
Field<Dim<3>, int>& extractFieldForNodeListInt3d(const FieldList<Dim<3>, int>* self,
                                                             const NodeSpace::NodeList<Dim<3> >& nodeList);
Field<Dim<3>, ULL>& extractFieldForNodeListULL3d(const FieldList<Dim<3>, ULL>* self,
                                                                const NodeSpace::NodeList<Dim<3> >& nodeList);
Field<Dim<3>, Dim<3>::Scalar>& extractFieldForNodeListScalar3d(const FieldList<Dim<3>, Dim<3>::Scalar>* self,
                                                                           const NodeSpace::NodeList<Dim<3> >& nodeList);
Field<Dim<3>, Dim<3>::Vector>& extractFieldForNodeListVector3d(const FieldList<Dim<3>, Dim<3>::Vector>* self,
                                                                           const NodeSpace::NodeList<Dim<3> >& nodeList);
Field<Dim<3>, Dim<3>::Vector3d>& extractFieldForNodeListVector3d3d(const FieldList<Dim<3>, Dim<3>::Vector3d>* self,
                                                                               const NodeSpace::NodeList<Dim<3> >& nodeList);
Field<Dim<3>, Dim<3>::Tensor>& extractFieldForNodeListTensor3d(const FieldList<Dim<3>, Dim<3>::Tensor>* self,
                                                                           const NodeSpace::NodeList<Dim<3> >& nodeList);
Field<Dim<3>, Dim<3>::SymTensor>& extractFieldForNodeListSymTensor3d(const FieldList<Dim<3>, Dim<3>::SymTensor>* self,
                                                                                 const NodeSpace::NodeList<Dim<3> >& nodeList);
Field<Dim<3>, Dim<3>::ThirdRankTensor>& extractFieldForNodeListThirdRankTensor3d(const FieldList<Dim<3>, Dim<3>::ThirdRankTensor>* self,
                                                                                 const NodeSpace::NodeList<Dim<3> >& nodeList);

int sampleValueFieldListInt1d(const FieldList<Dim<1>, int>* self, const Dim<1>::Vector& r, const KernelSpace::TableKernel<Dim<1> >& W);
ULL sampleValueFieldListULL1d(const FieldList<Dim<1>, ULL>* self, const Dim<1>::Vector& r, const KernelSpace::TableKernel<Dim<1> >& W);
Dim<1>::Scalar sampleValueFieldListScalar1d(const FieldList<Dim<1>, Dim<1>::Scalar>* self, const Dim<1>::Vector& r, const KernelSpace::TableKernel<Dim<1> >& W);
Dim<1>::Vector sampleValueFieldListVector1d(const FieldList<Dim<1>, Dim<1>::Vector>* self, const Dim<1>::Vector& r, const KernelSpace::TableKernel<Dim<1> >& W);
Dim<1>::Vector3d sampleValueFieldListVector3d1d(const FieldList<Dim<1>, Dim<1>::Vector3d>* self, const Dim<1>::Vector& r, const KernelSpace::TableKernel<Dim<1> >& W);
Dim<1>::Tensor sampleValueFieldListTensor1d(const FieldList<Dim<1>, Dim<1>::Tensor>* self, const Dim<1>::Vector& r, const KernelSpace::TableKernel<Dim<1> >& W);
Dim<1>::SymTensor sampleValueFieldListSymTensor1d(const FieldList<Dim<1>, Dim<1>::SymTensor>* self, const Dim<1>::Vector& r, const KernelSpace::TableKernel<Dim<1> >& W);
Dim<1>::ThirdRankTensor sampleValueFieldListThirdRankTensor1d(const FieldList<Dim<1>, Dim<1>::ThirdRankTensor>* self, const Dim<1>::Vector& r, const KernelSpace::TableKernel<Dim<1> >& W);

int sampleValueFieldListInt2d(const FieldList<Dim<2>, int>* self, const Dim<2>::Vector& r, const KernelSpace::TableKernel<Dim<2> >& W);
ULL sampleValueFieldListULL2d(const FieldList<Dim<2>, ULL>* self, const Dim<2>::Vector& r, const KernelSpace::TableKernel<Dim<2> >& W);
Dim<2>::Scalar sampleValueFieldListScalar2d(const FieldList<Dim<2>, Dim<2>::Scalar>* self, const Dim<2>::Vector& r, const KernelSpace::TableKernel<Dim<2> >& W);
Dim<2>::Vector sampleValueFieldListVector2d(const FieldList<Dim<2>, Dim<2>::Vector>* self, const Dim<2>::Vector& r, const KernelSpace::TableKernel<Dim<2> >& W);
Dim<2>::Vector3d sampleValueFieldListVector3d2d(const FieldList<Dim<2>, Dim<2>::Vector3d>* self, const Dim<2>::Vector& r, const KernelSpace::TableKernel<Dim<2> >& W);
Dim<2>::Tensor sampleValueFieldListTensor2d(const FieldList<Dim<2>, Dim<2>::Tensor>* self, const Dim<2>::Vector& r, const KernelSpace::TableKernel<Dim<2> >& W);
Dim<2>::SymTensor sampleValueFieldListSymTensor2d(const FieldList<Dim<2>, Dim<2>::SymTensor>* self, const Dim<2>::Vector& r, const KernelSpace::TableKernel<Dim<2> >& W);
Dim<2>::ThirdRankTensor sampleValueFieldListThirdRankTensor2d(const FieldList<Dim<2>, Dim<2>::ThirdRankTensor>* self, const Dim<2>::Vector& r, const KernelSpace::TableKernel<Dim<2> >& W);

int sampleValueFieldListInt3d(const FieldList<Dim<3>, int>* self, const Dim<3>::Vector& r, const KernelSpace::TableKernel<Dim<3> >& W);
ULL sampleValueFieldListULL3d(const FieldList<Dim<3>, ULL>* self, const Dim<3>::Vector& r, const KernelSpace::TableKernel<Dim<3> >& W);
Dim<3>::Scalar sampleValueFieldListScalar3d(const FieldList<Dim<3>, Dim<3>::Scalar>* self, const Dim<3>::Vector& r, const KernelSpace::TableKernel<Dim<3> >& W);
Dim<3>::Vector sampleValueFieldListVector3d(const FieldList<Dim<3>, Dim<3>::Vector>* self, const Dim<3>::Vector& r, const KernelSpace::TableKernel<Dim<3> >& W);
Dim<3>::Vector3d sampleValueFieldListVector3d3d(const FieldList<Dim<3>, Dim<3>::Vector3d>* self, const Dim<3>::Vector& r, const KernelSpace::TableKernel<Dim<3> >& W);
Dim<3>::Tensor sampleValueFieldListTensor3d(const FieldList<Dim<3>, Dim<3>::Tensor>* self, const Dim<3>::Vector& r, const KernelSpace::TableKernel<Dim<3> >& W);
Dim<3>::SymTensor sampleValueFieldListSymTensor3d(const FieldList<Dim<3>, Dim<3>::SymTensor>* self, const Dim<3>::Vector& r, const KernelSpace::TableKernel<Dim<3> >& W);
Dim<3>::ThirdRankTensor sampleValueFieldListThirdRankTensor3d(const FieldList<Dim<3>, Dim<3>::ThirdRankTensor>* self, const Dim<3>::Vector& r, const KernelSpace::TableKernel<Dim<3> >& W);

}
}

#endif

#ifndef __GCCXML__

#ifndef __ExtendFieldListMethods_hh__
#define __ExtendFieldListMethods_hh__

namespace Spheral {
namespace FieldSpace {

//------------------------------------------------------------------------------
// Index into a FieldList (1-D).
//------------------------------------------------------------------------------
inline
IntField1d*
getItemIntFieldList1d(IntFieldList1d& self,
                      const unsigned int index) {
  return getItemByValue<IntFieldList1d, IntField1d*>(self, index);
}

inline
ULLField1d*
getItemULLFieldList1d(ULLFieldList1d& self,
                      const unsigned int index) {
  return getItemByValue<ULLFieldList1d, ULLField1d*>(self, index);
}

inline
ScalarField1d*
getItemScalarFieldList1d(ScalarFieldList1d& self,
                         const unsigned int index) {
  return getItemByValue<ScalarFieldList1d, ScalarField1d*>(self, index);
}

inline
Field<Dim<1>, Dim<1>::Vector>*
getItemVectorFieldList1d(FieldList<Dim<1>, Dim<1>::Vector>& self,
                         const unsigned int index) {
  return getItemByValue<VectorFieldList1d, VectorField1d*>(self, index);
}

inline
Field<Dim<1>, Dim<1>::Vector3d>*
getItemVector3dFieldList1d(FieldList<Dim<1>, Dim<1>::Vector3d>& self,
                           const unsigned int index) {
  return getItemByValue<Vector3dFieldList1d, Vector3dField1d*>(self, index);
}

inline
Field<Dim<1>, Dim<1>::Tensor>*
getItemTensorFieldList1d(FieldList<Dim<1>, Dim<1>::Tensor>& self,
                         const unsigned int index) {
  return getItemByValue<TensorFieldList1d, TensorField1d*>(self, index);
}

inline
Field<Dim<1>, Dim<1>::SymTensor>*
getItemSymTensorFieldList1d(FieldList<Dim<1>, Dim<1>::SymTensor>& self,
                            const unsigned int index) {
  return getItemByValue<SymTensorFieldList1d, SymTensorField1d*>(self, index);
}

inline
Field<Dim<1>, Dim<1>::ThirdRankTensor>*
getItemThirdRankTensorFieldList1d(FieldList<Dim<1>, Dim<1>::ThirdRankTensor>& self,
                            const unsigned int index) {
  return getItemByValue<ThirdRankTensorFieldList1d, ThirdRankTensorField1d*>(self, index);
}

//------------------------------------------------------------------------------
// Index into a FieldList (2-D).
//------------------------------------------------------------------------------
inline
IntField2d*
getItemIntFieldList2d(IntFieldList2d& self,
                         const unsigned int index) {
  return getItemByValue<IntFieldList2d, IntField2d*>(self, index);
}

inline
ULLField2d*
getItemULLFieldList2d(ULLFieldList2d& self,
                      const unsigned int index) {
  return getItemByValue<ULLFieldList2d, ULLField2d*>(self, index);
}

inline
ScalarField2d*
getItemScalarFieldList2d(ScalarFieldList2d& self,
                         const unsigned int index) {
  return getItemByValue<ScalarFieldList2d, ScalarField2d*>(self, index);
}

inline
Field<Dim<2>, Dim<2>::Vector>*
getItemVectorFieldList2d(FieldList<Dim<2>, Dim<2>::Vector>& self,
                         const unsigned int index) {
  return getItemByValue<VectorFieldList2d, VectorField2d*>(self, index);
}

inline
Field<Dim<2>, Dim<2>::Vector3d>*
getItemVector3dFieldList2d(FieldList<Dim<2>, Dim<2>::Vector3d>& self,
                           const unsigned int index) {
  return getItemByValue<Vector3dFieldList2d, Vector3dField2d*>(self, index);
}

inline
Field<Dim<2>, Dim<2>::Tensor>*
getItemTensorFieldList2d(FieldList<Dim<2>, Dim<2>::Tensor>& self,
                         const unsigned int index) {
  return getItemByValue<TensorFieldList2d, TensorField2d*>(self, index);
}

inline
Field<Dim<2>, Dim<2>::SymTensor>*
getItemSymTensorFieldList2d(FieldList<Dim<2>, Dim<2>::SymTensor>& self,
                            const unsigned int index) {
  return getItemByValue<SymTensorFieldList2d, SymTensorField2d*>(self, index);
}

inline
Field<Dim<2>, Dim<2>::ThirdRankTensor>*
getItemThirdRankTensorFieldList2d(FieldList<Dim<2>, Dim<2>::ThirdRankTensor>& self,
                            const unsigned int index) {
  return getItemByValue<ThirdRankTensorFieldList2d, ThirdRankTensorField2d*>(self, index);
}

//------------------------------------------------------------------------------
// Index into a FieldList (3-D).
//------------------------------------------------------------------------------
inline
IntField3d*
getItemIntFieldList3d(IntFieldList3d& self,
                      const unsigned int index) {
  return getItemByValue<IntFieldList3d, IntField3d*>(self, index);
}

inline
ULLField3d*
getItemULLFieldList3d(ULLFieldList3d& self,
                      const unsigned int index) {
  return getItemByValue<ULLFieldList3d, ULLField3d*>(self, index);
}

inline
ScalarField3d*
getItemScalarFieldList3d(ScalarFieldList3d& self,
                         const unsigned int index) {
  return getItemByValue<ScalarFieldList3d, ScalarField3d*>(self, index);
}

inline
Field<Dim<3>, Dim<3>::Vector>*
getItemVectorFieldList3d(FieldList<Dim<3>, Dim<3>::Vector>& self,
                         const unsigned int index) {
  return getItemByValue<VectorFieldList3d, VectorField3d*>(self, index);
}

inline
Field<Dim<3>, Dim<3>::Vector3d>*
getItemVector3dFieldList3d(FieldList<Dim<3>, Dim<3>::Vector3d>& self,
                           const unsigned int index) {
  return getItemByValue<Vector3dFieldList3d, Vector3dField3d*>(self, index);
}

inline
Field<Dim<3>, Dim<3>::Tensor>*
getItemTensorFieldList3d(FieldList<Dim<3>, Dim<3>::Tensor>& self,
                         const unsigned int index) {
  return getItemByValue<TensorFieldList3d, TensorField3d*>(self, index);
}

inline
Field<Dim<3>, Dim<3>::SymTensor>*
getItemSymTensorFieldList3d(FieldList<Dim<3>, Dim<3>::SymTensor>& self,
                            const unsigned int index) {
  return getItemByValue<SymTensorFieldList3d, SymTensorField3d*>(self, index);
}

inline
Field<Dim<3>, Dim<3>::ThirdRankTensor>*
getItemThirdRankTensorFieldList3d(FieldList<Dim<3>, Dim<3>::ThirdRankTensor>& self,
                            const unsigned int index) {
  return getItemByValue<ThirdRankTensorFieldList3d, ThirdRankTensorField3d*>(self, index);
}

//------------------------------------------------------------------------------
// Provide a list interface to the field members (1-D).
//------------------------------------------------------------------------------
inline
boost::python::list
intFields1d(FieldList<Dim<1>, int>* self) {
  return iteratorsAsListByRef(self->begin(),
                                       self->end());
}

inline
boost::python::list
ULLFields1d(FieldList<Dim<1>, ULL>* self) {
  return iteratorsAsListByRef(self->begin(),
                              self->end());
}

inline
boost::python::list
scalarFields1d(FieldList<Dim<1>, Dim<1>::Scalar>* self) {
  return iteratorsAsListByRef(self->begin(),
                                       self->end());
}

inline
boost::python::list
vectorFields1d(FieldList<Dim<1>, Dim<1>::Vector>* self) {
  return iteratorsAsListByRef(self->begin(),
                                       self->end());
}

inline
boost::python::list
vector3dFields1d(FieldList<Dim<1>, Dim<1>::Vector3d>* self) {
  return iteratorsAsListByRef(self->begin(),
                                       self->end());
}

inline
boost::python::list
tensorFields1d(FieldList<Dim<1>, Dim<1>::Tensor>* self) {
  return iteratorsAsListByRef(self->begin(),
                                       self->end());
}

inline
boost::python::list
symTensorFields1d(FieldList<Dim<1>, Dim<1>::SymTensor>* self) {
  return iteratorsAsListByRef(self->begin(),
                                       self->end());
}

inline
boost::python::list
thirdRankTensorFields1d(FieldList<Dim<1>, Dim<1>::ThirdRankTensor>* self) {
  return iteratorsAsListByRef(self->begin(),
                                       self->end());
}

//------------------------------------------------------------------------------
// Provide a list interface to the field members (2-D).
//------------------------------------------------------------------------------
inline
boost::python::list
intFields2d(FieldList<Dim<2>, int>* self) {
  return iteratorsAsListByRef(self->begin(),
                                       self->end());
}

inline
boost::python::list
ULLFields2d(FieldList<Dim<2>, ULL>* self) {
  return iteratorsAsListByRef(self->begin(),
                              self->end());
}

inline
boost::python::list
scalarFields2d(FieldList<Dim<2>, Dim<2>::Scalar>* self) {
  return iteratorsAsListByRef(self->begin(),
                                       self->end());
}

inline
boost::python::list
vectorFields2d(FieldList<Dim<2>, Dim<2>::Vector>* self) {
  return iteratorsAsListByRef(self->begin(),
                                       self->end());
}

inline
boost::python::list
vector3dFields2d(FieldList<Dim<2>, Dim<2>::Vector3d>* self) {
  return iteratorsAsListByRef(self->begin(),
                                       self->end());
}

inline
boost::python::list
tensorFields2d(FieldList<Dim<2>, Dim<2>::Tensor>* self) {
  return iteratorsAsListByRef(self->begin(),
                                       self->end());
}

inline
boost::python::list
symTensorFields2d(FieldList<Dim<2>, Dim<2>::SymTensor>* self) {
  return iteratorsAsListByRef(self->begin(),
                                       self->end());
}

inline
boost::python::list
thirdRankTensorFields2d(FieldList<Dim<2>, Dim<2>::ThirdRankTensor>* self) {
  return iteratorsAsListByRef(self->begin(),
                                       self->end());
}

//------------------------------------------------------------------------------
// Provide a list interface to the field members (3-D).
//------------------------------------------------------------------------------
inline
boost::python::list
intFields3d(FieldList<Dim<3>, int>* self) {
  return iteratorsAsListByRef(self->begin(),
                                       self->end());
}

inline
boost::python::list
ULLFields3d(FieldList<Dim<3>, ULL>* self) {
  return iteratorsAsListByRef(self->begin(),
                              self->end());
}

inline
boost::python::list
scalarFields3d(FieldList<Dim<3>, Dim<3>::Scalar>* self) {
  return iteratorsAsListByRef(self->begin(),
                                       self->end());
}

inline
boost::python::list
vectorFields3d(FieldList<Dim<3>, Dim<3>::Vector>* self) {
  return iteratorsAsListByRef(self->begin(),
                                       self->end());
}

inline
boost::python::list
vector3dFields3d(FieldList<Dim<3>, Dim<3>::Vector3d>* self) {
  return iteratorsAsListByRef(self->begin(),
                                       self->end());
}

inline
boost::python::list
tensorFields3d(FieldList<Dim<3>, Dim<3>::Tensor>* self) {
  return iteratorsAsListByRef(self->begin(),
                                       self->end());
}

inline
boost::python::list
symTensorFields3d(FieldList<Dim<3>, Dim<3>::SymTensor>* self) {
  return iteratorsAsListByRef(self->begin(),
                                       self->end());
}

inline
boost::python::list
thirdRankTensorFields3d(FieldList<Dim<3>, Dim<3>::ThirdRankTensor>* self) {
  return iteratorsAsListByRef(self->begin(),
                                       self->end());
}

//------------------------------------------------------------------------------
// Extract the field corresponding to the given NodeList from a FieldList.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>&
extractFieldForNodeList(const FieldList<Dimension, DataType>* self,
                        const NodeSpace::NodeList<Dimension>& nodeList) {
  VERIFY(self->haveNodeList(nodeList));
  return **(self->fieldForNodeList(nodeList));
}

// 1-D Int
inline
Field<Dim<1>, int>&
extractFieldForNodeListInt1d(const FieldList<Dim<1>,
                                                                  int>* self,
                             const NodeSpace::NodeList<Dim<1> >& nodeList) {
  return extractFieldForNodeList(self, nodeList);
}

// 1-D ULL
inline
Field<Dim<1>, ULL>&
extractFieldForNodeListULL1d(const FieldList<Dim<1>, ULL>* self,
                             const NodeSpace::NodeList<Dim<1> >& nodeList) {
  return extractFieldForNodeList(self, nodeList);
}

// 1-D Scalar
inline
Field<Dim<1>, Dim<1>::Scalar>&
extractFieldForNodeListScalar1d(const FieldList<Dim<1>,
                                                                     Dim<1>::Scalar>* self,
                                const NodeSpace::NodeList<Dim<1> >& nodeList) {
  return extractFieldForNodeList(self, nodeList);
}

// 1-D Vector
inline
Field<Dim<1>, Dim<1>::Vector>&
extractFieldForNodeListVector1d(const FieldList<Dim<1>,
                                                                     Dim<1>::Vector>* self,
                                const NodeSpace::NodeList<Dim<1> >& nodeList) {
  return extractFieldForNodeList(self, nodeList);
}

// 1-D Vector3d
inline
Field<Dim<1>, Dim<1>::Vector3d>&
extractFieldForNodeListVector3d1d(const FieldList<Dim<1>,
                                                                     Dim<1>::Vector3d>* self,
                                const NodeSpace::NodeList<Dim<1> >& nodeList) {
  return extractFieldForNodeList(self, nodeList);
}

// 1-D Tensor
inline
Field<Dim<1>, Dim<1>::Tensor>&
extractFieldForNodeListTensor1d(const FieldList<Dim<1>,
                                                                     Dim<1>::Tensor>* self,
                                const NodeSpace::NodeList<Dim<1> >& nodeList) {
  return extractFieldForNodeList(self, nodeList);
}

// 1-D SymTensor
inline
Field<Dim<1>, Dim<1>::SymTensor>&
extractFieldForNodeListSymTensor1d(const FieldList<Dim<1>,
                                                                     Dim<1>::SymTensor>* self,
                                const NodeSpace::NodeList<Dim<1> >& nodeList) {
  return extractFieldForNodeList(self, nodeList);
}

// 1-D ThirdRankTensor
inline
Field<Dim<1>, Dim<1>::ThirdRankTensor>&
extractFieldForNodeListThirdRankTensor1d(const FieldList<Dim<1>,
                                                                     Dim<1>::ThirdRankTensor>* self,
                                const NodeSpace::NodeList<Dim<1> >& nodeList) {
  return extractFieldForNodeList(self, nodeList);
}

// 2-D Int
inline
Field<Dim<2>, int>&
extractFieldForNodeListInt2d(const FieldList<Dim<2>,
                                                                  int>* self,
                             const NodeSpace::NodeList<Dim<2> >& nodeList) {
  return extractFieldForNodeList(self, nodeList);
}

// 2-D ULL
inline
Field<Dim<2>, ULL>&
extractFieldForNodeListULL2d(const FieldList<Dim<2>, ULL>* self,
                             const NodeSpace::NodeList<Dim<2> >& nodeList) {
  return extractFieldForNodeList(self, nodeList);
}

// 2-D Scalar
inline
Field<Dim<2>, Dim<2>::Scalar>&
extractFieldForNodeListScalar2d(const FieldList<Dim<2>,
                                                                     Dim<2>::Scalar>* self,
                                const NodeSpace::NodeList<Dim<2> >& nodeList) {
  return extractFieldForNodeList(self, nodeList);
}

// 2-D Vector
inline
Field<Dim<2>, Dim<2>::Vector>&
extractFieldForNodeListVector2d(const FieldList<Dim<2>,
                                                                     Dim<2>::Vector>* self,
                                const NodeSpace::NodeList<Dim<2> >& nodeList) {
  return extractFieldForNodeList(self, nodeList);
}

// 2-D Vector3d
inline
Field<Dim<2>, Dim<2>::Vector3d>&
extractFieldForNodeListVector3d2d(const FieldList<Dim<2>,
                                                                     Dim<2>::Vector3d>* self,
                                const NodeSpace::NodeList<Dim<2> >& nodeList) {
  return extractFieldForNodeList(self, nodeList);
}

// 2-D Tensor
inline
Field<Dim<2>, Dim<2>::Tensor>&
extractFieldForNodeListTensor2d(const FieldList<Dim<2>,
                                                                     Dim<2>::Tensor>* self,
                                const NodeSpace::NodeList<Dim<2> >& nodeList) {
  return extractFieldForNodeList(self, nodeList);
}

// 2-D SymTensor
inline
Field<Dim<2>, Dim<2>::SymTensor>&
extractFieldForNodeListSymTensor2d(const FieldList<Dim<2>,
                                                                     Dim<2>::SymTensor>* self,
                                const NodeSpace::NodeList<Dim<2> >& nodeList) {
  return extractFieldForNodeList(self, nodeList);
}

// 2-D ThirdRankTensor
inline
Field<Dim<2>, Dim<2>::ThirdRankTensor>&
extractFieldForNodeListThirdRankTensor2d(const FieldList<Dim<2>,
                                                                     Dim<2>::ThirdRankTensor>* self,
                                const NodeSpace::NodeList<Dim<2> >& nodeList) {
  return extractFieldForNodeList(self, nodeList);
}

// 3-D Int
inline
Field<Dim<3>, int>&
extractFieldForNodeListInt3d(const FieldList<Dim<3>,
                                                                  int>* self,
                             const NodeSpace::NodeList<Dim<3> >& nodeList) {
  return extractFieldForNodeList(self, nodeList);
}

// 3-D ULL
inline
Field<Dim<3>, ULL>&
extractFieldForNodeListULL3d(const FieldList<Dim<3>, ULL>* self,
                             const NodeSpace::NodeList<Dim<3> >& nodeList) {
  return extractFieldForNodeList(self, nodeList);
}

// 3-D Scalar
inline
Field<Dim<3>, Dim<3>::Scalar>&
extractFieldForNodeListScalar3d(const FieldList<Dim<3>,
                                                                     Dim<3>::Scalar>* self,
                                const NodeSpace::NodeList<Dim<3> >& nodeList) {
  return extractFieldForNodeList(self, nodeList);
}

// 3-D Vector
inline
Field<Dim<3>, Dim<3>::Vector>&
extractFieldForNodeListVector3d(const FieldList<Dim<3>,
                                                                     Dim<3>::Vector>* self,
                                const NodeSpace::NodeList<Dim<3> >& nodeList) {
  return extractFieldForNodeList(self, nodeList);
}

// 3-D Vector3d
inline
Field<Dim<3>, Dim<3>::Vector3d>&
extractFieldForNodeListVector3d3d(const FieldList<Dim<3>,
                                                                     Dim<3>::Vector3d>* self,
                                const NodeSpace::NodeList<Dim<3> >& nodeList) {
  return extractFieldForNodeList(self, nodeList);
}

// 3-D Tensor
inline
Field<Dim<3>, Dim<3>::Tensor>&
extractFieldForNodeListTensor3d(const FieldList<Dim<3>, Dim<3>::Tensor>* self,
                                const NodeSpace::NodeList<Dim<3> >& nodeList) {
  return extractFieldForNodeList(self, nodeList);
}

// 3-D SymTensor
inline
Field<Dim<3>, Dim<3>::SymTensor>&
extractFieldForNodeListSymTensor3d(const FieldList<Dim<3>,
                                   Dim<3>::SymTensor>* self,
                                   const NodeSpace::NodeList<Dim<3> >& nodeList) {
  return extractFieldForNodeList(self, nodeList);
}

// 3-D ThirdRankTensor
inline
Field<Dim<3>, Dim<3>::ThirdRankTensor>&
extractFieldForNodeListThirdRankTensor3d(const FieldList<Dim<3>,
                                   Dim<3>::ThirdRankTensor>* self,
                                   const NodeSpace::NodeList<Dim<3> >& nodeList) {
  return extractFieldForNodeList(self, nodeList);
}

//------------------------------------------------------------------------------
// Kernel sampling of a FieldList.
//------------------------------------------------------------------------------
// 1-D Int
inline
int
sampleValueFieldListInt1d(const FieldList<Dim<1>, int>* self,
                          const Dim<1>::Vector& r,
                          const KernelSpace::TableKernel<Dim<1> >& W) {
  return (*self)(r, W);
}

// 1-D ULL
inline
ULL
sampleValueFieldListULL1d(const FieldList<Dim<1>, ULL>* self,
                          const Dim<1>::Vector& r,
                          const KernelSpace::TableKernel<Dim<1> >& W) {
  return (*self)(r, W);
}

// 1-D Scalar
inline
Dim<1>::Scalar
sampleValueFieldListScalar1d(const FieldList<Dim<1>, Dim<1>::Scalar>* self,
                          const Dim<1>::Vector& r,
                          const KernelSpace::TableKernel<Dim<1> >& W) {
  return (*self)(r, W);
}

// 1-D Vector
inline
Dim<1>::Vector
sampleValueFieldListVector1d(const FieldList<Dim<1>, Dim<1>::Vector>* self,
                          const Dim<1>::Vector& r,
                          const KernelSpace::TableKernel<Dim<1> >& W) {
  return (*self)(r, W);
}

// 1-D Vector3d
inline
Dim<1>::Vector3d
sampleValueFieldListVector3d1d(const FieldList<Dim<1>, Dim<1>::Vector3d>* self,
                          const Dim<1>::Vector& r,
                          const KernelSpace::TableKernel<Dim<1> >& W) {
  return (*self)(r, W);
}

// 1-D Tensor
inline
Dim<1>::Tensor
sampleValueFieldListTensor1d(const FieldList<Dim<1>, Dim<1>::Tensor>* self,
                          const Dim<1>::Vector& r,
                          const KernelSpace::TableKernel<Dim<1> >& W) {
  return (*self)(r, W);
}

// 1-D SymTensor
inline
Dim<1>::SymTensor
sampleValueFieldListSymTensor1d(const FieldList<Dim<1>, Dim<1>::SymTensor>* self,
                          const Dim<1>::Vector& r,
                          const KernelSpace::TableKernel<Dim<1> >& W) {
  return (*self)(r, W);
}

// 1-D ThirdRankTensor
inline
Dim<1>::ThirdRankTensor
sampleValueFieldListThirdRankTensor1d(const FieldList<Dim<1>, Dim<1>::ThirdRankTensor>* self,
                          const Dim<1>::Vector& r,
                          const KernelSpace::TableKernel<Dim<1> >& W) {
  return (*self)(r, W);
}

// 2-D Int
inline
int
sampleValueFieldListInt2d(const FieldList<Dim<2>, int>* self,
                          const Dim<2>::Vector& r,
                          const KernelSpace::TableKernel<Dim<2> >& W) {
  return (*self)(r, W);
}

// 2-D ULL
inline
ULL
sampleValueFieldListULL2d(const FieldList<Dim<2>, ULL>* self,
                          const Dim<2>::Vector& r,
                          const KernelSpace::TableKernel<Dim<2> >& W) {
  return (*self)(r, W);
}

// 2-D Scalar
inline
Dim<2>::Scalar
sampleValueFieldListScalar2d(const FieldList<Dim<2>, Dim<2>::Scalar>* self,
                          const Dim<2>::Vector& r,
                          const KernelSpace::TableKernel<Dim<2> >& W) {
  return (*self)(r, W);
}

// 2-D Vector
inline
Dim<2>::Vector
sampleValueFieldListVector2d(const FieldList<Dim<2>, Dim<2>::Vector>* self,
                          const Dim<2>::Vector& r,
                          const KernelSpace::TableKernel<Dim<2> >& W) {
  return (*self)(r, W);
}

// 2-D Vector3d
inline
Dim<2>::Vector3d
sampleValueFieldListVector3d2d(const FieldList<Dim<2>, Dim<2>::Vector3d>* self,
                          const Dim<2>::Vector& r,
                          const KernelSpace::TableKernel<Dim<2> >& W) {
  return (*self)(r, W);
}

// 2-D Tensor
inline
Dim<2>::Tensor
sampleValueFieldListTensor2d(const FieldList<Dim<2>, Dim<2>::Tensor>* self,
                          const Dim<2>::Vector& r,
                          const KernelSpace::TableKernel<Dim<2> >& W) {
  return (*self)(r, W);
}

// 2-D SymTensor
inline
Dim<2>::SymTensor
sampleValueFieldListSymTensor2d(const FieldList<Dim<2>, Dim<2>::SymTensor>* self,
                          const Dim<2>::Vector& r,
                          const KernelSpace::TableKernel<Dim<2> >& W) {
  return (*self)(r, W);
}

// 2-D ThirdRankTensor
inline
Dim<2>::ThirdRankTensor
sampleValueFieldListThirdRankTensor2d(const FieldList<Dim<2>, Dim<2>::ThirdRankTensor>* self,
                          const Dim<2>::Vector& r,
                          const KernelSpace::TableKernel<Dim<2> >& W) {
  return (*self)(r, W);
}

// 3-D Int
inline
int
sampleValueFieldListInt3d(const FieldList<Dim<3>, int>* self,
                          const Dim<3>::Vector& r,
                          const KernelSpace::TableKernel<Dim<3> >& W) {
  return (*self)(r, W);
}

// 3-D ULL
inline
ULL
sampleValueFieldListULL3d(const FieldList<Dim<3>, ULL>* self,
                          const Dim<3>::Vector& r,
                          const KernelSpace::TableKernel<Dim<3> >& W) {
  return (*self)(r, W);
}

// 3-D Scalar
inline
Dim<3>::Scalar
sampleValueFieldListScalar3d(const FieldList<Dim<3>, Dim<3>::Scalar>* self,
                          const Dim<3>::Vector& r,
                          const KernelSpace::TableKernel<Dim<3> >& W) {
  return (*self)(r, W);
}

// 3-D Vector
inline
Dim<3>::Vector
sampleValueFieldListVector3d(const FieldList<Dim<3>, Dim<3>::Vector>* self,
                          const Dim<3>::Vector& r,
                          const KernelSpace::TableKernel<Dim<3> >& W) {
  return (*self)(r, W);
}

// 3-D Vector3d
inline
Dim<3>::Vector3d
sampleValueFieldListVector3d3d(const FieldList<Dim<3>, Dim<3>::Vector3d>* self,
                          const Dim<3>::Vector& r,
                          const KernelSpace::TableKernel<Dim<3> >& W) {
  return (*self)(r, W);
}

// 3-D Tensor
inline
Dim<3>::Tensor
sampleValueFieldListTensor3d(const FieldList<Dim<3>, Dim<3>::Tensor>* self,
                          const Dim<3>::Vector& r,
                          const KernelSpace::TableKernel<Dim<3> >& W) {
  return (*self)(r, W);
}

// 3-D SymTensor
inline
Dim<3>::SymTensor
sampleValueFieldListSymTensor3d(const FieldList<Dim<3>, Dim<3>::SymTensor>* self,
                          const Dim<3>::Vector& r,
                          const KernelSpace::TableKernel<Dim<3> >& W) {
  return (*self)(r, W);
}

// 3-D ThirdRankTensor
inline
Dim<3>::ThirdRankTensor
sampleValueFieldListThirdRankTensor3d(const FieldList<Dim<3>, Dim<3>::ThirdRankTensor>* self,
                          const Dim<3>::Vector& r,
                          const KernelSpace::TableKernel<Dim<3> >& W) {
  return (*self)(r, W);
}

}
}

#endif
#endif
