#ifndef __PBGWRAPS_FIELDTYPES__
#define __PBGWRAPS_FIELDTYPES__

#include <vector>
#include "Geometry/Dimension.hh"
#include "Field/FieldBase.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Field/FieldListSet.hh"
#include "Utilities/FieldDataTypeTraits.hh"
#include "PBGWraps/CXXTypes/CXXTypes.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Names!
//------------------------------------------------------------------------------
// This is kind of stupid, but it's convenient with our naming scheme.
typedef int Int;

typedef FieldBase<Dim<1> > FieldBase1d;
typedef FieldBase<Dim<2> > FieldBase2d;
typedef FieldBase<Dim<3> > FieldBase3d;

typedef FieldListBase<Dim<1> > FieldListBase1d;
typedef FieldListBase<Dim<2> > FieldListBase2d;
typedef FieldListBase<Dim<3> > FieldListBase3d;

typedef Field<Dim<1>, int> IntField1d;
typedef Field<Dim<1>, uint64_t> ULLField1d;
typedef Field<Dim<1>, Dim<1>::Scalar > ScalarField1d;
typedef Field<Dim<1>, Dim<1>::Vector > VectorField1d;
typedef Field<Dim<1>, Dim<1>::Vector3d > Vector3dField1d;
typedef Field<Dim<1>, Dim<1>::Tensor > TensorField1d;
typedef Field<Dim<1>, Dim<1>::SymTensor > SymTensorField1d;
typedef Field<Dim<1>, Dim<1>::ThirdRankTensor > ThirdRankTensorField1d;
typedef Field<Dim<1>, Dim<1>::FourthRankTensor > FourthRankTensorField1d;
typedef Field<Dim<1>, Dim<1>::FifthRankTensor > FifthRankTensorField1d;
typedef Field<Dim<1>, Dim<1>::FacetedVolume > FacetedVolumeField1d;
typedef Field<Dim<1>, std::vector<double> > VectorDoubleField1d;
typedef Field<Dim<1>, std::vector<Dim<1>::Vector> > VectorVectorField1d;
typedef Field<Dim<1>, std::vector<Dim<1>::Tensor> > VectorTensorField1d;
typedef Field<Dim<1>, std::vector<Dim<1>::SymTensor> > VectorSymTensorField1d;

typedef Field<Dim<2>, int> IntField2d;
typedef Field<Dim<2>, uint64_t> ULLField2d;
typedef Field<Dim<2>, Dim<2>::Scalar > ScalarField2d;
typedef Field<Dim<2>, Dim<2>::Vector3d > Vector3dField2d;
typedef Field<Dim<2>, Dim<2>::Vector > VectorField2d;
typedef Field<Dim<2>, Dim<2>::Tensor > TensorField2d;
typedef Field<Dim<2>, Dim<2>::SymTensor > SymTensorField2d;
typedef Field<Dim<2>, Dim<2>::ThirdRankTensor > ThirdRankTensorField2d;
typedef Field<Dim<2>, Dim<2>::FourthRankTensor > FourthRankTensorField2d;
typedef Field<Dim<2>, Dim<2>::FifthRankTensor > FifthRankTensorField2d;
typedef Field<Dim<2>, Dim<2>::FacetedVolume > FacetedVolumeField2d;
typedef Field<Dim<2>, std::vector<double> > VectorDoubleField2d;
typedef Field<Dim<2>, std::vector<Dim<2>::Vector> > VectorVectorField2d;
typedef Field<Dim<2>, std::vector<Dim<2>::Tensor> > VectorTensorField2d;
typedef Field<Dim<2>, std::vector<Dim<2>::SymTensor> > VectorSymTensorField2d;

typedef Field<Dim<3>, int> IntField3d;
typedef Field<Dim<3>, uint64_t> ULLField3d;
typedef Field<Dim<3>, Dim<3>::Scalar > ScalarField3d;
typedef Field<Dim<3>, Dim<3>::Vector > VectorField3d;
typedef Field<Dim<3>, Dim<3>::Vector3d > Vector3dField3d;
typedef Field<Dim<3>, Dim<3>::Tensor > TensorField3d;
typedef Field<Dim<3>, Dim<3>::SymTensor > SymTensorField3d;
typedef Field<Dim<3>, Dim<3>::ThirdRankTensor > ThirdRankTensorField3d;
typedef Field<Dim<3>, Dim<3>::FourthRankTensor > FourthRankTensorField3d;
typedef Field<Dim<3>, Dim<3>::FifthRankTensor > FifthRankTensorField3d;
typedef Field<Dim<3>, Dim<3>::FacetedVolume > FacetedVolumeField3d;
typedef Field<Dim<3>, std::vector<double> > VectorDoubleField3d;
typedef Field<Dim<3>, std::vector<Dim<3>::Vector> > VectorVectorField3d;
typedef Field<Dim<3>, std::vector<Dim<3>::Tensor> > VectorTensorField3d;
typedef Field<Dim<3>, std::vector<Dim<3>::SymTensor> > VectorSymTensorField3d;

typedef FieldList<Dim<1>, int> IntFieldList1d;
typedef FieldList<Dim<1>, uint64_t> ULLFieldList1d;
typedef FieldList<Dim<1>, Dim<1>::Scalar > ScalarFieldList1d;
typedef FieldList<Dim<1>, Dim<1>::Vector > VectorFieldList1d;
typedef FieldList<Dim<1>, Dim<1>::Vector3d > Vector3dFieldList1d;
typedef FieldList<Dim<1>, Dim<1>::Tensor > TensorFieldList1d;
typedef FieldList<Dim<1>, Dim<1>::SymTensor > SymTensorFieldList1d;
typedef FieldList<Dim<1>, Dim<1>::ThirdRankTensor > ThirdRankTensorFieldList1d;
typedef FieldList<Dim<1>, Dim<1>::FourthRankTensor > FourthRankTensorFieldList1d;
typedef FieldList<Dim<1>, Dim<1>::FifthRankTensor > FifthRankTensorFieldList1d;
typedef FieldList<Dim<1>, Dim<1>::FacetedVolume > FacetedVolumeFieldList1d;
typedef FieldList<Dim<1>, std::vector<double> > VectorDoubleFieldList1d;
typedef FieldList<Dim<1>, std::vector<Dim<1>::Vector> > VectorVectorFieldList1d;
typedef FieldList<Dim<1>, std::vector<Dim<1>::Tensor> > VectorTensorFieldList1d;
typedef FieldList<Dim<1>, std::vector<Dim<1>::SymTensor> > VectorSymTensorFieldList1d;

typedef FieldList<Dim<2>, int> IntFieldList2d;
typedef FieldList<Dim<2>, uint64_t> ULLFieldList2d;
typedef FieldList<Dim<2>, Dim<2>::Scalar > ScalarFieldList2d;
typedef FieldList<Dim<2>, Dim<2>::Vector > VectorFieldList2d;
typedef FieldList<Dim<2>, Dim<2>::Vector3d > Vector3dFieldList2d;
typedef FieldList<Dim<2>, Dim<2>::Tensor > TensorFieldList2d;
typedef FieldList<Dim<2>, Dim<2>::SymTensor > SymTensorFieldList2d;
typedef FieldList<Dim<2>, Dim<2>::ThirdRankTensor > ThirdRankTensorFieldList2d;
typedef FieldList<Dim<2>, Dim<2>::FourthRankTensor > FourthRankTensorFieldList2d;
typedef FieldList<Dim<2>, Dim<2>::FifthRankTensor > FifthRankTensorFieldList2d;
typedef FieldList<Dim<2>, Dim<2>::FacetedVolume > FacetedVolumeFieldList2d;
typedef FieldList<Dim<2>, std::vector<double> > VectorDoubleFieldList2d;
typedef FieldList<Dim<2>, std::vector<Dim<2>::Vector> > VectorVectorFieldList2d;
typedef FieldList<Dim<2>, std::vector<Dim<2>::Tensor> > VectorTensorFieldList2d;
typedef FieldList<Dim<2>, std::vector<Dim<2>::SymTensor> > VectorSymTensorFieldList2d;

typedef FieldList<Dim<3>, int> IntFieldList3d;
typedef FieldList<Dim<3>, uint64_t> ULLFieldList3d;
typedef FieldList<Dim<3>, Dim<3>::Scalar > ScalarFieldList3d;
typedef FieldList<Dim<3>, Dim<3>::Vector > VectorFieldList3d;
typedef FieldList<Dim<3>, Dim<3>::Vector3d > Vector3dFieldList3d;
typedef FieldList<Dim<3>, Dim<3>::Tensor > TensorFieldList3d;
typedef FieldList<Dim<3>, Dim<3>::SymTensor > SymTensorFieldList3d;
typedef FieldList<Dim<3>, Dim<3>::ThirdRankTensor > ThirdRankTensorFieldList3d;
typedef FieldList<Dim<3>, Dim<3>::FourthRankTensor > FourthRankTensorFieldList3d;
typedef FieldList<Dim<3>, Dim<3>::FifthRankTensor > FifthRankTensorFieldList3d;
typedef FieldList<Dim<3>, Dim<3>::FacetedVolume > FacetedVolumeFieldList3d;
typedef FieldList<Dim<3>, std::vector<double> > VectorDoubleFieldList3d;
typedef FieldList<Dim<3>, std::vector<Dim<3>::Vector> > VectorVectorFieldList3d;
typedef FieldList<Dim<3>, std::vector<Dim<3>::Tensor> > VectorTensorFieldList3d;
typedef FieldList<Dim<3>, std::vector<Dim<3>::SymTensor> > VectorSymTensorFieldList3d;

typedef FieldListSet<Dim<1> > FieldListSet1d;
typedef FieldListSet<Dim<2> > FieldListSet2d;
typedef FieldListSet<Dim<3> > FieldListSet3d;

}

// std::vectors of Field*'s
typedef std::vector<Spheral::Field<Dim<1>, int>*> vector_of_IntFieldPtr1d;
typedef std::vector<Spheral::Field<Dim<1>, Dim<1>::Scalar>*> vector_of_ScalarFieldPtr1d;
typedef std::vector<Spheral::Field<Dim<1>, Dim<1>::Vector>*> vector_of_VectorFieldPtr1d;
typedef std::vector<Spheral::Field<Dim<1>, Dim<1>::Tensor>*> vector_of_TensorFieldPtr1d;
typedef std::vector<Spheral::Field<Dim<1>, Dim<1>::SymTensor>*> vector_of_SymTensorFieldPtr1d;

typedef std::vector<Spheral::Field<Dim<2>, int>*> vector_of_IntFieldPtr2d;
typedef std::vector<Spheral::Field<Dim<2>, Dim<2>::Scalar>*> vector_of_ScalarFieldPtr2d;
typedef std::vector<Spheral::Field<Dim<2>, Dim<2>::Vector>*> vector_of_VectorFieldPtr2d;
typedef std::vector<Spheral::Field<Dim<2>, Dim<2>::Tensor>*> vector_of_TensorFieldPtr2d;
typedef std::vector<Spheral::Field<Dim<2>, Dim<2>::SymTensor>*> vector_of_SymTensorFieldPtr2d;

typedef std::vector<Spheral::Field<Dim<3>, int>*> vector_of_IntFieldPtr3d;
typedef std::vector<Spheral::Field<Dim<3>, Dim<3>::Scalar>*> vector_of_ScalarFieldPtr3d;
typedef std::vector<Spheral::Field<Dim<3>, Dim<3>::Vector>*> vector_of_VectorFieldPtr3d;
typedef std::vector<Spheral::Field<Dim<3>, Dim<3>::Tensor>*> vector_of_TensorFieldPtr3d;
typedef std::vector<Spheral::Field<Dim<3>, Dim<3>::SymTensor>*> vector_of_SymTensorFieldPtr3d;

// std::vectors of FieldLists
typedef std::vector<Spheral::FieldList<Dim<1>, int> > vector_of_IntFieldList1d;
typedef std::vector<Spheral::FieldList<Dim<1>, Dim<1>::Scalar> > vector_of_ScalarFieldList1d;
typedef std::vector<Spheral::FieldList<Dim<1>, Dim<1>::Vector> > vector_of_VectorFieldList1d;
typedef std::vector<Spheral::FieldList<Dim<1>, Dim<1>::Tensor> > vector_of_TensorFieldList1d;
typedef std::vector<Spheral::FieldList<Dim<1>, Dim<1>::SymTensor> > vector_of_SymTensorFieldList1d;

typedef std::vector<Spheral::FieldList<Dim<2>, int> > vector_of_IntFieldList2d;
typedef std::vector<Spheral::FieldList<Dim<2>, Dim<2>::Scalar> > vector_of_ScalarFieldList2d;
typedef std::vector<Spheral::FieldList<Dim<2>, Dim<2>::Vector> > vector_of_VectorFieldList2d;
typedef std::vector<Spheral::FieldList<Dim<2>, Dim<2>::Tensor> > vector_of_TensorFieldList2d;
typedef std::vector<Spheral::FieldList<Dim<2>, Dim<2>::SymTensor> > vector_of_SymTensorFieldList2d;

typedef std::vector<Spheral::FieldList<Dim<3>, int> > vector_of_IntFieldList3d;
typedef std::vector<Spheral::FieldList<Dim<3>, Dim<3>::Scalar> > vector_of_ScalarFieldList3d;
typedef std::vector<Spheral::FieldList<Dim<3>, Dim<3>::Vector> > vector_of_VectorFieldList3d;
typedef std::vector<Spheral::FieldList<Dim<3>, Dim<3>::Tensor> > vector_of_TensorFieldList3d;
typedef std::vector<Spheral::FieldList<Dim<3>, Dim<3>::SymTensor> > vector_of_SymTensorFieldList3d;

namespace Spheral {

// //------------------------------------------------------------------------------
// // Index into a Field.
// //------------------------------------------------------------------------------
// template<typename FieldType>
// inline
// typename FieldType::FieldDataType
// indexField(FieldType& container, 
//            int index) {
//   if (index < 0) index += container.size();
//   if (index < container.size()) {
//     return container[index];
//   } else {
//     PyErr_SetString(PyExc_IndexError, "Field index out of range");
//     return typename FieldType::FieldDataType();
//   }
// }

// template<typename FieldType>
// inline
// typename FieldType::FieldDataType*
// indexFieldAsPointer(FieldType& container, 
//                     int index) {
//   if (index < 0) index += container.size();
//   if (index < container.size()) {
//     return &(container[index]);
//   } else {
//     PyErr_SetString(PyExc_IndexError, "Field index out of range");
//     return NULL;
//   }
// }

// //------------------------------------------------------------------------------
// // Assign to a postion in a Field.
// //------------------------------------------------------------------------------
// template<typename FieldType>
// inline
// int
// assignToFieldIndex(FieldType& container, 
//                    int index,
//                    const typename FieldType::FieldDataType& value) {
//   if (index < 0) index += container.size();
//   if (index >= container.size()) {
//     PyErr_SetString(PyExc_IndexError, "Field index out of range");
//   } else {
//     container[index] = value;
//   }
// }

// //------------------------------------------------------------------------------
// // Index into a FieldList (extracting Fields).
// //------------------------------------------------------------------------------
// template<typename FieldListType>
// inline
// typename FieldListType::ElementType
// indexFieldList(FieldListType& container, 
//                const size_t index) {
//   if (index < container.size()) {
//     return container[index];
//   } else {
//     PyErr_SetString(PyExc_IndexError, "FieldList index out of range");
//     return NULL;
//   }
// }

template<typename Dimension, typename FieldListType>
inline
typename FieldListType::ElementType
fieldForNodeList(FieldListType& container, 
                 const NodeList<Dimension>& nodeList) {
  return *(container.fieldForNodeList(nodeList));
}

//------------------------------------------------------------------------------
// Index into a FieldList (extracting individual values).
//------------------------------------------------------------------------------
template<typename FieldListType>
inline
typename FieldListType::FieldDataType
indexFieldListForValue(FieldListType& container, 
                       const size_t fieldIndex,
                       const size_t nodeIndex) {
  if (fieldIndex < container.size() and
      nodeIndex < container[fieldIndex]->size()) {
    return container(fieldIndex, nodeIndex);
  } else {
    PyErr_SetString(PyExc_IndexError, "FieldList (fieldIndex, nodeIndex) indicies out of range");
    return typename FieldListType::FieldDataType();
  }
}

template<typename FieldListType>
inline
typename FieldListType::FieldDataType*
indexFieldListForValuePointer(FieldListType& container, 
                              const size_t fieldIndex,
                              const size_t nodeIndex) {
  if (fieldIndex < container.size() and
      nodeIndex < container[fieldIndex]->size()) {
    return &(container(fieldIndex, nodeIndex));
  } else {
    PyErr_SetString(PyExc_IndexError, "FieldList (fieldIndex, nodeIndex) indicies out of range");
    return NULL;
  }
}

}

#endif
