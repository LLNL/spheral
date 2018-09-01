namespace Spheral {

//------------------------------------------------------------------------------
// Access the member data as pointers.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
std::vector< FieldList<Dimension, typename Dimension::Scalar> >*
FieldListSet<Dimension>::
ScalarFieldListPtrs() {
  return &ScalarFieldLists;
}

template<typename Dimension>
inline
std::vector< FieldList<Dimension, typename Dimension::Vector> >*
FieldListSet<Dimension>::
VectorFieldListPtrs() {
  return &VectorFieldLists;
}

template<typename Dimension>
inline
std::vector< FieldList<Dimension, typename Dimension::Tensor> >*
FieldListSet<Dimension>::
TensorFieldListPtrs() {
  return &TensorFieldLists;
}

template<typename Dimension>
inline
std::vector< FieldList<Dimension, typename Dimension::SymTensor> >*
FieldListSet<Dimension>::
SymTensorFieldListPtrs() {
  return &SymTensorFieldLists;
}

}

