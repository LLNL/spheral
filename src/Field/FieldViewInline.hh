
// Inlined methods.
namespace Spheral {

template<typename Dimension, typename DataType>
inline
FieldView<Dimension, DataType>::
FieldView() {}

template<typename Dimension, typename DataType>
inline
FieldView<Dimension, DataType>::
FieldView(const FieldView::FIELD_TYPE& field) :
  mDataView{field.mDataArray},
  mDataPoolView{field.mDataArray},
  mFieldPtr{const_cast<FIELD_TYPE*>(&field)}
{}

template<typename Dimension, typename DataType>
inline
FieldView<Dimension, DataType>::
FieldView(const FieldView::FIELD_TYPE& field,
          const FieldView::FIELD_TYPE& pool) :
  mDataView{field.mDataArray},
  mDataPoolView{pool.mDataArray},
  mFieldPtr{const_cast<FIELD_TYPE*>(&field)}
{}

template<typename Dimension, typename DataType>
inline
void
FieldView<Dimension, DataType>::
move(LvArray::MemorySpace const& space, bool touch) const {
  mDataView.move(space,touch);
}

template<typename Dimension, typename DataType>
inline
typename FieldView<Dimension, DataType>::ARRAY_VIEW_TYPE&
FieldView<Dimension, DataType>::
getView() {
  return mDataView;
}

template<typename Dimension, typename DataType>
inline
typename FieldView<Dimension, DataType>::FIELD_TYPE&
FieldView<Dimension, DataType>::
operator*() const{
  return *mFieldPtr;
}

template<typename Dimension, typename DataType>
inline
typename FieldView<Dimension, DataType>::FIELD_TYPE&
FieldView<Dimension, DataType>::
get() const {
  return *mFieldPtr;
}

template<typename Dimension, typename DataType>
inline
typename FieldView<Dimension, DataType>::FIELD_TYPE*
FieldView<Dimension, DataType>::
operator->() const {
  return mFieldPtr;
}

template<typename Dimension, typename DataType>
inline
DataType&
FieldView<Dimension, DataType>::
operator[](const unsigned int index) {
  return mDataView(index);
}

template<typename Dimension, typename DataType>
inline
DataType&
FieldView<Dimension, DataType>::
operator[](const unsigned int index) const {
  return mDataView(index);
}

template<typename Dimension, typename DataType>
inline
DataType&
FieldView<Dimension, DataType>::
pool(const unsigned index) const {
  return mDataPoolView(index);
}

template<typename Dimension, typename DataType>
inline constexpr 
FieldView<Dimension, DataType>::
FieldView(const FieldView & source) noexcept :
  mDataView{source.mDataView},
  mDataPoolView{source.mDataPoolView},
  mFieldPtr{source.mFieldPtr}
{}

template<typename Dimension, typename DataType>
inline
FieldView<Dimension, DataType>&
FieldView<Dimension, DataType>::
operator=(const FieldView& rhs) {
  mDataView = rhs.mDataView;
  mDataPoolView = rhs.mDataPoolView;
  mFieldPtr = rhs.mFieldPtr;
  return *this;
}

template<typename Dimension, typename DataType>
inline
bool
FieldView<Dimension, DataType>::
operator==(const FIELD_TYPE& rhs) const { 
  const FieldView temp(const_cast<FIELD_TYPE&>(rhs));
  return *this == temp;
}

template<typename Dimension, typename DataType>
inline
bool
FieldView<Dimension, DataType>::
operator==(const FieldView& rhs) const { 
  return mFieldPtr == rhs.mFieldPtr;
}

} // namespace Spheral
