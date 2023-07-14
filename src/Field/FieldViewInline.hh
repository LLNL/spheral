
// Inlined methods.
namespace Spheral {

template<typename Dimension, typename DataType>
inline
FieldView<Dimension, DataType>::
FieldView() {}

template<typename Dimension, typename DataType>
inline
FieldView<Dimension, DataType>::
FieldView(const FieldView::FieldType& field) :
  mDataView{field.mDataArray},
  mDataPoolView{field.mDataArray},
  mFieldPtr{const_cast<FieldType*>(&field)}
{}

template<typename Dimension, typename DataType>
inline
FieldView<Dimension, DataType>::
FieldView(const FieldView::FieldType& field,
          const FieldView::FieldType& pool) :
  mDataView{field.mDataArray},
  mDataPoolView{pool.mDataArray},
  mFieldPtr{const_cast<FieldType*>(&field)}
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
typename FieldView<Dimension, DataType>::ArrayViewType&
FieldView<Dimension, DataType>::
getView() {
  return mDataView;
}

template<typename Dimension, typename DataType>
inline
typename FieldView<Dimension, DataType>::FieldType&
FieldView<Dimension, DataType>::
operator*() const{
  return *mFieldPtr;
}

template<typename Dimension, typename DataType>
inline
typename FieldView<Dimension, DataType>::FieldType&
FieldView<Dimension, DataType>::
get() const {
  return *mFieldPtr;
}

template<typename Dimension, typename DataType>
inline
typename FieldView<Dimension, DataType>::FieldType*
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
operator==(const FieldType& rhs) const { 
  const FieldView temp(const_cast<FieldType&>(rhs));
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
