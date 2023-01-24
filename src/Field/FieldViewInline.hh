
// Inlined methods.
namespace Spheral {

template<typename Dimension, typename DataType>
inline
FieldView<Dimension, DataType>::
FieldView(const FieldView::FIELD_TYPE& field) :
  mDataView{field.mDataArray},
  mDataPoolView{field.mDataArray}
{}

template<typename Dimension, typename DataType>
inline
FieldView<Dimension, DataType>::
FieldView(const FieldView::FIELD_TYPE& field,
          const FieldView::FIELD_TYPE& pool) :
  mDataView{field.mDataArray},
  mDataPoolView{pool.mDataArray}
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
FieldView( FieldView const & source) noexcept :
  mDataView{source.mDataView},
  mDataPoolView{source.mDataPoolView}
{}

} // namespace Spheral
