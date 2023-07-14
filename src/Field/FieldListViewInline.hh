
// Inlined methods.
namespace Spheral {

template<typename Dimension, typename DataType>
FieldListView<Dimension, DataType>::
FieldListView(const FieldList<Dimension, DataType>& field) :
  mFieldViewsView{field.mFieldPtrs.toView()}
{}


template<typename Dimension, typename DataType>
inline
void
FieldListView<Dimension, DataType>::
move(LvArray::MemorySpace const& space, bool touch) const {
  mFieldViewsView.move(space,touch);
}


template<typename Dimension, typename DataType>
inline
DataType&
FieldListView<Dimension, DataType>::
operator()(const unsigned int field_id, const unsigned int idx) const {
  return mFieldViewsView(field_id)[idx];
}

template<typename Dimension, typename DataType>
inline
typename FieldListView<Dimension, DataType>::FieldViewType&
FieldListView<Dimension, DataType>::
operator[](const unsigned int index) {
  return mFieldViewsView(index);
}

template<typename Dimension, typename DataType>
inline
typename FieldListView<Dimension, DataType>::FieldViewType&
FieldListView<Dimension, DataType>::
operator[](const unsigned int index) const {
  return mFieldViewsView(index);
}


template<typename Dimension, typename DataType>
inline
DataType&
FieldListView<Dimension, DataType>::
pool(const unsigned field_id, const unsigned idx) const {
  return mFieldViewsView(field_id).pool(idx);
}

template<typename Dimension, typename DataType>
inline
typename FieldListView<Dimension, DataType>::AtomicDataType&
FieldListView<Dimension, DataType>::
atomic(const unsigned field_id, const unsigned idx) const {
  return *reinterpret_cast<AtomicDataType*>(&mFieldViewsView(field_id)[idx]);
}

template<typename Dimension, typename DataType>
inline
typename FieldListView<Dimension, DataType>::AtomicDataType&
FieldListView<Dimension, DataType>::
pool_atomic(const unsigned field_id, const unsigned idx) const {
  return *reinterpret_cast<AtomicDataType*>(&mFieldViewsView(field_id).pool(idx));
}

template<typename Dimension, typename DataType>
inline
FieldListView<Dimension, DataType>::
FieldListView( FieldListView const & source) noexcept :
  mFieldViewsView{source.mFieldViewsView}
{}



} // namespace Spheral
