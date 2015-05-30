namespace Spheral {

//------------------------------------------------------------------------------
// Min value.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
inline
BoundValueType
ReplaceBoundedState<Dimension, ValueType, BoundValueType>::
minValue() const {
  return mMinValue;
}

//------------------------------------------------------------------------------
// Max value.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
inline
BoundValueType
ReplaceBoundedState<Dimension, ValueType, BoundValueType>::
maxValue() const {
  return mMaxValue;
}

}
