namespace Spheral {

//------------------------------------------------------------------------------
// Enroll the given policy.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
State<Dimension>::
enroll(const typename State<Dimension>::KeyType& key,
       typename State<Dimension>::PolicyPointer polptr) {
  KeyType fieldKey, nodeKey;
  this->splitFieldKey(key, fieldKey, nodeKey);
  mPolicyMap[fieldKey][key] = polptr;
}

//------------------------------------------------------------------------------
// Remove the policy associated with a FieldList.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
void
State<Dimension>::
removePolicy(FieldList<Dimension, DataType>& fieldList) {
  this->removePolicy(StateBase<Dimension>::key(fieldList));
}

//------------------------------------------------------------------------------
// Enroll the given field and policy.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
State<Dimension>::
enroll(FieldBase<Dimension>& field,
       typename State<Dimension>::PolicyPointer polptr) {
  this->enroll(field);
  this->enroll(this->key(field), polptr);
}

//------------------------------------------------------------------------------
// Enroll the given FieldList and policy.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
inline
void
State<Dimension>::
enroll(FieldList<Dimension, DataType>& fieldList,
       typename State<Dimension>::PolicyPointer polptr) {
  this->enroll(fieldList);
  this->enroll(this->key(fieldList), polptr);
}

//------------------------------------------------------------------------------
// Enroll the given field.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
State<Dimension>::
enroll(FieldBase<Dimension>& field) {
  StateBase<Dimension>::enroll(field);
}

//------------------------------------------------------------------------------
// Enroll the given field shared_pointer.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
State<Dimension>::
enroll(std::shared_ptr<FieldBase<Dimension>>& fieldPtr) {
  StateBase<Dimension>::enroll(fieldPtr);
}

//------------------------------------------------------------------------------
// Enroll the given field list.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
inline
void
State<Dimension>::
enroll(FieldList<Dimension, DataType>& fieldList) {
  StateBase<Dimension>::enroll(fieldList);
}

//------------------------------------------------------------------------------
// Return the policy for the specified field.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename Value>
inline
typename State<Dimension>::PolicyPointer
State<Dimension>::
policy(const Field<Dimension, Value>& field) const {
  const KeyType key = StateBase<Dimension>::key(field);
  return this->policy(key);
}
  
//------------------------------------------------------------------------------
// Optionally trip a flag indicating policies should time advance only -- no replacing state!
// This is useful when you're trying to cheat and reuse derivatives from a prior advance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
State<Dimension>::
timeAdvanceOnly() const {
  return mTimeAdvanceOnly;
}

template<typename Dimension>
inline
void
State<Dimension>::
timeAdvanceOnly(const bool x) {
  mTimeAdvanceOnly = x;
}

}
