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
  mPolicyMap[key] = polptr;
}

//------------------------------------------------------------------------------
// Enroll the given field and policy.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename Value>
inline
void
State<Dimension>::
enroll(FieldSpace::Field<Dimension, Value>& field,
       typename State<Dimension>::PolicyPointer polptr) {
  this->enroll(field);
  this->enroll(this->key(field), polptr);
}

//------------------------------------------------------------------------------
// Enroll the given field.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
State<Dimension>::
enroll(FieldSpace::FieldBase<Dimension>& field) {
  StateBase<Dimension>::enroll(field);
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
