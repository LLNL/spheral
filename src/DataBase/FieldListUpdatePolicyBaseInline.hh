namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldListUpdatePolicyBase<Dimension, DataType>::
FieldListUpdatePolicyBase():
  UpdatePolicyBase<Dimension>() {
}

template<typename Dimension, typename DataType>
inline
FieldListUpdatePolicyBase<Dimension, DataType>::
FieldListUpdatePolicyBase(const std::string& depend0):
  UpdatePolicyBase<Dimension>(depend0) {
}

template<typename Dimension, typename DataType>
inline
FieldListUpdatePolicyBase<Dimension, DataType>::
FieldListUpdatePolicyBase(const std::string& depend0,
                          const std::string& depend1):
  UpdatePolicyBase<Dimension>(depend0, depend1) {
}

template<typename Dimension, typename DataType>
inline
FieldListUpdatePolicyBase<Dimension, DataType>::
FieldListUpdatePolicyBase(const std::string& depend0,
                          const std::string& depend1,
                          const std::string& depend2):
  UpdatePolicyBase<Dimension>(depend0, depend1, depend2) {
}

template<typename Dimension, typename DataType>
inline
FieldListUpdatePolicyBase<Dimension, DataType>::
FieldListUpdatePolicyBase(const std::string& depend0,
                          const std::string& depend1,
                          const std::string& depend2,
                          const std::string& depend3):
  UpdatePolicyBase<Dimension>(depend0, depend1, depend2, depend3) {
}

template<typename Dimension, typename DataType>
inline
FieldListUpdatePolicyBase<Dimension, DataType>::
FieldListUpdatePolicyBase(const std::string& depend0,
                          const std::string& depend1,
                          const std::string& depend2,
                          const std::string& depend3,
                          const std::string& depend4):
  UpdatePolicyBase<Dimension>(depend0, depend1, depend2, depend3, depend4) {
}

template<typename Dimension, typename DataType>
inline
FieldListUpdatePolicyBase<Dimension, DataType>::
FieldListUpdatePolicyBase(const std::string& depend0,
                          const std::string& depend1,
                          const std::string& depend2,
                          const std::string& depend3,
                          const std::string& depend4,
                          const std::string& depend5):
  UpdatePolicyBase<Dimension>(depend0, depend1, depend2, depend3, depend4, depend5) {
}

template<typename Dimension, typename DataType>
inline
FieldListUpdatePolicyBase<Dimension, DataType>::
FieldListUpdatePolicyBase(const std::string& depend0,
                          const std::string& depend1,
                          const std::string& depend2,
                          const std::string& depend3,
                          const std::string& depend4,
                          const std::string& depend5,
                          const std::string& depend6):
  UpdatePolicyBase<Dimension>(depend0, depend1, depend2, depend3, depend4, depend5, depend6) {
}

//------------------------------------------------------------------------------
// Overload the methods describing how to update FieldLists.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldListUpdatePolicyBase<Dimension, DataType>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {
  
  // Get the field name portion of the key.
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  CHECK(nodeListKey == UpdatePolicyBase<Dimension>::wildcard());

  // Now get the Fields that match from the state
  FieldList<Dimension, ValueType> fieldList = state.fields(fieldKey, ValueType());

  // Walk each Field, and use the corresponding individual NodeList policies to update it.
  for (auto* fieldPtr: fieldList) {
    auto policy = this->policyForNodeList(fieldPtr->nodeList().name());
    auto fieldKey = StateBase<Dimension>::key(*fieldPtr);
    policy->update(fieldKey, state, derivs, multiplier, t, dt);
  }
}

//------------------------------------------------------------------------------
// An alternate method to be called when you want to specify that the derivative information
// should be assumed to not necessarily be properly time-centered, and therefore you should 
// only use time advancement ideas, no "replace" or more sophisticated approaches.
// Default to just calling the generic method.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldListUpdatePolicyBase<Dimension, DataType>::
updateAsIncrement(const KeyType& key,
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs,
                  const double multiplier,
                  const double t,
                  const double dt) {

  // Get the field name portion of the key.
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  CHECK(nodeListKey == UpdatePolicyBase<Dimension>::wildcard());

  // Now get the Fields that match from the state
  FieldList<Dimension, ValueType> fieldList = state.fields(fieldKey, ValueType());

  // Walk each Field, and use the corresponding individual NodeList policies to update it.
  for (auto* fieldPtr: fieldList) {
    auto policy = this->policyForNodeList(fieldPtr->nodeList().name());
    auto fieldKey = StateBase<Dimension>::key(*fieldPtr);
    policy->updateAsIncrement(fieldKey, state, derivs, multiplier, t, dt);
  }
}

//------------------------------------------------------------------------------
// Enroll a policy specialized for a given NodeList
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldListUpdatePolicyBase<Dimension, DataType>::
enroll(const std::string& nodeListName, std::shared_ptr<UpdatePolicyBase<Dimension>> policy) {
  mNodeListPolicies[nodeListName] = policy;
}
  
//------------------------------------------------------------------------------
// Check if there is a specialized policy for a NodeList
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldListUpdatePolicyBase<Dimension, DataType>::
haveSpecializedPolicy(const std::string& nodeListName) const {
  return mNodeListPolicies.find(nodeListName) == nodeListPolicies.end();
}
  
//------------------------------------------------------------------------------
// Return the specialized policy for a NodeList
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
std::shared_ptr<UpdatePolicyBase<Dimension>>
FieldListUpdatePolicyBase<Dimension, DataType>::
policyForNodeList(const std::string& nodeListName) const {
  auto result = mNodeListPolicies.find(nodeListName);
  CHECK(result != nodeListPolicies.end());
  return result;
}

}
