namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
FieldUpdatePolicy<Dimension>::
FieldUpdatePolicy(std::initializer_list<std::string> depends):
  UpdatePolicyBase<Dimension>(depends) {
}

}
