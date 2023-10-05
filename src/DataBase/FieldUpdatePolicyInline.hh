namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
FieldUpdatePolicy<Dimension>::
FieldUpdatePolicy():
  UpdatePolicyBase<Dimension>() {
}

template<typename Dimension>
inline
FieldUpdatePolicy<Dimension>::
FieldUpdatePolicy(const std::string& depend0):
  UpdatePolicyBase<Dimension>(depend0) {
}

template<typename Dimension>
inline
FieldUpdatePolicy<Dimension>::
FieldUpdatePolicy(const std::string& depend0,
                  const std::string& depend1):
  UpdatePolicyBase<Dimension>(depend0, depend1) {
}

template<typename Dimension>
inline
FieldUpdatePolicy<Dimension>::
FieldUpdatePolicy(const std::string& depend0,
                  const std::string& depend1,
                  const std::string& depend2):
  UpdatePolicyBase<Dimension>(depend0, depend1, depend2) {
}

template<typename Dimension>
inline
FieldUpdatePolicy<Dimension>::
FieldUpdatePolicy(const std::string& depend0,
                  const std::string& depend1,
                  const std::string& depend2,
                  const std::string& depend3):
  UpdatePolicyBase<Dimension>(depend0, depend1, depend2, depend3) {
}

template<typename Dimension>
inline
FieldUpdatePolicy<Dimension>::
FieldUpdatePolicy(const std::string& depend0,
                  const std::string& depend1,
                  const std::string& depend2,
                  const std::string& depend3,
                  const std::string& depend4):
  UpdatePolicyBase<Dimension>(depend0, depend1, depend2, depend3, depend4) {
}

template<typename Dimension>
inline
FieldUpdatePolicy<Dimension>::
FieldUpdatePolicy(const std::string& depend0,
                  const std::string& depend1,
                  const std::string& depend2,
                  const std::string& depend3,
                  const std::string& depend4,
                  const std::string& depend5):
  UpdatePolicyBase<Dimension>(depend0, depend1, depend2, depend3, depend4, depend5) {
}

template<typename Dimension>
inline
FieldUpdatePolicy<Dimension>::
FieldUpdatePolicy(const std::string& depend0,
                  const std::string& depend1,
                  const std::string& depend2,
                  const std::string& depend3,
                  const std::string& depend4,
                  const std::string& depend5,
                  const std::string& depend6):
  UpdatePolicyBase<Dimension>(depend0, depend1, depend2, depend3, depend4, depend5, depend6) {
}

}
