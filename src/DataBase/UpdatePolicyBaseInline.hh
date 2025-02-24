#include "Utilities/DBC.hh"

#include <algorithm>

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
UpdatePolicyBase<Dimension>::
UpdatePolicyBase(std::initializer_list<std::string> depends):
  mDependencies(depends) {
  std::sort(mDependencies.begin(), mDependencies.end());
}

//------------------------------------------------------------------------------
// The != operator.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
UpdatePolicyBase<Dimension>::
operator!=(const UpdatePolicyBase& rhs) const {
  return !(*this == rhs);
}

//------------------------------------------------------------------------------
// The set of field names this state is dependent on.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector<std::string>&
UpdatePolicyBase<Dimension>::
dependencies() const {
  return mDependencies;
}

//------------------------------------------------------------------------------
// Allow the addition of dependencies.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
UpdatePolicyBase<Dimension>::
addDependency(const std::string& depend) {
  mDependencies.push_back(depend);
  std::sort(mDependencies.begin(), mDependencies.end());
}

//------------------------------------------------------------------------------
// Is this state dependent?
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
UpdatePolicyBase<Dimension>::
independent() const {
  return mDependencies.size() == 0;
}

//------------------------------------------------------------------------------
// Serialize the underlying data.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
UpdatePolicyBase<Dimension>::
serializeData(const KeyType& key,
              const State<Dimension>& state,
              std::vector<double>& buf) const {
  VERIFY2(false, "UpdatePolicyBase ERROR: attempt to call base serialize method on " + key);
}

//------------------------------------------------------------------------------
// Deserialize the underlying data.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
size_t
UpdatePolicyBase<Dimension>::
deserializeData(const KeyType& key,
                const State<Dimension>& state,
                const std::vector<double>& buf,
                const size_t offset) const {
  VERIFY2(false, "UpdatePolicyBase ERROR: attempt to call base deserialize method on " + key);
  return offset;
}

}
