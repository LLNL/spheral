#include <algorithm>

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
UpdatePolicyBase<Dimension>::
UpdatePolicyBase():
  mDependencies() {
}

template<typename Dimension>
inline
UpdatePolicyBase<Dimension>::
UpdatePolicyBase(const std::string& depend0):
  mDependencies() {
    mDependencies.push_back(depend0);
}

template<typename Dimension>
inline
UpdatePolicyBase<Dimension>::
UpdatePolicyBase(const std::string& depend0,
                 const std::string& depend1):
  mDependencies() {
    mDependencies.push_back(depend0);
    mDependencies.push_back(depend1);
    std::sort(mDependencies.begin(), mDependencies.end());
}

template<typename Dimension>
inline
UpdatePolicyBase<Dimension>::
UpdatePolicyBase(const std::string& depend0,
                 const std::string& depend1,
                 const std::string& depend2):
  mDependencies() {
    mDependencies.push_back(depend0);
    mDependencies.push_back(depend1);
    mDependencies.push_back(depend2);
    std::sort(mDependencies.begin(), mDependencies.end());
}

template<typename Dimension>
inline
UpdatePolicyBase<Dimension>::
UpdatePolicyBase(const std::string& depend0,
                 const std::string& depend1,
                 const std::string& depend2,
                 const std::string& depend3):
  mDependencies() {
    mDependencies.push_back(depend0);
    mDependencies.push_back(depend1);
    mDependencies.push_back(depend2);
    mDependencies.push_back(depend3);
    std::sort(mDependencies.begin(), mDependencies.end());
}

template<typename Dimension>
inline
UpdatePolicyBase<Dimension>::
UpdatePolicyBase(const std::string& depend0,
                 const std::string& depend1,
                 const std::string& depend2,
                 const std::string& depend3,
                 const std::string& depend4):
  mDependencies() {
    mDependencies.push_back(depend0);
    mDependencies.push_back(depend1);
    mDependencies.push_back(depend2);
    mDependencies.push_back(depend3);
    mDependencies.push_back(depend4);
    std::sort(mDependencies.begin(), mDependencies.end());
}

template<typename Dimension>
inline
UpdatePolicyBase<Dimension>::
UpdatePolicyBase(const std::string& depend0,
                 const std::string& depend1,
                 const std::string& depend2,
                 const std::string& depend3,
                 const std::string& depend4,
                 const std::string& depend5):
  mDependencies() {
    mDependencies.push_back(depend0);
    mDependencies.push_back(depend1);
    mDependencies.push_back(depend2);
    mDependencies.push_back(depend3);
    mDependencies.push_back(depend4);
    mDependencies.push_back(depend5);
    std::sort(mDependencies.begin(), mDependencies.end());
}

template<typename Dimension>
inline
UpdatePolicyBase<Dimension>::
UpdatePolicyBase(const std::string& depend0,
                 const std::string& depend1,
                 const std::string& depend2,
                 const std::string& depend3,
                 const std::string& depend4,
                 const std::string& depend5,
                 const std::string& depend6):
  mDependencies() {
    mDependencies.push_back(depend0);
    mDependencies.push_back(depend1);
    mDependencies.push_back(depend2);
    mDependencies.push_back(depend3);
    mDependencies.push_back(depend4);
    mDependencies.push_back(depend5);
    mDependencies.push_back(depend6);
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
// Is this state dependent?
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
UpdatePolicyBase<Dimension>::
independent() const {
  return mDependencies.size() == 0;
}

template<typename Dimension>
inline
bool
UpdatePolicyBase<Dimension>::
dependent() const {
  return mDependencies.size() != 0;
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

}
