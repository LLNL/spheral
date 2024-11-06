#include "FieldBase.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
FieldListBase<Dimension>::
FieldListBase():
  mNewCoarseNodes(true),
  mNewRefineNodes(true) {
}

//------------------------------------------------------------------------------
// Copy constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
FieldListBase<Dimension>::
FieldListBase(const FieldListBase<Dimension>& /*fieldListBase*/):
  mNewCoarseNodes(true),
  mNewRefineNodes(true) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
FieldListBase<Dimension>::
~FieldListBase() {
}

//------------------------------------------------------------------------------
// Assignment operator.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
FieldListBase<Dimension>&
FieldListBase<Dimension>::
operator=(const FieldListBase<Dimension>& rhs) {
  if (this != &rhs) {
    mNewCoarseNodes = true;
    mNewRefineNodes = true;
  }
  return *this;
}

////------------------------------------------------------------------------------
//// Register this FieldList with the given Field.
////------------------------------------------------------------------------------
//template<typename Dimension>
//inline
//void
//FieldListBase<Dimension>::
//registerWithField(const FieldBaseView<Dimension>& fieldBase) const {
//  fieldBase.registerFieldList(*this);
//}
//
////------------------------------------------------------------------------------
//// Unregister this FieldList from the given Field.
////------------------------------------------------------------------------------
//template<typename Dimension>
//inline
//void
//FieldListBase<Dimension>::
//unregisterFromField(const FieldBaseView<Dimension>& fieldBase) const {
//  fieldBase.unregisterFieldList(*this);
//}

}
