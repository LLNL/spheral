#include "FieldBase.hh"

namespace Spheral {
namespace FieldSpace {

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
inline
FieldListBase::
FieldListBase():
  mNewCoarseNodes(true),
  mNewRefineNodes(true) {
}

//------------------------------------------------------------------------------
// Copy constructor.
//------------------------------------------------------------------------------
inline
FieldListBase::
FieldListBase(const FieldListBase& fieldListBase):
  mNewCoarseNodes(true),
  mNewRefineNodes(true) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
inline
FieldListBase::
~FieldListBase() {
}

//------------------------------------------------------------------------------
// Assignment operator.
//------------------------------------------------------------------------------
inline
FieldListBase&
FieldListBase::
operator=(const FieldListBase& rhs) {
  if (this != &rhs) {
    mNewCoarseNodes = true;
    mNewRefineNodes = true;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Register this FieldList with the given Field.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
FieldListBase::
registerWithField(const FieldBase<Dimension>& fieldBase) const {
  fieldBase.registerFieldList(*this);
}

//------------------------------------------------------------------------------
// Unregister this FieldList from the given Field.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
FieldListBase::
unregisterFromField(const FieldBase<Dimension>& fieldBase) const {
  fieldBase.unregisterFieldList(*this);
}

}
}
