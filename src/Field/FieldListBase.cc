//---------------------------------Spheral++----------------------------------//
// FieldListBase -- A base class to provide a generic handle on FieldLists.
//
// Created by JMO, Sun Apr 11 14:28:59 2004
//----------------------------------------------------------------------------//
#include "FieldListBase.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
FieldListBase::
FieldListBase():
  mNewCoarseNodes(true),
  mNewRefineNodes(true) {
}

//------------------------------------------------------------------------------
// Copy constructor.
//------------------------------------------------------------------------------
FieldListBase::
FieldListBase(const FieldListBase& fieldListBase):
  mNewCoarseNodes(true),
  mNewRefineNodes(true) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
FieldListBase::
~FieldListBase() {
}

//------------------------------------------------------------------------------
// Assignment operator.
//------------------------------------------------------------------------------
FieldListBase&
FieldListBase::
operator=(const FieldListBase& rhs) {
  if (this != &rhs) {
    mNewCoarseNodes = true;
    mNewRefineNodes = true;
  }
  return *this;
}

// //------------------------------------------------------------------------------
// // Trip the flag indicating there are new coarse nodes.
// //------------------------------------------------------------------------------
// void
// FieldListBase::
// notifyNewCoarseNodes() const {
//   mNewCoarseNodes = true;
// }

// //------------------------------------------------------------------------------
// // Trip the flag indicating there are new refine nodes.
// //------------------------------------------------------------------------------
// void
// FieldListBase::
// notifyNewRefineNodes() const {
//   mNewRefineNodes = true;
// }

}
