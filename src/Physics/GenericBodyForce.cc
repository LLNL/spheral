//---------------------------------Spheral++----------------------------------//
// GenericBodyForce -- The base class for all Spheral++ body force
// implementations.
//
// Created by JMO, Wed May 24 14:23:10 PDT 2000
//----------------------------------------------------------------------------//
#include <string>

#include "GenericBodyForce.hh"
#include "DataBase/IncrementState.hh"
#include "Field/Field.hh"
#include "Hydro/HydroFieldNames.hh"

namespace Spheral {
namespace PhysicsSpace {

using namespace std;
using DataBaseSpace::DataBase;
using FieldSpace::Field;

//------------------------------------------------------------------------------
// Default constructor
//------------------------------------------------------------------------------
template<typename Dimension>
GenericBodyForce<Dimension>::
GenericBodyForce():
  Physics<Dimension>(),
  mDxDt(FieldList<Dimension, Vector>::Copy),
  mDvDt(FieldList<Dimension, Vector>::Copy) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
GenericBodyForce<Dimension>::~GenericBodyForce() {
}

//------------------------------------------------------------------------------
// Create and register an acceleration source.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GenericBodyForce<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  typedef typename State<Dimension>::PolicyPointer PolicyPointer;

  // These state fields may be registered by other physics as well, but
  // since they are shared it is harmless.
  for (typename DataBase<Dimension>::NodeListIterator itr = dataBase.nodeListBegin();
       itr != dataBase.nodeListEnd();
       ++itr) {
    PolicyPointer positionPolicy(new IncrementState<Dimension, Vector>);
    PolicyPointer velocityPolicy(new IncrementState<Dimension, Vector>);
    state.enroll((*itr)->mass());
    state.enroll((*itr)->positions(), positionPolicy);
    state.enroll((*itr)->velocity(), velocityPolicy);
    state.enroll((*itr)->Hfield());
  }
}

//------------------------------------------------------------------------------
// Create and register an acceleration source.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GenericBodyForce<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
  typedef typename StateDerivatives<Dimension>::KeyType Key;
  const Key DxDtName = IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::position;
  const Key DvDtName = IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::velocity;

  // These derivative fields *may* conflict with other physics packages.
  // Since we create local storage for derivative fields we have to be
  // careful here, and only create storage for those NodeLists that are
  // not already registered by someone else.  We also deliberately do not
  // zero out the fields at this stage!
  dataBase.resizeGlobalFieldList(mDxDt, Vector::zero, IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::position, false);
  dataBase.resizeGlobalFieldList(mDvDt, Vector::zero, IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::velocity, false);
  size_t nodeListi = 0;
  for (typename DataBase<Dimension>::NodeListIterator itr = dataBase.nodeListBegin();
       itr != dataBase.nodeListEnd();
       ++itr, ++nodeListi) {
    const Key DxDtKey = StateBase<Dimension>::key(*mDxDt[nodeListi]);
    const Key DvDtKey = StateBase<Dimension>::key(*mDvDt[nodeListi]);
    if (not derivs.registered(DxDtKey)) derivs.enroll(*mDxDt[nodeListi]);
    if (not derivs.registered(DvDtKey)) derivs.enroll(*mDvDt[nodeListi]);
  }
}

//------------------------------------------------------------------------------
// Our derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
const FieldSpace::FieldList<Dimension, typename Dimension::Vector>&
GenericBodyForce<Dimension>::
DxDt() const {
  return mDxDt;
}

template<typename Dimension>
const FieldSpace::FieldList<Dimension, typename Dimension::Vector>&
GenericBodyForce<Dimension>::
DvDt() const {
  return mDvDt;
}

}
}

