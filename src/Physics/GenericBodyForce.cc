//---------------------------------Spheral++----------------------------------//
// GenericBodyForce -- The base class for all Spheral++ body force
// implementations.
//
// Created by JMO, Wed May 24 14:23:10 PDT 2000
//----------------------------------------------------------------------------//
#include <string>

#include "GenericBodyForce.hh"
#include "DataBase/IncrementFieldList.hh"
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
  mDxDt(FieldSpace::Copy),
  mDvDt(FieldSpace::Copy) {
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
  FieldList<Dimension, Scalar> mass = dataBase.globalMass();
  FieldList<Dimension, Vector> position = dataBase.globalPosition();
  FieldList<Dimension, Vector> velocity = dataBase.globalVelocity();
  FieldList<Dimension, SymTensor> Hfield = dataBase.globalHfield();

  PolicyPointer positionPolicy(new IncrementFieldList<Dimension, Vector>());
  PolicyPointer velocityPolicy(new IncrementFieldList<Dimension, Vector>());
  state.enroll(mass);
  state.enroll(position, positionPolicy);
  state.enroll(velocity, velocityPolicy);
  state.enroll(Hfield);
}

//------------------------------------------------------------------------------
// Create and register an acceleration source.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GenericBodyForce<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
  // These derivative fields *may* conflict with other physics packages.
  // Since we create local storage for derivative fields we have to be
  // careful here, and only create storage for those NodeLists that are
  // not already registered by someone else.  We also deliberately do not
  // zero out the fields at this stage!
  dataBase.resizeGlobalFieldList(mDxDt, Vector::zero, IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, false);
  dataBase.resizeGlobalFieldList(mDvDt, Vector::zero, IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::velocity, false);
  if (not derivs.registered(mDxDt)) derivs.enroll(mDxDt);
  if (not derivs.registered(mDvDt)) derivs.enroll(mDvDt);
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

