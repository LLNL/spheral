//---------------------------------Spheral++----------------------------------//
// MassFluxPolicy -- This is basically a direct copy of the standard 
//                      position policy but instead we're substituting in 
//                      the nodal velocity as the derivative.
//
// J. M. Pearl 2023
//----------------------------------------------------------------------------//

#include "MassFluxPolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "NodeList/FluidNodeList.hh"
#include "DataBase/FieldListUpdatePolicyBase.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {


//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
MassFluxPolicy<Dimension>::
MassFluxPolicy():
  IncrementFieldList<Dimension, typename Dimension::Scalar>() {
}

template<typename Dimension>
MassFluxPolicy<Dimension>::
MassFluxPolicy(const std::string& depend0):
  IncrementFieldList<Dimension, typename Dimension::Scalar>(depend0) {
}

template<typename Dimension>
MassFluxPolicy<Dimension>::
MassFluxPolicy(const std::string& depend0,const std::string& depend1):
  IncrementFieldList<Dimension, typename Dimension::Scalar>(depend0,depend1) {
}

template<typename Dimension>
MassFluxPolicy<Dimension>::
MassFluxPolicy(const std::string& depend0,const std::string& depend1,const std::string& depend2):
  IncrementFieldList<Dimension, typename Dimension::Scalar>(depend0,depend1,depend2) {
}

template<typename Dimension>
MassFluxPolicy<Dimension>::
MassFluxPolicy(const std::string& depend0,const std::string& depend1,const std::string& depend2,const std::string& depend3):
  IncrementFieldList<Dimension, typename Dimension::Scalar>(depend0,depend1,depend2,depend3) {
}

template<typename Dimension>
MassFluxPolicy<Dimension>::
MassFluxPolicy(const std::string& depend0,const std::string& depend1,const std::string& depend2,const std::string& depend3,const std::string& depend4):
  IncrementFieldList<Dimension, typename Dimension::Scalar>(depend0,depend1,depend2,depend3,depend4) {
}

template<typename Dimension>
MassFluxPolicy<Dimension>::
MassFluxPolicy(const std::string& depend0,const std::string& depend1,const std::string& depend2,const std::string& depend3,const std::string& depend4,const std::string& depend5):
  IncrementFieldList<Dimension, typename Dimension::Scalar>(depend0,depend1,depend2,depend3,depend4,depend5) {
}
  
//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
MassFluxPolicy<Dimension>::
~MassFluxPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MassFluxPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double /*t*/,
       const double /*dt*/) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::mass and 
          nodeListKey == UpdatePolicyBase<Dimension>::wildcard());

  // state
  auto m = state.fields(fieldKey, 0.0);
  const auto numFields = m.numFields();

  // deriv
  const auto dmdt = derivs.fields(IncrementFieldList<Dimension, Scalar>::prefix() + fieldKey, 0.0);

  // Walk the fields.
  for (auto i = 0u; i != numFields; ++i) {
    const auto n = m[i]->numInternalElements();
    for (auto j = 0u; j < n; ++j) {
      m(i,j) += std::max(multiplier*dmdt(i,j),-m(i,j));
    }
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
MassFluxPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const MassFluxPolicy<Dimension>* rhsPtr = dynamic_cast<const MassFluxPolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

