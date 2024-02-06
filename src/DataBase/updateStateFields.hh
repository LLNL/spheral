//---------------------------------Spheral++----------------------------------//
// Convenience function to use policies to update a FieldList
//
// Created by JMO, Wed Jan 31 13:52:31 PST 2024
//----------------------------------------------------------------------------//
#ifndef __Spheral_updateStateFields_hh__
#define __Spheral_updateStateFields_hh__

#include "Field/FieldBase.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"

namespace Spheral {

template<typename Dimension>
inline
void updateStateFields(const typename FieldBase<Dimension>::FieldName& fieldName,
                       State<Dimension>& state,
                       StateDerivatives<Dimension>& derivs,
                       const double multiplier = 1.0,
                       const double t = 0.0,
                       const double dt = 0.0) {
  const auto policies = state.policies(fieldName);
  CHECK(not policies.empty());
  for (auto& [key, policy]: policies) {
    policy->update(key, state, derivs, multiplier, t, dt);
  }
}

}

#endif

