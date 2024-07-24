//---------------------------------Spheral++----------------------------------//
// Convenience function to apply a policy to all the Fields in a FieldList
//
// Created by JMO, Thu Aug 26 14:28:07 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_applyPolicyToFieldList_hh__
#define __Spheral_applyPolicyToFieldList_hh__

#include "Field/FieldList.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"

namespace Spheral {

template<typename FieldListType>
inline
void applyPolicyToFieldList(FieldListType& fieldList,
                            std::shared_ptr<UpdatePolicyBase<typename FieldListType::FieldDimension>> policyPtr,
                            State<typename FieldListType::FieldDimension>& state,
                            StateDerivatives<typename FieldListType::FieldDimension>& derivs,
                            const double multiplier = 1.0,
                            const double t = 0.0,
                            const double dt = 0.0) {
  using Dimension = typename FieldListType::FieldDimension;
  if (policyPtr->clonePerField()) {
    for (auto* fieldPtr: fieldList) {
      const auto key = StateBase<Dimension>::key(*fieldPtr);
      policyPtr->update(key, state, derivs, multiplier, t, dt);
    }
  } else {
    if (fieldList.numFields() > 0) {
      const auto key = StateBase<Dimension>::buildFieldKey(fieldList[0]->name(), policyPtr->wildcard());
      policyPtr->update(key, state, derivs, multiplier, t, dt);
    }
  }
}

}

#endif

