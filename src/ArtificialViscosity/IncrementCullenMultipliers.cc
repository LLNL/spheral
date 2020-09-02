//---------------------------------Spheral++----------------------------------//
// IncrementCullenMultipliers
//
// Created by JMO, Mon Dec 28 22:18:44 PST 2015
//----------------------------------------------------------------------------//
#include "IncrementCullenMultipliers.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"
#include "Hydro/HydroFieldNames.hh"

using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension>
IncrementCullenMultipliers<Dimension>::
IncrementCullenMultipliers(const typename Dimension::Scalar minValue,
                           const typename Dimension::Scalar maxValue,
                           const bool hopkinsForm):
  IncrementBoundedFieldList<Dimension, typename Dimension::Scalar, typename Dimension::Scalar>(minValue, 
                                                                                               maxValue),
  mHopkinsForm(hopkinsForm) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
IncrementCullenMultipliers<Dimension>::
~IncrementCullenMultipliers() {
}

//------------------------------------------------------------------------------
// Update the Cullen Q multiplier factors.
//------------------------------------------------------------------------------
template<typename Dimension>
void
IncrementCullenMultipliers<Dimension>::
update(const KeyType&,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double /*t*/,
       const double /*dt*/) {

  // Get the state we're advancing and the needed derivatives.
  FieldList<Dimension, Scalar> rvQ = state.fields(HydroFieldNames::ArtificialViscousCqMultiplier, 0.0);
  FieldList<Dimension, Scalar> rvL = state.fields(HydroFieldNames::ArtificialViscousClMultiplier, 0.0);
  FieldList<Dimension, Scalar> alpha0 = state.fields("mCullAlpha", 0.0);
  const FieldList<Dimension, Scalar> alpha_local = derivs.fields("Cullen alpha local", 0.0);
  const FieldList<Dimension, Scalar> DalphaDt = derivs.fields("Cullen alpha delta", 0.0);
  const FieldList<Dimension, Scalar> alpha_tmp = derivs.fields("mCullAlpha2", 0.0);

  const unsigned numNodeLists = rvQ.size();
  CHECK(rvL.size() == numNodeLists);
  CHECK(alpha0.size() == numNodeLists);
  CHECK(alpha_local.size() == numNodeLists);
  CHECK(DalphaDt.size() == numNodeLists);

  // Set the new values.
  const Scalar alphaMin = this->minValue();
  const Scalar alphaMax = this->maxValue();
  for (unsigned k = 0; k != numNodeLists; ++k) {
    const unsigned n = rvQ[k]->numInternalElements();
    for (unsigned i = 0; i < n; ++i) {
      if (mHopkinsForm) {
        // Hopkins 2014
        const Scalar alphai = std::max(alphaMin, std::min(alphaMax, alpha_local(k, i)));
        rvQ(k, i) = alphai;
        rvL(k, i) = alphai;
        alpha0(k, i) = DalphaDt(k, i);
      } else {
        // Cullen & Dehnen 2010
        const Scalar alphai = std::max(alphaMin, std::min(alphaMax, max(alpha_local(k, i), rvQ(k, i)) + multiplier*DalphaDt(k, i)));
        rvQ(k, i) = alphai;
        rvL(k, i) = alphai;
      }
    }
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
IncrementCullenMultipliers<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const IncrementCullenMultipliers<Dimension>* rhsPtr = dynamic_cast<const IncrementCullenMultipliers<Dimension>*>(&rhs);
  if (rhsPtr == 0) return false;

  // Ok, now do we agree on min & max?
  return (this->minValue() == rhsPtr->minValue()) && (this->maxValue() == rhsPtr->maxValue());
}

}

