//---------------------------------Spheral++----------------------------------//
// IncrementASPHHtensor
//
// Specialized version of FieldUpdatePolicy for time integrating the H tensor.
//
// Created by JMO, Mon Oct  7 13:31:02 PDT 2024
//----------------------------------------------------------------------------//
#include "IncrementASPHHtensor.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Utilities/rotationMatrix.hh"
#include "Utilities/GeometricUtilities.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
IncrementASPHHtensor<Dimension>::
IncrementASPHHtensor(const Scalar hmin,
                     const Scalar hmax,
                     const Scalar hminratio,
                     const bool fixShape,
                     const bool radialOnly):
  FieldUpdatePolicy<Dimension>(),
  mhmin(hmin),
  mhmax(hmax),
  mhminratio(hminratio),
  mFixShape(fixShape),
  mRadialOnly(radialOnly) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
IncrementASPHHtensor<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  // Get the field name portion of the key.
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  CHECK(fieldKey == HydroFieldNames::H);

  const auto hminInv = 1.0/mhmin;
  const auto hmaxInv = 1.0/mhmax;

  // Get the state we're updating.
  auto&       H = state.field(key, SymTensor::zero);
  const auto& DHDt = derivs.field(prefix() + StateBase<Dimension>::buildFieldKey(HydroFieldNames::H, nodeListKey), SymTensor::zero);
  const auto& pos = state.field(StateBase<Dimension>::buildFieldKey(HydroFieldNames::position, nodeListKey), Vector::zero);  // Only needed if we're using radial scaling

  // Walk the nodes and update H (with limiting)
  const auto n = H.numInternalElements();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    
    // Check for special update rules
    if (mFixShape) {

      // Fix the shape (only volume scaling allowed)

      auto fi = Dimension::rootnu((H(i) + multiplier*DHDt(i)).Determinant()/H(i).Determinant());
      H(i) *= fi;

    } else if (mRadialOnly) {

      // Force only the radial component of H to be scaled
      const auto nhat = pos(i).unitVector();
      const auto T = rotationMatrix(nhat);
      H(i).rotationalTransform(T);          // Should have one eigenvector aligned with the x' axis in this frame
      auto DHDti = DHDt(i);
      DHDti.rotationalTransform(T);
      H(i)[0] += multiplier * DHDti[0];
      H(i).rotationalTransform(T.Transpose());

    } else {

      H(i) += multiplier * DHDt(i);

    }

    // Apply limiting
    const auto hev = H(i).eigenVectors();
    const auto hminEffInv = min(hminInv, max(hmaxInv, hev.eigenValues.minElement())/mhminratio);
    H(i) = constructSymTensorWithBoundedDiagonal(hev.eigenValues, hmaxInv, hminEffInv);
    H(i).rotationalTransform(hev.eigenVectors);
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
IncrementASPHHtensor<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {
  const auto rhsPtr = dynamic_cast<const IncrementASPHHtensor<Dimension>*>(&rhs);
  if (rhsPtr == nullptr) return false;

  // Ok, now do we agree on min & max?
  return (hmin() == rhsPtr->hmin() and
          hmax() == rhsPtr->hmax() and
          hminratio() == rhsPtr->hminratio() and
          fixShape() == rhsPtr->fixShape() and
          radialOnly() == rhsPtr->radialOnly());
}

}
