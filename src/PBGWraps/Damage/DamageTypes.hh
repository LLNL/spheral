#ifndef __PBGWRAPS_DAMAGETYPES__
#define __PBGWRAPS_DAMAGETYPES__

#include "Geometry/Dimension.hh"
#include "NodeList/SolidNodeList.hh"
#include "Strength/SolidFieldNames.hh"
#include "Damage/DamageModel.hh"
#include "Damage/TensorDamageModel.hh"
#include "Damage/JohnsonCookDamage.hh"
#include "Damage/weibullFlawDistribution.hh"
#include "Damage/computeFragmentField.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Names!
//------------------------------------------------------------------------------
typedef DamageModel<Dim<1> > DamageModel1d;
typedef DamageModel<Dim<2> > DamageModel2d;
typedef DamageModel<Dim<3> > DamageModel3d;

typedef TensorDamageModel<Dim<1> > TensorDamageModel1d;
typedef TensorDamageModel<Dim<2> > TensorDamageModel2d;
typedef TensorDamageModel<Dim<3> > TensorDamageModel3d;

typedef JohnsonCookDamage<Dim<1> > JohnsonCookDamage1d;
typedef JohnsonCookDamage<Dim<2> > JohnsonCookDamage2d;
typedef JohnsonCookDamage<Dim<3> > JohnsonCookDamage3d;

//------------------------------------------------------------------------------
// Extract fields from DamageModel.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
Field<Dimension, typename Dimension::Scalar>*
youngsModulusFromDamageModel(const DamageModel<Dimension>& self) {
  return &const_cast<Field<Dimension, typename Dimension::Scalar>&>(self.youngsModulus());
}

template<typename Dimension>
inline
Field<Dimension, typename Dimension::Scalar>*
longitudinalSoundSpeedFromDamageModel(const DamageModel<Dimension>& self) {
  return &const_cast<Field<Dimension, typename Dimension::Scalar>&>(self.longitudinalSoundSpeed());
}

}

#endif
