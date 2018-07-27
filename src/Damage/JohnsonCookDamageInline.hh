namespace Spheral {
namespace PhysicsSpace {

//------------------------------------------------------------------------------
// Access the attributes
//------------------------------------------------------------------------------
template<typename Dimension>
const NodeSpace::SolidNodeList<Dimension>&
JohnsonCookDamageBase<Dimension>::
nodeList() const {
  return mNodeList;
}

template<typename Dimension>
const FieldSpace::Field<Dimension, typename Dimension::Scalar>&
JohnsonCookDamageBase<Dimension>::
failureStrain() const {
  return mFailureStrain;
}

template<typename Dimension>
const FieldSpace::Field<Dimension, typename Dimension::Scalar>&
JohnsonCookDamageBase<Dimension>::
meltSpecificEnergy() const {
  return mMeltSpecificEnergy;
}

template<typename Dimension>
const FieldSpace::Field<Dimension, typename Dimension::SymTensor>&
JohnsonCookDamageBase<Dimension>::
newEffectiveDamage() const {
  return mNewEffectiveDamage;
}

template<typename Dimension>
const FieldSpace::Field<Dimension, typename Dimension::Scalar>&
JohnsonCookDamageBase<Dimension>::
D1() const {
  return mD1;
}

template<typename Dimension>
const FieldSpace::Field<Dimension, typename Dimension::Scalar>&
JohnsonCookDamageBase<Dimension>::
D2() const {
  return mD2;
}

template<typename Dimension>
double
JohnsonCookDamageBase<Dimension>::
D3() const {
  return mD3;
}

template<typename Dimension>
double
JohnsonCookDamageBase<Dimension>::
D4() const {
  return mD4;
}

template<typename Dimension>
double
JohnsonCookDamageBase<Dimension>::
D5() const {
  return mD5;
}

template<typename Dimension>
double
JohnsonCookDamageBase<Dimension>::
epsilondot0() const {
  return mepsilondot0;
}

template<typename Dimension>
double
JohnsonCookDamageBase<Dimension>::
Tcrit() const {
  return mTcrit;
}

template<typename Dimension>
double
JohnsonCookDamageBase<Dimension>::
sigmamax() const {
  return msigmamax;
}

template<typename Dimension>
double
JohnsonCookDamageBase<Dimension>::
efailmin() const {
  return mefailmin;
}

}
}
