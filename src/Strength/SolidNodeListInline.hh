namespace Spheral {
namespace SolidMaterial {

//------------------------------------------------------------------------------
// Access the deviatoric stress field.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
FieldSpace::Field<Dimension, typename Dimension::SymTensor>&
SolidNodeList<Dimension>::
deviatoricStress() {
  return mDeviatoricStress;
}

template<typename Dimension>
inline
const FieldSpace::Field<Dimension, typename Dimension::SymTensor>&
SolidNodeList<Dimension>::
deviatoricStress() const {
  return mDeviatoricStress;
}

//------------------------------------------------------------------------------
// Access the plastic strain field.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
FieldSpace::Field<Dimension, typename Dimension::Scalar>&
SolidNodeList<Dimension>::
plasticStrain() {
  return mPlasticStrain;
}

template<typename Dimension>
inline
const FieldSpace::Field<Dimension, typename Dimension::Scalar>&
SolidNodeList<Dimension>::
plasticStrain() const {
  return mPlasticStrain;
}

//------------------------------------------------------------------------------
// Access the plastic strain rate field.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
FieldSpace::Field<Dimension, typename Dimension::Scalar>&
SolidNodeList<Dimension>::
plasticStrainRate() {
  return mPlasticStrainRate;
}

template<typename Dimension>
inline
const FieldSpace::Field<Dimension, typename Dimension::Scalar>&
SolidNodeList<Dimension>::
plasticStrainRate() const {
  return mPlasticStrainRate;
}

//------------------------------------------------------------------------------
// Access the damage field.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
FieldSpace::Field<Dimension, typename Dimension::SymTensor>&
SolidNodeList<Dimension>::
damage() {
  return mDamage;
}

template<typename Dimension>
inline
const FieldSpace::Field<Dimension, typename Dimension::SymTensor>&
SolidNodeList<Dimension>::
damage() const {
  return mDamage;
}

//------------------------------------------------------------------------------
// Access the effective damage field.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
FieldSpace::Field<Dimension, typename Dimension::SymTensor>&
SolidNodeList<Dimension>::
effectiveDamage() {
  return mEffectiveDamage;
}

template<typename Dimension>
inline
const FieldSpace::Field<Dimension, typename Dimension::SymTensor>&
SolidNodeList<Dimension>::
effectiveDamage() const {
  return mEffectiveDamage;
}

//------------------------------------------------------------------------------
// Access the damage gradient field.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
FieldSpace::Field<Dimension, typename Dimension::Vector>&
SolidNodeList<Dimension>::
damageGradient() {
  return mDamageGradient;
}

template<typename Dimension>
inline
const FieldSpace::Field<Dimension, typename Dimension::Vector>&
SolidNodeList<Dimension>::
damageGradient() const {
  return mDamageGradient;
}

//------------------------------------------------------------------------------
// Access the strength model this solid node list is using.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const SolidMaterial::StrengthModel<Dimension>&
SolidNodeList<Dimension>::
strengthModel() const {
  return mStrength;
}

}
}
