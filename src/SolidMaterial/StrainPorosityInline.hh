namespace Spheral {

//------------------------------------------------------------------------------
// Access various local state variables.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
StrainPorosity<Dimension>::
epsE() const {
  return mEpsE;
}

template<typename Dimension>
inline
double
StrainPorosity<Dimension>::
epsX() const {
  return mEpsX;
}

template<typename Dimension>
inline
double
StrainPorosity<Dimension>::
kappa() const {
  return mKappa;
}

template<typename Dimension>
inline
double
StrainPorosity<Dimension>::
gammaS0() const {
  return mGammaS0;
}

template<typename Dimension>
inline
double
StrainPorosity<Dimension>::
cS0() const {
  return mcS0;
}

template<typename Dimension>
inline
const PorousEquationOfState<Dimension>&
StrainPorosity<Dimension>::
porousEOS() const {
  return mPorousEOS;
}

template<typename Dimension>
inline
const PorousStrengthModel<Dimension>&
StrainPorosity<Dimension>::
porousStrength() const {
  return mPorousStrength;
}

template<typename Dimension>
inline
const NodeList<Dimension>&
StrainPorosity<Dimension>::
nodeList() const {
  return mNodeList;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Scalar>&
StrainPorosity<Dimension>::
c0() const {
  return mc0;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Scalar>&
StrainPorosity<Dimension>::
alpha0() const {
  return mAlpha0;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Scalar>&
StrainPorosity<Dimension>::
alpha() const {
  return mAlpha;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Scalar>&
StrainPorosity<Dimension>::
DalphaDt() const {
  return mDalphaDt;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Scalar>&
StrainPorosity<Dimension>::
strain() const {
  return mStrain;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Scalar>&
StrainPorosity<Dimension>::
DstrainDt() const {
  return mDstrainDt;
}

}
