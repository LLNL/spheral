namespace Spheral {

//------------------------------------------------------------------------------
// Access various local state variables.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
PalphaPorosity<Dimension>::
Pe() const {
  return mPe;
}

template<typename Dimension>
inline
double
PalphaPorosity<Dimension>::
Pt() const {
  return mPt;
}

template<typename Dimension>
inline
double
PalphaPorosity<Dimension>::
Ps() const {
  return mPs;
}

template<typename Dimension>
inline
double
PalphaPorosity<Dimension>::
alphae() const {
  return mAlphae;
}

template<typename Dimension>
inline
double
PalphaPorosity<Dimension>::
alphat() const {
  return mAlphat;
}

template<typename Dimension>
inline
double
PalphaPorosity<Dimension>::
n1() const {
  return mn1;
}

template<typename Dimension>
inline
double
PalphaPorosity<Dimension>::
n2() const {
  return mn2;
}

template<typename Dimension>
inline
double
PalphaPorosity<Dimension>::
cS0() const {
  return mcS0;
}

template<typename Dimension>
inline
const PorousEquationOfState<Dimension>&
PalphaPorosity<Dimension>::
porousEOS() const {
  return mPorousEOS;
}

template<typename Dimension>
inline
const PorousStrengthModel<Dimension>&
PalphaPorosity<Dimension>::
porousStrength() const {
  return mPorousStrength;
}

template<typename Dimension>
inline
const NodeList<Dimension>&
PalphaPorosity<Dimension>::
nodeList() const {
  return mNodeList;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Scalar>&
PalphaPorosity<Dimension>::
c0() const {
  return mc0;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Scalar>&
PalphaPorosity<Dimension>::
alpha0() const {
  return mAlpha0;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Scalar>&
PalphaPorosity<Dimension>::
alpha() const {
  return mAlpha;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Scalar>&
PalphaPorosity<Dimension>::
DalphaDt() const {
  return mDalphaDt;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Scalar>&
PalphaPorosity<Dimension>::
partialPpartialEps() const {
  return mdPdU;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Scalar>&
PalphaPorosity<Dimension>::
partialPpartialRho() const {
  return mdPdR;
}

}
