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
