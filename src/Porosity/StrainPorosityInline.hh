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
