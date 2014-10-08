namespace Spheral {
namespace SolidMaterial {

//------------------------------------------------------------------------------
// Access various local state variables.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
StrainPorosity<Dimension>::
phi0() const {
  return 1.0 - 1.0/mAlpha0;
}

template<typename Dimension>
inline
double
StrainPorosity<Dimension>::
alpha0() const {
  return mAlpha0;
}

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
epsC() const {
  return mEpsC;
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
const PorousEquationOfState<Dimension>&
StrainPorosity<Dimension>::
porousEOS() const {
  return mPorousEOS;
}

template<typename Dimension>
inline
const NodeSpace::NodeList<Dimension>&
StrainPorosity<Dimension>::
nodeList() const {
  return mNodeList;
}

template<typename Dimension>
inline
const FieldSpace::Field<Dimension, typename Dimension::Scalar>&
StrainPorosity<Dimension>::
alpha() const {
  return mAlpha;
}

template<typename Dimension>
inline
const FieldSpace::Field<Dimension, typename Dimension::Scalar>&
StrainPorosity<Dimension>::
DalphaDt() const {
  return mDalphaDt;
}

template<typename Dimension>
inline
const FieldSpace::Field<Dimension, typename Dimension::Scalar>&
StrainPorosity<Dimension>::
strain() const {
  return mStrain;
}

template<typename Dimension>
inline
const FieldSpace::Field<Dimension, typename Dimension::Scalar>&
StrainPorosity<Dimension>::
DstrainDt() const {
  return mDstrainDt;
}

}
}
