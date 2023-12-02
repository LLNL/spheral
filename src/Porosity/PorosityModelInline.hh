namespace Spheral {

//------------------------------------------------------------------------------
// Access various local state variables.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
PorosityModel<Dimension>::
jutziStateUpdate() const {
  return mJutziStateUpdate;
}

template<typename Dimension>
inline
double
PorosityModel<Dimension>::
rhoS0() const {
  return mRhoS0;
}

template<typename Dimension>
inline
double
PorosityModel<Dimension>::
cS0() const {
  return mcS0;
}

template<typename Dimension>
inline
double
PorosityModel<Dimension>::
KS0() const {
  return mKS0;
}

template<typename Dimension>
inline
double
PorosityModel<Dimension>::
fdt() const {
  return mfdt;
}

template<typename Dimension>
inline
double
PorosityModel<Dimension>::
maxAbsDalphaDt() const {
  return mMaxAbsDalphaDt;
}

template<typename Dimension>
inline
const SolidNodeList<Dimension>&
PorosityModel<Dimension>::
nodeList() const {
  return mNodeList;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Scalar>&
PorosityModel<Dimension>::
alpha0() const {
  return mAlpha0;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Scalar>&
PorosityModel<Dimension>::
alpha() const {
  return mAlpha;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Scalar>&
PorosityModel<Dimension>::
DalphaDt() const {
  return mDalphaDt;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Scalar>&
PorosityModel<Dimension>::
solidMassDensity() const {
  return mSolidMassDensity;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Scalar>&
PorosityModel<Dimension>::
c0() const {
  return mc0;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Scalar>&
PorosityModel<Dimension>::
fDS() const {
  return mfDS;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Scalar>&
PorosityModel<Dimension>::
fDSnew() const {
  return mfDSnew;
}

//------------------------------------------------------------------------------
// setters
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
PorosityModel<Dimension>::
fdt(const double x) {
  mfdt = x;
}

}
