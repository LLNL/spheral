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
  VERIFY2(mJutziStateUpdate and mfDSptr, "PorosityModel: fDS is not available since jutziStateUpdate is not set");
  return *mfDSptr;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Scalar>&
PorosityModel<Dimension>::
fDS_new() const {
  VERIFY2(mJutziStateUpdate and mfDSnewPtr, "PorosityModel: fDS_new is not available since jutziStateUpdate is not set");
  return *mfDSnewPtr;
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
