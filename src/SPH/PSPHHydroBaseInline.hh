namespace Spheral {
namespace SPHSpace {

//------------------------------------------------------------------------------
// Access the flag determining if we're evolving total or specific energy
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
PSPHHydroBase<Dimension>::evolveTotalEnergy() const {
  return mEvolveTotalEnergy;
}

template<typename Dimension>
inline
void
PSPHHydroBase<Dimension>::evolveTotalEnergy(const bool val) {
  mEvolveTotalEnergy = val;
}

//------------------------------------------------------------------------------
// Access the flag determining if we're applying Hopkins 2014 artificial
// conductivity.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
PSPHHydroBase<Dimension>::HopkinsConductivity() const {
  return mHopkinsConductivity;
}

template<typename Dimension>
inline
void
PSPHHydroBase<Dimension>::HopkinsConductivity(const bool val) {
  mHopkinsConductivity = val;
}

//------------------------------------------------------------------------------
// The internal state field lists.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
PSPHHydroBase<Dimension>::
energy() const {
  return mEnergy;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
PSPHHydroBase<Dimension>::
gamma() const {
  return mGamma;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
PSPHHydroBase<Dimension>::
PSPHcorrection() const {
  return mPSPHcorrection;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
PSPHHydroBase<Dimension>::
DenergyDt() const {
  return mDenergyDt;
}

}
}
