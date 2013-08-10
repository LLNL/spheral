namespace Spheral {
namespace SVPHSpace {

//------------------------------------------------------------------------------
// Choose whether we want to sum for mass density, or integrate the continuity
// equation.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
PhysicsSpace::MassDensityType
SVPHHydroBase<Dimension>::densityUpdate() const {
  return mDensityUpdate;
}

template<typename Dimension>
inline
void
SVPHHydroBase<Dimension>::
densityUpdate(const PhysicsSpace::MassDensityType type) {
  mDensityUpdate = type;
}

//------------------------------------------------------------------------------
// Choose how we want to update the H tensor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
PhysicsSpace::HEvolutionType
SVPHHydroBase<Dimension>::HEvolution() const {
  return mHEvolution;
}

template<typename Dimension>
inline
void
SVPHHydroBase<Dimension>::
HEvolution(const PhysicsSpace::HEvolutionType type) {
  mHEvolution = type;
}

//------------------------------------------------------------------------------
// Access the flag determining if we're using the compatible energy evolution 
// algorithm.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
SVPHHydroBase<Dimension>::compatibleEnergyEvolution() const {
  return mCompatibleEnergyEvolution;
}

template<typename Dimension>
inline
void
SVPHHydroBase<Dimension>::compatibleEnergyEvolution(const bool val) {
  mCompatibleEnergyEvolution = val;
}

//------------------------------------------------------------------------------
// Access the flag determining if we're using the XSVPH algorithm.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
SVPHHydroBase<Dimension>::XSVPH() const {
  return mXSVPH;
}

template<typename Dimension>
inline
void
SVPHHydroBase<Dimension>::XSVPH(const bool val) {
  mXSVPH = val;
}

//------------------------------------------------------------------------------
// Access the optional min & max bounds for generating meshes.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename Dimension::Vector&
SVPHHydroBase<Dimension>::
xmin() const {
  return mXmin;
}

template<typename Dimension>
inline
const typename Dimension::Vector&
SVPHHydroBase<Dimension>::
xmax() const {
  return mXmax;
}

template<typename Dimension>
inline
void
SVPHHydroBase<Dimension>::
xmin(const typename Dimension::Vector& x) {
  mXmin = x;
}

template<typename Dimension>
inline
void
SVPHHydroBase<Dimension>::
xmax(const typename Dimension::Vector& x) {
  mXmax = x;
}

//------------------------------------------------------------------------------
// The object defining how smoothing scales are evolved.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const NodeSpace::SmoothingScaleBase<Dimension>&
SVPHHydroBase<Dimension>::
smoothingScaleMethod() const {
  return mSmoothingScaleMethod;
}

//------------------------------------------------------------------------------
// The mesh.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const MeshSpace::Mesh<Dimension>&
SVPHHydroBase<Dimension>::
mesh() const {
  return mMesh;
}

//------------------------------------------------------------------------------
// The internal state field lists.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SVPHHydroBase<Dimension>::
A() const {
  return mA;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Vector>&
SVPHHydroBase<Dimension>::
B() const {
  return mB;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Tensor>&
SVPHHydroBase<Dimension>::
gradB() const {
  return mGradB;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, int>&
SVPHHydroBase<Dimension>::
timeStepMask() const {
  return mTimeStepMask;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SVPHHydroBase<Dimension>::
pressure() const {
  return mPressure;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SVPHHydroBase<Dimension>::
soundSpeed() const {
  return mSoundSpeed;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SVPHHydroBase<Dimension>::
volume() const {
  return mVolume;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SVPHHydroBase<Dimension>::
specificThermalEnergy0() const {
  return mSpecificThermalEnergy0;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>&
SVPHHydroBase<Dimension>::
Hideal() const {
  return mHideal;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SVPHHydroBase<Dimension>::
maxViscousPressure() const {
  return mMaxViscousPressure;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SVPHHydroBase<Dimension>::
massDensitySum() const {
  return mMassDensitySum;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SVPHHydroBase<Dimension>::
weightedNeighborSum() const {
  return mWeightedNeighborSum;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>&
SVPHHydroBase<Dimension>::
massSecondMoment() const {
  return mMassSecondMoment;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Vector>&
SVPHHydroBase<Dimension>::
XSVPHDeltaV() const {
  return mXSVPHDeltaV;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Vector>&
SVPHHydroBase<Dimension>::
DxDt() const {
  return mDxDt;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Vector>&
SVPHHydroBase<Dimension>::
DvDt() const {
  return mDvDt;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SVPHHydroBase<Dimension>::
DmassDensityDt() const {
  return mDmassDensityDt;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SVPHHydroBase<Dimension>::
DspecificThermalEnergyDt() const {
  return mDspecificThermalEnergyDt;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>&
SVPHHydroBase<Dimension>::
DHDt() const {
  return mDHDt;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Tensor>&
SVPHHydroBase<Dimension>::
DvDx() const {
  return mDvDx;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Tensor>&
SVPHHydroBase<Dimension>::
internalDvDx() const {
  return mInternalDvDx;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, std::vector<typename Dimension::Vector> >&
SVPHHydroBase<Dimension>::
pairAccelerations() const {
  return mPairAccelerations;
}

}
}
