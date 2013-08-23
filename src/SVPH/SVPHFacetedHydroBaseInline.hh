namespace Spheral {
namespace SVPHSpace {

//------------------------------------------------------------------------------
// Choose whether we want to sum for mass density, or integrate the continuity
// equation.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
PhysicsSpace::MassDensityType
SVPHFacetedHydroBase<Dimension>::densityUpdate() const {
  return mDensityUpdate;
}

template<typename Dimension>
inline
void
SVPHFacetedHydroBase<Dimension>::
densityUpdate(const PhysicsSpace::MassDensityType type) {
  mDensityUpdate = type;
}

//------------------------------------------------------------------------------
// Choose how we want to update the H tensor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
PhysicsSpace::HEvolutionType
SVPHFacetedHydroBase<Dimension>::HEvolution() const {
  return mHEvolution;
}

template<typename Dimension>
inline
void
SVPHFacetedHydroBase<Dimension>::
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
SVPHFacetedHydroBase<Dimension>::compatibleEnergyEvolution() const {
  return mCompatibleEnergyEvolution;
}

template<typename Dimension>
inline
void
SVPHFacetedHydroBase<Dimension>::compatibleEnergyEvolution(const bool val) {
  mCompatibleEnergyEvolution = val;
}

//------------------------------------------------------------------------------
// Access the flag determining if we're using the XSVPH algorithm.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
SVPHFacetedHydroBase<Dimension>::XSVPH() const {
  return mXSVPH;
}

template<typename Dimension>
inline
void
SVPHFacetedHydroBase<Dimension>::XSVPH(const bool val) {
  mXSVPH = val;
}

//------------------------------------------------------------------------------
// Access the flag determining if we're using the linear consistency corretions.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
SVPHFacetedHydroBase<Dimension>::linearConsistent() const {
  return mLinearConsistent;
}

template<typename Dimension>
inline
void
SVPHFacetedHydroBase<Dimension>::linearConsistent(const bool val) {
  mLinearConsistent = val;
}

//------------------------------------------------------------------------------
// Access the flag determining if we're generating void points.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
SVPHFacetedHydroBase<Dimension>::generateVoid() const {
  return mGenerateVoid;
}

template<typename Dimension>
inline
void
SVPHFacetedHydroBase<Dimension>::generateVoid(const bool val) {
  mGenerateVoid = val;
}

//------------------------------------------------------------------------------
// Access the optional min & max bounds for generating meshes.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename Dimension::Vector&
SVPHFacetedHydroBase<Dimension>::
xmin() const {
  return mXmin;
}

template<typename Dimension>
inline
const typename Dimension::Vector&
SVPHFacetedHydroBase<Dimension>::
xmax() const {
  return mXmax;
}

template<typename Dimension>
inline
void
SVPHFacetedHydroBase<Dimension>::
xmin(const typename Dimension::Vector& x) {
  mXmin = x;
}

template<typename Dimension>
inline
void
SVPHFacetedHydroBase<Dimension>::
xmax(const typename Dimension::Vector& x) {
  mXmax = x;
}

//------------------------------------------------------------------------------
// The object defining how smoothing scales are evolved.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const NodeSpace::SmoothingScaleBase<Dimension>&
SVPHFacetedHydroBase<Dimension>::
smoothingScaleMethod() const {
  return mSmoothingScaleMethod;
}

//------------------------------------------------------------------------------
// The mesh.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const MeshSpace::Mesh<Dimension>&
SVPHFacetedHydroBase<Dimension>::
mesh() const {
  return *mMeshPtr;
}

//------------------------------------------------------------------------------
// The internal state field lists.
//------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
// SVPHFacetedHydroBase<Dimension>::
// A() const {
//   return mA;
// }

// template<typename Dimension>
// inline
// const FieldSpace::FieldList<Dimension, typename Dimension::Vector>&
// SVPHFacetedHydroBase<Dimension>::
// B() const {
//   return mB;
// }

// template<typename Dimension>
// inline
// const FieldSpace::FieldList<Dimension, typename Dimension::Tensor>&
// SVPHFacetedHydroBase<Dimension>::
// gradB() const {
//   return mGradB;
// }

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, int>&
SVPHFacetedHydroBase<Dimension>::
timeStepMask() const {
  return mTimeStepMask;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SVPHFacetedHydroBase<Dimension>::
pressure() const {
  return mPressure;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SVPHFacetedHydroBase<Dimension>::
soundSpeed() const {
  return mSoundSpeed;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SVPHFacetedHydroBase<Dimension>::
volume() const {
  return mVolume;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SVPHFacetedHydroBase<Dimension>::
specificThermalEnergy0() const {
  return mSpecificThermalEnergy0;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>&
SVPHFacetedHydroBase<Dimension>::
Hideal() const {
  return mHideal;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SVPHFacetedHydroBase<Dimension>::
maxViscousPressure() const {
  return mMaxViscousPressure;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SVPHFacetedHydroBase<Dimension>::
massDensitySum() const {
  return mMassDensitySum;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SVPHFacetedHydroBase<Dimension>::
weightedNeighborSum() const {
  return mWeightedNeighborSum;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>&
SVPHFacetedHydroBase<Dimension>::
massSecondMoment() const {
  return mMassSecondMoment;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Vector>&
SVPHFacetedHydroBase<Dimension>::
XSVPHDeltaV() const {
  return mXSVPHDeltaV;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Vector>&
SVPHFacetedHydroBase<Dimension>::
DxDt() const {
  return mDxDt;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Vector>&
SVPHFacetedHydroBase<Dimension>::
DvDt() const {
  return mDvDt;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SVPHFacetedHydroBase<Dimension>::
DmassDensityDt() const {
  return mDmassDensityDt;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SVPHFacetedHydroBase<Dimension>::
DspecificThermalEnergyDt() const {
  return mDspecificThermalEnergyDt;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>&
SVPHFacetedHydroBase<Dimension>::
DHDt() const {
  return mDHDt;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Tensor>&
SVPHFacetedHydroBase<Dimension>::
DvDx() const {
  return mDvDx;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Tensor>&
SVPHFacetedHydroBase<Dimension>::
internalDvDx() const {
  return mInternalDvDx;
}

// template<typename Dimension>
// inline
// const FieldSpace::FieldList<Dimension, std::vector<typename Dimension::Scalar> >&
// SVPHFacetedHydroBase<Dimension>::
// faceMass() const {
//   return mFaceMass;
// }

// template<typename Dimension>
// inline
// const FieldSpace::FieldList<Dimension, std::vector<typename Dimension::Vector> >&
// SVPHFacetedHydroBase<Dimension>::
// faceVelocity() const {
//   return mFaceVelocity;
// }

// template<typename Dimension>
// inline
// const FieldSpace::FieldList<Dimension, std::vector<typename Dimension::Vector> >&
// SVPHFacetedHydroBase<Dimension>::
// faceAcceleration() const {
//   return mFaceAcceleration;
// }

// template<typename Dimension>
// inline
// const FieldSpace::FieldList<Dimension, std::vector<typename Dimension::Scalar> >&
// SVPHFacetedHydroBase<Dimension>::
// faceSpecificThermalEnergy0() const {
//   return mFaceSpecificThermalEnergy0;
// }

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, std::vector<typename Dimension::Vector> >&
SVPHFacetedHydroBase<Dimension>::
faceForce() const {
  return mFaceForce;
}

}
}
