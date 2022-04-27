namespace Spheral{
//------------------------------------------------------------------------------
// Return the surface normal field list ref
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Vector>&
SlideSurface<Dimension>::
surfaceNormals() const {
  return mSurfaceNormals;
}

//------------------------------------------------------------------------------
// Return the surface fraction
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Scalar>&
SlideSurface<Dimension>::
surfaceFraction() const {
  return mSurfaceFraction;
}

//------------------------------------------------------------------------------
// smoothness metric for mixing interfaces
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Scalar>&
SlideSurface<Dimension>::
surfaceSmoothness() const {
  return mSurfaceSmoothness;
}

//------------------------------------------------------------------------------
// next time step  surface normal field list ref
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Vector>&
SlideSurface<Dimension>::
newSurfaceNormals() const {
  return mNewSurfaceNormals;
}

//------------------------------------------------------------------------------
// next time step  surface normal field list ref
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Vector>&
SlideSurface<Dimension>::
newSmoothedSurfaceNormals() const {
  return mNewSmoothedSurfaceNormals;
}

//------------------------------------------------------------------------------
// next time step  surface fraction
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Scalar>&
SlideSurface<Dimension>::
newSurfaceFraction() const {
  return mNewSurfaceFraction;
}

//------------------------------------------------------------------------------
// next time step smoothness metric for mixing interfaces
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Scalar>&
SlideSurface<Dimension>::
newSurfaceSmoothness() const {
  return mNewSurfaceSmoothness;
}

//------------------------------------------------------------------------------
// next time step smoothness metric for mixing interfaces
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Scalar>&
SlideSurface<Dimension>::
smoothnessNormalization() const {
  return mSmoothnessNormalization;
}

//------------------------------------------------------------------------------
// set/get bool list of interactions 
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
SlideSurface<Dimension>::
isSlideSurface(const std::vector<bool> x) {
  mIsSlideSurface = x;
}
template<typename Dimension>
inline
std::vector<bool>
SlideSurface<Dimension>::
isSlideSurface() const {
  return mIsSlideSurface;
}

//------------------------------------------------------------------------------
// set/get bool if any slides exist
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
SlideSurface<Dimension>::
isActive(const bool x) {
  mIsActive = x;
}
template<typename Dimension>
inline
bool
SlideSurface<Dimension>::
isActive() const {
  return mIsActive;
}

//------------------------------------------------------------------------------
// set/get bool to turn normal smoothing on or off
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
SlideSurface<Dimension>::
normalsAreSmoothed(const bool x) {
  mNormalsAreSmoothed = x;
}
template<typename Dimension>
inline
bool
SlideSurface<Dimension>::
normalsAreSmoothed() const {
  return mNormalsAreSmoothed;
}

//------------------------------------------------------------------------------
// set/get bool to turn off gradient correction for normal calculation
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
SlideSurface<Dimension>::
gradientsAreCorrected(const bool x) {
  mGradientsAreCorrected = x;
}
template<typename Dimension>
inline
bool
SlideSurface<Dimension>::
gradientsAreCorrected() const {
  return mGradientsAreCorrected;
}


//------------------------------------------------------------------------------
// set/get number of node lists
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
SlideSurface<Dimension>::
numNodeLists(const int x) {
  mNumNodeLists = x;
}
template<typename Dimension>
inline
int
SlideSurface<Dimension>::
numNodeLists() const {
  return mNumNodeLists;
}


//------------------------------------------------------------------------------
// set/get our method of calculating surface normals
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
SlideSurface<Dimension>::
surfaceNormalMethod(const SurfaceNormalMethod x) {
  mSurfaceNormalMethod = x;
}
template<typename Dimension>
inline
SurfaceNormalMethod
SlideSurface<Dimension>::
surfaceNormalMethod() const {
  return mSurfaceNormalMethod;
}

}