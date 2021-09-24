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
}