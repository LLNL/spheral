namespace Spheral {

//------------------------------------------------------------------------------
// The FacetedVolume
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename Dimension::FacetedVolume&
FacetedVolumeBoundary<Dimension>::polyVolume() const {
  return mPoly;
}

//------------------------------------------------------------------------------
// interior flag
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
FacetedVolumeBoundary<Dimension>::interiorBoundary() const {
  return mInteriorBoundary;
}

//------------------------------------------------------------------------------
// Create ghosts or not
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
FacetedVolumeBoundary<Dimension>::useGhosts() const {
  return mUseGhosts;
}

//------------------------------------------------------------------------------
// The reflection operator per facet of the FacetedVolume
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename Dimension::Tensor&
FacetedVolumeBoundary<Dimension>::reflectOperator(const unsigned facetID) const {
  REQUIRE(facetID < mReflectOperators.size());
  return mReflectOperators[facetID];
}

}
