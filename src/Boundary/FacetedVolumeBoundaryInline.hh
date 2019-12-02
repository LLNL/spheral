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
// The reflection operator per facet of the FacetedVolume
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename Dimension::Tensor&
FacetedVolumeBoundary<Dimension>::reflectOperator(const unsigned facetID) const {
  REQUIRE(facetID < mReflectOperators.size());
  return mReflectOperator[facetID];
}

}
