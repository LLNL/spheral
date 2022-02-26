namespace Spheral {


template<typename Dimension>
inline
const FieldList<Dimension, int>&
PairFields<Dimension>::
uniqueIndices() const {
  return mUniqueIndices;
}

template<typename Dimension>
inline
const FieldList<Dimension, vector<int>>&
PairFields<Dimension>::
neighborIndices() const {
  return mNeighborIndices;
}

template<typename Dimension>
inline
const FieldList<Dimension, vector<typename Dimension::Vector>>&
PairFields<Dimension>::
shearDisplacement() const {
  return mShearDisplacement;
}

template<typename Dimension>
inline
const FieldList<Dimension, vector<typename Dimension::Vector>>&
PairFields<Dimension>::
DDtShearDisplacement() const {
  return mDDtShearDisplacement;
}

template<typename Dimension>
inline
const FieldList<Dimension, vector<typename Dimension::Scalar>>&
PairFields<Dimension>::
equilibriumOverlap() const {
  return mEquilibriumOverlap;
}

}