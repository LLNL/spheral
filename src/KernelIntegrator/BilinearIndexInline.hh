#include <limits>
#include <tuple>

namespace Spheral {

//------------------------------------------------------------------------------
// Get specialized hash functions
//------------------------------------------------------------------------------
// For a pair, assume the two numbers are small and hash equally
template<>
inline
int
BilinearHash<std::pair<int, int>>::
operator()(const std::pair<int, int>& v) const {
  return (v.first << std::numeric_limits<int>::digits / 2) ^ v.second;
}

// For an array, we are hashing a normal, and the final element is either 0 or 1
template<>
inline
int
BilinearHash<std::array<int, 1>>::
operator()(const std::array<int, 1>& v) const {
  return v[0];
}
template<>
inline
int
BilinearHash<std::array<int, 2>>::
operator()(const std::array<int, 2>& v) const {
  return (v[0] << 1) ^ v[1];
}
template<>
inline
int
BilinearHash<std::array<int, 3>>::
operator()(const std::array<int, 3>& v) const {
  return (v[0] << std::numeric_limits<int>::digits / 2) ^ (v[1] << 1) ^ v[2];
}

//------------------------------------------------------------------------------
// Return the number of nonzero indices for this row
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
BilinearIndex<Dimension>::
numElementsForRow(const std::pair<int, int>& ni) const {
  CHECK(mConnectivityComputed);
  return mBilinearNodeIndex(ni.first, ni.second).size();
}

//------------------------------------------------------------------------------
// Return the flat index j, given the row ni and the column nj
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
BilinearIndex<Dimension>::
flatIndex(const std::pair<int, int>& ni,
          const std::pair<int, int>& nj) const {
  CHECK(mConnectivityComputed);
  CHECK(mBilinearFlatIndex(ni.first, ni.second).count(nj) > 0);
  return mBilinearFlatIndex(ni.first, ni.second).at(nj);
}

//------------------------------------------------------------------------------
// Return the node index nj, given the row ni and the index j
//------------------------------------------------------------------------------
template<typename Dimension>
inline
std::pair<int, int>
BilinearIndex<Dimension>::
nodeIndex(const std::pair<int, int>& ni,
          const int& j) const {
  CHECK(mConnectivityComputed);
  CHECK(size_t(j) < mBilinearNodeIndex(ni.first, ni.second).size());
  return mBilinearNodeIndex(ni.first, ni.second)[j];
}

//------------------------------------------------------------------------------
// Return the node indices for the row ni
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector<std::pair<int, int>>& 
BilinearIndex<Dimension>::
nodeIndices(const std::pair<int, int>& ni) const {
  CHECK(mConnectivityComputed);
  return mBilinearNodeIndex(ni.first, ni.second);
}

//------------------------------------------------------------------------------
// Return the number of unique surfaces for row ni
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
BilinearIndex<Dimension>::
numSurfacesForRow(const std::pair<int, int>& ni) const {
  CHECK(mSurfaceIndexingComputed);
  return mSurfaceNormal(ni.first, ni.second).size();
}

//------------------------------------------------------------------------------
// Convert from the normal to rounded integer values
// We can convert one dimension to only sign dependence due to normalization
//------------------------------------------------------------------------------
template<>
inline
void
BilinearIndex<Dim<1>>::
normalToArray(const Vector& normal,
              ArrayDim& values) const {
  CHECK(fuzzyEqual(normal.magnitude(), 1.0));
  values[0] = normal[0] > 0;
}

template<>
inline
void
BilinearIndex<Dim<2>>::
normalToArray(const Vector& normal,
              ArrayDim& values) const {
  CHECK(fuzzyEqual(normal.magnitude(), 1.0));
  values[0] = static_cast<int>(normal[0] * mRoundValue + 0.5);
  values[1] = static_cast<int>(normal[1] * mRoundValue + 0.5 > 0);
}

template<>
inline
void
BilinearIndex<Dim<3>>::
normalToArray(const Vector& normal,
              ArrayDim& values) const {
  CHECK(fuzzyEqual(normal.magnitude(), 1.0));
  values[0] = static_cast<int>(normal[0] * mRoundValue + 0.5);
  values[1] = static_cast<int>(normal[1] * mRoundValue + 0.5);
  values[2] = static_cast<int>(normal[2] * mRoundValue + 0.5 > 0);
}

//------------------------------------------------------------------------------
// Convert from the rounded integer values to the normal
//------------------------------------------------------------------------------
template<>
inline
void
BilinearIndex<Dim<1>>::
arrayToNormal(const ArrayDim& values,
              Vector& normal) const {
  CHECK(values[0] == 0 || values[0] == 1);
  normal[0] = values[0] ? 1 : -1;
  ENSURE(-1 <= normal[0] && normal[0] <= 1);
}

template<>
inline
void
BilinearIndex<Dim<2>>::
arrayToNormal(const ArrayDim& values,
              Vector& normal) const {
  normal[0] = values[0] * mRoundValueInv;
  const auto normal1Mag = sqrt(1 - normal[0] * normal[0]);
  const auto normal1Sign = values[0] ? 1 : -1;
  normal[1] = normal1Mag * normal1Sign;
  ENSURE(fuzzyEqual(normal.magnitude(), 1.0));
}

template<>
inline
void
BilinearIndex<Dim<3>>::
arrayToNormal(const ArrayDim& values,
              Vector& normal) const {
  normal[0] = values[0] * mRoundValueInv;
  normal[1] = values[1] * mRoundValueInv;
  const auto normal2Mag = sqrt(1 - normal[0] * normal[0] - normal[1] * normal[1]);
  const auto normal2Sign = values[0] ? 1 : -1;
  normal[2] = normal2Mag * normal2Sign;
  ENSURE(fuzzyEqual(normal.magnitude(), 1.0));
}

//------------------------------------------------------------------------------
// Return the surface index for the given i and j
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
BilinearIndex<Dimension>::
surfaceIndex(const std::pair<int, int>& ni,
             const Vector& normal) const {
  normalToArray(normal,
                mScratchArray);
  CHECK(mSurfaceIndexingComputed);
  CHECK(mSurfaceFlatIndex(ni.first, ni.second).count(mScratchArray) > 0);
  return mSurfaceFlatIndex(ni.first, ni.second).at(mScratchArray);
}

//------------------------------------------------------------------------------
// Return the surface index for the given i and j
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
BilinearIndex<Dimension>::
surfaceIndex(const std::pair<int, int>& ni,
             const ArrayDim& values) const {
  CHECK(mSurfaceIndexingComputed);
  CHECK(mSurfaceFlatIndex(ni.first, ni.second).count(values) > 0);
  return mSurfaceFlatIndex(ni.first, ni.second).at(values);
}

//------------------------------------------------------------------------------
// Return the normal values
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector<typename Dimension::Vector>&
BilinearIndex<Dimension>::
normal(const std::pair<int, int>& ni) const {
  CHECK(mSurfaceIndexingComputed);
  return mSurfaceNormal(ni.first, ni.second);
}

//------------------------------------------------------------------------------
// Return the surface index for the given i and j
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector<int>&
BilinearIndex<Dimension>::
surfaceFlags(const std::pair<int, int>& ni) const {
  return mSurfaceFlags(ni.first, ni.second);
}

} // end namespace Spheral
