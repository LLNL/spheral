namespace Spheral {

//------------------------------------------------------------------------------
// Return empty statistics
//------------------------------------------------------------------------------
inline
std::vector<std::shared_ptr<IncrementalStatistic<double>>>
LinearSolver::
statistics() const {
  return std::vector<std::shared_ptr<IncrementalStatistic<double>>>();
}

//------------------------------------------------------------------------------
// Return whether map and data have been set
//------------------------------------------------------------------------------
inline
bool
LinearSolver::
mapSet() const {
  return static_cast<bool>(mMap);
}

inline
bool
LinearSolver::
dataSet() const {
  return static_cast<bool>(mData);
}

//------------------------------------------------------------------------------
// Return the matrix map and data
//------------------------------------------------------------------------------
inline
const std::shared_ptr<MatrixMap>
LinearSolver::
map() const {
  return mMap;
}

inline
const std::shared_ptr<MatrixData>
LinearSolver::
data() const {
  return mData;
}

} // end namespace Spheral
