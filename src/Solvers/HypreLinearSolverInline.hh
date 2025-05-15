namespace Spheral {

//------------------------------------------------------------------------------
// Return statistics for number of iterations and tolerance reached
//------------------------------------------------------------------------------
inline
std::vector<std::shared_ptr<IncrementalStatistic<double>>>
HypreLinearSolver::
statistics() const {
  return {mIterationStatistics, mFinalResidualStatistics};
}

//------------------------------------------------------------------------------
// Return whether we are ready to call solve
//------------------------------------------------------------------------------
inline
bool
HypreLinearSolver::
readyToSolve() const {
  return (!mGraphChangedSinceFill &&
          this->mapSet() && this->dataSet());
}

} // end namespace Spheral

