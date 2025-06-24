namespace Spheral {

//------------------------------------------------------------------------------
// Return whether we are ready to call solve
//------------------------------------------------------------------------------
inline
bool
EigenLinearSolver::
readyToSolve() const {
  return (!mGraphChangedSinceFill && !mMatrixChangedSinceFactorization &&
          this->mapSet() && this->dataSet());
}

} // end namespace Spheral

