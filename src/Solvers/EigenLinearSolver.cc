//---------------------------------Spheral++----------------------------------//
// EigenLinearSolver
//
// Represents a solver for a linear system of equations
//----------------------------------------------------------------------------//

#include "EigenLinearSolver.hh"

#include <map>
#include <numeric> // for std::accumulate
#include "Utilities/DBC.hh"
#include "Distributed/Process.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
EigenLinearSolver::
EigenLinearSolver(std::shared_ptr<EigenOptions> options):
  mOptions(options),
  mGraphChangedSinceFill(true),
  mMatrixChangedSinceFactorization(true) {
}

//------------------------------------------------------------------------------
// Initialize the graph
//------------------------------------------------------------------------------
void
EigenLinearSolver::
initializeGraph() {
  VERIFY(this->map());
  VERIFY(Process::getTotalNumberOfProcesses() == 1);

  // Optionally do not reinitialize
  if (mOptions->reuseSolver && this->readyToSolve()) {
    return;
  }
  
  const auto map = this->map();
  const auto numNodes = map->numLocalElements();
  
  mMatrix.resize(numNodes, numNodes);
  mRhs.resize(numNodes);
  mLhs.resize(numNodes);
  
  mGraphChangedSinceFill = true;
}

//------------------------------------------------------------------------------
// Initialize the matrix
//------------------------------------------------------------------------------
void
EigenLinearSolver::
initializeMatrix() {
  VERIFY(this->map());
  VERIFY(this->data());
  
  // Optionally do not reinitialize
  if (mOptions->reuseSolver && this->readyToSolve()) {
    return;
  }
  
  const auto map = this->map();
  const auto data = this->data();
  
  const auto numElements = map->numLocalElements();
  const auto numElementsPerRow = map->numElementsPerRow();
  const auto numEntries = std::accumulate(numElementsPerRow.begin(), numElementsPerRow.end(), 0);
  
  std::vector<DataType> triplets;
  triplets.reserve(numEntries);
  std::vector<int> indices;
  std::vector<double> values;
  for (auto i = 0; i < numElements; ++i) {
    data->getRowValues(i,
                       indices,
                       values);
    CHECK(indices.size() == values.size());
    const auto numVals = indices.size();
    for (auto j = 0u; j < numVals; ++j) {
      triplets.emplace_back(i, indices[j], values[j]);
    }
  }
  
  mMatrix.setFromTriplets(triplets.begin(), triplets.end());
  mMatrix.makeCompressed();
  
  mGraphChangedSinceFill = false;
  mMatrixChangedSinceFactorization = true;
}

//------------------------------------------------------------------------------
// Initialize the solver
//------------------------------------------------------------------------------
void
EigenLinearSolver::
initializeSolver() {
  VERIFY(!mGraphChangedSinceFill);

  // Optionally do not reinitialize
  if (mOptions->reuseSolver && this->readyToSolve()) {
    return;
  }
  
  if (mOptions->qr) {
    mQRSolver.analyzePattern(mMatrix);
    mQRSolver.factorize(mMatrix);
  }
  else {
    mLUSolver.analyzePattern(mMatrix);
    mLUSolver.factorize(mMatrix);
  }

  mMatrixChangedSinceFactorization = false;
}

//------------------------------------------------------------------------------
// Solve the system
//------------------------------------------------------------------------------
void
EigenLinearSolver::
solve(const std::vector<double>& input,
      std::vector<double>& output) {
  VERIFY(this->readyToSolve());
  
  const auto map = this->map();
  const auto numElements = map->numLocalElements();

  for (auto i = 0; i < numElements; ++i) {
    mRhs(i) = input[i];
  }

  if (mOptions->qr) {
    mLhs = mQRSolver.solve(mRhs);
  }
  else {
    mLhs = mLUSolver.solve(mRhs);
  }

  for (auto i = 0; i < numElements; ++i) {
    output[i] = mLhs(i);
  }
}

//------------------------------------------------------------------------------
// Multiply matrix by a vector
//------------------------------------------------------------------------------
void
EigenLinearSolver::
multiply(const std::vector<double>& input,
         std::vector<double>& output) {
  VERIFY2(false, "Eigen sparse multiplication not implemented");
  VERIFY(this->readyToSolve());

  const auto map = this->map();
  const auto numElements = map->numLocalElements();

  for (auto i = 0; i < numElements; ++i) {
    mRhs(i) = input[i];
  }

  mLhs = mMatrix * mRhs;

  for (auto i = 0; i < numElements; ++i) {
    output[i] = mLhs(i);
  }
}

} // end namespace Spheral
