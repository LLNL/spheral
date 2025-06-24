//---------------------------------Spheral++----------------------------------//
// LinearSolver
//
// Represents a solver for a linear system of equations
//----------------------------------------------------------------------------//

#include "LinearSolver.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
LinearSolver::
LinearSolver() {
}

//------------------------------------------------------------------------------
// Initialize the map data
//------------------------------------------------------------------------------
void
LinearSolver::
setMap(std::shared_ptr<MatrixMap> map) {
  // Set map
  mMap = map;

  // Remove reference to data, as if we change the map the data must also change
  mData = std::shared_ptr<MatrixData>();
  
  // Initialize graph
  initializeGraph();
}

//------------------------------------------------------------------------------
// Initialize the map data
//------------------------------------------------------------------------------
void
LinearSolver::
setData(std::shared_ptr<MatrixData> data) {
  // Make sure map has been set
  VERIFY(mapSet());
  
  // Set data
  mData = data;

  // Initialize matrix and solver
  initializeMatrix();
  initializeSolver();
}

} // end namespace Spheral
