//---------------------------------Spheral++----------------------------------//
// LinearSolver
//
// Represents a solver for a linear system of equations
// A given type only needs to fill in the pure virtual methods, and the rest
// is done through the setMap and setData functions. 
//----------------------------------------------------------------------------//
#ifndef __Spheral_LinearSolver_hh__
#define __Spheral_LinearSolver_hh__

#include <memory>
#include <vector>

#include "MatrixMap.hh"
#include "MatrixData.hh"

namespace Spheral {

template<class DataType> class IncrementalStatistic;

class LinearSolver {
public:
  // Constructor
  LinearSolver();

  // Check whether map and data have been set
  virtual bool mapSet() const;
  virtual bool dataSet() const;

  // Check whether we are ready to solve
  virtual bool readyToSolve() const = 0;
  
  // Set map, initialize graph if applicable
  virtual void setMap(std::shared_ptr<MatrixMap> map);
  
  // Set data, initialize matrix and solver if applicable
  virtual void setData(std::shared_ptr<MatrixData> data);

  // Solve the system of equations
  virtual void solve(const std::vector<double>& input,
                     std::vector<double>& output) = 0;

  // Multiply the matrix by a vector
  virtual void multiply(const std::vector<double>& input,
                        std::vector<double>& output) = 0;

  // Return the map and data
  virtual const std::shared_ptr<MatrixMap> map() const;
  virtual const std::shared_ptr<MatrixData> data() const;
  
  // Return statistics if available, empty vector if not
  virtual std::vector<std::shared_ptr<IncrementalStatistic<double>>> statistics() const;
  
protected:
  // Initialize the graph, matrix, and solver
  virtual void initializeGraph() = 0;
  virtual void initializeMatrix() = 0;
  virtual void initializeSolver() = 0;
  
private:
  // Map and data
  std::shared_ptr<MatrixMap> mMap;
  std::shared_ptr<MatrixData> mData;
}; // end class LinearSolver

} // end namespace Spheral

#include "LinearSolverInline.hh"

#endif
