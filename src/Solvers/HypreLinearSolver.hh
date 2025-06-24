//---------------------------------Spheral++----------------------------------//
// HypreLinearSolver
//
// Represents a solver for a linear system of equations
//----------------------------------------------------------------------------//
#ifndef __Spheral_HypreLinearSolver_hh__
#define __Spheral_HypreLinearSolver_hh__

#include <memory>
#include <vector>

#include "HypreOptions.hh"
#include "LinearSolver.hh"

struct hypre_IJMatrix_struct;
struct hypre_IJVector_struct;
struct hypre_Solver_struct;

namespace Spheral {

template<typename DataType> class IncrementalStatistic;

class HypreLinearSolver : public LinearSolver {
public:
  typedef struct hypre_IJMatrix_struct MatrixType;
  typedef struct hypre_IJVector_struct VectorType;
  typedef struct hypre_Solver_struct SolverType;
  
  // Constructor
  HypreLinearSolver(std::shared_ptr<HypreOptions> options);

  // Check whether we are ready to solve
  virtual bool readyToSolve() const override;
  
  // Solve the system of equations
  virtual void solve(const std::vector<double>& input,
                     std::vector<double>& output) override;

  // Multiply the matrix by a vector
  virtual void multiply(const std::vector<double>& input,
                        std::vector<double>& output) override;
  
  // Return statistics
  virtual std::vector<std::shared_ptr<IncrementalStatistic<double>>> statistics() const override;
  
protected:
  // Initialize the graph, matrix, and solver
  virtual void initializeGraph() override;
  virtual void initializeMatrix() override;
  virtual void initializeSolver() override;
  
private:

  // Get Hypre objects that self-destruct
  virtual std::shared_ptr<VectorType> createVector() const;
  virtual std::shared_ptr<MatrixType> createMatrix() const;
  virtual std::shared_ptr<SolverType> createSolver() const;
  virtual std::shared_ptr<SolverType> createPreconditioner(std::shared_ptr<SolverType> solver) const;

  // Fill given matrix with coefficients from MatrixData
  virtual void fillMatrix(std::shared_ptr<MatrixType> hypreMatrix) const;

  // Set or get the values of a Hypre vector
  virtual void setVectorValues(std::vector<double> const& x,
                               std::shared_ptr<VectorType> hypreVector) const;
  virtual void getVectorValues(std::shared_ptr<VectorType> hypreVector,
                               std::vector<double> &x) const;

  // Solve system that is stored, return whether solve was successful
  bool solve();

  // Multiply system that is stored
  void multiply();

  // Has graph changed since fill?
  bool mGraphChangedSinceFill;
  
  // Data
  std::shared_ptr<HypreOptions> mHypreOptions;
  std::shared_ptr<MatrixType> mHypreMatrix;
  std::shared_ptr<VectorType> mHypreLHS;
  std::shared_ptr<VectorType> mHypreRHS;
  std::shared_ptr<VectorType> mHypreResidual;
  std::shared_ptr<SolverType> mHypreSolver;
  std::shared_ptr<SolverType> mHyprePreconditioner;

  // Iteration statistics
  std::shared_ptr<IncrementalStatistic<double>> mIterationStatistics;
  std::shared_ptr<IncrementalStatistic<double>> mFinalResidualStatistics;
}; // end class HypreLinearSolver

} // end namespace Spheral

#include "HypreLinearSolverInline.hh"

#endif
