//---------------------------------Spheral++----------------------------------//
// EigenLinearSolver
//
// Solver linear system using Eigen
// Only works on one processor
//----------------------------------------------------------------------------//
#ifndef __Spheral_EigenLinearSolver_hh__
#define __Spheral_EigenLinearSolver_hh__

#include <memory>
#include <vector>

#include "Eigen/OrderingMethods"
#include "Eigen/Sparse"
#include "Eigen/SparseLU"
#include "Eigen/SparseQR"

#include "EigenOptions.hh"
#include "LinearSolver.hh"

namespace Spheral {

class EigenLinearSolver : public LinearSolver {
public:
  typedef Eigen::SparseMatrix<double> MatrixType;
  typedef Eigen::Triplet<double> DataType;
  typedef Eigen::VectorXd VectorType;
  typedef Eigen::SparseLU<MatrixType, Eigen::COLAMDOrdering<int>> LUSolverType;
  typedef Eigen::SparseQR<MatrixType, Eigen::COLAMDOrdering<int>> QRSolverType;
  
  // Constructor
  EigenLinearSolver(std::shared_ptr<EigenOptions> options);

  // Check whether we are ready to solve
  virtual bool readyToSolve() const override;
  
  // Solve the system of equations
  virtual void solve(const std::vector<double>& input,
                     std::vector<double>& output) override;

  // Multiply the matrix by a vector
  virtual void multiply(const std::vector<double>& input,
                        std::vector<double>& output) override;

protected:
  // Initialize the graph, matrix, and solver
  virtual void initializeGraph() override;
  virtual void initializeMatrix() override;
  virtual void initializeSolver() override;
  
private:
  // Input data
  std::shared_ptr<EigenOptions> mOptions;
  
  // Have we initialized things correctly?
  bool mGraphChangedSinceFill;
  bool mMatrixChangedSinceFactorization;
  
  // Eigen data
  MatrixType mMatrix;
  LUSolverType mLUSolver;
  QRSolverType mQRSolver;
  VectorType mRhs;
  VectorType mLhs;
  
}; // end class EigenLinearSolver

} // end namespace Spheral

#include "EigenLinearSolverInline.hh"

#endif
