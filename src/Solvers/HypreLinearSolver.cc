//---------------------------------Spheral++----------------------------------//
// HypreLinearSolver
//
// Represents a solver for a linear system of equations
//----------------------------------------------------------------------------//

#include "HypreLinearSolver.hh"

#include <map>
#include <numeric> // for std::iota
#include "_hypre_parcsr_ls.h"
#include "HYPRE.h"
#include "HYPRE_IJ_mv.h"
#include "Distributed/Communicator.hh"
#include "Solvers/IncrementalStatistic.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
HypreLinearSolver::
HypreLinearSolver(std::shared_ptr<HypreOptions> options) :
  mGraphChangedSinceFill(true),
  mHypreOptions(options),
  mIterationStatistics(std::make_shared<IncrementalStatistic<double>>(options->meanIterationsGuess,
                                                                      "HypreLinearSolver iterations",
                                                                      options->printIterations)), // print
  mFinalResidualStatistics(std::make_shared<IncrementalStatistic<double>>(options->toleranceL2,
                                                                          "HypreLinearSolver residual",
                                                                          false)) { // print
}

//------------------------------------------------------------------------------
// Initialize the graph
//------------------------------------------------------------------------------
void
HypreLinearSolver::
initializeGraph() {
  // Create vectors and matrix
  mHypreLHS = createVector();
  mHypreRHS = createVector();
  mHypreResidual = createVector();
  mHypreMatrix = createMatrix();
  mGraphChangedSinceFill = true;
}

//------------------------------------------------------------------------------
// Initialize the matrix
//------------------------------------------------------------------------------
void
HypreLinearSolver::
initializeMatrix() {
  // Make sure things are initialized
  VERIFY(mHypreMatrix);
  
  // Fill matrix
  fillMatrix(mHypreMatrix);
  mGraphChangedSinceFill = false;

  // Set statistic description to match data
  mIterationStatistics->setName("HypreLinearSolver (" + this->data()->description() + ") iterations");
  mFinalResidualStatistics->setName("HypreLinearSolver (" + this->data()->description() + ") residual");
}

//------------------------------------------------------------------------------
// Initialize the solver
//------------------------------------------------------------------------------
void
HypreLinearSolver::
initializeSolver() {
  // Create solver and preconditioner
  mHypreSolver = createSolver();
  mHyprePreconditioner = createPreconditioner(mHypreSolver);
}

//------------------------------------------------------------------------------
// Solve the system
//------------------------------------------------------------------------------
void
HypreLinearSolver::
solve(const std::vector<double>& input,
      std::vector<double>& output) {
  // Make sure data is initialized
  VERIFY(this->mapSet() && this->dataSet());
  VERIFY(!mGraphChangedSinceFill);
  
  // Set the source
  setVectorValues(input, mHypreRHS); // include override values
  
  // Set the initial guess to also be the source, unless otherwise specified
  const auto data = this->data();
  if (data->initialGuessAvailable()) {
    data->getInitialGuess(output);
    setVectorValues(output, mHypreLHS);
  }
  else {
    setVectorValues(input, mHypreLHS);
  }

  // Solve the system
  const auto solved = solve();
  if (!solved) {
    if (mHypreOptions->quitIfDiverged) {
      VERIFY2(false, "Hypre solver (" + this->data()->description() + ") did not converge");
    }
    else {
      if (mHypreOptions->warnIfDiverged && Process::getRank() == 0) {
        std::cerr << "Hypre solver (" + this->data()->description() + ") did not converge" << std::endl;
      }
    }
  }
  
  // Get the result
  getVectorValues(mHypreLHS, output);
}

//------------------------------------------------------------------------------
// Multiply matrix by a vector
//------------------------------------------------------------------------------
void
HypreLinearSolver::
multiply(const std::vector<double>& input,
         std::vector<double>& output) {
  // Make sure data is initialized
  VERIFY(this->mapSet() && this->dataSet());
  VERIFY(!mGraphChangedSinceFill);

  // Set the source
  setVectorValues(input, mHypreLHS); // include override values
  
  // Solve the system
  multiply();
  
  // Get the result
  getVectorValues(mHypreRHS, output);
}
//------------------------------------------------------------------------------
// Apply inverse to given vector and return result in place
//------------------------------------------------------------------------------
static int solveCall = 0;
bool
HypreLinearSolver::
solve() {
  VERIFY(this->readyToSolve());
  
  // Get storage in CRS format
  HYPRE_ParCSRMatrix hypreParMatrix;
  HYPRE_ParVector hypreParLHS;
  HYPRE_ParVector hypreParRHS;

  int hypreStatus;
  hypreStatus = HYPRE_IJMatrixGetObject(mHypreMatrix.get(),
                                        (void **)&hypreParMatrix);
  VERIFY(hypreStatus == 0);

  hypreStatus = HYPRE_IJVectorGetObject(mHypreLHS.get(),
                                        (void **)&hypreParLHS);
  VERIFY(hypreStatus == 0);

  hypreStatus = HYPRE_IJVectorGetObject(mHypreRHS.get(),
                                        (void **)&hypreParRHS);
  VERIFY(hypreStatus == 0);
  
  if (mHypreOptions->saveLinearSystem) {
    const auto name = this->data()->description();
    const auto matName = "HypreMatrix_" + name + "_" + std::to_string(solveCall);
    const auto rhsName = "HypreRHS_" + name + "_" + std::to_string(solveCall);
    const auto initName = "HypreInit_" + name + "_" + std::to_string(solveCall);
    HYPRE_IJMatrixPrint(mHypreMatrix.get(), matName.c_str());
    HYPRE_IJVectorPrint(mHypreRHS.get(), rhsName.c_str());
    HYPRE_IJVectorPrint(mHypreLHS.get(), initName.c_str());
  }
  
  hypreStatus = HYPRE_ParCSRGMRESSetLogging(mHypreSolver.get(), mHypreOptions->logLevel);
  VERIFY(hypreStatus == 0);

  hypreStatus = HYPRE_ParCSRGMRESSetPrintLevel(mHypreSolver.get(), mHypreOptions->printLevel);
  VERIFY(hypreStatus == 0);

  hypreStatus
    = HYPRE_ParCSRGMRESSetMaxIter(mHypreSolver.get(), mHypreOptions->maxNumberOfIterations);
  VERIFY(hypreStatus == 0);

  hypreStatus = HYPRE_ParCSRGMRESSetTol(mHypreSolver.get(), mHypreOptions->toleranceL2);
  VERIFY(hypreStatus == 0);

  if (mHypreOptions->useRobustTolerance) {
    hypreStatus = HYPRE_GMRESSetSkipRealResidualCheck(mHypreSolver.get(), 1);
    VERIFY(hypreStatus == 0);
  }

  hypreStatus
    = HYPRE_ParCSRGMRESSetAbsoluteTol(mHypreSolver.get(), mHypreOptions->absoluteTolerance);
  VERIFY(hypreStatus == 0);

  hypreStatus = HYPRE_ParCSRGMRESSetup(mHypreSolver.get(),
                                       hypreParMatrix,
                                       hypreParRHS,
                                       hypreParLHS);
  VERIFY2(hypreStatus == 0, "ParCSRGMRESSetup");

  int solveStatus = HYPRE_ParCSRGMRESSolve(mHypreSolver.get(),
                                           hypreParMatrix,
                                           hypreParRHS,
                                           hypreParLHS);
  int numIterations = 0;
  if (solveStatus) {
    VERIFY2(!HYPRE_CheckError(solveStatus, HYPRE_ERROR_GENERIC),
            "INF or NaN in Matrix, RHS, or initial guess");
    VERIFY2(!HYPRE_CheckError(solveStatus, HYPRE_ERROR_MEMORY),
            "HYPRE cannot allocate memory");
    VERIFY2(HYPRE_CheckError(solveStatus, HYPRE_ERROR_CONV),
            "GMRES failed with (unknown) solver status: " << solveStatus);
    HYPRE_ClearError(HYPRE_ERROR_CONV);
    numIterations = mHypreOptions->maxNumberOfIterations;
  }
  else {
    hypreStatus
      = HYPRE_ParCSRGMRESGetNumIterations(mHypreSolver.get(), &numIterations);
    VERIFY(hypreStatus == 0);
  }
  mIterationStatistics->add(numIterations);

  double finalResNorm;
  hypreStatus = HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm(mHypreSolver.get(),
                                                               &finalResNorm);
  VERIFY(hypreStatus == 0);
  mFinalResidualStatistics->add(finalResNorm);

  if (mHypreOptions->saveLinearSystem) {
    const auto name = this->data()->description();
    const auto lhsName = "HypreLHS_" + name + "_" + std::to_string(solveCall);
    HYPRE_IJVectorPrint(mHypreLHS.get(), lhsName.c_str());
    ++solveCall;
  }
  
  return numIterations < mHypreOptions->maxNumberOfIterations;
}

//------------------------------------------------------------------------------
// Apply matrix-vector product
//------------------------------------------------------------------------------
// static int solve_call = 0;
void
HypreLinearSolver::
multiply() {
  // Get storage in CRS format
  HYPRE_ParCSRMatrix hypreParMatrix;
  HYPRE_ParVector hypreParLHS;
  HYPRE_ParVector hypreParRHS;

  int hypreStatus;
  hypreStatus = HYPRE_IJMatrixGetObject(mHypreMatrix.get(),
                                        (void **)&hypreParMatrix);
  VERIFY(hypreStatus == 0);

  hypreStatus = HYPRE_IJVectorGetObject(mHypreLHS.get(),
                                        (void **)&hypreParLHS);
  VERIFY(hypreStatus == 0);

  hypreStatus = HYPRE_IJVectorGetObject(mHypreRHS.get(),
                                        (void **)&hypreParRHS);
  VERIFY(hypreStatus == 0);

  // Perform matrix multiplication, y = \alpha * A x + \beta y
  hypreStatus = HYPRE_ParCSRMatrixMatvec(1.0, // \alpha
                                         hypreParMatrix, // A
                                         hypreParLHS, // x
                                         0.0, // \beta
                                         hypreParRHS); // y
  VERIFY(hypreStatus == 0);
}

//------------------------------------------------------------------------------
// Put values into hypre vector
//------------------------------------------------------------------------------
void
HypreLinearSolver::
setVectorValues(std::vector<double> const& x,
                std::shared_ptr<VectorType> hypreVector) const {
  // Get map data
  const auto map = this->map();
  const auto firstGlobalIndex = map->firstGlobalIndex();
  const auto numLocalElements = map->numLocalElements();
  std::vector<int> globalIndices(numLocalElements);
  std::iota(globalIndices.begin(), globalIndices.end(), firstGlobalIndex);
  VERIFY(x.size() == size_t(numLocalElements));
  
  // Check for NaN
  BEGIN_CONTRACT_SCOPE
  {
    for (auto i = 0; i < numLocalElements; ++i) {
      CHECK2(!std::isnan(x[i]), "NaN found in input to Hypre at index " + std::to_string(i));
    }
  }
  END_CONTRACT_SCOPE
  
  // Set vector values
  int hypreStatus;
  hypreStatus = HYPRE_IJVectorSetValues(hypreVector.get(),
                                        numLocalElements,
                                        &globalIndices[0],
                                        &x[0]);
  VERIFY(hypreStatus == 0);
  
}

//------------------------------------------------------------------------------
// Get values from hypre vector
//------------------------------------------------------------------------------
void
HypreLinearSolver::
getVectorValues(std::shared_ptr<VectorType> hypreVector,
                std::vector<double> & x) const {
  // Get map data
  const auto map = this->map();
  const auto firstGlobalIndex = map->firstGlobalIndex();
  const auto numLocalElements = map->numLocalElements();
  std::vector<int> globalIndices(numLocalElements);
  std::iota(globalIndices.begin(), globalIndices.end(), firstGlobalIndex);
  
  // Get vector values
  x.resize(numLocalElements);
  int hypreStatus;
  hypreStatus = HYPRE_IJVectorGetValues(hypreVector.get(),
                                        numLocalElements,
                                        &globalIndices[0],
                                        &x[0]);
  VERIFY(hypreStatus == 0);

  // Check for NaN
  BEGIN_CONTRACT_SCOPE
  {
    for (auto i = 0; i < numLocalElements; ++i) {
      CHECK2(!std::isnan(x[i]), "NaN found in Hypre output at index " + std::to_string(i));
    }
  }
  END_CONTRACT_SCOPE
}

//------------------------------------------------------------------------------
// Create Hypre vector that should self-destruct when it goes out of scope.
//------------------------------------------------------------------------------
std::shared_ptr<HypreLinearSolver::VectorType>
HypreLinearSolver::
createVector() const {
  // Get map data
  const auto map = this->map();
  const auto firstGlobalIndex = map->firstGlobalIndex();
  const auto lastGlobalIndex = map->lastGlobalIndex();

  // For now, require at least one element per processor
  VERIFY(map->numLocalElements() > 0);
  
  int hypreStatus;
  VectorType* tempVector;
  
  hypreStatus = HYPRE_IJVectorCreate(Communicator::communicator(),
                                     firstGlobalIndex,
                                     lastGlobalIndex,
                                     &tempVector);
  VERIFY(hypreStatus == 0);
  
  hypreStatus = HYPRE_IJVectorSetObjectType(tempVector, HYPRE_PARCSR);
  VERIFY(hypreStatus == 0);

  hypreStatus = HYPRE_IJVectorInitialize(tempVector);
  VERIFY(hypreStatus == 0);
  
  return std::shared_ptr<VectorType>(tempVector, HYPRE_IJVectorDestroy);
}

//------------------------------------------------------------------------------
// Create Hypre matrix that should self-destruct when it goes out of scope.
//------------------------------------------------------------------------------
std::shared_ptr<HypreLinearSolver::MatrixType>
HypreLinearSolver::
createMatrix() const {
  // Get map data
  const auto map = this->map();
  const auto firstGlobalIndex = map->firstGlobalIndex();
  const auto lastGlobalIndex = map->lastGlobalIndex();
  const std::vector<int> numElementsPerRow = map->numElementsPerRow();

  int hypreStatus;
  MatrixType* tempMatrix;
  hypreStatus = HYPRE_IJMatrixCreate(Communicator::communicator(),
                                     firstGlobalIndex,
                                     lastGlobalIndex,
                                     firstGlobalIndex,
                                     lastGlobalIndex,
                                     &tempMatrix);
  VERIFY(hypreStatus == 0);

  hypreStatus = HYPRE_IJMatrixSetObjectType(tempMatrix, HYPRE_PARCSR);
  VERIFY(hypreStatus == 0);
  
  hypreStatus = HYPRE_IJMatrixSetRowSizes(tempMatrix, &numElementsPerRow[0]);
  VERIFY(hypreStatus == 0);
  
  hypreStatus = HYPRE_IJMatrixInitialize(tempMatrix);
  VERIFY(hypreStatus == 0);
  
  return std::shared_ptr<MatrixType>(tempMatrix, HYPRE_IJMatrixDestroy);
}

//------------------------------------------------------------------------------
// Create Hypre solver that should self-destruct when it goes out of scope.
// Right now, the constants and type (GMRES) are hardcoded. At the least, the
// kDim should be changed to a variable. 
//------------------------------------------------------------------------------
std::shared_ptr<HypreLinearSolver::SolverType>
HypreLinearSolver::
createSolver() const {
  int hypreStatus;
  SolverType* tempSolver;
  hypreStatus = HYPRE_ParCSRGMRESCreate(Communicator::communicator(), &tempSolver);
  VERIFY(hypreStatus == 0);

  hypreStatus = HYPRE_ParCSRGMRESSetKDim(tempSolver, mHypreOptions->kDim);
  VERIFY(hypreStatus == 0);

  hypreStatus = HYPRE_ParCSRGMRESSetMinIter(tempSolver, mHypreOptions->minIters);
  VERIFY(hypreStatus == 0);

  hypreStatus = HYPRE_ParCSRGMRESSetLogging(tempSolver, mHypreOptions->logLevel);
  VERIFY(hypreStatus == 0);

  hypreStatus = HYPRE_ParCSRGMRESSetPrintLevel(tempSolver, mHypreOptions->printLevel);
  VERIFY(hypreStatus == 0);

  return std::shared_ptr<SolverType>(tempSolver, HYPRE_ParCSRGMRESDestroy);
}

//------------------------------------------------------------------------------
// Create Hypre preconditioner that should self-destruct when it goes out of
// scope.
// The constants and type (AMG) are hardcoded to start out with.
//------------------------------------------------------------------------------
std::shared_ptr<HypreLinearSolver::SolverType>
HypreLinearSolver::
createPreconditioner(std::shared_ptr<SolverType> solver) const {
  CHECK(solver);

  int hypreStatus;
  SolverType* tempPreconditioner;

  switch (mHypreOptions->preconditionerType) {
  case HypreOptions::HyprePreconditionerType::NoPreconditioner:
    return std::shared_ptr<SolverType>();
  case HypreOptions::HyprePreconditionerType::AMGPreconditioner:
    // Begin copied code
    hypreStatus = HYPRE_BoomerAMGCreate(&tempPreconditioner);
    VERIFY(hypreStatus == 0);
    
    hypreStatus = HYPRE_BoomerAMGSetCoarsenType(tempPreconditioner, mHypreOptions->coarsenTypeAMG);
    VERIFY(hypreStatus == 0);

    hypreStatus = HYPRE_BoomerAMGSetMeasureType(tempPreconditioner, mHypreOptions->measure_type);
    VERIFY(hypreStatus == 0);

    hypreStatus = HYPRE_BoomerAMGSetTol(tempPreconditioner, mHypreOptions->pcTol);
    VERIFY(hypreStatus == 0);

    hypreStatus
      = HYPRE_BoomerAMGSetStrongThreshold(tempPreconditioner, mHypreOptions->strongThresholdAMG);
    VERIFY(hypreStatus == 0);

    // Determine when we're a diagonally dominant row. We've modified it.
    hypreStatus = HYPRE_BoomerAMGSetMaxRowSum(tempPreconditioner, mHypreOptions->maxRowSumAMG);
    VERIFY(hypreStatus == 0);

    if (mHypreOptions->interpTypeAMG >= 0) {
      hypreStatus = HYPRE_BoomerAMGSetInterpType(tempPreconditioner, mHypreOptions->interpTypeAMG);
      VERIFY(hypreStatus == 0);
    }

    if (mHypreOptions->aggNumLevelsAMG >= 0) {
      hypreStatus
        = HYPRE_BoomerAMGSetAggNumLevels(tempPreconditioner, mHypreOptions->aggNumLevelsAMG);
      VERIFY(hypreStatus == 0);
    }
    if (mHypreOptions->aggInterpTypeAMG >= 0) {
      hypreStatus
        = HYPRE_BoomerAMGSetAggInterpType(tempPreconditioner, mHypreOptions->aggInterpTypeAMG);
      VERIFY(hypreStatus == 0);
    }
    if (mHypreOptions->pMaxElmtsAMG >= 0){
      hypreStatus = HYPRE_BoomerAMGSetPMaxElmts(tempPreconditioner, mHypreOptions->pMaxElmtsAMG);
      VERIFY(hypreStatus == 0);
    }

    hypreStatus = HYPRE_BoomerAMGSetTruncFactor(tempPreconditioner, mHypreOptions->truncFactorAMG);
    VERIFY(hypreStatus == 0);

    // Go with global value, unless we specified something else.
    if (mHypreOptions->printLevelAMG == -1) {
      mHypreOptions->printLevelAMG = mHypreOptions->printLevel;
    }
    hypreStatus = HYPRE_BoomerAMGSetPrintLevel(tempPreconditioner, mHypreOptions->printLevelAMG);
    VERIFY(hypreStatus == 0);

    hypreStatus = HYPRE_BoomerAMGSetLogging(tempPreconditioner, mHypreOptions->logLevelAMG);
    VERIFY(hypreStatus == 0);

    hypreStatus = HYPRE_BoomerAMGSetMinIter(tempPreconditioner, mHypreOptions->minItersAMG);
    VERIFY(hypreStatus == 0);

    hypreStatus = HYPRE_BoomerAMGSetMaxIter(tempPreconditioner, mHypreOptions->maxItersAMG);
    VERIFY(hypreStatus == 0);

    hypreStatus = HYPRE_BoomerAMGSetCycleType(tempPreconditioner, mHypreOptions->cycleType);
    VERIFY(hypreStatus == 0);

    hypreStatus = HYPRE_BoomerAMGSetRelaxWt(tempPreconditioner, mHypreOptions->relaxWeightAMG);
    VERIFY(hypreStatus == 0);

    if (mHypreOptions->relaxTypeCoarseAMG == -1) {
      mHypreOptions->relaxTypeCoarseAMG = mHypreOptions->relaxTypeAMG;
    }

    // Up sweeps
    hypreStatus
      = HYPRE_BoomerAMGSetCycleRelaxType(tempPreconditioner,
                                         mHypreOptions->relaxTypeAMG, 1);
    VERIFY(hypreStatus == 0);
    // Down sweeps
    hypreStatus
      = HYPRE_BoomerAMGSetCycleRelaxType(tempPreconditioner,
                                         mHypreOptions->relaxTypeAMG, 2);
    VERIFY(hypreStatus == 0);
    // Coarsest level
    hypreStatus
      = HYPRE_BoomerAMGSetCycleRelaxType(tempPreconditioner,
                                         mHypreOptions->relaxTypeCoarseAMG, 3);
    VERIFY(hypreStatus == 0);

    // Up sweeps
    hypreStatus
      = HYPRE_BoomerAMGSetCycleNumSweeps(tempPreconditioner,
                                         mHypreOptions->cycleNumSweepsAMG, 1);
    VERIFY(hypreStatus == 0);
    // down sweeps
    hypreStatus
      = HYPRE_BoomerAMGSetCycleNumSweeps(tempPreconditioner,
                                         mHypreOptions->cycleNumSweepsAMG, 2);
    VERIFY(hypreStatus == 0);
    // coarsest
    if (mHypreOptions->relaxTypeCoarseAMG == 19
        || mHypreOptions->relaxTypeCoarseAMG == 29
        || mHypreOptions->relaxTypeCoarseAMG == 9) {
      mHypreOptions->cycleNumSweepsCoarseAMG = 1;
    }
    else if (mHypreOptions->cycleNumSweepsCoarseAMG == -1) {
      mHypreOptions->cycleNumSweepsCoarseAMG = mHypreOptions->cycleNumSweepsAMG;
    }
    hypreStatus = HYPRE_BoomerAMGSetCycleNumSweeps(tempPreconditioner,
                                                   mHypreOptions->cycleNumSweepsCoarseAMG, 3);
    VERIFY(hypreStatus == 0);

    hypreStatus = HYPRE_BoomerAMGSetMaxLevels(tempPreconditioner, mHypreOptions->maxLevelsAMG);
    VERIFY(hypreStatus == 0);
    // End copied code

    hypreStatus = HYPRE_ParCSRGMRESSetPrecond(solver.get(),
                                              HYPRE_BoomerAMGSolve,
                                              HYPRE_BoomerAMGSetup,
                                              tempPreconditioner);
    VERIFY(hypreStatus == 0);
      
    // Set up preconditioner
    hypreStatus = HYPRE_BoomerAMGSetSetupType(tempPreconditioner, 1);
    VERIFY(hypreStatus == 0);
  
    return std::shared_ptr<SolverType>(tempPreconditioner, HYPRE_BoomerAMGDestroy);
  case HypreOptions::HyprePreconditionerType::ILUPreconditioner:
    hypreStatus = HYPRE_EuclidCreate(Communicator::communicator(), &tempPreconditioner);
    VERIFY(hypreStatus == 0);
    
    hypreStatus = HYPRE_EuclidSetLevel(tempPreconditioner, mHypreOptions->factorLevelILU);
    VERIFY(hypreStatus == 0);

    if (mHypreOptions->printLevelILU == -1) {
      mHypreOptions->printLevelILU = mHypreOptions->printLevel;
    }
    hypreStatus = HYPRE_EuclidSetStats(tempPreconditioner, mHypreOptions->printLevelILU);
    VERIFY(hypreStatus == 0);

    if (mHypreOptions->useILUT) {
      hypreStatus = HYPRE_EuclidSetILUT(tempPreconditioner, mHypreOptions->dropToleranceILU);
    }
    else {
      hypreStatus = HYPRE_EuclidSetSparseA(tempPreconditioner, mHypreOptions->dropToleranceILU);
      VERIFY(hypreStatus == 0);
    }

    hypreStatus = HYPRE_EuclidSetRowScale(tempPreconditioner, mHypreOptions->rowScaleILU);
    VERIFY(hypreStatus == 0);
    
    hypreStatus = HYPRE_ParCSRGMRESSetPrecond(solver.get(),
                                              HYPRE_EuclidSolve,
                                              HYPRE_EuclidSetup,
                                              tempPreconditioner);
    VERIFY(hypreStatus == 0);
    
    return std::shared_ptr<SolverType>(tempPreconditioner, HYPRE_EuclidDestroy);
  default:
    VERIFY2(false, "incorrect preconditioner type");
    return std::shared_ptr<SolverType>(); // So compiler doesn't complain
  }
}

//------------------------------------------------------------------------------
// Fill the given matrix with values from MatrixData.
// Assumes that matrix was created using createMatrix (and is thus properly
// initialized).
//------------------------------------------------------------------------------
void
HypreLinearSolver::
fillMatrix(std::shared_ptr<MatrixType> hypreMatrix) const {
  CHECK(hypreMatrix);
  
  // Get map data
  const auto map = this->map();
  const auto data = this->data();
  const auto numLocalElements = map->numLocalElements();

  std::vector<int> globalColumnIndices;
  std::vector<double> globalColumnValues;
  
  // Fill matrix
  for (int localRowIndex = 0; localRowIndex < numLocalElements; ++localRowIndex) {
    const auto globalRowIndex = map->getGlobalIndex(localRowIndex);
    data->getRowValues(localRowIndex,
                       globalColumnIndices,
                       globalColumnValues);
    int numValues = globalColumnIndices.size();
  
    // Sum duplicate values
    if (mHypreOptions->sumDuplicates) {
      // Get map with all values
      std::map<int, double> values;
      for (auto i = 0; i < numValues; ++i) {
        values[globalColumnIndices[i]] += globalColumnValues[i];
      }
    
      // Convert map back to indices and values
      globalColumnIndices.clear();
      globalColumnValues.clear();
      for (std::map<int, double>::iterator it = values.begin(); it != values.end(); ++it) {
        globalColumnIndices.push_back(it->first);
        globalColumnValues.push_back(it->second);
      }
      numValues = globalColumnIndices.size();

      // Make sure that after summing, there are fewer values than initialized columns
      VERIFY(numValues <= map->numElementsPerRow()[localRowIndex]);
    }
    else {
      // If not summing values, the number of columns specified should be the number we get
      VERIFY(numValues == map->numElementsPerRow()[localRowIndex]);
    }
  
    // Make sure the values and indices have the same size
    VERIFY(globalColumnValues.size() == size_t(numValues));
  
    // Check for NaN values
    BEGIN_CONTRACT_SCOPE
      {
        for (auto i = 0; i < numValues; ++i) {
          CHECK2(!(std::isnan(globalColumnIndices[i]) || std::isnan(globalColumnValues[i])),
                 "NaN found in input to Hypre matrix for (row, col) = (" + std::to_string(globalRowIndex) + std::to_string(globalColumnIndices[i]) + ")");
        }
      }
    END_CONTRACT_SCOPE

      // Set the data
      int hypreStatus;
    if (mHypreOptions->addToValues) {
      hypreStatus = HYPRE_IJMatrixAddToValues(hypreMatrix.get(),
                                              1, // number of rows
                                              &numValues,
                                              &globalRowIndex,
                                              &globalColumnIndices[0],
                                              &globalColumnValues[0]);
    }
    else {
      hypreStatus = HYPRE_IJMatrixSetValues(hypreMatrix.get(),
                                            1, // number of rows
                                            &numValues,
                                            &globalRowIndex,
                                            &globalColumnIndices[0],
                                            &globalColumnValues[0]);
    }
    VERIFY(hypreStatus == 0);
  }
  
  // Assemble matrix
  int hypreStatus;
  hypreStatus = HYPRE_IJMatrixAssemble(hypreMatrix.get());
  VERIFY(hypreStatus == 0);
}

} // end namespace Spheral
