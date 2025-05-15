include_directories(.)
set(Solvers_inst
  NodeMap
  OverlapNodeMap
  )

set(Solvers_sources
  EigenLinearSolver.cc
  HypreLinearSolver.cc
  KINSOL.cc
  LinearSolver.cc
  MatrixMap.cc
  SimpleMatrixData.cc
  SimpleMatrixMap.cc
)

instantiate(Solvers_inst Solvers_sources)

set(Solvers_headers
  containsConstantNodes.hh
  containsConstantNodesInline.hh
  EigenLinearSolver.hh
  EigenLinearSolverInline.hh
  EigenOptions.hh
  HypreLinearSolver.hh
  HypreLinearSolverInline.hh
  HypreOptions.hh
  IncrementalStatistic.hh
  IncrementalStatisticInline.hh
  KINSOL.hh
  LinearSolver.hh
  LinearSolverInline.hh
  MatrixData.hh
  MatrixDataInline.hh
  MatrixMap.hh
  NodeMap.hh
  OverlapNodeMap.hh
  SimpleMatrixData.hh
  SimpleMatrixMap.hh
  SolverFunction.hh
)

spheral_add_obj_library(Solvers SPHERAL_OBJ_LIBS)
