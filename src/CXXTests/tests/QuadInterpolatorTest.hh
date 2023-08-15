#ifndef SPHERAL_QUAD_INTERP_TEST
#define SPHERAL_QUAD_INTERP_TEST

#define CHAI_DEBUG

#include "RAJA/RAJA.hpp"
#include "RAJA/util/Timer.hpp"
#include "chai/ManagedArray.hpp"
#include "umpire/ResourceManager.hpp"
#include "umpire/strategy/QuickPool.hpp"

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <omp.h>

#include "NodeList/FluidNodeList.hh"
#include "Utilities/QuadraticInterpolator.hh"
#include "Kernel/GaussianKernel.hh"

//*****************************************************************************
// Initialize execution platform
//*****************************************************************************

#define USE_DEVICE 1

#include "ExecStrategy.hh"
#include "memoryManager.hh"
#include "LvField.hh"

#include "SphArray.hh"

//*****************************************************************************
// Set up problem size
//*****************************************************************************

#if 0 // X Large Problem
#define N_PAIRS 500000000
#define DATA_SZ 1000000

#elif 1 // Large Problem
#define N_PAIRS  5000000
#define DATA_SZ  50000

#elif 1 // Medium Problem
#define N_PAIRS  1000000
#define DATA_SZ  10000

#else // Small Problem
#define N_PAIRS 2000
#define DATA_SZ 200
#endif


//*****************************************************************************
#define PRINT_DATA(X, SZ) \
  std::cout << #X << "\t = "; \
  for (int i = 0; i < (SZ < 50 ? SZ : 50); i++) std::cout << X[i] << ", "; \
  std::cout << std::endl;

//*****************************************************************************

//#define VEC3 1
//#include "vec.hh"

void myUserCallback(const chai::PointerRecord*, chai::Action action, chai::ExecutionSpace) {
  switch(action){
    case chai::Action::ACTION_ALLOC:
      printf("Alloc\n");
      break;
    case chai::Action::ACTION_FREE:
      printf("Free\n");
      break;
    case chai::Action::ACTION_MOVE:
      printf("Move\n");
      break;
  }

}

void QuadInterpolatorTest()
{
  
  // Setup Timers
  srand(3);
  RAJA::Timer seq_timer;
  RAJA::Timer launch_timer;
  RAJA::Timer timer_pair;
  RAJA::Timer timer_red;

  using DIM = Spheral::Dim<3>;

  using DATA_TYPE = double;
  //using DATA_TYPE = Spheral::GeomVector<3>;
  using TRS_UINT = RAJA::TypedRangeSegment<unsigned>;
  using FIELD_TYPE = Spheral::Field<DIM, DATA_TYPE>;

  //---------------------------------------------------------------------------
  //
  // Setup Eval Derivs Execution.
  //
  //---------------------------------------------------------------------------
  unsigned int n_pairs = N_PAIRS;
  unsigned int data_sz = DATA_SZ;

  // Initialize an execution strategy for EvalDerivs problem
  // TODO: Make this a runtime switch.
#if USE_DEVICE && defined(RAJA_ENABLE_CUDA)
  ExecutionStrategy strat(n_pairs, data_sz, RAJA::Platform::cuda);
  using PAIR_EXEC_POL = RAJA::cuda_exec<1024>;
  using DATA_EXEC_POL = RAJA::loop_exec;
#else
  ExecutionStrategy strat(n_pairs, data_sz, RAJA::Platform::host);
  using PAIR_EXEC_POL = RAJA::loop_exec;
  using DATA_EXEC_POL = RAJA::loop_exec;
#endif
  strat.print();
  
  //---------------------------------------------------------------------------
  //
  // Initialize Prototype data and allocate temp memory.
  //
  //---------------------------------------------------------------------------
  
  chai::ArrayManager::getInstance()->setGlobalUserCallback(myUserCallback);

#if 1
  auto kernel = Spheral::GaussianKernel<DIM>(5.0);
  auto interp = Spheral::QuadraticInterpolator(0.0, kernel.kernelExtent(), 20, [&](const double x){return kernel(x, 1.0);});
#else
  std::vector<double> vals = {0, 2, 3, 5, 8, 11, 19, 30, 49, 79, 128, 209};
  auto interp = ArrayContainerType(vals);
#endif

  data_sz = interp.size();

  Spheral::NodeList<DIM> data_node_list("DataNodeList", data_sz, 0);

  FIELD_TYPE One("One", data_node_list);
  for (size_t i = 0; i < data_sz; i++) One[i] = DATA_TYPE(1.0);
  auto test_ar_v = One.toView(); 

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  //
  // Eval Deriv Problem...
  //
  //---------------------------------------------------------------------------
#if 1
  std::cout << "RAJA Teams Implementation Idx.\n";
  launch_timer.start();
  timer_pair.start();

  bool result = true;


  RAJA::forall<PAIR_EXEC_POL>(TRS_UINT(0, data_sz),
    [=] RAJA_HOST_DEVICE (unsigned t_idx) {
      test_ar_v[t_idx] = interp.coeffs()[t_idx] * 2;
      //test_ar_v[t_idx] = interp[t_idx] * 2;
    });
  timer_pair.stop();

  timer_red.start();
  RAJA::forall<DATA_EXEC_POL>(TRS_UINT(0, data_sz),
    [=, &result] (unsigned t_idx) {
      if (test_ar_v[t_idx] != interp.coeffs()[t_idx]*2) result = false;
      std::cout << t_idx <<" : "<< test_ar_v[t_idx] << "\n";
    });

  timer_red.stop();
  launch_timer.stop();

#endif
  //---------------------------------------------------------------------------

  if (result) std::cout << "\n-- PASSED --\n\n";
  else std::cout << "\n-- FAILED --\n\n";

}

#endif //  SPHERAL_QUAD_INTERP_TEST
