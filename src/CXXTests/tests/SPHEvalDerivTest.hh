#ifndef SPHERAL_SPHEVALDERIVTEST
#define SPHERAL_SPHEVALDERIVTEST

#include "RAJA/RAJA.hpp"
#include "RAJA/util/Timer.hpp"
#include "chai/ManagedArray.hpp"
#include "umpire/ResourceManager.hpp"
#include "umpire/strategy/QuickPool.hpp"

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <omp.h>

#include "memoryManager.hh"
#include "ExecStrategy.hh"

//*****************************************************************************
// Set up problem size
//*****************************************************************************

#if 0 // X Large Problem
#define N_PAIRS 500000000
#define DATA_SZ 1000000

#elif 0 // Large Problem
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
// Initialize execution platform
//*****************************************************************************
#define USE_DEVICE 0

//*****************************************************************************
#define PRINT_DATA(X, SZ) \
  std::cout << #X << "\t = "; \
  for (int i = 0; i < (SZ < 50 ? SZ : 50); i++) std::cout << X[i] << ", "; \
  std::cout << std::endl;

//*****************************************************************************
void SpheralEvalDerivTest()
{
  // Setup Timers
  srand(3);
  RAJA::Timer seq_timer;
  RAJA::Timer launch_timer;
  RAJA::Timer timer_pair;
  RAJA::Timer timer_red;

  // Data Types to Use
  //using DATA_TYPE = Spheral::GeomVector<3>;
  using DATA_TYPE = double;
  using TRS_UINT = RAJA::TypedRangeSegment<unsigned>;

  //---------------------------------------------------------------------------
  //
  // Setup Eval Derivs Execution.
  //
  //---------------------------------------------------------------------------

  // These are going to be runtime mutable irl...
  unsigned int n_pairs = N_PAIRS;
  unsigned int data_sz = DATA_SZ;

  // Initialize an execution strategy for EvalDerivs problem
  // We need to: 
  //  - define how many copies of temp memory we need
  //  - calculate the blocks / block sizes depeding on execution space.
  //  - define execution policies for CPU and GPU execution.
  // TODO: Make this a runtime switch.
#if USE_DEVICE
  ExecutionStrategy strat(n_pairs, data_sz, RAJA::ExecPlace::DEVICE);
  strat.block_sz = 1024;
  strat.n_blocks = RAJA_DIVIDE_CEILING_INT(strat.n_pairs, strat.block_sz);
  using PAIR_EXEC_POL = RAJA::cuda_exec<1024>;
  using DATA_EXEC_POL = RAJA::cuda_exec<1024>;
#else
  ExecutionStrategy strat(n_pairs, data_sz, RAJA::ExecPlace::HOST);
  strat.block_sz = n_pairs / omp_get_max_threads();
  strat.n_blocks = RAJA_DIVIDE_CEILING_INT(strat.n_pairs, strat.block_sz);
  using PAIR_EXEC_POL = RAJA::omp_parallel_for_exec;
  using DATA_EXEC_POL = RAJA::omp_parallel_for_exec;
#endif
  strat.print();

  //---------------------------------------------------------------------------
  //
  // Initialize Prototype data and allocate temp memory.
  //
  //---------------------------------------------------------------------------
  
  chai::ManagedArray<unsigned> pairs(n_pairs);
  
  for (unsigned int i = 0; i < n_pairs; i++) pairs[i] = rand() % DATA_SZ;
  pairs.registerTouch(chai::CPU);
  PRINT_DATA(pairs, N_PAIRS)

  chai::ManagedArray<DATA_TYPE> A(data_sz);
  chai::ManagedArray<DATA_TYPE> B(data_sz);
  chai::ManagedArray<DATA_TYPE> C(data_sz);

#if USE_DEVICE
  chai::ManagedArray<DATA_TYPE> g_A(strat.n_blocks * data_sz, chai::GPU);
  chai::ManagedArray<DATA_TYPE> g_B(strat.n_blocks * data_sz, chai::GPU);
  chai::ManagedArray<DATA_TYPE> g_C(strat.n_blocks * data_sz, chai::GPU);
  
  pairs.move(chai::GPU);
  A.move(chai::GPU);
  B.move(chai::GPU);
  C.move(chai::GPU);
#else
  chai::ManagedArray<DATA_TYPE> g_A(strat.n_blocks * data_sz);
  chai::ManagedArray<DATA_TYPE> g_B(strat.n_blocks * data_sz);
  chai::ManagedArray<DATA_TYPE> g_C(strat.n_blocks * data_sz);
#endif

  std::cout << "Test\n";

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

  // Initial Pair loop: We are using RAJA::forall over teams for chai::ManagedArray support.
  RAJA::forall<PAIR_EXEC_POL>(TRS_UINT(0, strat.n_pairs),
    [=] RAJA_HOST_DEVICE (unsigned t_idx) {

      auto pair_idx = pairs[t_idx];
      auto b_idx = t_idx / strat.block_sz;
      auto g_idx = b_idx*data_sz + pair_idx;

#if USE_DEVICE
      // We use atomics for Device code as blocks per memory pool is greater than 1.
      // g_X arrays are implicitly copied to the device with chai.
      ATOMIC_ADD(&g_A[g_idx], DATA_TYPE(1.0));
      ATOMIC_ADD(&g_B[g_idx], DATA_TYPE(1.0));
      ATOMIC_ADD(&g_C[g_idx], DATA_TYPE(1.0));
#else
      // When executing on host we create one memory pool per omp thread, atomics are not needed.
      g_A[g_idx] += DATA_TYPE(1.0);
      g_B[g_idx] += DATA_TYPE(1.0);
      g_C[g_idx] += DATA_TYPE(1.0);
#endif
    });
  timer_pair.stop();

  timer_red.start();
  // We need to perform an array reduction accross the memory pools. This is also performed on the 
  // device to reduce when applicable to minize data movement back to the CPU.
  RAJA::forall<DATA_EXEC_POL>(TRS_UINT(0, data_sz),
    [=] RAJA_HOST_DEVICE (unsigned t_idx) {
      for (int b_idx = 0; b_idx < strat.n_blocks; b_idx++) {
        auto g_idx = b_idx*data_sz + t_idx;
        //printf("%d, ", g_idx);
        A[t_idx] += g_A[g_idx];
        B[t_idx] += g_B[g_idx];
        C[t_idx] += g_C[g_idx];
      }
    });

  A.move(chai::CPU);
  B.move(chai::CPU);
  C.move(chai::CPU);

  timer_red.stop();
  launch_timer.stop();

  PRINT_DATA(A, DATA_SZ)
  PRINT_DATA(B, DATA_SZ)
  PRINT_DATA(C, DATA_SZ)

#endif

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  pairs.move(chai::CPU);

  DATA_TYPE check_A[data_sz] = {0};
  DATA_TYPE check_B[data_sz] = {0};
  DATA_TYPE check_C[data_sz] = {0};

  seq_timer.start();
  std::cout << "C++ Sequential Implementation.\n";
  for(int i = 0; i < N_PAIRS; i++){
    check_A[pairs[i]] += DATA_TYPE(1.0);  
    check_B[pairs[i]] += DATA_TYPE(1.0);  
    check_C[pairs[i]] += DATA_TYPE(1.0);  
  }
  seq_timer.stop();

  PRINT_DATA(check_A, DATA_SZ)
  PRINT_DATA(check_B, DATA_SZ)
  PRINT_DATA(check_C, DATA_SZ)

  //---------------------------------------------------------------------------
  
  bool result = true;
  for (int i = 0; i<DATA_SZ; i++) {
    if (A[i] != check_A[i]) { std::cout << "Failed at A[" << i << "] " << A[i] << ", " << check_A[i] << "\n"; result = false; }
    if (B[i] != check_B[i]) { std::cout << "Failed at B[" << i << "] " << B[i] << ", " << check_B[i] << "\n"; result = false; }
    if (C[i] != check_C[i]) { std::cout << "Failed at C[" << i << "] " << C[i] << ", " << check_C[i] << "\n"; result = false; }
    if (!result) break;
  }

  if (result) std::cout << "\n-- PASSED --\n\n";
  else std::cout << "\n-- FAILED --\n\n";

  std::cout << "   Seq Impl : " << seq_timer.elapsed() << " seconds\n";
  std::cout << "Launch Impl : " << launch_timer.elapsed() << " seconds (" << seq_timer.elapsed() / launch_timer.elapsed() << "x) speedup\n";
  std::cout << "       pair : " << timer_pair.elapsed() << " seconds (" << (timer_pair.elapsed() / launch_timer.elapsed())*100 << "%)\n";
  std::cout << "        red : " << timer_red.elapsed() << " seconds (" << (timer_red.elapsed() / launch_timer.elapsed())*100 << "%)\n";


  g_A.free();
  g_B.free();
  g_C.free();

  A.free();
  B.free();
  C.free();

  pairs.free();
}

#endif //  SPHERAL_SPHEVALDERIVTEST
