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

#include "ExecStrategy.hh"
#include "NodeList/FluidNodeList.hh"

#include "memoryManager.hh"
#include "LvField.hh"

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
#define USE_DEVICE 1

//*****************************************************************************
#define PRINT_DATA(X, SZ) \
  std::cout << #X << "\t = "; \
  for (int i = 0; i < (SZ < 50 ? SZ : 50); i++) std::cout << X[i] << ", "; \
  std::cout << std::endl;

//*****************************************************************************

#define VEC3 1
#include "vec.hh"

void SpheralEvalDerivTest()
{
  
  // Setup Timers
  srand(3);
  RAJA::Timer seq_timer;
  RAJA::Timer launch_timer;
  RAJA::Timer timer_pair;
  RAJA::Timer timer_red;

  using DATA_TYPE = Spheral::GeomVector<3>;
  //using DATA_TYPE = vec;
  //using DATA_TYPE = double;
  using TRS_UINT = RAJA::TypedRangeSegment<unsigned>;

  using FIELD_TYPE = LvField<DATA_TYPE>;
  using VIEW_TYPE = FIELD_TYPE::VIEW_TYPE;

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
#define MEM_SPACE RAJA::Platform::cuda
  ExecutionStrategy strat(n_pairs, data_sz, RAJA::ExecPlace::DEVICE);
  strat.block_sz = 1024;
  strat.n_blocks = RAJA_DIVIDE_CEILING_INT(strat.n_pairs, strat.block_sz);
  using PAIR_EXEC_POL = RAJA::cuda_exec<1024>;
  using DATA_EXEC_POL = RAJA::cuda_exec<1024>;
#else
#define MEM_SPACE RAJA::Platform::host
  ExecutionStrategy strat(n_pairs, data_sz, RAJA::ExecPlace::HOST);
  strat.block_sz = n_pairs / omp_get_max_threads();
  strat.n_blocks = RAJA_DIVIDE_CEILING_INT(strat.n_pairs, strat.block_sz);
  using PAIR_EXEC_POL = RAJA::omp_parallel_for_exec;
  using DATA_EXEC_POL = RAJA::omp_parallel_for_exec;
#endif
#define HOST_SPACE RAJA::Platform::host
  strat.print();
  
  //---------------------------------------------------------------------------
  //
  // Initialize Prototype data and allocate temp memory.
  //
  //---------------------------------------------------------------------------
  
  // Generate pair data...
  LvField<unsigned> pair_data(n_pairs, "Pairs");
  for (unsigned int i = 0; i < n_pairs; i++) pair_data[i] = rand() % DATA_SZ;
  PRINT_DATA(pair_data, N_PAIRS)
  const LvFieldView<unsigned> pairs(pair_data);
  pairs.move(MEM_SPACE);

  // Setting up our "Field Data", this is done through simulation setup in spheral e.g. node generation.
  FIELD_TYPE A(data_sz, "A");
  FIELD_TYPE B(data_sz, "B");
  FIELD_TYPE C(data_sz, "C");

  FIELD_TYPE One(data_sz, "One");
  for (size_t i = 0; i < data_sz; i++) One[i] = DATA_TYPE(1.0);
  std::cout << "one : " << One[0] << std::endl;

  // Creating "FieldLists" In evalDerivs we call STATE_TYPE::fields(...) to return a fieldList.
  // In evalderivs we will want to return Something like LvFieldListView types for RAJA lambda 
  // capture and chai data migration. 
  LvFieldList<DATA_TYPE> fl("MyFirstFieldList");
  LvFieldList<DATA_TYPE> fl2("MySecondFieldList");

  // Setting up global device pool memory for each Field...
  auto g_A = A.make_pool_field(strat.n_blocks, MEM_SPACE);
  auto g_B = B.make_pool_field(strat.n_blocks, MEM_SPACE);
  auto g_C = C.make_pool_field(strat.n_blocks, MEM_SPACE);

  // Wrap the fields with their pools this isn't entirely necessary.
  auto Av = A.toViewWithPool(g_A);
  auto Bv = B.toViewWithPool(g_B);
  auto Cv = C.toViewWithPool(g_C);

  fl.appendField(Av);
  fl.appendField(Bv);
  fl2.appendField(Cv);
  fl2.appendField(Av);

  LvFieldList<DATA_TYPE> flo("const FL One");
  flo.appendField(One);

  // The FieldList types used in evalderivs.
  LvFieldListView<DATA_TYPE> flv(fl);
  LvFieldListView<DATA_TYPE> flv2(fl2);
  
  const LvFieldListView<DATA_TYPE> fl_one(flo);

  flv.move(MEM_SPACE);
  flv2.move(MEM_SPACE);
  fl_one.move(MEM_SPACE);

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

  RAJA::forall<PAIR_EXEC_POL>(TRS_UINT(0, strat.n_pairs),
    [=] RAJA_HOST_DEVICE (unsigned t_idx) {

      auto pair_idx = pairs[t_idx];
      auto b_idx = t_idx / strat.block_sz;
      auto g_pair_idx = b_idx*data_sz + pair_idx;

      const auto& one = fl_one(0, pair_idx);

      auto& a =  flv.pool_atomic(0, g_pair_idx);
      auto& b =  flv.pool_atomic(1, g_pair_idx);
      auto& c = flv2.pool_atomic(0, g_pair_idx);

      a += one;
      b += one;
      c += one;

    });
  timer_pair.stop();

  timer_red.start();
  // We need to perform an array reduction accross the memory pools. This is also performed on the 
  // device to reduce when applicable to minize data movement back to the CPU.
  RAJA::forall<DATA_EXEC_POL>(TRS_UINT(0, data_sz),
    [=] RAJA_HOST_DEVICE (unsigned t_idx) {
      for (int b_idx = 0; b_idx < strat.n_blocks; b_idx++) {
        auto g_idx = b_idx*data_sz + t_idx;

        auto& a = flv(0, t_idx);
        auto& b = flv(1, t_idx);
        auto& c = flv2(0, t_idx);

        a +=  flv.pool(0, g_idx);
        b +=  flv.pool(1, g_idx);
        c += flv2.pool(0, g_idx);
      }
    });

  flv.move(HOST_SPACE);
  flv2.move(HOST_SPACE);
  pairs.move(HOST_SPACE);

  timer_red.stop();
  launch_timer.stop();

  PRINT_DATA(A, DATA_SZ)
  PRINT_DATA(B, DATA_SZ)
  PRINT_DATA(C, DATA_SZ)

#endif

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  DATA_TYPE check_A[data_sz];
  DATA_TYPE check_B[data_sz];
  DATA_TYPE check_C[data_sz];

  seq_timer.start();
  std::cout << "C++ Sequential Implementation.\n";
  for(int i = 0; i < N_PAIRS; i++){
    DATA_TYPE inc(1.0);
    if(i == 0) std::cout << "inc : " << inc << std::endl;
    check_A[pairs[i]] += inc;  
    check_B[pairs[i]] += inc;  
    check_C[pairs[i]] += inc;  
    //check_A[pairs[i]] += DATA_TYPE(1.0);  
    //check_B[pairs[i]] += DATA_TYPE(1.0);  
    //check_C[pairs[i]] += DATA_TYPE(1.0);  
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

}

#endif //  SPHERAL_SPHEVALDERIVTEST
