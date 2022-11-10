#ifndef SPHERAL_SPHEVALDERIVTEST
#define SPHERAL_SPHEVALDERIVTEST

#include "RAJA/RAJA.hpp"
#include "RAJA/util/Timer.hpp"
#include "chai/ManagedArray.hpp"

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <omp.h>

#include "memoryManager.hh"
#include "ExecStrategy.hh"

//*****************************************************************************
//
//*****************************************************************************
#if 0
#define N_PAIRS 500000000
#define DATA_SZ 1000000
//#define BLOCK_SZ 200000
#elif 1
#define N_PAIRS  1000000
#define DATA_SZ  10000
#else
#define N_PAIRS 22
#define DATA_SZ 20
#endif

//*****************************************************************************
#define USE_DEVICE 0
#define USE_CHAI 1

#if USE_DEVICE
#define USE_UNIFIED_MEM
#endif
//*****************************************************************************

#define PRINT_DATA(X, SZ) \
  std::cout << #X << "\t = "; \
  for (int i = 0; i < (SZ < 50 ? SZ : 50); i++) std::cout << X[i] << ", "; \
  std::cout << std::endl;

//*****************************************************************************
void SpheralEvalDerivTest()
{
  srand(3);
  RAJA::Timer seq_timer;
  RAJA::Timer launch_timer;
  RAJA::Timer timer_pair;
  RAJA::Timer timer_red;

  using DATA_TYPE = double;
  using TRS_UINT = RAJA::TypedRangeSegment<unsigned>;

  //---------------------------------------------------------------------------
  unsigned int n_pairs = N_PAIRS;
  unsigned int data_sz = DATA_SZ;

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
  unsigned *pairs = (unsigned*) malloc(sizeof(unsigned) * n_pairs);
  for (unsigned int i = 0; i < n_pairs; i++) pairs[i] = rand() % DATA_SZ;
  //PRINT_DATA(pairs, N_PAIRS)

#if USE_CHAI
  chai::ManagedArray<DATA_TYPE> A(data_sz, chai::CPU);
  chai::ManagedArray<DATA_TYPE> B(data_sz, chai::CPU);
  chai::ManagedArray<DATA_TYPE> C(data_sz, chai::CPU);
#else
  DATA_TYPE A[data_sz] = {0};
  DATA_TYPE B[data_sz] = {0};
  DATA_TYPE C[data_sz] = {0};
#endif

  // Switch between Unified Memory Allocation and Malloc
#ifdef USE_UNIFIED_MEM
  #if USE_CHAI
  chai::ManagedArray<DATA_TYPE> g_A(strat.n_blocks * data_sz, chai::GPU);
  chai::ManagedArray<DATA_TYPE> g_B(strat.n_blocks * data_sz, chai::GPU);
  chai::ManagedArray<DATA_TYPE> g_C(strat.n_blocks * data_sz, chai::GPU);
  #else
  DATA_TYPE* g_A = memoryManager::allocate<DATA_TYPE>(strat.n_blocks * data_sz);
  DATA_TYPE* g_B = memoryManager::allocate<DATA_TYPE>(strat.n_blocks * data_sz);
  DATA_TYPE* g_C = memoryManager::allocate<DATA_TYPE>(strat.n_blocks * data_sz);
  #endif
#else
  DATA_TYPE *g_A = (DATA_TYPE*) malloc(sizeof(DATA_TYPE) * strat.n_blocks * data_sz);
  DATA_TYPE *g_B = (DATA_TYPE*) malloc(sizeof(DATA_TYPE) * strat.n_blocks * data_sz);
  DATA_TYPE *g_C = (DATA_TYPE*) malloc(sizeof(DATA_TYPE) * strat.n_blocks * data_sz);
#endif

  //---------------------------------------------------------------------------
//#ifdef ENABLE_CUDA 
#if 1
  using pair_pol = RAJA::LoopPolicy<RAJA::omp_for_exec, RAJA::cuda_global_thread_x>;
  using data_pol = RAJA::LoopPolicy<RAJA::omp_for_exec, RAJA::cuda_global_thread_x>;
  using loop_pol = RAJA::LoopPolicy<RAJA::loop_exec, RAJA::cuda_global_thread_x>;
  using launch_policy = RAJA::LaunchPolicy<RAJA::omp_launch_t, RAJA::cuda_launch_t<false>>;
  auto lp = RAJA::LaunchParams(RAJA::Teams(strat.n_blocks), RAJA::Threads(strat.block_sz));
#else
  using pair_pol     = RAJA::LoopPolicy<RAJA::omp_for_exec>;
  using data_pol      = RAJA::LoopPolicy<RAJA::omp_for_exec>;
  using loop_pol      = RAJA::LoopPolicy<RAJA::loop_exec>;
  using launch_policy = RAJA::LaunchPolicy<RAJA::omp_launch_t>;
  auto lp = RAJA::LaunchParams();
#endif

  RAJA::View<DATA_TYPE, RAJA::Layout<1>> Av(A, data_sz);
  RAJA::View<DATA_TYPE, RAJA::Layout<1>> Bv(B, data_sz);
  RAJA::View<DATA_TYPE, RAJA::Layout<1>> Cv(C, data_sz);

  RAJA::View<DATA_TYPE, RAJA::Layout<2>> g_Av(g_A, strat.n_blocks, data_sz);
  RAJA::View<DATA_TYPE, RAJA::Layout<2>> g_Bv(g_B, strat.n_blocks, data_sz);
  RAJA::View<DATA_TYPE, RAJA::Layout<2>> g_Cv(g_C, strat.n_blocks, data_sz);

#if 1
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  std::cout << "RAJA Teams Implementation Idx.\n";
  if (strat.sequential) std::cout << "Running Sequential\n";
  launch_timer.start();

  if (strat.sequential)
  {

    RAJA::launch<launch_policy>(strat.exec_context, lp,
      [=] RAJA_HOST_DEVICE (RAJA::LaunchContext ctx) {
      RAJA::loop<pair_pol>(ctx, TRS_UINT(0, strat.n_pairs), [&] (unsigned t_idx) {
        auto pair_idx = pairs[t_idx];
        Av(pair_idx) += 1.0;
        Bv(pair_idx) += 1.0;
        Cv(pair_idx) += 1.0;
      });
    });

  }else{
    timer_pair.start();
    //RAJA::launch<launch_policy>(strat.exec_context, lp,
    //  [=] RAJA_HOST_DEVICE (RAJA::LaunchContext ctx) {
    //  RAJA::loop<pair_pol>(ctx, TRS_UINT(0, strat.n_pairs), [&] (unsigned t_idx) {
    RAJA::forall<PAIR_EXEC_POL>(TRS_UINT(0, strat.n_pairs), [=] RAJA_HOST_DEVICE (unsigned t_idx) {
        auto pair_idx = pairs[t_idx];
        unsigned b_idx = t_idx / strat.block_sz;
        auto g_idx = b_idx*data_sz + pair_idx;

        // Device
#if USE_CHAI
  #if USE_DEVICE
        ATOMIC_ADD(&g_A[g_idx], 1.0);
        ATOMIC_ADD(&g_B[g_idx], 1.0);
        ATOMIC_ADD(&g_C[g_idx], 1.0);
  #else
        g_A[g_idx] += 1.0;
        g_B[g_idx] += 1.0;
        g_C[g_idx] += 1.0;
  #endif
#else
        ATOMIC_ADD(&g_Av(b_idx, pair_idx), 1.0);
        ATOMIC_ADD(&g_Bv(b_idx, pair_idx), 1.0);
        ATOMIC_ADD(&g_Cv(b_idx, pair_idx), 1.0);
#endif
    //  });
    });
    timer_pair.stop();

    timer_red.start();
    //RAJA::launch<launch_policy>(lp,
    //  [=] RAJA_HOST_DEVICE (RAJA::LaunchContext ctx) {
    //  RAJA::loop<data_pol>(ctx, TRS_UINT(0, data_sz), [&] (unsigned t_idx) {
    //    RAJA::loop<loop_pol>(ctx, TRS_UINT(0, strat.n_blocks), [&] (unsigned b_idx) {
    RAJA::forall<DATA_EXEC_POL>(TRS_UINT(0, data_sz), [=] RAJA_HOST_DEVICE (unsigned t_idx) {
        for (int b_idx = 0; b_idx < strat.n_blocks; b_idx++) {
          auto g_idx = b_idx*data_sz + t_idx;
          //printf("%d, ", g_idx);
  #if USE_CHAI
          A[t_idx] += g_A[g_idx];
          B[t_idx] += g_B[g_idx];
          C[t_idx] += g_C[g_idx];
  #else
          Av(t_idx) += g_Av(b_idx, t_idx);
          Bv(t_idx) += g_Bv(b_idx, t_idx);
          Cv(t_idx) += g_Cv(b_idx, t_idx);
  #endif
        }

    //    });
    //  });
    });
    timer_red.stop();

  }

  launch_timer.stop();

  PRINT_DATA(g_A, DATA_SZ * 2);
  PRINT_DATA(g_B, DATA_SZ * 2);
  PRINT_DATA(g_C, DATA_SZ * 2);

  PRINT_DATA(A, DATA_SZ)
  PRINT_DATA(B, DATA_SZ)
  PRINT_DATA(C, DATA_SZ)

#endif

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  DATA_TYPE check_A[data_sz] = {0};
  DATA_TYPE check_B[data_sz] = {0};
  DATA_TYPE check_C[data_sz] = {0};

  seq_timer.start();
  std::cout << "C++ Sequential Implementation.\n";
  for(int i = 0; i < N_PAIRS; i++){
    check_A[pairs[i]] += 1.0;  
    check_B[pairs[i]] += 1.0;  
    check_C[pairs[i]] += 1.0;  
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


#ifdef USE_UNIFIED_MEM
  #if USE_CHAI
  g_A.free();
  g_B.free();
  g_C.free();
  #else
  memoryManager::deallocate(g_A);
  memoryManager::deallocate(g_B);
  memoryManager::deallocate(g_C);
  #endif
#else
  free(g_A);
  free(g_B);
  free(g_C);
#endif


}

#endif //  SPHERAL_SPHEVALDERIVTEST
