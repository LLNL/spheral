#ifndef SPHERAL_SPHVECTORTEST
#define SPHERAL_SPHVECTORTEST

#include "RAJA/RAJA.hpp"
#include "RAJA/util/Timer.hpp"
#include "chai/ManagedArray.hpp"
#include "umpire/ResourceManager.hpp"
#include "umpire/strategy/QuickPool.hpp"


#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <omp.h>

//#include "NodeList/FluidNodeList.hh"
//#include "Utilities/QuadraticInterpolator.hh"

//*****************************************************************************
// Initialize execution platform
//*****************************************************************************

#define USE_DEVICE 1

#include "ExecStrategy.hh"
#include "memoryManager.hh"
#include "SphVectorField.hh"

//*****************************************************************************
// Set up problem size
//*****************************************************************************

#if 0 // X Large Problem
#define N_PAIRS 500000000
#define DATA_SZ 1000000

#elif 0 // Large Problem
#define N_PAIRS  5000000
#define DATA_SZ  50000

#elif 0 // Medium Problem
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

void SpheralVectorTest()
{
  // Setup Timers
  srand(3);
  RAJA::Timer seq_timer;
  RAJA::Timer launch_timer;
  RAJA::Timer timer_pair;
  RAJA::Timer timer_red;

  using DIM = Spheral::Dim<3>;

  using DATA_TYPE = Spheral::GeomVector<3>;
  //using DATA_TYPE = vec;
  //using DATA_TYPE = double;
  using TRS_UINT = RAJA::TypedRangeSegment<unsigned>;

  using FIELD_TYPE = Spheral::Field<DIM, DATA_TYPE>;
  //using FIELD_TYPE = LvField<DATA_TYPE>;

  using FIELDLIST_TYPE = Spheral::FieldList<DIM, DATA_TYPE>;
  using FIELDLISTVIEW_TYPE = Spheral::FieldListView<DIM, DATA_TYPE>;

  //using VIEW_TYPE = FIELD_TYPE::view_type;

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
  std::cout << "Using Cuda exec policies...\n";
  ExecutionStrategy strat(n_pairs, data_sz, RAJA::Platform::cuda);
  using PAIR_EXEC_POL = RAJA::cuda_exec<1024>;
  using DATA_EXEC_POL = RAJA::loop_exec;
#else
  ExecutionStrategy strat(n_pairs, data_sz, RAJA::Platform::host);
  using PAIR_EXEC_POL = RAJA::loop_exec;
  using DATA_EXEC_POL = RAJA::loop_exec;
  //using PAIR_EXEC_POL = RAJA::omp_parallel_for_exec;
  //using DATA_EXEC_POL = RAJA::omp_parallel_for_exec;
#endif
  strat.print();
  
  //---------------------------------------------------------------------------
  //
  // Initialize Prototype data and allocate temp memory.
  //
  //---------------------------------------------------------------------------

  //SphVector<DATA_TYPE> a(10);
  //a.push_back(DATA_TYPE(5,5,5));
  //a.emplace_back(4,4,4);
  //a.push_back(DATA_TYPE(6,6,6));
  //chai::ManagedArray<DATA_TYPE> b(10, chai::CPU);

  //Spheral::NodeList<DIM> pair_node_list("PairNodeList", n_pairs, 0);
  //Spheral::Field<DIM, unsigned> pair_data("Pairs", pair_node_list);
  //for (unsigned int i = 0; i < n_pairs; i++) pair_data[i] = rand() % DATA_SZ;
  //PRINT_DATA(pair_data, N_PAIRS)

  //SphVector<DATA_TYPE> test_field(4);
  //DATA_TYPE foo(5,5,5);
  //test_field.push_back(foo);
  //test_field.emplace_back(4,4,4);
  //test_field.push_back(DATA_TYPE(6,6,6));

  Spheral::NodeList<DIM> data_node_list_a("DataNodeListA", data_sz, 0);
  //Spheral::NodeList<DIM> data_node_list_b("DataNodeListB", data_sz, 0);
  //Spheral::NodeList<DIM> data_node_list_c("DataNodeListC", data_sz, 0);
  Spheral::NodeList<DIM> data_node_list_bc("DataNodeListBC", data_sz, 0);
  
  FIELD_TYPE a_field("MyFieldA", data_node_list_a);
  //FIELD_TYPE b_field("MyFieldB", data_node_list_b);
  //FIELD_TYPE c_field("MyFieldC", data_node_list_c);
  FIELD_TYPE b_field("MyFieldB", data_node_list_bc);
  FIELD_TYPE c_field("MyFieldC", data_node_list_bc);

  FIELDLIST_TYPE fieldlist;

  fieldlist.appendField(a_field);
  fieldlist.appendField(b_field);
  //fieldlist.appendField(c_field);

  FIELDLIST_TYPE fieldlist2;

  fieldlist2.appendField(c_field);
  fieldlist2.appendField(a_field);
  //fieldlist2.appendField(b_field);

  FIELDLISTVIEW_TYPE fieldlist2_v = fieldlist2.toView();
  FIELDLISTVIEW_TYPE fieldlist_v = fieldlist.toView();

  for(auto& elem : fieldlist_v) std::cout << elem->name() << std::endl;
  for(auto& elem : fieldlist2_v) std::cout << elem->name() << std::endl;
  
  std::cout << "RAJA Teams Implementation Idx.\n";
  launch_timer.start();
  timer_pair.start();

  RAJA::forall<DATA_EXEC_POL>(TRS_UINT(0, data_sz),
    [=] RAJA_HOST (unsigned t_idx) {
      //a_field[t_idx] = DATA_TYPE(1, 1, 1);
      auto& x = fieldlist_v(0, t_idx);
      auto& x2= fieldlist2_v(1, t_idx);
      x = DATA_TYPE(1, 1, 1);
      x2= DATA_TYPE(1, 1, 1);
    });

  std::cout << "check2\n";
  //std::cout << "CHECK : " << fieldlist.name() << "\n";
#if 1
  RAJA::forall<PAIR_EXEC_POL>(TRS_UINT(0, data_sz),
    [=] RAJA_HOST_DEVICE (unsigned t_idx) {

      //a_field[t_idx] += DATA_TYPE(2, 2, 1);
      auto& x = fieldlist_v(0, t_idx);
      auto& x2= fieldlist2_v(1, t_idx);
      x += DATA_TYPE(t_idx, 2, 1);
      x2+= DATA_TYPE(t_idx, 2, 1);

    });
  timer_pair.stop();
  std::cout << "check3\n";

  timer_red.start();
  // We need to perform an array reduction accross the memory pools. This is also performed on the 
  // device to reduce when applicable to minize data movement back to the CPU.
  RAJA::forall<DATA_EXEC_POL>(TRS_UINT(0, data_sz),
    [=] RAJA_HOST (unsigned t_idx) {
      //a_field[t_idx] += DATA_TYPE(1, 1, 1);
      auto& x = fieldlist_v(0, t_idx);
      auto& x2= fieldlist2_v(1, t_idx);
      x += DATA_TYPE(1, 1, 1);
      x2+= DATA_TYPE(1, 1, 1);
      //b[t_idx] += DATA_TYPE(1, 1, 1);
    });


  timer_red.stop();
  launch_timer.stop();


  PRINT_DATA(a_field, a_field.size())
  PRINT_DATA(b_field, b_field.size())
  PRINT_DATA(c_field, c_field.size())
  //PRINT_DATA(a, a.size())
  //PRINT_DATA(b, DATA_SZ)

#endif

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  //DATA_TYPE check_A[DATA_SZ];
  //DATA_TYPE check_B[DATA_SZ];
  //DATA_TYPE check_C[DATA_SZ];

  //seq_timer.start();
  //std::cout << "C++ Sequential Implementation.\n";
  //for(int i = 0; i < N_PAIRS; i++){
  //  DATA_TYPE inc(1.0);
  //  check_A[pairs[i]] += inc;  
  //  check_B[pairs[i]] += inc;  
  //  check_C[pairs[i]] += inc;  
  //}
  //seq_timer.stop();

  //PRINT_DATA(check_A, DATA_SZ)
  //PRINT_DATA(check_B, DATA_SZ)
  //PRINT_DATA(check_C, DATA_SZ)

  ////---------------------------------------------------------------------------
  //
  //bool result = true;
  //for (int i = 0; i<DATA_SZ; i++) {
  //  if (A[i] != check_A[i]) { std::cout << "Failed at A[" << i << "] " << A[i] << ", " << check_A[i] << "\n"; result = false; }
  //  if (B[i] != check_B[i]) { std::cout << "Failed at B[" << i << "] " << B[i] << ", " << check_B[i] << "\n"; result = false; }
  //  if (C[i] != check_C[i]) { std::cout << "Failed at C[" << i << "] " << C[i] << ", " << check_C[i] << "\n"; result = false; }
  //  if (!result) break;
  //}

  //if (result) std::cout << "\n-- PASSED --\n\n";
  //else std::cout << "\n-- FAILED --\n\n";

  //std::cout << "   Seq Impl : " << seq_timer.elapsed() << " seconds\n";
  //std::cout << "Launch Impl : " << launch_timer.elapsed() << " seconds (" << seq_timer.elapsed() / launch_timer.elapsed() << "x) speedup\n";
  //std::cout << "       pair : " << timer_pair.elapsed() << " seconds (" << (timer_pair.elapsed() / launch_timer.elapsed())*100 << "%)\n";
  //std::cout << "        red : " << timer_red.elapsed() << " seconds (" << (timer_red.elapsed() / launch_timer.elapsed())*100 << "%)\n";



}

#endif //  SPHERAL_SPHVECTORTEST
