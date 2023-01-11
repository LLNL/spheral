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
#include "NodeList/FluidNodeList.hh"
#include "Material/GammaLawGas.hh"

#include "LvArray/Array.hpp"
#include "LvArray/ArrayOfArrays.hpp"
#include "LvArray/ChaiBuffer.hpp"
#include "LvArray/bufferManipulation.hpp"


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
// Initialize execution platform
//*****************************************************************************
#define USE_DEVICE 1

//*****************************************************************************
#define PRINT_DATA(X, SZ) \
  std::cout << #X << "\t = "; \
  for (int i = 0; i < (SZ < 50 ? SZ : 50); i++) std::cout << X[i] << ", "; \
  std::cout << std::endl;

//*****************************************************************************


class pdouble {
public:
  using DATA_TYPE = double;

  pdouble() : data(0) {};
  pdouble(DATA_TYPE d) : data(d) {};

  std::string string() const {return std::to_string(data);}

  DATA_TYPE get_data() const {return data;}


  pdouble& operator+=(const pdouble& rhs) {data += rhs.get_data(); return *this;}
  friend pdouble operator+(pdouble lhs, const pdouble& rhs) {lhs+=rhs; return lhs;}
  bool operator!=(const pdouble& rhs) const {return !(data == rhs.get_data());}

private:
    DATA_TYPE data = 0;
};




std::ostream& operator<<(std::ostream& out, const pdouble& d) {
  return out << d.string();
}

template<typename DATA_TYPE>
using LV_ARRAY_CHAI_1D = LvArray::Array<DATA_TYPE, 1, camp::idx_seq<0>, std::ptrdiff_t, LvArray::ChaiBuffer>;
template<typename DATA_TYPE>
using LV_ARRAY_VIEW_CHAI_1D = LvArray::ArrayView<DATA_TYPE, 1, 0, std::ptrdiff_t, LvArray::ChaiBuffer>;

template<typename DATA_TYPE>
class LvFieldBaseView;

template<typename DATA_TYPE>
class LvFieldListView;

template<typename DATA_TYPE>
class LvFieldBase {
public:

  friend class LvFieldBaseView<DATA_TYPE>;

  using ARRAY_TYPE = LV_ARRAY_CHAI_1D<DATA_TYPE>;
  using VIEW_TYPE = LvFieldBaseView<DATA_TYPE>;

  LvFieldBase(size_t elems) { mDataArray = ARRAY_TYPE(elems); }

  LvFieldBase(size_t elems, RAJA::Platform platform) { 
    mDataArray = ARRAY_TYPE(elems);
    mDataArray.move(platform);
  }

  LvFieldBase(size_t elems, std::string name) { 
    mDataArray = ARRAY_TYPE(elems);
    mDataArray.setName(name);
  }

  LvFieldBase(size_t elems, std::string name, RAJA::Platform platform) { 
    mDataArray = ARRAY_TYPE(elems);
    mDataArray.setName(name);
    mDataArray.move(platform);
  }

  void setName(std::string name) {mDataArray.setName(name);}
  void move (RAJA::Platform platform) {mDataArray.move(platform);}

  // Index operator.
  RAJA_HOST_DEVICE
  DATA_TYPE& operator[](const unsigned int index) {return mDataArray[index];}
  RAJA_HOST_DEVICE
  DATA_TYPE& operator[](const unsigned int index) const {return mDataArray[index];}

  ARRAY_TYPE& getArray() {return mDataArray;}

private:
  ARRAY_TYPE mDataArray;

};

template<typename DATA_TYPE>
class LvFieldBaseView {
public:
  using ARRAY_VIEW_TYPE = LV_ARRAY_VIEW_CHAI_1D<DATA_TYPE>;
  using FIELD_TYPE = LvFieldBase<DATA_TYPE>;

  friend class LvFieldListView<DATA_TYPE>;

  LvFieldBaseView(const FIELD_TYPE& field) : mDataView{field.mDataArray} {}

  void move(LvArray::MemorySpace const& space, bool touch = true) const {
    mDataView.move(space,touch);
  }

  ARRAY_VIEW_TYPE& getView() {return mDataView;}

  RAJA_HOST_DEVICE
  DATA_TYPE& operator[](const unsigned int index) {return mDataView(index);}

  RAJA_HOST_DEVICE
  DATA_TYPE& operator[](const unsigned int index) const {return mDataView(index);}

  RAJA_HOST_DEVICE
  inline constexpr LvFieldBaseView( LvFieldBaseView const & source) noexcept : mDataView{source.mDataView} {}

private:
  ARRAY_VIEW_TYPE mDataView;
};


template<typename DATA_TYPE>
class LvFieldList{
public:

  friend class LvFieldListView<DATA_TYPE>;

  using FIELD_VIEW_TYPE = LvFieldBaseView<DATA_TYPE>;
  using ARRAY_TYPE = LV_ARRAY_CHAI_1D<FIELD_VIEW_TYPE>;

  LvFieldList() {}

  LvFieldList(std::string name) {
    mFieldArray.setName(name);
  }

  void appendField(const FIELD_VIEW_TYPE& field) { 
    mFieldArray.emplace_back(field);
  }

private:
  ARRAY_TYPE mFieldArray;
};


template<typename DATA_TYPE>
class LvFieldListView{
public:

  using FIELD_VIEW_TYPE = LvFieldBaseView<DATA_TYPE>;
  using ARRAY_VIEW_TYPE = LV_ARRAY_VIEW_CHAI_1D<FIELD_VIEW_TYPE>;

  RAJA_HOST_DEVICE
  LvFieldListView(const LvFieldList<DATA_TYPE>& field) : mFieldView{field.mFieldArray.toView()} {}

  void move(LvArray::MemorySpace const& space, bool touch = true) const {
    mFieldView.move(space,touch);
    //using namespace LvArray;
    //auto prev_space = mFieldView.getPreviousSpace();
    //std::cout<< prev_space << " : " << space << std::endl;

    //if(space != LvArray::MemorySpace::undefined && space != prev_space) {
    //  std::cout<<"FieldView"<< &(mFieldView.dataBuffer()) <<std::endl;
    //  for (size_t i = 0; i < mFieldView.size(); i++) {
    //    std::cout<<"DataView "<<i<< " : " << &(mFieldView[i].mDataView.dataBuffer()) <<std::endl;
    //    mFieldView[i].mDataView.dataBuffer().move(space,touch);
    //  }
    //  mFieldView.dataBuffer().move(space,touch);
    //}

  }

  RAJA_HOST_DEVICE
  FIELD_VIEW_TYPE& operator[](const unsigned int index) {return mFieldView(index);}

  RAJA_HOST_DEVICE
  FIELD_VIEW_TYPE& operator[](const unsigned int index) const {return mFieldView(index);}

  RAJA_HOST_DEVICE
  inline LvFieldListView( LvFieldListView const & source) noexcept : mFieldView{source.mFieldView} {}
//  inline LvFieldListView( LvFieldListView const & source) noexcept {
//#if defined(LVARRAY_USE_CUDA) && !defined(__CUDA_ARCH__)
//    auto space = LvArray::internal::toMemorySpace(
//            LvArray::internal::getArrayManager().getExecutionSpace());
//    source.move(space);
//#endif
//    mFieldView = source.mFieldView;
//  }

private:
  ARRAY_VIEW_TYPE mFieldView;
};



void SpheralEvalDerivTest()
{

  
  // Setup Timers
  srand(3);
  RAJA::Timer seq_timer;
  RAJA::Timer launch_timer;
  RAJA::Timer timer_pair;
  RAJA::Timer timer_red;

  // Data Types to Use
  using DIM = Spheral::Dim<1>;

  //using DATA_TYPE = Spheral::GeomVector<3>;
  using DATA_TYPE = pdouble;
  //using DATA_TYPE = double;
  using TRS_UINT = RAJA::TypedRangeSegment<unsigned>;

  using FIELD_TYPE = LvFieldBase<DATA_TYPE>;
  using VIEW_TYPE = FIELD_TYPE::VIEW_TYPE;

  using FIELDLIST_TYPE = LvFieldList<DATA_TYPE>;

#define MEM_SPACE RAJA::Platform::cuda
#define HOST_SPACE RAJA::Platform::host

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


  // Setting up our "Field Data", this is done through simulation setup in spheral e.g. node generation.
  FIELD_TYPE A(data_sz, "A");
  FIELD_TYPE B(data_sz, "B");
  FIELD_TYPE C(data_sz, "C");

  VIEW_TYPE Av(A);
  VIEW_TYPE Bv(B);
  VIEW_TYPE Cv(C);

#if USE_DEVICE
  FIELD_TYPE g_A(strat.n_blocks * data_sz, "g_A", MEM_SPACE);
  FIELD_TYPE g_B(strat.n_blocks * data_sz, "g_B", MEM_SPACE);
  FIELD_TYPE g_C(strat.n_blocks * data_sz, "g_C", MEM_SPACE);

  pairs.move(chai::GPU);
#else
  FIELD_TYPE g_A(strat.n_blocks * data_sz, "g_A");
  FIELD_TYPE g_B(strat.n_blocks * data_sz, "g_B");
  FIELD_TYPE g_C(strat.n_blocks * data_sz, "g_C");
#endif

  VIEW_TYPE g_Av(g_A);
  VIEW_TYPE g_Bv(g_B);
  VIEW_TYPE g_Cv(g_C);

  // Creating "FieldLists" In evalDerivs we call STATE_TYPE::fields(...) to return a fieldList.
  // In evalderivs we will want to return Something like LvFieldListView types for RAJA lambda 
  // capture and chai data migration. 
  LvFieldList<DATA_TYPE> fl("MyFirstFieldList");
  LvFieldList<DATA_TYPE> fl2("MySecondFieldList");
  fl.appendField(Av);
  fl.appendField(Bv);
  fl2.appendField(Cv);
  fl2.appendField(Av);

  // The FieldList types used in evalderivs.
  LvFieldListView<DATA_TYPE> flv(fl);
  LvFieldListView<DATA_TYPE> flv2(fl2);

  flv.move(MEM_SPACE);
  flv2.move(MEM_SPACE);

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
      ATOMIC_ADD(&g_Av[g_idx], DATA_TYPE(1.0));
      ATOMIC_ADD(&g_Bv[g_idx], DATA_TYPE(1.0));
      ATOMIC_ADD(&g_Cv[g_idx], DATA_TYPE(1.0));
#else
      // When executing on host we create one memory pool per omp thread, atomics are not needed.
      g_Av[g_idx] += DATA_TYPE(1.0);
      g_Bv[g_idx] += DATA_TYPE(1.0);
      g_Cv[g_idx] += DATA_TYPE(1.0);
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
        flv[0][t_idx] += g_Av[g_idx];
        flv[1][t_idx] += g_Bv[g_idx];
        flv2[0][t_idx] += g_Cv[g_idx];
      }
    });

  std::cout<<"Calc done\n";
  flv.move(HOST_SPACE);
  flv2.move(HOST_SPACE);

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

  DATA_TYPE check_A[data_sz];
  DATA_TYPE check_B[data_sz];
  DATA_TYPE check_C[data_sz];

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

  pairs.free();
}

#endif //  SPHERAL_SPHEVALDERIVTEST
