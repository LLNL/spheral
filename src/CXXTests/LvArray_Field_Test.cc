#include "Field/Field.hh"

#include "LvArray/Array.hpp"
#include "LvArray/ChaiBuffer.hpp"

using EXEC_POL = RAJA::cuda_exec<256>;
using HOST_POL = RAJA::seq_exec;

template<typename FieldList>
struct GPUFieldList {

  GPUFieldList()
  {
#if !defined(CHAI_DEVICE_COMPILE)
    m_array_manager = chai::ArrayManager::getInstance();
    m_prev_space = new chai::ExecutionSpace{chai::CPU};
#endif
  }

  GPUFieldList(FieldList* flp) : GPUFieldList() { m_fieldptrs = flp; };

  RAJA_HOST_DEVICE GPUFieldList( const GPUFieldList& other)
    : m_array_manager(other.m_array_manager),
      m_prev_space(other.m_prev_space),
      m_fieldptrs(other.m_fieldptrs)
  {
#if !defined(CHAI_DEVICE_COMPILE) 
    move(m_array_manager->getExecutionSpace());
#endif
  }

  void move(chai::ExecutionSpace space)
  {
    if (space != chai::NONE /*&& *m_prev_space != space*/)
    {
      for(size_t i = 0; i < m_fieldptrs->size(); i++) {
        (*m_fieldptrs)[i]->mDataArray.move(LvArray::internal::toMemorySpace(space));
      }
      *m_prev_space = space;
    }
  }

  using Field = typename std::remove_pointer<typename FieldList::value_type>::type;
  using FieldView = typename Field::ContainerTypeView;
  using ValueType = typename Field::ValueType;

  RAJA_HOST_DEVICE FieldView operator[](const unsigned int index){ return (*m_fieldptrs)[index]->getAccessorView(); }
  RAJA_HOST_DEVICE const FieldView operator[](const unsigned int index) const { return (*m_fieldptrs)[index]->getAccessorView(); }

  RAJA_HOST_DEVICE ValueType& operator()(const unsigned int fieldIdx,
                                         const unsigned int nodeIdx)
  {
    return (*m_fieldptrs)[fieldIdx]->getAccessorView()[nodeIdx];
  }

  FieldList* m_fieldptrs = nullptr;
  chai::ArrayManager* m_array_manager = nullptr;
  chai::ExecutionSpace* m_prev_space;
};


  RAJA_HOST_DEVICE
  void atomicAdd(Spheral::GeomVector<3>* addr, Spheral::GeomVector<3> rhs) {
    RAJA::atomicAdd<RAJA::cuda_atomic>(&addr[0][0], rhs[0]);
    RAJA::atomicAdd<RAJA::cuda_atomic>(&addr[0][1], rhs[1]);
    RAJA::atomicAdd<RAJA::cuda_atomic>(&addr[0][2], rhs[2]);
  }


int main() {

  constexpr int N = 1000000;

#if 1
  {

    std::cout << "Prototype FieldList data manipulation w/ atomics.\n";

    using Dim = Spheral::Dim<3>;
    Spheral::NodeList<Dim> node_list("example node list", N, 0);
    auto n_pos = node_list.positions();

    using FLCT = LvArray::Array< decltype(n_pos)*, 1, camp::idx_seq<0>, std::ptrdiff_t, LvArray::ChaiBuffer >;
    FLCT fake_field_list;
    fake_field_list.emplace_back(&n_pos);

    GPUFieldList <decltype(fake_field_list)>gpufl(&fake_field_list);

    auto field_view = n_pos.getAccessorView();

    RAJA::RangeSegment range(0, field_view.size());
    RAJA::forall<HOST_POL>(range,
      [=](unsigned int kk) mutable {
        auto& v = gpufl(0, kk);
        v = Spheral::GeomVector<3>(1,1,1);
    });

    RAJA::forall<EXEC_POL>(range, 
      [=] RAJA_HOST_DEVICE (int kk) mutable {
        if (kk < 50) {
          auto& v0 = gpufl(0,0);
          ::atomicAdd(&v0, Spheral::GeomVector<3>(1,1,1));
        }
        if (kk != 0) {
          auto& v = gpufl(0, kk);
          v += Spheral::GeomVector<3>(1,1,1);
        }
    });

    bool correctness = true;
    RAJA::forall<HOST_POL>(range,
      [=, &correctness] (int kk) {
        if (kk == 0) {
          double ans = N > 50 ? 50 : N; ans++;
          if (gpufl[0][kk] != Spheral::GeomVector<3>(ans,ans,ans)) {
            correctness = false;
            std::cout << gpufl[0][kk] << "\n";
          }
        }else{
          if (gpufl[0][kk] != Spheral::GeomVector<3>(2,2,2)) correctness = false;
        }
    });

    if (correctness)
      std::cout << "PASSED\n\n";
    else
      std::cout << "FAILED\n\n";
  }
#endif

#if 1
  {

    std::cout << "Basic field data manipulation.\n";

    using Dim = Spheral::Dim<3>;
    Spheral::NodeList<Dim> node_list("example node list", N, 0);
    auto n_pos = node_list.positions();

    auto field_view = n_pos.getAccessorView();

    RAJA::RangeSegment range(0, field_view.size());
    RAJA::forall<HOST_POL>(range,
      [=](unsigned int kk) {
        field_view[kk][0]++;
        field_view[kk][1]++;
        field_view[kk][2]++;
    });

    RAJA::forall<EXEC_POL>(range, 
      [=] RAJA_HOST_DEVICE (int kk) {
        field_view[kk][0]++;
        field_view[kk][1]++;
        field_view[kk][2]++;
    });

    bool correctness = true;
    RAJA::forall<HOST_POL>(range,
      [=, &correctness] (int kk) {
        if (field_view[kk] != Spheral::GeomVector<3>(2,2,2)) correctness = false;
    });

    if (correctness)
      std::cout << "PASSED\n\n";
    else
      std::cout << "FAILED\n\n";
  }
#endif

#if 1
  {
    std::cout << "Atomic field data manipulation.\n";

    using Dim = Spheral::Dim<3>;
    Spheral::NodeList<Dim> node_list("example node list", N, 0);
    auto n_pos = node_list.positions();

    auto field_view = n_pos.getAccessorView();

    RAJA::RangeSegment range(0, field_view.size());
    RAJA::forall<HOST_POL>(range,
      [=](unsigned int kk) {
        field_view[kk][0] = 0;
        field_view[kk][1] = 0;
        field_view[kk][2] = 0;
    });

    RAJA::forall<EXEC_POL>(range, 
      [=] RAJA_HOST_DEVICE (int kk) {
        if (kk < 50) {
          RAJA::atomicAdd<RAJA::cuda_atomic>(&field_view[0][0], 1.0);
          RAJA::atomicAdd<RAJA::cuda_atomic>(&field_view[0][1], 1.0);
          RAJA::atomicAdd<RAJA::cuda_atomic>(&field_view[0][2], 1.0);
        }
    });

    bool correctness = true;
    RAJA::forall<HOST_POL>(range,
      [=, &correctness] (int kk) {
        
        if (kk == 0) {
          double ans = N > 50 ? 50 : N;
          if (field_view[kk] != Spheral::GeomVector<3>(ans,ans,ans)) {
            correctness = false;
          }
        }
    });

    if (correctness)
      std::cout << "PASSED\n\n";
    else
      std::cout << "FAILED\n\n";
  }
#endif

  return EXIT_SUCCESS;
}
