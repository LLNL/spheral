#include "Field/Field.hh"

#include "LvArray/Array.hpp"
#include "LvArray/ChaiBuffer.hpp"

using EXEC_POL = RAJA::cuda_exec<256>;
using HOST_POL = RAJA::seq_exec;

template<typename T>
using Array1D = LvArray::Array< T, 1, camp::idx_seq<0>, std::ptrdiff_t, LvArray::ChaiBuffer >;

template<typename T>
using Array1DView = LvArray::ArrayView< T, 1, 0, std::ptrdiff_t, LvArray::ChaiBuffer >;


namespace Spheral{
  namespace detail{

#define DEVICE_ACCESSOR(VAL) DeviceAccessor<decltype(VAL)>(VAL)

  template<typename T>
  class DeviceAccessor {

    using value_type = typename T::ValueType;
    using array_type = typename T::ContainerType;
    using view_type  = typename T::ContainerTypeView;

  public:
    DeviceAccessor(const T& val) : array_parent(val.mDataArray), view(val.mDataArray) {}

    unsigned size() const { return view.size(); }

    template<typename IDX_TYPE>
    RAJA_HOST_DEVICE value_type& operator[](const IDX_TYPE idx) const { return view[idx]; }

    void move( const LvArray::MemorySpace space ) {
    #if defined(RAJA_ENABLE_CUDA)
      array_parent.move(space);
    #else
      RAJA_UNUSED_VAR(space);
    #endif
    }
                                                                       
  private:
    const array_type& array_parent;
    const view_type& view;
  };

  } // namespace detail
} // namespace Spheral


int main() {

  constexpr int N = 10;


  using Dim = Spheral::Dim<3>;
  Spheral::NodeList<Dim> node_list("example node list", N, 0);
  auto n_pos = node_list.positions();

  auto field_view = Spheral::detail::DEVICE_ACCESSOR(n_pos);

  RAJA::RangeSegment range(0, field_view.size());
  RAJA::forall<HOST_POL>(range,
    [=](unsigned int kk) {
      field_view[kk][0]++;
      field_view[kk][1]++;
      field_view[kk][2]++;
  });

  field_view.move(LvArray::MemorySpace::cuda);

  RAJA::forall<EXEC_POL>(range, 
    [=] RAJA_HOST_DEVICE (int kk) {
      field_view[kk][0]++;
      field_view[kk][1]++;
      field_view[kk][2]++;
  });

  field_view.move(LvArray::MemorySpace::host);

  bool correctness = true;
  RAJA::forall<HOST_POL>(range,
    [=, &correctness] (int kk) {
      if (field_view[kk] != Spheral::GeomVector<3>(2,2,2)) correctness = false;
  });

  if (correctness)
    std::cout << "PASSED\n";
  else
    std::cout << "FAILED\n";

  return EXIT_SUCCESS;
}
