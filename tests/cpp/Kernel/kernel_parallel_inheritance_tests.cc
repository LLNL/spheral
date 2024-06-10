#include "test-utilities.hh"
#include "test-basic-exec-policies.hh"

#include "Field/SphArray.hh"
#include "Utilities/ValueViewInterface.hh"

//--------------------------------
// Impl Interface

namespace impl {

template<typename Dim, typename Desc>
class Kernel : public Spheral::SPHERALCopyable<Kernel<Dim, Desc>> {

  SPHERAL_HOST_DEVICE Desc& asDesc() {
    return static_cast<Desc&>(const_cast<Kernel<Dim, Desc>&>(*this));
  }
public:
	SPHERAL_HOST_DEVICE void doSomething() { printf("Ki HD doSomething()\n"); asDesc().doSomething(); }

  friend Kernel deepCopy(Kernel const& rhs) {
    Kernel result(rhs);
    return result;
  }

};

template<typename Dim>
class TableKernel : public Kernel<Dim, TableKernel<Dim>> {
public:
  TableKernel() = default;
	SPHERAL_HOST_DEVICE void doSomething() { printf("TKi HD doSomething()\n"); printf("TableKernel doSomething\n"); }
};

template<typename Dim>
class OtherKernel : public Kernel<Dim, OtherKernel<Dim>> {
public:
  OtherKernel() = default;
	SPHERAL_HOST_DEVICE void doSomething() { printf("OKi HD doSomething()\n"); printf("OtherKernel doSomething\n"); }
};

} // namespace impl

//--------------------------------
// View Interface

template<typename Dim, typename Desc>
class Kernel;
template<typename Dim>
class TableKernel;
template<typename Dim>
class OtherKernel;

template<typename Dim, typename Desc>
VIEW_INTERFACE_METACLASS_DECLARATION( (Kernel<Dim, Desc>), (KernelView), (impl::Kernel<Dim, typename Desc::ImplType>) )

template<typename Dim>
VIEW_INTERFACE_METACLASS_DECLARATION( (TableKernel<Dim>), (TableKernelView), (impl::TableKernel<Dim>))

template<typename Dim>
VIEW_INTERFACE_METACLASS_DECLARATION( (OtherKernel<Dim>), (OtherKernelView), (impl::OtherKernel<Dim>))

//--------------------------------
// Value Interface

template<typename Dim, typename Desc>
class Kernel : public Spheral::SpheralValueInterface<KernelView<Dim, Desc>>
{
  VALUE_TYPE_ALIASES((KernelView<Dim, Desc>))
public:
  VALUE_DEF_CTOR(Kernel)
  VALUE_COPY_CTOR(Kernel)
  VALUE_ASSIGNEMT_OP()
  VALUE_TOVIEW_OP()
  ViewType operator&() { return toView(); }

	SPHERAL_HOST void doSomething() { printf("K H doSomething()\n"); this->sptr_data().doSomething(); }
};

template<typename Dim>
class TableKernel : public Spheral::SpheralValueInterface<TableKernelView<Dim>>
{
  VALUE_TYPE_ALIASES((TableKernelView<Dim>))
public:
  VALUE_DEF_CTOR(TableKernel)
  VALUE_COPY_CTOR(TableKernel)
  VALUE_ASSIGNEMT_OP()
  VALUE_TOVIEW_OP()
  ViewType operator&() { return toView(); }

	SPHERAL_HOST void doSomething() { printf("TK H doSomething()\n"); this->sptr_data().doSomething(); }
};

template<typename Dim>
class OtherKernel : public Spheral::SpheralValueInterface<OtherKernelView<Dim>>
{
  VALUE_TYPE_ALIASES((OtherKernelView<Dim>))
public:
  VALUE_DEF_CTOR(OtherKernel)
  VALUE_COPY_CTOR(OtherKernel)
  VALUE_ASSIGNEMT_OP()
  VALUE_TOVIEW_OP()
  ViewType operator&() { return toView(); }

	SPHERAL_HOST void doSomething() { printf("OK H doSomething()\n"); this->sptr_data().doSomething(); }
};

class Dim1 {};

//
// Setting up G Test for QuadraticInterpolator
template<typename T>
class KernelParallelInterface : public::testing::Test {};

// All QuadraticInterpolatorTets cases will run over each type in EXEC_TYPES.
TYPED_TEST_CASE(KernelParallelInterface, EXEC_TYPES);



TEST(KernelParallelImpl, Impl)
{

  impl::TableKernel<Dim1> tk_impl;
  
  tk_impl.doSomething();

  impl::Kernel<Dim1, impl::TableKernel<Dim1>> k_impl;
  
  k_impl.doSomething();

}

GPU_TYPED_TEST(KernelParallelInterface, TableKernelInterface)
{
  using WORK_EXEC_POLICY = TypeParam;

  TableKernel<Dim1> tk;
  
  tk.doSomething();

  TableKernelView<Dim1> tkv = &tk;

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    tkv->doSomething();
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST(KernelParallelInterface, KernelInterface)
{
  using WORK_EXEC_POLICY = TypeParam;

  Kernel<Dim1, TableKernel<Dim1>> k;
  
  k.doSomething();

  KernelView<Dim1, TableKernel<Dim1>> kv = &k;

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
  kv->doSomething();
  EXEC_IN_SPACE_END()

  KernelView<Dim1, TableKernelView<Dim1>> kv_tkv = &k;

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
  kv_tkv->doSomething();
  EXEC_IN_SPACE_END()

  Kernel<Dim1, TableKernelView<Dim1>> k_tkv;

  k_tkv->doSomething();

  KernelView<Dim1, TableKernelView<Dim1>> kv_tkv2 = &k_tkv;

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
  kv_tkv2->doSomething();
  EXEC_IN_SPACE_END()
  
}

GPU_TYPED_TEST(KernelParallelInterface, OtherKernelInterface)
{
  using WORK_EXEC_POLICY = TypeParam;

  Kernel<Dim1, OtherKernel<Dim1>> k;
  
  k->doSomething();

  KernelView<Dim1, OtherKernel<Dim1>> kv = &k;

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
  kv->doSomething();
  EXEC_IN_SPACE_END()

  KernelView<Dim1, OtherKernelView<Dim1>> kv_tkv = &k;

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
  kv_tkv->doSomething();
  EXEC_IN_SPACE_END()

  Kernel<Dim1, OtherKernelView<Dim1>> k_tkv;

  k_tkv->doSomething();

  KernelView<Dim1, OtherKernelView<Dim1>> kv_tkv2 = &k_tkv;

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
  kv_tkv2->doSomething();
  EXEC_IN_SPACE_END()
  
}
