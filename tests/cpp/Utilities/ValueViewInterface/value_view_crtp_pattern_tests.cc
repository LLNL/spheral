#include "test-utilities.hh"
#include "test-basic-exec-policies.hh"

#include "Utilities/ValueViewInterface.hh"
#include "Utilities/ManagedVector.hh"

//--------------------------------
// Impl Interface

VVI_IMPL_BEGIN

template<typename Dim, typename Desc>
class Kernel : public Spheral::SPHERALCopyable {

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

VVI_IMPL_END

//--------------------------------
// View Interface

template<typename Dim, typename Desc>
class Kernel;
template<typename Dim>
class TableKernel;
template<typename Dim>
class OtherKernel;

#define KernelView__(code) PTR_VIEW_METACLASS_DECL( (Kernel<Dim, Desc>), (KernelView), (vvimpl::Kernel<Dim, typename Desc::ImplType>), (code))
#define TableKernelView__(code) PTR_VIEW_METACLASS_DECL( (TableKernel<Dim>), (TableKernelView), (vvimpl::TableKernel<Dim>), (code))
#define OtherKernelView__(code) PTR_VIEW_METACLASS_DECL( (OtherKernel<Dim>), (OtherKernelView), (vvimpl::OtherKernel<Dim>), (code))

#define Kernel__(code) PTR_VALUE_METACLASS_DECL((Kernel), (KernelView<Dim, Desc>), (code))
#define TableKernel__(code) PTR_VALUE_METACLASS_DECL((TableKernel), (TableKernelView<Dim>), (code))
#define OtherKernel__(code) PTR_VALUE_METACLASS_DECL((OtherKernel), (OtherKernelView<Dim>), (code))

template<typename Dim, typename Desc>
class KernelView__( VVI_VIEW_DEFAULT(KernelView) );

template<typename Dim>
class TableKernelView__( VVI_VIEW_DEFAULT(TableKernelView) );

template<typename Dim>
class OtherKernelView__( VVI_VIEW_DEFAULT(OtherKernelView) );

//--------------------------------
// Value Interface

template<typename Dim, typename Desc>
class Kernel__(
public:
  VVI_VALUE_DEF_CTOR(Kernel)
	void doSomething() { printf("K H doSomething()\n"); VVI_IMPL_INST().doSomething(); }
);

template<typename Dim>
class TableKernel__(
public:
  VVI_VALUE_DEF_CTOR(TableKernel)
	void doSomething() { printf("TK H doSomething()\n"); VVI_IMPL_INST().doSomething(); }
);

template<typename Dim>
class OtherKernel__(
public:
  VVI_VALUE_DEF_CTOR(OtherKernel)
	void doSomething() { printf("OK H doSomething()\n"); VVI_IMPL_INST().doSomething(); }
);

class Dim1 {};

//
// Setting up G Test for QuadraticInterpolator
template<typename T>
class KernelParallelInterface : public::testing::Test {};

// All QuadraticInterpolatorTets cases will run over each type in EXEC_TYPES.
TYPED_TEST_CASE(KernelParallelInterface, EXEC_TYPES);



TEST(KernelParallelImpl, Impl)
{

  vvimpl::TableKernel<Dim1> tk_impl;
  
  tk_impl.doSomething();

  vvimpl::Kernel<Dim1, vvimpl::TableKernel<Dim1>> k_impl;
  
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
