#include "test-utilities.hh"

#include "Field/SphArray.hh"
#include "Utilities/ValueViewInterface.hh"

// New macro for this implementation.
#define SPTR_FWD_CTOR(type) \
  type(SmartPtrType&& rhs) : Base(std::forward<SmartPtrType>(rhs)) {}

//--------------------------------
// Impl Interface

template<typename Dim, typename Desc>
class KernelImpl : public Spheral::SPHERALCopyable<KernelImpl<Dim, Desc>> {

  Desc& asDesc() {
    return static_cast<Desc&>(const_cast<KernelImpl<Dim, Desc>&>(*this));
  }
public:
	SPHERAL_HOST_DEVICE void doSomething() { std::cout << "Ki HD doSomething()\n"; asDesc().doSomething(); }

  void free() {}
  KernelImpl& operator=(std::nullptr_t) { return *this; }
  void shallowCopy(KernelImpl const& rhs) {*this = rhs;}

  friend KernelImpl deepCopy(KernelImpl const& rhs) {
    KernelImpl result(rhs);
    return result;
  }

};

template<typename Dim>
class TableKernelImpl : public KernelImpl<Dim, TableKernelImpl<Dim>> {
public:
  TableKernelImpl() = default;
	SPHERAL_HOST_DEVICE void doSomething() { std::cout << "TKi HD doSomething()\n"; std::cout << "TableKernel doSomething\n"; }
};

//--------------------------------
// View Interface

template<typename Dim, typename Desc>
class KernelView :
	public Spheral::SpheralViewInterface<
		KernelView<Dim, Desc>,
		KernelImpl<Dim, typename Desc::ImplType>
	>
{
  using DescImplType = typename Desc::ImplType;
  UNPACK_VIEW_TYPE_ALIASES( (KernelView<Dim, Desc>), (KernelImpl<Dim, DescImplType>) )
  UNPACK_VIEW_DEFINE_ALLOC_CTOR(KernelView)

public:
  SPTR_FWD_CTOR(KernelView)

	SPHERAL_HOST_DEVICE void doSomething() { std::cout << "Kv HD doSomething()\n"; this->sptr_data().doSomething(); }
};

template<typename Dim>
class TableKernelView : 
	public Spheral::SpheralViewInterface<
		TableKernelView<Dim>,
		TableKernelImpl<Dim>
	>
{
  UNPACK_VIEW_TYPE_ALIASES((TableKernelView), (TableKernelImpl<Dim>))
  UNPACK_VIEW_DEFINE_ALLOC_CTOR(TableKernelView)
public:
  SPTR_FWD_CTOR(TableKernelView)

	SPHERAL_HOST_DEVICE void doSomething() { std::cout << "TKv HD doSomething()\n"; this->sptr_data().doSomething(); }
};

//--------------------------------
// Value Interface

template<typename Dim, typename Desc>
class Kernel :
	public Spheral::SpheralValueInterface<
		KernelView<Dim, Desc>,
		KernelImpl<Dim, typename Desc::ImplType>
	>
{
  using DescImplType = typename Desc::ImplType;
  UNPACK_VALUE_TYPE_ALIASES( (Kernel<Dim, Desc>), (KernelView<Dim, Desc>), (KernelImpl<Dim, DescImplType>) )

public:
  UNPACK_VALUE_DEF_CTOR(Kernel)
  UNPACK_VALUE_COPY_CTOR(Kernel)
  UNPACK_VALUE_ASSIGNEMT_OP()
  VALUE_TOVIEW_OP()

	SPHERAL_HOST_DEVICE void doSomething() { std::cout << "K H doSomething()\n"; this->sptr_data().doSomething(); }
};

template<typename Dim>
class TableKernel :
	public Spheral::SpheralValueInterface<
		TableKernelView<Dim>,
		TableKernelImpl<Dim>
	>
{
  UNPACK_VALUE_TYPE_ALIASES( (TableKernel<Dim>), (TableKernelView<Dim>), (TableKernelImpl<Dim>) )

public:
  UNPACK_VALUE_DEF_CTOR(TableKernel)
  UNPACK_VALUE_COPY_CTOR(TableKernel)
  UNPACK_VALUE_ASSIGNEMT_OP()
  VALUE_TOVIEW_OP()

	SPHERAL_HOST_DEVICE void doSomething() { std::cout << "TK H doSomething()\n"; this->sptr_data().doSomething(); }
};

class Dim1 {};


TEST(KernelParallelInterface, Impl)
{

  TableKernelImpl<Dim1> tk_impl;
  
  tk_impl.doSomething();

  KernelImpl<Dim1, TableKernelImpl<Dim1>> k_impl;
  
  k_impl.doSomething();

}

TEST(KernelParallelInterface, TableKernelInterface)
{

  TableKernel<Dim1> tk;
  
  tk.doSomething();

  TableKernelView<Dim1> tkv = tk.toView();

  tkv.doSomething();
}

TEST(KernelParallelInterface, KernelInterface)
{

  Kernel<Dim1, TableKernel<Dim1>> k;
  
  k.doSomething();

  KernelView<Dim1, TableKernel<Dim1>> kv = k.toView();

  kv.doSomething();
}
