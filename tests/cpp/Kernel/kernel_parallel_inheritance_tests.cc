#include "test-utilities.hh"

#include "Field/SphArray.hh"
#include "Utilities/ValueViewInterface.hh"

// New macro for this implementation.
#define SPTR_MOVE_CTOR(type) \
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
class Kernel;
template<typename Dim>
class TableKernel;

template<typename Dim, typename Desc>
class KernelView :
	public Spheral::SpheralViewInterface< KernelView<Dim, Desc>, KernelImpl<Dim, typename Desc::ImplType> >
{
  VIEW_TYPE_ALIASES( (Kernel<Dim, Desc>), (KernelView<Dim, Desc>), (KernelImpl<Dim, typename Desc::ImplType>) )
  VIEW_DEFINE_ALLOC_CTOR(KernelView)
public:
  SPTR_MOVE_CTOR(KernelView)

	SPHERAL_HOST_DEVICE void doSomething() { std::cout << "Kv HD doSomething()\n"; this->sptr_data().doSomething(); }
};

template<typename Dim>
class TableKernelView : 
	public Spheral::SpheralViewInterface<TableKernelView<Dim>,TableKernelImpl<Dim>>
{
  VIEW_TYPE_ALIASES( (TableKernel<Dim>), (TableKernelView), (TableKernelImpl<Dim>))
  VIEW_DEFINE_ALLOC_CTOR(TableKernelView)
public:
  SPTR_MOVE_CTOR(TableKernelView)

	SPHERAL_HOST_DEVICE void doSomething() { std::cout << "TKv HD doSomething()\n"; this->sptr_data().doSomething(); }
};

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

	SPHERAL_HOST void doSomething() { std::cout << "K H doSomething()\n"; this->sptr_data().doSomething(); }
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

	SPHERAL_HOST void doSomething() { std::cout << "TK H doSomething()\n"; this->sptr_data().doSomething(); }
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

  KernelView<Dim1, TableKernelView<Dim1>> kv_tkv = k.toView();

  kv_tkv.doSomething();

  Kernel<Dim1, TableKernelView<Dim1>> k_tkv;

  k_tkv.doSomething();

  KernelView<Dim1, TableKernelView<Dim1>> kv_tkv2 = k_tkv.toView();

  kv_tkv2.doSomething();
  
}
