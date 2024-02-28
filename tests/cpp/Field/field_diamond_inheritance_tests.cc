#include "test-utilities.hh"

#include "Field/SphArray.hh"
#include "Utilities/ValueViewInterface.hh"

namespace Spheral {


// Some macro redefinitions that can be reworked to be variabdic and replaced in ValueViewInterface.
#define VIEW_DEFINE_ALLOC_CTOR_DESCENDENT(view_t, impl_t, parent_t) \
public: \
  view_t(impl_t* rhs) : Base(SmartPtrType(rhs, [](impl_t *p) { p->free(); } )), parent_t(SPTR_REF()){}

#define VALUE_DEF_CTOR_DESCENDENT(type, impl_t, parent_t) \
  type() : Base(new impl_t()), parent_t(SPTR_REF()) {}

#define VALUE_COPY_CTOR_DESCENDENT(type, impl_t, parent_t) \
  type(type const& rhs) : Base(new impl_t( deepCopy( rhs.SPTR_DATA_REF() ) )), parent_t(SPTR_REF()) {}

// New macro for this implementation.
#define PARENT_SPTR_FWD_CTOR(type) \
  type(SmartPtrType&& rhs) : Base(std::forward<SmartPtrType>(rhs)) {}


class FBi {
  public:
    FBi() {}
    FBi(size_t h) : hash(h) {}

    void free() {}

    size_t getHash() const { return hash; }

    size_t hash = 0;

    virtual void resize(size_t sz) = 0;
    virtual size_t size() = 0;
};

template<typename T>
class Fi : public FBi, Spheral::SPHERALCopyable<Fi<T>>{
public:
  Spheral::ManagedVector<T> m_data;
  Fi() : FBi(), m_data(0) {}
  Fi(size_t h, size_t sz) : FBi(h), m_data(sz) {}

  T* data() { return &m_data[0]; }

  void free() {m_data.free();}
  Fi& operator=(std::nullptr_t) {m_data = nullptr; return *this; }
  void shallowCopy(Fi const& rhs) {*this = rhs;}

  friend Fi deepCopy(Fi const& rhs) {
    Fi result(rhs);
    result.m_data = Spheral::deepCopy(rhs.m_data);
    return result;
  }

  virtual void resize(size_t sz) override {std::cout << "Resize : " << sz << std::endl; m_data.resize(sz); } 
  virtual size_t size() override { return m_data.size(); }
};


class FBV : protected Spheral::SpheralViewInterface<FBV, FBi> {
  VIEW_TYPE_ALIASES(FBV, FBi)
  VIEW_DEFINE_ALLOC_CTOR(FBV, FBi)
protected:
  PARENT_SPTR_FWD_CTOR(FBV)
public:
  SPHERAL_HOST_DEVICE virtual size_t getHash() const { return sptr_data().getHash(); }
  virtual size_t size() = 0;
};


class FB : public Spheral::SpheralValueInterface<FBV, FBi> {
  VALUE_TYPE_ALIASES(FB, FBV, FBi)
protected:
  PARENT_SPTR_FWD_CTOR(FB)
public:
  virtual void resize(size_t sz) = 0;
};



template<typename T>
class FV : protected Spheral::SpheralViewInterface<FV<T>, Fi<T>>, public FBV {
  VIEW_TYPE_ALIASES(FV, Fi<T>)
public:
  VIEW_DEFINE_ALLOC_CTOR_DESCENDENT(FV, Fi<T>, FBV)

  SPHERAL_HOST_DEVICE T* data() {return SPTR_DATA_REF().data(); }
  virtual size_t size() override { return SPTR_DATA_REF().size(); }
};

  // We need to pass the shared_ptr created by FV<T> to the Ctor of 
  // FB so we are calling a copy on the smart pointer rather than a 
  // new allocation. 
  //
  // We use macros for accessing sptr() and sptr_data() so we can be explicit 
  // about using the View version of those funciton calls in the given class.
template<typename T>
class F : public Spheral::SpheralValueInterface<FV<T>, Fi<T>>, public FB {
  VALUE_TYPE_ALIASES(F, FV<T>, Fi<T>)

public:
  // Def Ctor, Copy Ctor, Assign Op,toView
  VALUE_DEF_CTOR_DESCENDENT(F, Fi<T>, FB)
  VALUE_COPY_CTOR_DESCENDENT(F, Fi<T>, FB)
  VALUE_ASSIGNEMT_OP(F, Fi<T>)
  VALUE_TOVIEW_OP()

  // Ctor 
  F(size_t h, size_t sz) : Base(new Fi<T>(h, sz)), FB(SPTR_REF()) {}

  virtual void resize(size_t sz) override { SPTR_DATA_REF().resize(sz); } 

  // Methods defined in base implementation need to be explicilty
  // defined in the top level class to remove ambiguous function calls.
  size_t getHash() const { return FB::getHash(); }
  size_t size() { return FV<T>::size(); }
};

}// namespace Spheral


TEST(FieldDiamondInheritance, AccessPattern)
{

  Spheral::F<double> f(2, 200);
  Spheral::FB* fbptr = &f;
  auto f_v = f.toView();
  Spheral::FBV* fb_v = &f_v;
  Spheral::FV<double>* f_v2 = dynamic_cast<Spheral::FV<double>*>(fb_v);

  std::cout << "fb_v : " << fb_v->getHash() << " , " << std::endl;
  std::cout << "fbptr: " << fbptr->getHash()<< " , " << std::endl;
  std::cout << "f_v  : " << f_v.getHash()   << " , " << f_v.data()   << " , " << f_v.size() << std::endl;
  std::cout << "f.FB : " << f.FB::getHash() << " , " << std::endl;
  std::cout << "f    : " << f.getHash()     << " , " << f.data()     << " , " << f.size()   << std::endl;
  std::cout << "f_v2 : " << f_v2->getHash() << " , " << f_v2->data() << std::endl;
  

  fbptr->resize(1123);


  std::cout << "f    : " << f.getHash()     << " , " << f.data()     << " , " << f.size()   << std::endl;
  std::cout << "f_v  : " << f_v.getHash()   << " , " << f_v.data()   << " , " << f_v.size() << std::endl;

  std::cout << fbptr->size() << std::endl;
  std::cout << f_v2->size() << std::endl;


  Spheral::ManagedVector<Spheral::FBV*> vec_fbv(5);
  std::cout << vec_fbv.size() << std::endl;
  vec_fbv.free();


}
