#include "test-utilities.hh"

#include "Field/SphArray.hh"
#include "Utilities/ValueViewInterface.hh"
#include <memory>

namespace Spheral {

// New macro for this implementation.
#define SPTR_FWD_CTOR(type) \
  type(SmartPtrType&& rhs) : Base(std::forward<SmartPtrType>(rhs)) {}


class FBi {
//class FBi : public Spheral::SPHERALCopyable<FBi> {
  public:
    FBi() {}
    FBi(size_t h) : hash(h) {}

    void free() {}
    //FBi& operator=(std::nullptr_t) {return *this; }
    //void shallowCopy(FBi const& rhs) {*this = rhs;}

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

template<typename T>
class FV;
template<typename T>
class F;


class FBV : protected Spheral::SpheralViewInterface<FBV, FBi> {
  VIEW_TYPE_ALIASES(FBV, FBi)
public:
  VIEW_DEFINE_ALLOC_CTOR(FBV, FBi)
  VIEW_COPY_CTOR(FBV)
  SPHERAL_HOST_DEVICE virtual size_t getHash() const { return sptr_data().getHash(); }
  size_t size() const {return sptr_data().size(); }

  SPTR_FWD_CTOR(FBV)

  template<typename T>
  explicit operator FV<T>() const {return FV<T>(std::dynamic_pointer_cast<Fi<T>>(this->sptr()));}
};


class FB : public Spheral::SpheralValueInterface<FBV, FBi> {
  VALUE_TYPE_ALIASES(FB, FBV, FBi)
protected:
public:
  //VALUE_DEF_CTOR(F, Fi<T>)
  //VALUE_COPY_CTOR(FB, FBi)
  //VALUE_ASSIGNEMT_OP(FB, FBi)
  VALUE_TOVIEW_OP()

  // Can no longer be pure
  void resize(size_t sz) {SPTR_DATA_REF().resize(sz);}

  // Underlying smart pointer conversion happens in here...
  SPTR_FWD_CTOR(FB)
};



template<typename T>
class FV : public Spheral::SpheralViewInterface<FV<T>, Fi<T>> {
  VIEW_TYPE_ALIASES(FV, Fi<T>)
public:
  VIEW_DEFINE_ALLOC_CTOR(FV, Fi<T>)

  SPHERAL_HOST_DEVICE T* data() {return SPTR_DATA_REF().data(); }
  size_t size() { return SPTR_DATA_REF().size(); }

  explicit operator FBV() const {return FBV(this->sptr());}

  size_t getHash() const { return SPTR_DATA_REF().getHash(); }
  SPTR_FWD_CTOR(FV)
};

template<typename T>
class F : public Spheral::SpheralValueInterface<FV<T>, Fi<T>> {
  VALUE_TYPE_ALIASES(F, FV<T>, Fi<T>)
public:
  // Def Ctor, Copy Ctor, Assign Op,toView
  VALUE_DEF_CTOR(F, Fi<T>)
  VALUE_COPY_CTOR(F, Fi<T>)
  VALUE_ASSIGNEMT_OP(F, Fi<T>)
  VALUE_TOVIEW_OP()

  explicit operator FB() const {return FB(this->sptr());}

  // Ctor 
  F(size_t h, size_t sz) : Base(new Fi<T>(h, sz)) {}

  void resize(size_t sz) { SPTR_DATA_REF().resize(sz); } 

  // Methods defined in base implementation need to be explicilty
  // defined in the top level class to remove ambiguous function calls.
  size_t size() { return FV<T>::size(); }
};

}// namespace Spheral


TEST(FieldOrthInheritance, AccessPattern)
{

  Spheral::F<double> f(2, 200);
  //Spheral::FB* fbptr = &f; // Can not do this.
  Spheral::FB fb = (Spheral::FB)f; // Can not do this.
  //Spheral::FB fb2;
  auto f_v = f.toView();
  Spheral::FBV fb_v = (Spheral::FBV)f_v;

  Spheral::FV<double> f_v2 = (Spheral::FV<double>)(fb_v);

  //std::cout << "fb_v : " << fb_v.getHash() << " , " << std::endl;
  std::cout << "fb   : " << fb.getHash()<< " , " << std::endl;
  //std::cout << "f_v  : " << f_v.getHash()   << " , " << f_v.data()   << " , " << f_v.size() << std::endl;
  //std::cout << "f    : " << f.getHash() << " , " << std::endl;
  std::cout << "f    : " << f.getHash()     << " , " << f.data()     << " , " << f.size()   << std::endl;
  std::cout << "f_v2 : " << f_v2.getHash() << " , " << f_v2.data() << std::endl;
  //

  ////Spheral::FB* fbptr = &f;
  fb.resize(1123);
  //fbptr->resize(1123);


  
  std::cout << "f    : " << f.getHash()     << " , " << f.data()     << " , " << f.size()   << std::endl;
  std::cout << "f_v  : " << f_v.getHash()   << " , " << f_v.data()   << " , " << f_v.size() << std::endl;

  std::cout << fb.size() << std::endl;
}
