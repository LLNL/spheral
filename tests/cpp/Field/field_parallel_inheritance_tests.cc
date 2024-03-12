#include "test-utilities.hh"

#include "Field/SphArray.hh"
#include "Utilities/ValueViewInterface.hh"
#include <memory>

namespace Spheral {

// New macro for this implementation.
#define SPTR_FWD_CTOR(type) \
  type(SmartPtrType&& rhs) : Base(std::forward<SmartPtrType>(rhs)) {}

#define UPCAST_CONVERSION_OP(parent_t) \
  explicit operator parent_t() const {return parent_t(this->sptr());}

#define DOWNCAST_CONVERSION_OP(child_t, impl_t) \
  explicit operator child_t() const {return child_t(std::dynamic_pointer_cast<impl_t>(this->sptr()));}

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


// We need to forward declare child classes for use in conversion functions.
template<typename T>
class FV;
template<typename T>
class F;

class FB;

class FBV : protected Spheral::SpheralViewInterface<FBV, FBi> {
  VIEW_TYPE_ALIASES((FB), (FBV), (FBi))
  VIEW_DEFINE_ALLOC_CTOR(FBV)

public:
  VIEW_DEF_CTOR(FBV)
  VIEW_COPY_CTOR(FBV)
  //VIEW_ASSIGNEMT_OP()
  FBV& operator=(FBV const&) = default;

  SPTR_FWD_CTOR(FBV)

  template<typename T>
  DOWNCAST_CONVERSION_OP(FV<T>, Fi<T>)

  SPHERAL_HOST_DEVICE size_t getHash() const { return sptr_data().getHash(); }

  // We want to forward the pure virtual interface here.
  SPHERAL_HOST_DEVICE size_t size() const {return sptr_data().size(); }
  SPHERAL_HOST_DEVICE void resize(size_t sz) { sptr_data().resize(sz); }
};

// Because the underlying impl type is pure virtual we can not allow
// construction of FB class with default Ctor, Copy Ctor or assignment Op. 
// We can only construct a FB object from an existing smart_ptr type so 
// they must be constructed from F type objects directly...
class FB : public Spheral::SpheralValueInterface<FBV, FBi> {
  VALUE_TYPE_ALIASES((FBV))
  SPTR_FWD_CTOR(FB)
public:
  FB() = delete;
  FB(FB const&) = delete;
  FB& operator=(FB const&) = delete;
};



template<typename T>
class FV : public Spheral::SpheralViewInterface<FV<T>, Fi<T>> {
  VIEW_TYPE_ALIASES((F<T>), (FV), (Fi<T>))
  VIEW_DEFINE_ALLOC_CTOR(FV)
public:
  SPTR_FWD_CTOR(FV)

  SPHERAL_HOST_DEVICE T* data() {return SPTR_DATA_REF().data(); }
  SPHERAL_HOST_DEVICE size_t size() const { return this->sptr_data().size(); }

  SPHERAL_HOST_DEVICE size_t getHash() const { return this->sptr_data().getHash(); }
};

template<typename T>
class F : public Spheral::SpheralValueInterface<FV<T>, Fi<T>> {
  VALUE_TYPE_ALIASES((FV<T>))

public:
  VALUE_DEF_CTOR(F)
  VALUE_COPY_CTOR(F)
  VALUE_ASSIGNEMT_OP()
  VALUE_TOVIEW_OP()

  // Ctor 
  F(size_t h, size_t sz) : Base(new Fi<T>(h, sz)) {}

  void resize(size_t sz) { this->sptr_data().resize(sz); } 
};

}// namespace Spheral


TEST(FieldParallelInheritance, AccessPattern)
{

  Spheral::F<double> f(2, 200);
  //Spheral::FB fb = (Spheral::FB)f; // I don't think you should be able to do this in this case...
  //Spheral::FB* fbptr = &f; // Can not do this...
  //Spheral::FB fb2; //Can not do this...
  //Spheral::FB fb2 = fb; // Can not do this...
  auto f_v = f.toView();
  Spheral::FBV fb_v = (Spheral::FBV)f_v;

  Spheral::FV<double> f_v2 = (Spheral::FV<double>)(fb_v);

  std::cout << "fb_v : " << fb_v.getHash() << " , " << std::endl;
  //std::cout << "fb   : " << fb.getHash()<< " , " << std::endl;
  std::cout << "f_v  : " << f_v.getHash()   << " , " << f_v.data()   << " , " << f_v.size() << std::endl;
  std::cout << "f    : " << f.getHash()     << " , " << f.data()     << " , " << f.size()   << std::endl;
  std::cout << "f_v2 : " << f_v2.getHash() << " , " << f_v2.data() << std::endl;
  
  //fb.resize(1123);
  fb_v.resize(1123);
  
  std::cout << "f    : " << f.getHash()     << " , " << f.data()     << " , " << f.size()   << std::endl;
  std::cout << "f_v  : " << f_v.getHash()   << " , " << f_v.data()   << " , " << f_v.size() << std::endl;

  //std::cout << fb.size() << std::endl;
  std::cout << f_v2.size() << std::endl;

  
  Spheral::F<double> f0(0, 0);
  Spheral::F<double> f1(1, 1);
  Spheral::F<double> f2(2, 2);
  Spheral::F<double> f3(3, 3);
  Spheral::F<double> f4(4, 4);

  Spheral::ManagedVector<Spheral::FBV> vec_fbv(0);
  std::cout << "chekc\n";
  vec_fbv.push_back(  (Spheral::FBV)f0.toView()  );
  vec_fbv.push_back(  (Spheral::FBV)f1.toView()  );
  vec_fbv.push_back(  (Spheral::FBV)f2.toView()  );
  vec_fbv.push_back(  (Spheral::FBV)f3.toView()  );
  vec_fbv.push_back(  (Spheral::FBV)f4.toView()  );

  for(auto elem : vec_fbv) { std::cout << elem.size() << std::endl; elem.resize(elem.size()*2); }

  vec_fbv.free();

  std::cout << f0.size() << std::endl;
  std::cout << f1.size() << std::endl;
  std::cout << f2.size() << std::endl;
  std::cout << f3.size() << std::endl;
  std::cout << f4.size() << std::endl;

}
