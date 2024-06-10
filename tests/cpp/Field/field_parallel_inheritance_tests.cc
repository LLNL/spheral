#include "test-utilities.hh"
#include "test-basic-exec-policies.hh"


#include "Field/SphArray.hh"
#include "Utilities/ValueViewInterface.hh"
#include <memory>

namespace Spheral {

// New macro for this implementation.
#define SPTR_FWD_CTOR(type) \
  type(SmartPtrType&& rhs) : Base(std::forward<SmartPtrType>(rhs)) {}

#define UPCAST_CONVERSION_OP(parent_t) \
  explicit operator parent_t() const {return parent_t(this->sptr());}


class FBi : chai::CHAIPoly{
  public:
    SPHERAL_HOST_DEVICE FBi() {}
    FBi(size_t h) : hash(h) {}


    SPHERAL_HOST_DEVICE size_t getHash() const { return hash; }

    size_t hash = 0;

    virtual void resize(size_t sz) = 0;
    SPHERAL_HOST_DEVICE virtual size_t size() = 0;
};

template<typename T>
class Fi : public FBi, Spheral::SPHERALCopyable<Fi<T>>{
public:
  Spheral::ManagedVector<T> m_data;

  SPHERAL_HOST_DEVICE Fi() : FBi() {}
  Fi(size_t h, size_t sz) : FBi(h), m_data(sz) {}

  SPHERAL_HOST_DEVICE T* data() { return &m_data[0]; }
  SPHERAL_HOST_DEVICE T& operator()(size_t idx) { return m_data[idx]; }

  friend Fi deepCopy(Fi const& rhs) {
    Fi result(rhs);
    result.m_data = Spheral::deepCopy(rhs.m_data);
    return result;
  }

  virtual void resize(size_t sz) override {std::cout << "Resize : " << sz << std::endl; m_data.resize(sz); } 
  SPHERAL_HOST_DEVICE virtual size_t size() override { return m_data.size(); }
};


// We need to forward declare child classes for use in conversion functions.
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
  SPHERAL_HOST_DEVICE FBV& operator=(FBV const&) = default;

  SPHERAL_HOST_DEVICE ImplType& operator*() const { return SPTR_DATA_REF(); }
  SPHERAL_HOST_DEVICE ImplType* operator->() const { return &SPTR_DATA_REF(); }

  void shallowCopy(FBV const& rhs) {*this = rhs;}

  //SPHERAL_HOST_DEVICE size_t getHash() const { return sptr_data().getHash(); }

  //// We want to forward the pure virtual interface here.
  //SPHERAL_HOST_DEVICE size_t size() const {return sptr_data().size(); }
  //void resize(size_t sz) { sptr_data().resize(sz); }
};

// Because the underlying impl type is pure virtual we can not allow
// construction of FB class with default Ctor, Copy Ctor or assignment Op. 
// We can only construct a FB object from an existing smart_ptr type so 
// they must be constructed from F type objects directly...
class FB : public Spheral::SpheralValueInterface<FBV> {
  VALUE_TYPE_ALIASES((FBV))
  //SPTR_FWD_CTOR(FB)
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
  VIEW_DEF_CTOR(FV)

  SPHERAL_HOST_DEVICE ImplType& operator*() { return SPTR_DATA_REF(); }
  SPHERAL_HOST_DEVICE ImplType* operator->() { return &SPTR_DATA_REF(); }

  void shallowCopy(FV const& rhs) {*this = rhs;}
};


template<typename T>
class F : public Spheral::SpheralValueInterface<FV<T>> {
  VALUE_TYPE_ALIASES((FV<T>))

public:
  VALUE_DEF_CTOR(F)
  VALUE_COPY_CTOR(F)
  VALUE_ASSIGNEMT_OP()
  VALUE_TOVIEW_OP()

  ViewType operator&() { return toView(); }

  // Ctor 
  F(size_t h, size_t sz) : Base(chai::make_shared<Fi<T>>(h, sz)) {}
  ~F() { this->sptr()->m_data.free(); }

  void resize(size_t sz) { this->sptr_data().resize(sz); } 

  // Moved from View Interface
  SPHERAL_HOST_DEVICE T* data() {return SPTR_DATA_REF().data(); }
  SPHERAL_HOST_DEVICE T& operator()(size_t idx) { return this->sptr_data().operator()(idx); }
  SPHERAL_HOST_DEVICE size_t size() const { return this->sptr_data().size(); }

  SPHERAL_HOST_DEVICE size_t getHash() const { return this->sptr_data().getHash(); }
};

}// namespace Spheral

// Setting up G Test for QuadraticInterpolator
template<typename T>
class FieldParallelInheritanceTypedTest : public::testing::Test {};

// All QuadraticInterpolatorTets cases will run over each type in EXEC_TYPES.
TYPED_TEST_CASE(FieldParallelInheritanceTypedTest, EXEC_TYPES);


//TEST(FieldParallelInheritance, AccessPattern)
GPU_TYPED_TEST(FieldParallelInheritanceTypedTest, AccessPattern)
{
  {

  using WORK_EXEC_POLICY = TypeParam;

  Spheral::F<double> f(2, 200);

  auto f_v = &f;
  Spheral::FBV fb_v = (Spheral::FBV)f_v;

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    printf("--- GPU BEGIN ---\n");
    printf("%ld, %ld\n", fb_v->getHash(), fb_v->size());
    printf("--- GPU END ---\n");
  EXEC_IN_SPACE_END()

  std::cout << "fb_v : " << fb_v->getHash() << " , " << std::endl;
  std::cout << "f_v  : " << f_v->getHash()   << " , " << f_v->data()   << " , " << (*f_v).size() << std::endl;
  std::cout << "f    : " << f.getHash()     << " , " << f.data()     << " , " << f.size()   << std::endl;
  
  fb_v->resize(1123);
  
  std::cout << "f    : " << f.getHash()     << " , " << f.data()     << " , " << f.size()   << std::endl;
  std::cout << "f_v  : " << f_v->getHash()   << " , " << f_v->data()   << " , " << f_v->size() << std::endl;

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    printf("--- GPU BEGIN ---\n");
    printf("%ld, %ld\n", fb_v->getHash(), fb_v->size());
    printf("--- GPU END ---\n");
  EXEC_IN_SPACE_END()
  
  Spheral::F<double> f0(0, 0);
  Spheral::F<double> f1(1, 1);
  Spheral::F<double> f2(2, 2);
  Spheral::F<double> f3(3, 3);
  Spheral::F<double> f4(4, 4);

  Spheral::ManagedVector<Spheral::FBV> vec_fbv;
  vec_fbv.reserve(5);
  Spheral::ManagedVector<Spheral::FV<double>>  vec_fv;
  vec_fv.reserve(5);

  vec_fbv.push_back(  (Spheral::FBV)(&f0)  );
  vec_fbv.push_back(  (Spheral::FBV)(&f1)  );
  vec_fbv.push_back(  (Spheral::FBV)(&f2)  );
  vec_fbv.push_back(  (Spheral::FBV)(&f3)  );
  vec_fbv.push_back(  (Spheral::FBV)(&f4)  );

  vec_fv.push_back( &f0 );
  vec_fv.push_back( &f1 );
  vec_fv.push_back( &f2 );
  vec_fv.push_back( &f3 );
  vec_fv.push_back( &f4 );

  SPHERAL_ASSERT_EQ(vec_fbv.size(), vec_fv.size());
  for(int i = 0; i < vec_fbv.size(); i++) SPHERAL_ASSERT_EQ(vec_fbv[i]->size(), (*vec_fv[i]).size());
  std::cout << "Arr  Map Sz : " << chai::ArrayManager::getInstance()->getPointerMap().size() << std::endl;

  for(int i = 0; i < vec_fbv.size(); i++)
  {
    auto& elem = *vec_fbv[i];

    printf("%ld, %ld\n", elem.getHash(), elem.size());
    size_t sz = elem.size();
    elem.resize(sz*2);
  }
  std::cout << "Arr  Map Sz : " << chai::ArrayManager::getInstance()->getPointerMap().size() << std::endl;

  for(size_t i = 0; i < f0.size(); i++) { f0(i) = f0.getHash(); }
  for(size_t i = 0; i < f1.size(); i++) { f1(i) = f1.getHash(); }
  for(size_t i = 0; i < f2.size(); i++) { f2(i) = f2.getHash(); }
  for(size_t i = 0; i < f3.size(); i++) { f3(i) = f3.getHash(); }
  for(size_t i = 0; i < f4.size(); i++) { f4(i) = f4.getHash(); }


  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    printf("--- GPU BEGIN ---\n");
    printf("%ld\n", vec_fv.size());
    for(int i = 0; i < vec_fbv.size(); i++)
    {
      auto& elem_b = *vec_fbv[i];
      auto& elem_v = *vec_fv[i];
      //elem_v.resize(120);

      printf("%ld, %ld\n", elem_b.getHash(), elem_b.size());
      printf("%p\n", &elem_v(0));
      for(size_t j = 0; j < elem_v.size(); j++) { printf("%f, ",vec_fv[i]->operator()(j));} 
      printf("\n");
    }
    printf("--- GPU END ---\n");
  EXEC_IN_SPACE_END()

  vec_fbv.free();
  vec_fv.free();

  std::cout << f0.size() << std::endl;
  std::cout << f1.size() << std::endl;
  std::cout << f2.size() << std::endl;
  std::cout << f3.size() << std::endl;
  std::cout << f4.size() << std::endl;

  }
  std::cout << "Sptr Map Sz : " << chai::SharedPtrManager::getInstance()->getPointerMap().size() << std::endl;
  std::cout << "Arr  Map Sz : " << chai::ArrayManager::getInstance()->getPointerMap().size() << std::endl;
}
