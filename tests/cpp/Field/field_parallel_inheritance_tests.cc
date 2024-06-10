#include "test-utilities.hh"
#include "test-basic-exec-policies.hh"


#include "Field/SphArray.hh"
#include "Utilities/ValueViewInterface.hh"
#include <memory>

namespace Spheral {


namespace impl {

class FB : chai::CHAIPoly{
  public:
    SPHERAL_HOST_DEVICE FB() {}
    FB(size_t h) : hash(h) {}

    SPHERAL_HOST_DEVICE size_t getHash() const { return hash; }

    size_t hash = 0;

    virtual void resize(size_t sz) = 0;
    SPHERAL_HOST_DEVICE virtual size_t size() = 0;
};

template<typename T>
class F : public FB, SPHERALCopyable<F<T>>{
public:
  ManagedVector<T> m_data;

  SPHERAL_HOST_DEVICE F() : FB() {}
  F(size_t h, size_t sz) : FB(h), m_data(sz) {}

  SPHERAL_HOST_DEVICE T* data() { return &m_data[0]; }
  SPHERAL_HOST_DEVICE T& operator()(size_t idx) { return m_data[idx]; }

  friend F deepCopy(F const& rhs) {
    F result(rhs);
    result.m_data = deepCopy(rhs.m_data);
    return result;
  }

  virtual void resize(size_t sz) override {std::cout << "Resize : " << sz << std::endl; m_data.resize(sz); } 
  SPHERAL_HOST_DEVICE virtual size_t size() override { return m_data.size(); }
};

} // namespace impl

// We need to forward declare child classes for use in conversion functions.
template<typename T>
class F;
class FB;

#define FBV__(code) VIEW_INTERFACE_METACLASS_DECLARATION((FB), (FBV), (impl::FB), (code))
#define FB__(code)  VALUE_INTERFACE_METACLASS_DECLARATION((FB), (FBV), (code))

#define FV__(code) VIEW_INTERFACE_METACLASS_DECLARATION((F<T>), (FV), (impl::F<T>), (code))
#define F__(code) VALUE_INTERFACE_METACLASS_DECLARATION((F), (FV<T>), (code))

class FBV__(
  DEFAULT()
);

// Because the underlying impl type is pure virtual we can not allow
// construction of FB class with default Ctor, Copy Ctor or assignment Op. 
// We can only construct a FB object from an existing smart_ptr type so 
// they must be constructed from F type objects directly...
class FB__(
  DELETED_INTERFACE(FB)
);


template<typename T>
class FV__(
  UPCAST_CONVERSION_OP(FBV)
);


template<typename T>
class F__(
public:
  VALUE_DEF_CTOR(F)
  VALUE_COPY_CTOR(F)
  VALUE_ASSIGNEMT_OP()

  // Ctor 
  F(size_t h, size_t sz) : Base(chai::make_shared<ImplType>(h, sz)) {}
  ~F() { this->sptr()->m_data.free(); }

  void resize(size_t sz) { this->sptr_data().resize(sz); } 

  // Moved from View Interface
  SPHERAL_HOST_DEVICE T* data() {return SPTR_DATA_REF().data(); }
  SPHERAL_HOST_DEVICE T& operator()(size_t idx) { return this->sptr_data().operator()(idx); }
  SPHERAL_HOST_DEVICE size_t size() const { return this->sptr_data().size(); }

  SPHERAL_HOST_DEVICE size_t getHash() const { return this->sptr_data().getHash(); }
);
  

}// namespace Spheral



//-----------------------------------------------------------------------------
// TEST SUITE
//-----------------------------------------------------------------------------

// Setting up G Test for QuadraticInterpolator
template<typename T>
class FieldParallelInheritanceTypedTest : public::testing::Test {};

// All QuadraticInterpolatorTets cases will run over each type in EXEC_TYPES.
TYPED_TEST_CASE(FieldParallelInheritanceTypedTest, EXEC_TYPES);

//TEST(FieldParallelInheritance, AccessPattern)
GPU_TYPED_TEST(FieldParallelInheritanceTypedTest, AccessPattern)
{
  {

  // --------------------------------------------------------------------------
  // Field Access 
  
  using WORK_EXEC_POLICY = TypeParam;

  Spheral::F<double> f(2, 200);

  auto f_v = &f;
  Spheral::FB::ViewType fb_v = f_v;

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
  // --------------------------------------------------------------------------
  
  // --------------------------------------------------------------------------
  // Arrays of Fields
  Spheral::F<double> f0(0, 0);
  Spheral::F<double> f1(1, 1);
  Spheral::F<double> f2(2, 2);
  Spheral::F<double> f3(3, 3);
  Spheral::F<double> f4(4, 4);

  Spheral::ManagedVector<Spheral::FB::ViewType> vec_fbv;
  vec_fbv.reserve(5);
  Spheral::ManagedVector<Spheral::F<double>::ViewType>  vec_fv;
  vec_fv.reserve(5);

  vec_fbv.push_back( &f0 );
  vec_fbv.push_back( &f1 );
  vec_fbv.push_back( &f2 );
  vec_fbv.push_back( &f3 );
  vec_fbv.push_back( &f4 );

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
  // --------------------------------------------------------------------------

  }
  std::cout << "Sptr Map Sz : " << chai::SharedPtrManager::getInstance()->getPointerMap().size() << std::endl;
  std::cout << "Arr  Map Sz : " << chai::ArrayManager::getInstance()->getPointerMap().size() << std::endl;
}
