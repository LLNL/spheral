#include "test-utilities.hh"
#include "test-basic-exec-policies.hh"


#include "Utilities/ValueViewInterface.hh"
#include "Utilities/ManagedVector.hh"
#include <memory>

namespace Spheral {


//-----------------------------------------------------------------------------
// Implementation of Existing Class in Spheral...
//-----------------------------------------------------------------------------

VVI_IMPL_BEGIN

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
class F : public FB {
public:
  vvi::vector<T> m_data;

  SPHERAL_HOST_DEVICE F() : FB() {}
  F(size_t h, size_t sz) : FB(h), m_data(sz) {}

  SPHERAL_HOST_DEVICE T* data() { return &m_data[0]; }
  SPHERAL_HOST_DEVICE T& operator()(size_t idx) { return m_data[idx]; }

  virtual void resize(size_t sz) override {std::cout << "Resize : " << sz << std::endl; m_data.resize(sz); } 
  SPHERAL_HOST_DEVICE virtual size_t size() override { return m_data.size(); }

  VVI_IMPL_DEEPCOPY(F, m_data)
  VVI_IMPL_COMPARE(F, m_data, hash)
};

VVI_IMPL_END



#ifdef VVI_ENABLED
//-----------------------------------------------------------------------------
// Interface to support porting to the GPU.
//-----------------------------------------------------------------------------

// We need to forward declare value classes for view interface definitions.
template<typename T>
class F;
class FB;

// Define Metaclass macros for Value/View relationships
class PTR_VIEW_METACLASS_DEFAULT((FB), (FBV), (vvimpl::FB))
class PTR_VALUE_METACLASS_DELETED((FB), (FBV), (vvimpl::FB))

#define FV__(code) PTR_VIEW_METACLASS_DECL((F<T>), (FV), (vvimpl::F<T>), (code))
#define F__(code) PTR_VALUE_METACLASS_DECL((F), (FV<T>), (code))


template<typename T>
class FV__(
  // Field inherits from FieldBase so we would like to be able to implicitly
  // upcast a fieldview object to a fieldbaseview.
  VVI_UPCAST_CONVERSION_OP(FBV)
);


template<typename T>
class F__(
public:
  // Custom Ctor, note we need to create the underlying implementation 
  // object on ctor of value interfaces.
  F(size_t h, size_t sz) : VVI_VALUE_CTOR_ARGS( (h, sz) ) {}

  // Value semantics dictate that we free underlying data upon destruction.
  ~F() { VVI_IMPL_INST().m_data.free(); }
  
  // HOST only interface
  void resize(size_t sz) { VVI_IMPL_INST().resize(sz); } 

  // Moved from old View Interface pattern.
  T* data() {return VVI_IMPL_INST().data(); }
  T& operator()(size_t idx) { return VVI_IMPL_INST().operator()(idx); }
  size_t size() const { return VVI_IMPL_INST().size(); }

  size_t getHash() const { return VVI_IMPL_INST().getHash(); }
);
  
#endif // !defined(SPHERAL_ENABLE_VVI)

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
  //auto f_2 = f;
  Spheral::F<double> f_2(2, 200);
  std::cout << (f == f_2) << std::endl;

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

  vvi::vector<Spheral::FB::ViewType> vec_fbv;
  vec_fbv.reserve(5);
  vvi::vector<Spheral::F<double>::ViewType>  vec_fv;
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
  for(size_t i = 0; i < vec_fbv.size(); i++) SPHERAL_ASSERT_EQ(vec_fbv[i]->size(), (*vec_fv[i]).size());
  std::cout << "Arr  Map Sz : " << chai::ArrayManager::getInstance()->getPointerMap().size() << std::endl;

  for(size_t i = 0; i < vec_fbv.size(); i++)
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
    for(size_t i = 0; i < vec_fbv.size(); i++)
    {
      auto& elem_b = *vec_fbv[i];
      auto& elem_v = *vec_fv[i];

      printf("%ld, %ld\n", elem_b.getHash(), elem_b.size());
      printf("%p\n", (void*)&elem_v(0));
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
