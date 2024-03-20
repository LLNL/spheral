
#include "test-utilities.hh"
#include "test-basic-exec-policies.hh"

#include <memory>
#include <stdio.h>

#include "Utilities/SharedPtr.hh"

template<typename T>
using shared_ptr_t = Spheral::shared_ptr<T>;
//using shared_ptr_t = Spheral::managed_shared_ptr<T>;
//using shared_ptr_t = std::shared_ptr<T>;


class Base {
public:
  Base() {printf("CTor Base\n");}
  void doSomething(){ printf("Base doSomething\n"); }
  ~Base() { printf("DTor Base\n"); }
};

class Derived : public Base {
public:
  Derived() {printf("CTor Derived\n");}
  void doSomething(){ printf("Derived doSomething\n"); }
  ~Derived() { printf("DTor Derived\n"); }
};


TEST(SharedPtrTest, BasicCount)
{
  shared_ptr_t< Derived > s_ptr = shared_ptr_t<Derived>(new Derived());

  s_ptr->doSomething();
  SPHERAL_ASSERT_EQ(s_ptr.use_count(), 1);

  shared_ptr_t< Derived > s_ptr2 = s_ptr;

  s_ptr2->doSomething();

  SPHERAL_ASSERT_EQ(s_ptr.use_count(), 2);
  SPHERAL_ASSERT_EQ(s_ptr2.use_count(), 2);

  SPHERAL_ASSERT_EQ(s_ptr.get(), s_ptr2.get());
}

TEST(SharedPtrTest, UpcastCtor)
{
  shared_ptr_t< Base > s_ptr = shared_ptr_t<Base>(new Derived());

  s_ptr->doSomething();
  SPHERAL_ASSERT_EQ(s_ptr.use_count(), 1);
}

TEST(SharedPtrTest, UpcastCopyCtor)
{
  shared_ptr_t< Derived > s_ptr = shared_ptr_t<Derived>(new Derived());

  s_ptr->doSomething();
  SPHERAL_ASSERT_EQ(s_ptr.use_count(), 1);

  shared_ptr_t< Base > s_ptr_b = s_ptr;

  s_ptr_b->doSomething();

  SPHERAL_ASSERT_EQ(s_ptr.use_count(), 2);
  SPHERAL_ASSERT_EQ(s_ptr_b.use_count(), 2);

  SPHERAL_ASSERT_EQ(s_ptr.get(), s_ptr_b.get());
}

TEST(SharedPtrTest, SharedCount)
{
  
  //Spheral::__shared_count(new Derived(), [](Derived* d){printf("Custom Deleter\n"); d->~Derived();});
  Spheral::__shared_ptr<Derived> sptr(new Derived(), [](Derived* d){printf("Custom Deleter\n"); d->~Derived();});



}


// Setting up G Test for SharedPtr
template<typename T>
class SharedPtrTypedTest : public::testing::Test {};

// All SharedPtrTets cases will run over each type in EXEC_TYPES.
TYPED_TEST_CASE(SharedPtrTypedTest, EXEC_TYPES);

