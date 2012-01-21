// A set of faked up Boost classes for use with GCCXML, since GCCXML is kind of 
// broken with some of the Boost code.

#ifndef __FAKEBOOST__
#define __FAKEBOOST__

namespace boost {
  template<typename T> class shared_ptr;
  template<typename T> class weak_ptr;
  template<typename T> T& ref(T& x);
  namespace python {
    struct list {
      list();
      template<typename T> list(T x);
      template<typename T> void append(T);
    };
    struct tuple { template<typename X> X operator[](unsigned index); };
    struct object {};
    template<typename T> unsigned len(T x);
    template<typename T> T extract(object);
    template<typename T1> tuple make_tuple(T1 x1);
    template<typename T1, typename T2> tuple make_tuple(T1 x1, T2 x2);
    template<typename T1, typename T2, typename T3> tuple make_tuple(T1 x1, T2 x2, T3 x3);
  }
}

#endif
