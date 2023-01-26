#ifndef __Spheral_fakestl__
#define __Spheral_fakestl__

// Some fake STL declarations.
namespace std {
  template<typename T> class vector {
  public:
    class iterator;
    class const_iterator;
    class reverse_iterator;
    class const_reverse_iterator;
  };
  template<typename T> class list {
  public:
    class iterator;
    class const_iterator;
    class reverse_iterator;
    class const_reverse_iterator;
  };
  template<typename T> class set {
  public:
    class iterator;
    class const_iterator;
    class reverse_iterator;
    class const_reverse_iterator;
  };
  template<typename T1, typename T2> class map {
  public:
    class iterator;
    class const_iterator;
    class reverse_iterator;
    class const_reverse_iterator;
  };
  template<typename T1, typename T2> class pair;
}

// A few fake boost types as well.
namespace boost {
  namespace numeric {
    namespace ublas {
      template<typename T> class vector {
      public:
        class iterator;
        class const_iterator;
      };
    }
  }
  namespace tuples {
    struct null_type {};
    inline const null_type cnull() { return null_type(); }
    template <
      class T0 = null_type, class T1 = null_type, class T2 = null_type,
      class T3 = null_type, class T4 = null_type, class T5 = null_type,
      class T6 = null_type, class T7 = null_type, class T8 = null_type,
      class T9 = null_type>
    class tuple;
  }
  namespace python {
    class object;
  }
}

#endif
