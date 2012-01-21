#include <map>
#include <string>
#include <sstream>

#include "NodeList/NodeList.hh"
#include "Geometry/Dimension.hh"

// BPL includes
#include "boost/python.hpp"
using namespace boost::python;

#if defined(__APPLE__) && defined(__MACH__) && defined(__GNUC__) && __GNUC__ == 3 && __GNUC_MINOR__ == 3
// Hack workaround for a bug in Apples gcc 3.3 Xcode compiler.
namespace boost
{
  template <class T, class U>
  struct is_polymorphic<std::pair<T,U> >
    : mpl::false_
  {};
}
#endif

namespace Spheral {

//------------------------------------------------------------------------------
// Print a pair.
//------------------------------------------------------------------------------
template<typename A, typename B>
inline
const std::string
printPairType(const std::pair<A, B>& self) {
  std::stringstream result;
  result << "(" << self.first << ", " << self.second << ")" << std::ends;
  return result.str();
}

//------------------------------------------------------------------------------
// Template for wrapping pairs.
//------------------------------------------------------------------------------
template<typename A, typename B>
inline
void
wrapPairType(const char* name) {
  class_<std::pair<A, B> >(name, init<>())
    .def(init<A, B>())
    .def("__str__", printPairType<A, B>)
    .def("__repr__", printPairType<A, B>)
    .def_readwrite("first", &std::pair<A, B>::first)
    .def_readwrite("second", &std::pair<A, B>::second)
    .def(self == self)
    .def(self != self)
    ;
}

//------------------------------------------------------------------------------
// Wrap various pairs.
//------------------------------------------------------------------------------
void
wrapPairs() {
  wrapPairType<double, double>("pair_double_double");
  wrapPairType<double, std::string>("pair_double_string");

  wrapPairType<const Spheral::NodeSpace::NodeList<Dim<1> >*, std::string>("pair_NodeList1dPtr_string");
  wrapPairType<const Spheral::NodeSpace::NodeList<Dim<2> >*, std::string>("pair_NodeList2dPtr_string");
  wrapPairType<const Spheral::NodeSpace::NodeList<Dim<3> >*, std::string>("pair_NodeList3dPtr_string");

  wrapPairType<Dim<1>::Vector, Dim<1>::Vector>("pair_of_Vector1d");
  wrapPairType<Dim<2>::Vector, Dim<2>::Vector>("pair_of_Vector2d");
  wrapPairType<Dim<3>::Vector, Dim<3>::Vector>("pair_of_Vector3d");

  wrapPairType<Dim<1>::Scalar, Dim<1>::Vector>("pair_of_Scalar_Vector1d");
  wrapPairType<Dim<2>::Scalar, Dim<2>::Vector>("pair_of_Scalar_Vector2d");
  wrapPairType<Dim<3>::Scalar, Dim<3>::Vector>("pair_of_Scalar_Vector3d");
}

}
