#include <boost/geometry.hpp>

#include "Geometry/Dimension.hh"

//------------------------------------------------------------------------------
// GeomVector<2> -> Boost.Geometry
//------------------------------------------------------------------------------
namespace boost {
namespace geometry {
namespace traits {

// Adapt Spheral::GeomVector<2> to Boost.Geometry

template<> struct tag<Spheral::GeomVector<2>>                             { using type = point_tag; };
template<> struct dimension<Spheral::GeomVector<2>> : boost::mpl::int_<2> {};
template<> struct coordinate_type<Spheral::GeomVector<2>>                 { using type = double; };
template<> struct coordinate_system<Spheral::GeomVector<2>>               { using type = cs::cartesian; };

template<std::size_t Index>
struct access<Spheral::GeomVector<2>, Index> {
  static_assert(Index < 2, "Index out of dimensional range");
  using Point = Spheral::GeomVector<2>;
  using CoordinateType = typename coordinate_type<Point>::type;
  static inline CoordinateType get(Point const& p) { return p[Index]; }
  static inline void set(Point& p, CoordinateType const& value) { p[Index] = value; }
};

//   static double get(Spheral::GeomVector<2> const& p) {
//    return p.x();
//  }

//  static void set(Spheral::GeomVector<2>& p, double const& value) {
//    p.x(value);
//  }
// };

// template<>
// struct access<Spheral::GeomVector<2>, 1> {
//  static double get(Spheral::GeomVector<2> const& p) {
//    return p.y();
//  }

//  static void set(Spheral::GeomVector<2>& p, double const& value)  {
//    p.y(value);
//  }
// };

}
}
} // namespace boost::geometry::traits

