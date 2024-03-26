#include <boost/geometry.hpp>

#include "Geometry/Dimension.hh"

//------------------------------------------------------------------------------
// GeomVector<2> -> Boost.Geometry
//------------------------------------------------------------------------------
namespace boost {
namespace geometry {
namespace traits {

// Adapt Spheral::GeomVector<2> to Boost.Geometry

template<> struct tag<Spheral::GeomVector<2>>                             { typedef point_tag type; };
template<> struct coordinate_type<Spheral::GeomVector<2>>                 { typedef double type; };
template<> struct coordinate_system<Spheral::GeomVector<2>>               { typedef cs::cartesian type; };
template<> struct dimension<Spheral::GeomVector<2>> : boost::mpl::int_<2> {};

template<>
struct access<Spheral::GeomVector<2>, 0> {
 static double get(Spheral::GeomVector<2> const& p) {
   return p.x();
 }

 static void set(Spheral::GeomVector<2>& p, double const& value) {
   p.x(value);
 }
};

template<>
struct access<Spheral::GeomVector<2>, 1> {
 static double get(Spheral::GeomVector<2> const& p) {
   return p.y();
 }

 static void set(Spheral::GeomVector<2>& p, double const& value)  {
   p.y(value);
 }
};

}
}
} // namespace boost::geometry::traits

