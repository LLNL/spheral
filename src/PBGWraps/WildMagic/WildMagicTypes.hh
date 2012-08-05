#ifndef __WM_WildMagicTypes__
#define __WM_WildMagicTypes__

#include <string>
#include <sstream>

#include "Wm5Vector2.h"
#include "Wm5Vector3.h"
#include "Wm5Box2.h"
#include "Wm5Box3.h"

namespace Wm5 {

typedef Vector2<double> WMVector2d;
typedef Vector3<double> WMVector3d;

typedef Box2<double> WMBox2d;
typedef Box3<double> WMBox3d;

//------------------------------------------------------------------------------
// Nice string representations (Vectors)
//------------------------------------------------------------------------------
inline
std::string
printReprWMVector2d(const WMVector2d& self) {
  std::stringstream s;
  s << "WMVector2d(" 
    << self.X() << " "
    << self.Y() << ")";
  return s.str();
}

inline
std::string
printReprWMVector3d(const WMVector3d& self) {
  std::stringstream s;
  s << "WMVector3d(" 
    << self.X() << " "
    << self.Y() << " "
    << self.Z() << ")";
  return s.str();
}

//------------------------------------------------------------------------------
// Nice string representations (Boxes)
//------------------------------------------------------------------------------
inline
std::string
printReprWMBox2d(const WMBox2d& self) {
  std::stringstream s;
  s << "WMBox2d(" 
    << "center=(" 
    << self.Center.X() << " "
    << self.Center.Y() << "), "
    << "axes=[" 
    << "(" << self.Axis[0].X() << " " << self.Axis[0].Y() << "), "
    << "(" << self.Axis[1].X() << " " << self.Axis[1].Y() << ")], "
    << "extents=[" << self.Extent[0] << " " << self.Extent[1] << "])";
  return s.str();
}

inline
std::string
printReprWMBox3d(const WMBox3d& self) {
  std::stringstream s;
  s << "WMBox3d(" 
    << "center=(" 
    << self.Center.X() << " "
    << self.Center.Y() << " "
    << self.Center.Z() << "), "
    << "axes=[" 
    << "(" << self.Axis[0].X() << " " << self.Axis[0].Y() << " " << self.Axis[0].Z() << "), "
    << "(" << self.Axis[1].X() << " " << self.Axis[1].Y() << " " << self.Axis[1].Z() << ")], "
    << "extents=[" << self.Extent[0] << " " << self.Extent[1] << " " << self.Extent[2] << "])";
  return s.str();
}

}

#endif
