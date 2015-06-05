//---------------------------------Spheral++----------------------------------//
// boundPointWithinBox
//
// Limit a position to be within a box.
//
// Created by JMO, Tue Mar 30 10:00:59 PDT 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_boundPointWithinBox__
#define __Spheral_boundPointWithinBox__

#include "Geometry/Dimension.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// 1-D
//------------------------------------------------------------------------------
inline
Dim<1>::Vector
boundPointWithinBox(const Dim<1>::Vector& point,
                    const Dim<1>::Vector& xmin, const Dim<1>::Vector& xmax) {
  return Dim<1>::Vector(std::max(xmin.x(), std::min(xmax.x(), point.x())));
}

//------------------------------------------------------------------------------
// 2-D
//------------------------------------------------------------------------------
inline
Dim<2>::Vector
boundPointWithinBox(const Dim<2>::Vector& point,
                    const Dim<2>::Vector& xmin, const Dim<2>::Vector& xmax) {
  return Dim<2>::Vector(std::max(xmin.x(), std::min(xmax.x(), point.x())),
                        std::max(xmin.y(), std::min(xmax.y(), point.y())));
}

//------------------------------------------------------------------------------
// 3-D
//------------------------------------------------------------------------------
inline
Dim<3>::Vector
boundPointWithinBox(const Dim<3>::Vector& point,
                    const Dim<3>::Vector& xmin, const Dim<3>::Vector& xmax) {
  return Dim<3>::Vector(std::max(xmin.x(), std::min(xmax.x(), point.x())),
                        std::max(xmin.y(), std::min(xmax.y(), point.y())),
                        std::max(xmin.z(), std::min(xmax.z(), point.z())));
}

}

#endif

