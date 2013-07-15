//---------------------------------Spheral++----------------------------------//
// testBoxIntersection
//
// Some simple methods specialized by dimension to rapidly test if two boxes
// intersect.
//
// Created by JMO, Sun Feb  6 13:44:49 PST 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral_testBoxIntersection__
#define __Spheral_testBoxIntersection__

#include "Geometry/Dimension.hh"
#include "spheralWildMagicConverters.hh"
#include "Wm5IntrBox2Box2.h"
#include "Wm5IntrBox3Box3.h"
#include "Wm5ContBox2.h"
#include "Wm5ContBox3.h"

namespace Spheral {

//------------------------------------------------------------------------------
// These methods specify the boxes by the their min and max coordinates, and
// assume that the boxes are aligned with the coordinate system.
//------------------------------------------------------------------------------
inline
bool
testBoxIntersection(const Dim<1>::Vector& xmin1, const Dim<1>::Vector& xmax1,
                    const Dim<1>::Vector& xmin2, const Dim<1>::Vector& xmax2,
                    const double tol = 1.0e-10) {
  if (xmax1.x() < (xmin2.x() - tol)) return false;
  if (xmax2.x() < (xmin1.x() - tol)) return false;
  return true;
}

inline
bool
testBoxIntersection(const Dim<2>::Vector& xmin1, const Dim<2>::Vector& xmax1,
                    const Dim<2>::Vector& xmin2, const Dim<2>::Vector& xmax2,
                    const double tol = 1.0e-10) {
  if (xmax1.x() < (xmin2.x() - tol)) return false;
  if (xmax2.x() < (xmin1.x() - tol)) return false;

  if (xmax1.y() < (xmin2.y() - tol)) return false;
  if (xmax2.y() < (xmin1.y() - tol)) return false;

  return true;
}

inline
bool
testBoxIntersection(const Dim<3>::Vector& xmin1, const Dim<3>::Vector& xmax1,
                    const Dim<3>::Vector& xmin2, const Dim<3>::Vector& xmax2,
                    const double tol = 1.0e-10) {
  if (xmax1.x() < (xmin2.x() - tol)) return false;
  if (xmax2.x() < (xmin1.x() - tol)) return false;

  if (xmax1.y() < (xmin2.y() - tol)) return false;
  if (xmax2.y() < (xmin1.y() - tol)) return false;

  if (xmax1.z() < (xmin2.z() - tol)) return false;
  if (xmax2.z() < (xmin1.z() - tol)) return false;

  return true;
}

//------------------------------------------------------------------------------
// Test if a point is contained in a box.
//------------------------------------------------------------------------------
inline
bool
testPointInBox(const Dim<1>::Vector& point,
               const Dim<1>::Vector& xmin, const Dim<1>::Vector& xmax,
               const double tol = 1.0e-10) {
  if (point.x() < (xmin.x() - tol)) return false;
  if (point.x() > (xmax.x() + tol)) return false;
  return true;
}

inline
bool
testPointInBox(const Dim<2>::Vector& point,
               const Dim<2>::Vector& xmin, const Dim<2>::Vector& xmax,
               const double tol = 1.0e-10) {
  if (point.x() < (xmin.x() - tol)) return false;
  if (point.x() > (xmax.x() + tol)) return false;

  if (point.y() < (xmin.y() - tol)) return false;
  if (point.y() > (xmax.y() + tol)) return false;

  return true;
}

inline
bool
testPointInBox(const Dim<3>::Vector& point,
               const Dim<3>::Vector& xmin, const Dim<3>::Vector& xmax,
               const double tol = 1.0e-10) {
  if (point.x() < (xmin.x() - tol)) return false;
  if (point.x() > (xmax.x() + tol)) return false;

  if (point.y() < (xmin.y() - tol)) return false;
  if (point.y() > (xmax.y() + tol)) return false;

  if (point.z() < (xmin.z() - tol)) return false;
  if (point.z() > (xmax.z() + tol)) return false;

  return true;
}

//------------------------------------------------------------------------------
// The more general case of arbitrary boxes not necessarily axis aligned.
//------------------------------------------------------------------------------
inline
bool
testBoxIntersection(const Dim<1>::Box& box1, const Dim<1>::Box& box2) {
  return testBoxIntersection(box1.center() - box1.extent(),
                             box1.center() + box1.extent(),
                             box2.center() - box2.extent(),
                             box2.center() + box2.extent());
}

inline
bool
testBoxIntersection(const Dim<2>::Box& box1, const Dim<2>::Box& box2) {
  Wm5::IntrBox2Box2<double> test(box1, box2);
  return test.Test();
}

inline
bool
testBoxIntersection(const Dim<3>::Box& box1, const Dim<3>::Box& box2) {
  Wm5::IntrBox3Box3<double> test(box1, box2);
  return test.Test();
}

//------------------------------------------------------------------------------
// A similar in spirit test is if a point is contained in a box.
//------------------------------------------------------------------------------
inline
bool
testPointInBox(const Dim<1>::Vector& point, const Dim<1>::Box& box) {
  return (point.x() >= box.center().x() - box.extent() and
          point.x() <= box.center().x() + box.extent());
}

inline
bool
testPointInBox(const Dim<2>::Vector& point, const Dim<2>::Box& box) {
  return Wm5::InBox(convertVectorToWMVector<Dim<2> >(point), box);
}

inline
bool
testPointInBox(const Dim<3>::Vector& point, const Dim<3>::Box& box) {
  return Wm5::InBox(convertVectorToWMVector<Dim<3> >(point), box);
}

}

#endif

