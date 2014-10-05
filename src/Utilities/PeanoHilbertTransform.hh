//---------------------------------Spheral++----------------------------------//
// PeanoHilbertTransform
//
// Helper classes to encapsulate the transforms for walking the Peano-Hilbert
// curves based on integer coordinates.
//
// Created by JMO, Sun Apr 13 23:17:29 PDT 2008
//----------------------------------------------------------------------------//
#ifndef __Spheral_PeanoHilbertTransforms__
#define __Spheral_PeanoHilbertTransforms__

#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// 2-D
//------------------------------------------------------------------------------
class PeanoHilbertTransform2d {
public:

  PeanoHilbertTransform2d(const int xx, const int xy,
                          const int yx, const int yy):
    mxx(xx), mxy(xy), 
    myx(yx), myy(yy) {
    REQUIRE(xx >= -1 and xx <= 1);
    REQUIRE(xy >= -1 and xy <= 1);
    REQUIRE(yx >= -1 and yx <= 1);
    REQUIRE(yy >= -1 and yy <= 1);
  }

  unsigned operator()(const int x,
                      const int y) {
    REQUIRE(x >= 0 and x <= 1);
    REQUIRE(y >= 0 and y <= 1);
    const int x0 = 2*x - 1;
    const int y0 = 2*y - 1;
    const int xt = mxx*x0 + mxy*y0;
    const int yt = myx*x0 + myy*y0;
    const int xi = (xt + 1)/2;
    const int yi = (yt + 1)/2;
    CHECK(xi >= 0 and xi <= 1);
    CHECK(yi >= 0 and yi <= 1);
    return morder[xi + 2*yi];
  }

  PeanoHilbertTransform2d operator()(const PeanoHilbertTransform2d& T0) {
    return PeanoHilbertTransform2d((mxx*T0.mxx + mxy*T0.myx),
                                   (mxx*T0.mxy + mxy*T0.myy),
                                   (myx*T0.mxx + myy*T0.myx),
                                   (myx*T0.mxy + myy*T0.myy));
  }

private:
  int mxx, mxy, myx, myy;
  static unsigned morder[];

  PeanoHilbertTransform2d();
};

//------------------------------------------------------------------------------
// 3-D
//------------------------------------------------------------------------------
class PeanoHilbertTransform3d {
public:

  PeanoHilbertTransform3d(const int xx, const int xy, const int xz,
                          const int yx, const int yy, const int yz,
                          const int zx, const int zy, const int zz):
    mxx(xx), mxy(xy), mxz(xz),
    myx(yx), myy(yy), myz(yz),
    mzx(zx), mzy(zy), mzz(zz) {
    REQUIRE(xx >= -1 and xx <= 1);
    REQUIRE(xy >= -1 and xy <= 1);
    REQUIRE(xz >= -1 and xz <= 1);
    REQUIRE(yx >= -1 and yx <= 1);
    REQUIRE(yy >= -1 and yy <= 1);
    REQUIRE(yz >= -1 and yz <= 1);
    REQUIRE(zx >= -1 and zx <= 1);
    REQUIRE(zy >= -1 and zy <= 1);
    REQUIRE(zz >= -1 and zz <= 1);
  }

  unsigned operator()(const int x,
                      const int y,
                      const int z) {
    REQUIRE(x >= 0 and x <= 1);
    REQUIRE(y >= 0 and y <= 1);
    REQUIRE(z >= 0 and z <= 1);
    const int x0 = 2*x - 1;
    const int y0 = 2*y - 1;
    const int z0 = 2*z - 1;
    const int xt = mxx*x0 + mxy*y0 + mxz*z0;
    const int yt = myx*x0 + myy*y0 + myz*z0;
    const int zt = mzx*x0 + mzy*y0 + mzz*z0;
    const int xi = (xt + 1)/2;
    const int yi = (yt + 1)/2;
    const int zi = (zt + 1)/2;
    CHECK(xi >= 0 and xi <= 1);
    CHECK(yi >= 0 and yi <= 1);
    CHECK(zi >= 0 and zi <= 1);
    return morder[xi + 2*yi + 4*zi];
  }

  PeanoHilbertTransform3d operator()(const PeanoHilbertTransform3d& T0) {
    return PeanoHilbertTransform3d(mxx*T0.mxx + mxy*T0.myx + mxz*T0.mzx,
                                   mxx*T0.mxy + mxy*T0.myy + mxz*T0.mzy,
                                   mxx*T0.mxz + mxy*T0.myz + mxz*T0.mzz,
                                   myx*T0.mxx + myy*T0.myx + myz*T0.mzx,
                                   myx*T0.mxy + myy*T0.myy + myz*T0.mzy,
                                   myx*T0.mxz + myy*T0.myz + myz*T0.mzz,
                                   mzx*T0.mxx + mzy*T0.myx + mzz*T0.mzx,
                                   mzx*T0.mxy + mzy*T0.myy + mzz*T0.mzy,
                                   mzx*T0.mxz + mzy*T0.myz + mzz*T0.mzz);
  }

private:
  int mxx, mxy, mxz;
  int myx, myy, myz;
  int mzx, mzy, mzz;
  static unsigned morder[];

  PeanoHilbertTransform3d();
};

}

#endif

