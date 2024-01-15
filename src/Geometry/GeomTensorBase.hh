//---------------------------------Spheral++----------------------------------//
// GeomTensorBase -- provides the dimension dependent storage for GeomTensor.
//
// Created by JMO, Wed Nov 10 16:40:29 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral__GeomTensorBase__
#define __Spheral__GeomTensorBase__

namespace Spheral {

template<int nDim> class GeomTensorBase {};

template<>
class GeomTensorBase<1> {
 public:
  GeomTensorBase(const double xx):
    mxx(xx) {}
 protected:
  double mxx;
 private:
  GeomTensorBase();
};

template<>
class GeomTensorBase<2> {
 public:
  GeomTensorBase(const double a):
    mxx(a),
    mxy(a),
    myx(a),
    myy(a) {}
  GeomTensorBase(const double xx, const double xy,
                 const double yx, const double yy):
    mxx(xx),
    mxy(xy),
    myx(yx),
    myy(yy) {}
 protected:
  double mxx, mxy,
         myx, myy;
 private:
  GeomTensorBase();
};

template<>
class GeomTensorBase<3> {
 public:
  GeomTensorBase(const double a):
    mxx(a),
    mxy(a),
    mxz(a),
    myx(a),
    myy(a),
    myz(a),
    mzx(a),
    mzy(a), 
    mzz(a) {}
  GeomTensorBase(const double xx, const double xy, const double xz,
                 const double yx, const double yy, const double yz,
                 const double zx, const double zy, const double zz):
    mxx(xx),
    mxy(xy),
    mxz(xz),
    myx(yx),
    myy(yy),
    myz(yz),
    mzx(zx),
    mzy(zy), 
    mzz(zz) {}
 protected:
  double mxx, mxy, mxz,
         myx, myy, myz,
         mzx, mzy, mzz;
 private:
  GeomTensorBase();
};

}

#endif
