//---------------------------------Spheral++----------------------------------//
// GeomSymmetricTensorBase -- provides the dimension dependent storage for
// GeomSymmetricTensor.
//
// Created by JMO, Wed Nov 10 21:33:17 PST 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral__GeomSymmetricTensorBase__
#define __Spheral__GeomSymmetricTensorBase__

namespace Spheral {

template<int nDim> class GeomSymmetricTensorBase {};

template<>
class GeomSymmetricTensorBase<1> {
 public:
  GeomSymmetricTensorBase(const double xx):
    mxx(xx) {}
 protected:
  double mxx;
 private:
  GeomSymmetricTensorBase();
};

template<>
class GeomSymmetricTensorBase<2> {
 public:
  GeomSymmetricTensorBase(const double a):
    mxx(a),
    mxy(a),
    myy(a) {}
  GeomSymmetricTensorBase(const double xx, const double xy,
                                           const double yy):
    mxx(xx),
    mxy(xy),
    myy(yy) {}
 protected:
  double mxx, mxy,
              myy;
 private:
  GeomSymmetricTensorBase();
};

template<>
class GeomSymmetricTensorBase<3> {
 public:
  GeomSymmetricTensorBase(const double a):
    mxx(a),
    mxy(a),
    mxz(a),
    myy(a),
    myz(a),
    mzz(a) {}
  GeomSymmetricTensorBase(const double xx, const double xy, const double xz,
                                           const double yy, const double yz,
                                                            const double zz):
    mxx(xx),
    mxy(xy),
    mxz(xz),
    myy(yy),
    myz(yz),
    mzz(zz) {}
 protected:
  double mxx, mxy, mxz,
              myy, myz,
                   mzz;
 private:
  GeomSymmetricTensorBase();
};

}

#endif
