#ifndef __PBGWRAP_UTILITIESTYPES__
#define __PBGWRAP_UTILITIESTYPES__

#include "boost/math/special_functions/legendre.hpp"
#include "Utilities/erff.hh"
#include "Utilities/newtonRaphson.hh"
#include "Utilities/simpsonsIntegration.hh"
#include "Utilities/globalNodeIDs.hh"
#include "Utilities/rotationMatrix.hh"
#include "Utilities/iterateIdealH.hh"
#include "Utilities/mortonOrderIndicies.hh"
#include "Utilities/peanoHilbertOrderIndicies.hh"
#include "Utilities/boundingBox.hh"
#include "Utilities/globalBoundingVolumes.hh"
#include "Utilities/testBoxIntersection.hh"
#include "Utilities/lineSegmentIntersections.hh"
#include "Utilities/pointDistances.hh"
#include "Utilities/segmentIntersectEdges.hh"
#include "Utilities/integrateThroughMeshAlongSegment.hh"
#include "Utilities/numberDensity.hh"
#include "Utilities/writeRectilinearMesh.hh"
#include "Utilities/planarReflectingOperator.hh"
#include "Utilities/KeyTraits.hh"
#include "Utilities/pointInPolygon.hh"
#include "Utilities/pointOnPolygon.hh"
#include "Utilities/pointInPolyhedron.hh"
#include "Utilities/pointOnPolyhedron.hh"
#include "Utilities/refinePolyhedron.hh"

typedef std::pair<double, double> pair_double_double;
using namespace Spheral::FieldSpace;

namespace Spheral {

// Expose some stuff from boost::math
using boost::math::legendre_p;

//------------------------------------------------------------------------------
// An overridable python class functor for use with the Newton-Raphson root 
// finder.
//------------------------------------------------------------------------------
class NewtonRaphsonFunction {
public:
  NewtonRaphsonFunction() {};
  virtual ~NewtonRaphsonFunction() {};
  virtual std::pair<double, double> operator()(const double x) const { return __call__(x); }
  virtual std::pair<double, double> __call__(const double x) const = 0;
};

// The method we expose to python.
inline
double newtonRaphsonFindRoot(const NewtonRaphsonFunction& functor,
                             float x1,
                             float x2,
                             const float xaccuracy = 1.0e-15,
                             const float yaccuracy = 1.0e-15,
                             const unsigned maxIterations = 100) {
  return Spheral::newtonRaphson(functor,
                                x1,
                                x2,
                                xaccuracy,
                                yaccuracy,
                                maxIterations);
}

//------------------------------------------------------------------------------
// An overridable python class functor for use with the Simpsons rule numerical
// integation function.
//------------------------------------------------------------------------------
class SimpsonsIntegrationDoubleFunction {
public:
  SimpsonsIntegrationDoubleFunction() {};
  virtual ~SimpsonsIntegrationDoubleFunction() {};
  virtual double operator()(const double x) const { return __call__(x); }
  virtual double __call__(const double x) const = 0;
};

// The method we expose to python.
inline
double simpsonsIntegrationDouble(const SimpsonsIntegrationDoubleFunction& functor,
                                 const double x0,
                                 const double x1,
                                 const unsigned numBins) {
  return Spheral::simpsonsIntegration<SimpsonsIntegrationDoubleFunction, double, double>(functor,
                                                                                         x0,
                                                                                         x1,
                                                                                         numBins);
}

//------------------------------------------------------------------------------
// Wrapper around our packElement method more convenient for exposing to Python.
//------------------------------------------------------------------------------
template<typename Element>
inline
string
convertElementToString(const Element& x) {
  vector<char> buffer;
  packElement(x, buffer);
  return string(buffer.begin(), buffer.end());
}

//------------------------------------------------------------------------------
// Wrapper around our unpackElement method more convenient for exposing to
//  Python.
//------------------------------------------------------------------------------
template<typename Element>
inline
Element
convertStringToElement(const string& x) {
  Element result;
  vector<char> buffer(x.begin(), x.end());
  vector<char>::const_iterator itr = buffer.begin();
  unpackElement(result, itr, buffer.end());
  return result;
}

//------------------------------------------------------------------------------
// Stupidly enough we have to disambiguate these functions ourselves.
//------------------------------------------------------------------------------
inline
Dim<1>::Tensor
rotationMatrix1d(const Dim<1>::Vector& runit) {
  return rotationMatrix(runit);
}

inline
Dim<2>::Tensor
rotationMatrix2d(const Dim<2>::Vector& runit) {
  return rotationMatrix(runit);
}

inline
Dim<3>::Tensor
rotationMatrix3d(const Dim<3>::Vector& runit) {
  return rotationMatrix(runit);
}

//------------------------------------------------------------------------------
// And these...
//------------------------------------------------------------------------------
inline
FieldSpace::FieldList<Dim<1>, uint64_t>
mortonOrderIndicies1d(const DataBaseSpace::DataBase<Dim<1> >& dataBase) {
  return mortonOrderIndicies(dataBase);
}

inline
FieldSpace::FieldList<Dim<2>, uint64_t>
mortonOrderIndicies2d(const DataBaseSpace::DataBase<Dim<2> >& dataBase) {
  return mortonOrderIndicies(dataBase);
}

inline
FieldSpace::FieldList<Dim<3>, uint64_t>
mortonOrderIndicies3d(const DataBaseSpace::DataBase<Dim<3> >& dataBase) {
  return mortonOrderIndicies(dataBase);
}

//------------------------------------------------------------------------------
// And these...
//------------------------------------------------------------------------------
inline
FieldSpace::FieldList<Dim<1>, uint64_t>
mortonOrderIndicies1d(const DataBaseSpace::DataBase<Dim<1> >& dataBase,
                      const FieldList<Dim<1>, int>& mask) {
  return mortonOrderIndicies(dataBase, mask);
}

inline
FieldSpace::FieldList<Dim<2>, uint64_t>
mortonOrderIndicies2d(const DataBaseSpace::DataBase<Dim<2> >& dataBase,
                      const FieldList<Dim<2>, int>& mask) {
  return mortonOrderIndicies(dataBase, mask);
}

inline
FieldSpace::FieldList<Dim<3>, uint64_t>
mortonOrderIndicies3d(const DataBaseSpace::DataBase<Dim<3> >& dataBase,
                      const FieldList<Dim<3>, int>& mask) {
  return mortonOrderIndicies(dataBase, mask);
}

//------------------------------------------------------------------------------
// And these...
//------------------------------------------------------------------------------
inline
FieldSpace::FieldList<Dim<1>, uint64_t>
peanoHilbertOrderIndicies1d(const DataBaseSpace::DataBase<Dim<1> >& dataBase) {
  return peanoHilbertOrderIndicies(dataBase);
}

inline
FieldSpace::FieldList<Dim<2>, uint64_t>
peanoHilbertOrderIndicies2d(const DataBaseSpace::DataBase<Dim<2> >& dataBase) {
  return peanoHilbertOrderIndicies(dataBase);
}

inline
FieldSpace::FieldList<Dim<3>, uint64_t>
peanoHilbertOrderIndicies3d(const DataBaseSpace::DataBase<Dim<3> >& dataBase) {
  return peanoHilbertOrderIndicies(dataBase);
}

//------------------------------------------------------------------------------
// And these...
//------------------------------------------------------------------------------
namespace NodeSpace {

template<typename Dimension>
inline
int
numGlobalNodesAll(const DataBaseSpace::DataBase<Dimension>& dataBase) {
  return numGlobalNodes(dataBase);
}

template<typename Dimension>
inline
FieldSpace::FieldList<Dimension, int>
globalNodeIDsAll(const DataBaseSpace::DataBase<Dimension>& dataBase) {
  return globalNodeIDs(dataBase);
}

}

//------------------------------------------------------------------------------
// testBoxIntersection
//------------------------------------------------------------------------------
inline
bool
testBoxIntersection1d(const Dim<1>::Vector& xmin1, const Dim<1>::Vector& xmax1,
                      const Dim<1>::Vector& xmin2, const Dim<1>::Vector& xmax2) {
  return testBoxIntersection(xmin1, xmax1, xmin2, xmax2);
}

inline
bool
testBoxIntersection2d(const Dim<2>::Vector& xmin1, const Dim<2>::Vector& xmax1,
                      const Dim<2>::Vector& xmin2, const Dim<2>::Vector& xmax2) {
  return testBoxIntersection(xmin1, xmax1, xmin2, xmax2);
}

inline
bool
testBoxIntersection3d(const Dim<3>::Vector& xmin1, const Dim<3>::Vector& xmax1,
                      const Dim<3>::Vector& xmin2, const Dim<3>::Vector& xmax2) {
  return testBoxIntersection(xmin1, xmax1, xmin2, xmax2);
}

//------------------------------------------------------------------------------
// testPointInBox
//------------------------------------------------------------------------------
inline
bool
testPointInBox1d(const Dim<1>::Vector& point, 
                 const Dim<1>::Vector& xmin, const Dim<1>::Vector& xmax) {
  return testPointInBox(point, xmin, xmax);
}

inline
bool
testPointInBox2d(const Dim<2>::Vector& point, 
                 const Dim<2>::Vector& xmin, const Dim<2>::Vector& xmax) {
  return testPointInBox(point, xmin, xmax);
}

inline
bool
testPointInBox3d(const Dim<3>::Vector& point, 
                 const Dim<3>::Vector& xmin, const Dim<3>::Vector& xmax) {
  return testPointInBox(point, xmin, xmax);
}

//------------------------------------------------------------------------------
// Disambiguate integrateThroughMeshAlongSegment
//------------------------------------------------------------------------------
inline
Dim<1>::Scalar
integrateThroughMeshAlongSegment1d(const std::vector<std::vector<Dim<1>::Scalar> >& values,
                                   const Dim<1>::Vector& xmin,
                                   const Dim<1>::Vector& xmax,
                                   const std::vector<unsigned>& ncells,
                                   const Dim<1>::Vector& s0,
                                   const Dim<1>::Vector& s1) {
  return integrateThroughMeshAlongSegment<Dim<1>, Dim<1>::Scalar>(values,
                                                                  xmin,
                                                                  xmax,
                                                                  ncells,
                                                                  s0,
                                                                  s1);
}

inline
Dim<2>::Scalar
integrateThroughMeshAlongSegment2d(const std::vector<std::vector<Dim<2>::Scalar> >& values,
                                   const Dim<2>::Vector& xmin,
                                   const Dim<2>::Vector& xmax,
                                   const std::vector<unsigned>& ncells,
                                   const Dim<2>::Vector& s0,
                                   const Dim<2>::Vector& s1) {
  return integrateThroughMeshAlongSegment<Dim<2>, Dim<2>::Scalar>(values,
                                                                  xmin,
                                                                  xmax,
                                                                  ncells,
                                                                  s0,
                                                                  s1);
}

inline
Dim<3>::Scalar
integrateThroughMeshAlongSegment3d(const std::vector<std::vector<Dim<3>::Scalar> >& values,
                                   const Dim<3>::Vector& xmin,
                                   const Dim<3>::Vector& xmax,
                                   const std::vector<unsigned>& ncells,
                                   const Dim<3>::Vector& s0,
                                   const Dim<3>::Vector& s1) {
  return integrateThroughMeshAlongSegment<Dim<3>, Dim<3>::Scalar>(values,
                                                                  xmin,
                                                                  xmax,
                                                                  ncells,
                                                                  s0,
                                                                  s1);
}

}

#endif
