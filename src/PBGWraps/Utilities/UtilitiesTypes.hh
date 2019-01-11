#ifndef __PBGWRAP_UTILITIESTYPES__
#define __PBGWRAP_UTILITIESTYPES__

#include <utility>

#include "Geometry/Dimension.hh"
#include "boost/math/special_functions/legendre.hpp"
#include "Utilities/Functors.hh"
#include "Utilities/erff.hh"
#include "Utilities/newtonRaphson.hh"
#include "Utilities/simpsonsIntegration.hh"
#include "Utilities/globalNodeIDs.hh"
#include "Utilities/rotationMatrix.hh"
#include "Utilities/iterateIdealH.hh"
#include "Utilities/nodeOrdering.hh"
#include "Utilities/mortonOrderIndices.hh"
#include "Utilities/peanoHilbertOrderIndices.hh"
#include "Utilities/boundingBox.hh"
#include "Utilities/globalBoundingVolumes.hh"
#include "Utilities/testBoxIntersection.hh"
#include "Utilities/lineSegmentIntersections.hh"
#include "Utilities/pointDistances.hh"
#include "Utilities/segmentIntersectEdges.hh"
#include "Utilities/integrateThroughMeshAlongSegment.hh"
#include "Utilities/numberDensity.hh"
//#include "Utilities/writeRectilinearMesh.hh"
#include "Utilities/planarReflectingOperator.hh"
#include "Utilities/KeyTraits.hh"
#include "Utilities/pointInPolygon.hh"
#include "Utilities/pointOnPolygon.hh"
#include "Utilities/pointInPolyhedron.hh"
#include "Utilities/pointOnPolyhedron.hh"
#include "Utilities/refinePolyhedron.hh"
#include "Utilities/overlayRemapFields.hh"
#include "Utilities/computeShepardsInterpolation.hh"
#include "Utilities/Timer.hh"

#ifndef NOR3D
#include "Utilities/r3d_utils.hh"
#else
//------------------------------------------------------------------------------
// Stub these methods out when we're not allowing R3D.
//------------------------------------------------------------------------------
namespace Spheral {
  inline Dim<2>::FacetedVolume clipFacetedVolume(const Dim<2>::FacetedVolume& poly,
                                                 const std::vector<GeomPlane<Dim<2> > >& planes) {
    VERIFY2(false, "ERROR: clipFacetedVolume unavailable without R3D.");
  }

  inline Dim<3>::FacetedVolume clipFacetedVolume(const Dim<3>::FacetedVolume& poly,
                                                 const std::vector<GeomPlane<Dim<3> > >& planes) {
    VERIFY2(false, "ERROR: clipFacetedVolume unavailable without R3D.");
  }
}
#endif

typedef std::pair<double, double> pair_double_double;

namespace Spheral {

// Expose some stuff from boost::math
using boost::math::legendre_p;

//------------------------------------------------------------------------------
// newtonRaphsonFindRoot
//------------------------------------------------------------------------------
inline
double newtonRaphsonFindRoot(const PythonBoundFunctors::SpheralFunctor<double, std::pair<double, double> >& functor,
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
// simpsonsIntegrationDouble
//------------------------------------------------------------------------------
inline
double simpsonsIntegrationDouble(const PythonBoundFunctors::SpheralFunctor<double, double>& functor,
                                 const double x0,
                                 const double x1,
                                 const unsigned numBins) {
  return Spheral::simpsonsIntegration<PythonBoundFunctors::SpheralFunctor<double, double>, double, double>(functor,
                                                                                                           x0,
                                                                                                           x1,
                                                                                                           numBins);
}

//------------------------------------------------------------------------------
// Wrapper around our packElement method more convenient for exposing to Python.
//------------------------------------------------------------------------------
template<typename Element>
inline
std::string
convertElementToString(const Element& x) {
  vector<char> buffer;
  packElement(x, buffer);
  return std::string(buffer.begin(), buffer.end());
}

//------------------------------------------------------------------------------
// Wrapper around our unpackElement method more convenient for exposing to
//  Python.
//------------------------------------------------------------------------------
template<typename Element>
inline
Element
convertStringToElement(const std::string& x) {
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
FieldList<Dim<1>, uint64_t>
mortonOrderIndices1d(const FieldList<Dim<1>, Dim<1>::Vector>& positions) {
  return mortonOrderIndices(positions);
}

inline
FieldList<Dim<2>, uint64_t>
mortonOrderIndices2d(const FieldList<Dim<2>, Dim<2>::Vector>& positions) {
  return mortonOrderIndices(positions);
}

inline
FieldList<Dim<3>, uint64_t>
mortonOrderIndices3d(const FieldList<Dim<3>, Dim<3>::Vector>& positions) {
  return mortonOrderIndices(positions);
}

//------------------------------------------------------------------------------
// And these...
//------------------------------------------------------------------------------
inline
FieldList<Dim<1>, uint64_t>
mortonOrderIndices1d(const DataBase<Dim<1> >& dataBase) {
  return mortonOrderIndices(dataBase);
}

inline
FieldList<Dim<2>, uint64_t>
mortonOrderIndices2d(const DataBase<Dim<2> >& dataBase) {
  return mortonOrderIndices(dataBase);
}

inline
FieldList<Dim<3>, uint64_t>
mortonOrderIndices3d(const DataBase<Dim<3> >& dataBase) {
  return mortonOrderIndices(dataBase);
}

//------------------------------------------------------------------------------
// And these...
//------------------------------------------------------------------------------
inline
FieldList<Dim<1>, uint64_t>
mortonOrderIndices1d(const DataBase<Dim<1> >& dataBase,
                      const FieldList<Dim<1>, int>& mask) {
  return mortonOrderIndices(dataBase, mask);
}

inline
FieldList<Dim<2>, uint64_t>
mortonOrderIndices2d(const DataBase<Dim<2> >& dataBase,
                      const FieldList<Dim<2>, int>& mask) {
  return mortonOrderIndices(dataBase, mask);
}

inline
FieldList<Dim<3>, uint64_t>
mortonOrderIndices3d(const DataBase<Dim<3> >& dataBase,
                      const FieldList<Dim<3>, int>& mask) {
  return mortonOrderIndices(dataBase, mask);
}

//------------------------------------------------------------------------------
// And these...
//------------------------------------------------------------------------------
inline
FieldList<Dim<1>, uint64_t>
peanoHilbertOrderIndices1d(const FieldList<Dim<1>, Dim<1>::Vector>& positions) {
  return peanoHilbertOrderIndices(positions);
}

inline
FieldList<Dim<2>, uint64_t>
peanoHilbertOrderIndices2d(const FieldList<Dim<2>, Dim<2>::Vector>& positions) {
  return peanoHilbertOrderIndices(positions);
}

inline
FieldList<Dim<3>, uint64_t>
peanoHilbertOrderIndices3d(const FieldList<Dim<3>, Dim<3>::Vector>& positions) {
  return peanoHilbertOrderIndices(positions);
}

//------------------------------------------------------------------------------
// And these...
//------------------------------------------------------------------------------
inline
FieldList<Dim<1>, uint64_t>
peanoHilbertOrderIndices1d(const DataBase<Dim<1> >& dataBase) {
  return peanoHilbertOrderIndices(dataBase);
}

inline
FieldList<Dim<2>, uint64_t>
peanoHilbertOrderIndices2d(const DataBase<Dim<2> >& dataBase) {
  return peanoHilbertOrderIndices(dataBase);
}

inline
FieldList<Dim<3>, uint64_t>
peanoHilbertOrderIndices3d(const DataBase<Dim<3> >& dataBase) {
  return peanoHilbertOrderIndices(dataBase);
}

//------------------------------------------------------------------------------
// And these...
//------------------------------------------------------------------------------
inline
FieldList<Dim<1>, int>
nodeOrdering1d(const FieldList<Dim<1>, uint64_t>& criteria) {
  return nodeOrdering(criteria);
}

inline
FieldList<Dim<2>, int>
nodeOrdering2d(const FieldList<Dim<2>, uint64_t>& criteria) {
  return nodeOrdering(criteria);
}

inline
FieldList<Dim<3>, int>
nodeOrdering3d(const FieldList<Dim<3>, uint64_t>& criteria) {
  return nodeOrdering(criteria);
}

//------------------------------------------------------------------------------
// And these...
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
numGlobalNodesAll(const DataBase<Dimension>& dataBase) {
  return numGlobalNodes(dataBase);
}

template<typename Dimension>
inline
FieldList<Dimension, int>
globalNodeIDsAll(const DataBase<Dimension>& dataBase) {
  return globalNodeIDs(dataBase);
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
