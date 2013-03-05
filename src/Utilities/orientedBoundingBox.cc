//---------------------------------Spheral++----------------------------------//
// orientedBoundingBox
//
// A collection of methods to compute oriented bounding boxes for collections
// of points, and test for containment in such boxes.
//
// Created by JMO, Sun Jan 24 16:13:16 PST 2010
//----------------------------------------------------------------------------//
#include <vector>
#include <algorithm>

#include "orientedBoundingBox.hh"
#include "spheralWildMagicConverters.hh"
#include "DataBase/DataBase.hh"
#include "Wm5ContBox2.h"
#include "Wm5ContBox3.h"

namespace Spheral {

using namespace std;
using DataBaseSpace::DataBase;
using FieldSpace::Field;
using FieldSpace::FieldList;

//------------------------------------------------------------------------------
// Compute the minimum volume box containing all the nodes in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
void
orientedBoundingBox(const DataBase<Dimension>& dataBase,
                    typename Dimension::Box& nodeBox,
                    typename Dimension::Box& sampleBox) {
  typedef typename Dimension::WMVector WMVector;
  const vector<WMVector> points = wildMagicPositions(dataBase);
  const vector<WMVector> samplePoints = wildMagicSamplingPositions(dataBase);
  nodeBox = orientedBoundingBox(points);
  sampleBox = orientedBoundingBox(samplePoints);
}

//------------------------------------------------------------------------------
// 1-D
//------------------------------------------------------------------------------
Dim<1>::Box
orientedBoundingBox(const vector<Dim<1>::WMVector>& points) {
  typedef Dim<1>::Vector Vector;

  if (points.size() > 0) {
    double 
      xminPoints = DBL_MAX,
      xmaxPoints = -DBL_MAX;
    for (int i = 0; i != points.size(); ++i) {
      xminPoints = min(xminPoints, points[i].x());
      xmaxPoints = max(xmaxPoints, points[i].x());
    }
    return Dim<1>::Box(0.5*(xminPoints + xmaxPoints), 0.5*(xmaxPoints - xminPoints));
  } else {
    return Dim<1>::Box(Vector(), 0.0);
  }
}

//------------------------------------------------------------------------------
// 2-D
//------------------------------------------------------------------------------
Dim<2>::Box
orientedBoundingBox(const vector<Dim<2>::WMVector>& points) {
  typedef Dim<2> Dimension;
  typedef Dimension::WMVector WMVector;

  if (points.size() > 0) {
    return Wm5::ContOrientedBox(points.size(),
                                &(points.front()));
  } else {
    return Dim<2>::Box();
  }
}

//------------------------------------------------------------------------------
// 3-D
//------------------------------------------------------------------------------
Dim<3>::Box
orientedBoundingBox(const vector<Dim<3>::WMVector>& points) {
  typedef Dim<3> Dimension;
  typedef Dimension::Vector Vector;

  if (points.size() > 0) {
    return Wm5::ContOrientedBox(points.size(), 
                                &(points.front()));
  } else {
    return Dim<3>::Box();
  }
}

}

