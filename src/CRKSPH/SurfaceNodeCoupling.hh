//---------------------------------Spheral++----------------------------------//
// SurfaceNodeCoupling
//
// A functor base class encapsulating how we couple pairs of nodes.
// This form is specialized to prevent coupling between materials, except for 
// points flagged as surface.
//
// Created by JMO, Wed Aug 16 16:21:38 PDT 2017
//----------------------------------------------------------------------------//
#ifndef __Spheral_SurfaceNodeCoupling__
#define __Spheral_SurfaceNodeCoupling__

#include "SPH/NodeCoupling.hh"
#include "Utilities/DBC.hh"
#include "Field/FieldList.hh"

namespace Spheral {

template<typename Dimension>
class SurfaceNodeCoupling: public NodeCoupling {
public:
  // Constructor, destructor.
  SurfaceNodeCoupling(const FieldList<Dimension, int>& surfacePoints):
    NodeCoupling(),
    mSurfacePoints(surfacePoints) {}

  virtual ~SurfaceNodeCoupling() {}

  // The coupling operator.
  virtual double operator()(const unsigned nodeListi, const unsigned i,
                            const unsigned nodeListj, const unsigned j) const {
    if (nodeListi == nodeListj) return 1.0;
    if (mSurfacePoints(nodeListi, i) != 0 or
        mSurfacePoints(nodeListj, j) != 0) return 1.0;
    return 0.0;
  }

private:
  const FieldList<Dimension, int>& mSurfacePoints;
};

}

#endif
