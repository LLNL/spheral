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

#include "Utilities/NodeCoupling.hh"
#include "Utilities/DBC.hh"
#include "Field/FieldList.hh"

namespace Spheral {

template<typename Dimension>
class SurfaceNodeCoupling: public NodeCoupling {
public:
  // Constructor, destructor.
  SurfaceNodeCoupling(const FieldList<Dimension, int>& surfacePoint):
    NodeCoupling(),
    mSurfacePoint(surfacePoint) {}

  virtual ~SurfaceNodeCoupling() {}

  // The coupling operator.
  virtual double operator()(const unsigned nodeListi, const unsigned i,
                            const unsigned nodeListj, const unsigned j) const {
    return (nodeListi == nodeListj          ? 1.0 :
            mSurfacePoint(nodeListj, j) > 0 ? 1.0 :
                                              0.0);
  }

private:
  const FieldList<Dimension, int>& mSurfacePoint;
};

}

#endif
