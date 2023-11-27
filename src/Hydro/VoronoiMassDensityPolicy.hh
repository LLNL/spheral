//---------------------------------Spheral++----------------------------------//
// VoronoiMassDensityPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the mass density according to the specific volume from the
// Voronoi tesselation.
//
// Created by JMO, Tue Mar  8 22:07:34 PST 2011
//----------------------------------------------------------------------------//
#ifndef __Spheral_VoronoiMassDensityPolicy_hh__
#define __Spheral_VoronoiMassDensityPolicy_hh__

#include "DataBase/ReplaceFieldList.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class VoronoiMassDensityPolicy: public ReplaceFieldList<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename ReplaceFieldList<Dimension, Scalar>::KeyType KeyType;

  // Constructors, destructor.
  VoronoiMassDensityPolicy(const double rhoMin,
                           const double rhoMax);
  virtual ~VoronoiMassDensityPolicy();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt);

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const;

private:
  //--------------------------- Private Interface ---------------------------//
  double mRhoMin, mRhoMax;
  VoronoiMassDensityPolicy(const VoronoiMassDensityPolicy& rhs);
  VoronoiMassDensityPolicy& operator=(const VoronoiMassDensityPolicy& rhs);
};

}

#endif
