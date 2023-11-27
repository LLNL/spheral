//---------------------------------Spheral++----------------------------------//
// SumVoronoiMassDensityPolicy -- An implementation of UpdatePolicyBase 
// specialized for the updating the mass density according to the specific 
// volume from the Voronoi tesselation.  This version uses a weighted sum of the
// local mass and volume to do the deed.
//
// Created by JMO, Mon Aug  1 10:48:03 PDT 2011
//----------------------------------------------------------------------------//
#ifndef __Spheral_SumVoronoiMassDensityPolicy_hh__
#define __Spheral_SumVoronoiMassDensityPolicy_hh__

#include <string>

#include "DataBase/ReplaceFieldList.hh"
#include "Kernel/TableKernel.hh"
#include "Physics/Physics.hh"

namespace Spheral {

template<typename Dimension>
class SumVoronoiMassDensityPolicy: public ReplaceFieldList<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename ReplaceFieldList<Dimension, Scalar>::KeyType KeyType;

  // Constructors, destructor.
  SumVoronoiMassDensityPolicy(const TableKernel<Dimension>& W,
                              const Physics<Dimension>& package,
                              const double rhoMin,
                              const double rhoMax);
  virtual ~SumVoronoiMassDensityPolicy();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt);

  virtual void updateAsIncrement(const KeyType& key,
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs,
                                 const double multiplier,
                                 const double t,
                                 const double dt);

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const;

private:
  //--------------------------- Private Interface ---------------------------//
  bool mEnforceBoundaries;
  const TableKernel<Dimension>& mW;
  const Physics<Dimension>& mPackage;
  double mRhoMin, mRhoMax;
  SumVoronoiMassDensityPolicy(const SumVoronoiMassDensityPolicy& rhs);
  SumVoronoiMassDensityPolicy& operator=(const SumVoronoiMassDensityPolicy& rhs);
};

}

#endif
