//---------------------------------Spheral++----------------------------------//
// DEMBoundaryPolicy -- general policy to allow solid boundaries to move in
//                      DEM. 
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#ifndef __Spheral_DEMBoundaryPolicy_hh__
#define __Spheral_DEMBoundaryPolicy_hh__

#include <float.h>
#include <vector>
#include "DataBase/UpdatePolicyBase.hh"

namespace Spheral {

template<typename Dimension> class SolidBoundary;

template<typename Dimension>
class DEMBoundaryPolicy: public UpdatePolicyBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename UpdatePolicyBase<Dimension>::KeyType KeyType;
  typedef typename std::vector<SolidBoundary<Dimension>*>::iterator SolidBoundaryIterator;
  typedef typename std::vector<SolidBoundary<Dimension>*>::const_iterator ConstSolidBoundaryIterator;

  // Constructors, destructor.
  DEMBoundaryPolicy();

  virtual ~DEMBoundaryPolicy();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt);



  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const;

  // Solid Bounderies 
  void appendSolidBoundary(SolidBoundary<Dimension>& boundary);
  void prependSolidBoundary(SolidBoundary<Dimension>& boundary);
  void clearSolidBoundaries();
  bool haveSolidBoundary(const SolidBoundary<Dimension>& boundary) const;
  
  const std::vector<SolidBoundary<Dimension>*>& solidBoundaryConditions() const;

private:
  //--------------------------- Private Interface ---------------------------//
  std::vector<SolidBoundary<Dimension>*> mSolidBoundaries;

  DEMBoundaryPolicy(const DEMBoundaryPolicy& rhs);
  DEMBoundaryPolicy& operator=(const DEMBoundaryPolicy& rhs);
};

}

#ifndef __GCCXML__
#include "DEMBoundaryPolicyInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class DEMBoundaryPolicy;
}

#endif
