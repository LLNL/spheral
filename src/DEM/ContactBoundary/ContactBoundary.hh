//---------------------------------Spheral++----------------------------------//
// DEMBase -- basic DEM package for Spheral++
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//
#ifndef __Spheral_ContactBoundary_hh__
#define __Spheral_ContactBoundary_hh__

#include <string>

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class DataBase;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;
class FileIO;
class RedistributionNotificationHandle;
struct ContactIndex;

template<typename Dimension>
class ContactBoundary {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename DEMDimension<Dimension>::AngularVector RotationType;

  typedef typename Physics<Dimension>::TimeStepType TimeStepType;
  typedef typename Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;
  typedef std::shared_ptr<RedistributionNotificationHandle> RedistributionRegistrationType;
  
  // Constructors.
  ContactBoundary(const DataBase<Dimension>& dataBase,
                  const Scalar stepsPerCollision,
          const Scalar neighborSearchBuffer,
          const Vector& xmin,
          const Vector& xmax);

  // Destructor.
  virtual ~ContactBoundary();


private:
  //--------------------------- Private Interface ---------------------------//
  // No default constructor, copying, or assignment.
  ContactBoundary();
  ContactBoundary(const ContactBoundary&);
  ContactBoundary& operator=(const ContactBoundary&);
};

}

#include "ContactBoundaryInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class ContactBoundary;
}

#endif
