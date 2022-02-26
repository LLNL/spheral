//---------------------------------Spheral++----------------------------------//
// PairFields -- basic DEM package for Spheral++.
//----------------------------------------------------------------------------//
#ifndef __Spheral_PairFields_hh__
#define __Spheral_PairFields_hh__

#include <string>
#include "DEM/DEMDimension.hh"

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class DataBase;
template<typename Dimension> class ContactModelBase;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;
class FileIO;
class RedistributionNotificationHandle;

template<typename Dimension>
class PairFields {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename DEMDimension<Dimension>::AngularVector RotationType;
  typedef std::shared_ptr<RedistributionNotificationHandle> RedistributionRegistrationType;
  
  // Constructors.
  PairFields(DataBase<Dimension>& dataBase);

  // Destructor.
  ~PairFields();

  // Tasks we do once on problem startup.
  void initialize();

  // we'll need to do some organizing 
  void initializeBeforeRedistribution();
  void finalizeAfterRedistribution();

  // access for fieldLists
  const FieldList<Dimension, int>& uniqueIndices() const;
  const FieldList<Dimension, std::vector<int>>& neighborIndices() const;
  const FieldList<Dimension, std::vector<Vector>>& shearDisplacement() const;
  const FieldList<Dimension, std::vector<Vector>>& DDtShearDisplacement() const;
  const FieldList<Dimension, std::vector<Scalar>>& equilibriumOverlap() const;

  virtual std::string label() const  { return "PairFields" ; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);

protected:
  //---------------------------  Protected Interface ---------------------------//
  FieldList<Dimension,int> mUniqueIndices;
  FieldList<Dimension,std::vector<int>> mNeighborIndices;
  FieldList<Dimension,std::vector<Vector>> mShearDisplacement;
  FieldList<Dimension,std::vector<Vector>> mDDtShearDisplacement;
  FieldList<Dimension,std::vector<Scalar>> mEquilibriumOverlap;

  RestartRegistrationType mRestart;
  RedistributionRegistrationType mRedistribute;
private:
  //--------------------------- Private Interface ---------------------------//
  // No default constructor, copying, or assignment.
  PairFields();
  PairFields(const PairFields&);
  PairFields& operator=(const PairFields&);
};

}

#include "PairFieldsInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class PairFields;
}

#endif
