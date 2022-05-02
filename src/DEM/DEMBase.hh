//---------------------------------Spheral++----------------------------------//
// DEMBase -- basic DEM package for Spheral++.
//----------------------------------------------------------------------------//
#ifndef __Spheral_DEMBase_hh__
#define __Spheral_DEMBase_hh__

#include <string>
#include "Physics/Physics.hh"
#include "DEM/DEMDimension.hh"

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
class DEMBase: public Physics<Dimension> {

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
  DEMBase(const DataBase<Dimension>& dataBase,
          const Scalar stepsPerCollision,
          const Vector& xmin,
          const Vector& xmax);

  // Destructor.
  virtual ~DEMBase();

  // Tasks we do once on problem startup.
  virtual
  void initializeProblemStartup(DataBase<Dimension>& dataBase) override;

  // Register the state Hydro expects to use and evolve.
  virtual 
  void registerState(DataBase<Dimension>& dataBase,
                     State<Dimension>& state) override;

  // Register the derivatives/change fields for updating state.
  virtual
  void registerDerivatives(DataBase<Dimension>& dataBase,
                           StateDerivatives<Dimension>& derivs) override;

  // This method is called once at the beginning of a timestep, after all state registration.
  virtual void preStepInitialize(const DataBase<Dimension>& dataBase, 
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs) override;

  // Initialize the Hydro before we start a derivative evaluation.
  virtual
  void initialize(const Scalar time,
                  const Scalar dt,
                  const DataBase<Dimension>& dataBase,
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) override;
                       

  // Finalize the derivatives.
  virtual
  void finalizeDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivs) const override;

  // Apply boundary conditions to the physics specific fields.
  virtual
  void applyGhostBoundaries(State<Dimension>& state,
                            StateDerivatives<Dimension>& derivs) override;

  // Enforce boundary conditions for the physics specific fields.
  virtual
  void enforceBoundaries(State<Dimension>& state,
                         StateDerivatives<Dimension>& derivs) override;
  
  // methods needed to square our pair fieldlists before redistribution
  void initializeBeforeRedistribution();
  void finalizeAfterRedistribution();

  // templated methods
  template<typename Value1, typename Value2>
  void addContactsToPairFieldList(Value1& pairFieldList, const Value2& newValue) const;

  template<typename Value>
  void kullInactiveContactsFromPairFieldList(Value& pairFieldList) const;
  
  //#############################################################################
  // special methods for the pair fields if you add pairFieldLists 
  // make sure to follow the pattern in these methods.
  virtual
  void resizeDerivativePairFieldLists(StateDerivatives<Dimension>& derivs) const;

  virtual
  void resizeStatePairFieldLists(State<Dimension>& state) const;

  virtual 
  void kullInactiveContactsFromStatePairFieldLists(State<Dimension>& state) const;
  //#############################################################################

  void addContacts(const DataBase<Dimension>& dataBase,
                         State<Dimension>& state,
                         StateDerivatives<Dimension>& derivs);

  void kullInactiveContacts(const DataBase<Dimension>& dataBase,
                                  State<Dimension>& state,
                                  StateDerivatives<Dimension>& derivs);

  // methods to find storage location for pairwise fields
  std::vector<int> storageNodeSelection(int nodeListi,
                                        int i,
                                        int nodeListj,
                                        int j) const;

  std::vector<int> findContactIndex(int nodeListi,
                                    int i,
                                    int nodeListj,
                                    int j) const;

  // Optionally we can provide a bounding box for use generating the mesh
  // for the Voronoi mass density update.
  const Vector& xmin() const;
  const Vector& xmax() const;
  void xmin(const Vector& x);
  void xmax(const Vector& x);

  Scalar stepsPerCollision() const;
  void   stepsPerCollision(Scalar x);

  // access for node fieldLists
  const FieldList<Dimension, int>&    timeStepMask() const;
  const FieldList<Dimension, Vector>& DxDt() const;
  const FieldList<Dimension, Vector>& DvDt() const;
  const FieldList<Dimension, RotationType>& omega() const;
  const FieldList<Dimension, RotationType>& DomegaDt() const;

  // access for pair fieldLists
  const FieldList<Dimension, int>& uniqueIndices() const;
  const FieldList<Dimension, std::vector<int>>& isActiveContact() const;
  const FieldList<Dimension, std::vector<int>>& neighborIndices() const;
  const FieldList<Dimension, std::vector<Vector>>& shearDisplacement() const;
  const FieldList<Dimension, std::vector<Vector>>& DDtShearDisplacement() const;
  const FieldList<Dimension, std::vector<Scalar>>& equilibriumOverlap() const;
  
  const std::vector<ContactIndex>& contactStorageIndices() const;

  // access a single element in our state pair fields
  const Vector& shearDisplacement(const int nodeListi,
                                  const int nodei,
                                  const int contacti) const;
  const Scalar  equilibriumOverlap(const int nodeListi,
                                   const int nodei,
                                   const int contacti) const;

  // inlined and specialized for different dimensions
  Scalar momentOfInertia(const Scalar massi,
                         const Scalar particleRadiusi) const;

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "DEMBase" ; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //****************************************************************************

protected:
  //---------------------------  Protected Interface ---------------------------//

  const DataBase<Dimension>& mDataBase;

  int mCyclesSinceLastKulling;
  int mKullFrequency;
  
  // number of steps per collision time-scale
  Scalar mStepsPerCollision;              

  // Optional bounding box for generating the mesh.
  Vector mxmin, mxmax;

  // fields attached to the nodes
  FieldList<Dimension, int>      mTimeStepMask;
  FieldList<Dimension, Vector>   mDxDt;
  FieldList<Dimension, Vector>   mDvDt;
  FieldList<Dimension, RotationType> mOmega;
  FieldList<Dimension, RotationType> mDomegaDt;

  // fields attached to the pair interactions
  FieldList<Dimension,int> mUniqueIndices;                         // each nodes global unqiue index
  FieldList<Dimension,std::vector<int>> mNeighborIndices;          // tracks unique indices of contacts-we upate these 
  FieldList<Dimension,std::vector<Scalar>> mEquilibriumOverlap;    // nonzero values for composite particles
  FieldList<Dimension,std::vector<Vector>> mShearDisplacement;     // displacement for friction spring
  FieldList<Dimension,std::vector<int>> mIsActiveContact;          // tracks if a interfaction is still active
  FieldList<Dimension,std::vector<Vector>> mDDtShearDisplacement;  // derivative to evolve frictional spring displacement

  std::vector<ContactIndex> mContactStorageIndices;

  // The restart registration.
  RestartRegistrationType mRestart;
  RedistributionRegistrationType mRedistribute;

private:
  //--------------------------- Private Interface ---------------------------//
  // No default constructor, copying, or assignment.
  DEMBase();
  DEMBase(const DEMBase&);
  DEMBase& operator=(const DEMBase&);
};

}

#include "DEMBaseInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class DEMBase;
}

#endif
