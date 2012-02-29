//---------------------------------Spheral++----------------------------------//
// OctTreeGravity -- An implementation of the tree n-body gravity solver.
//
// Created by JMO, 2012-02-28
//----------------------------------------------------------------------------//
#ifndef __Spheral_OctTreeGravity__
#define __Spheral_OctTreeGravity__

#include <stdint.h>
#include "boost/unordered_map.hpp"

#include "Geometry/Dimension.hh"
#include "Physics/GenericBodyForce.hh"
#include "Field/FieldList.hh"

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;

namespace GravitySpace {

class OctTreeGravity: public PhysicsSpace::GenericBodyForce<Dim<3> > {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef Dim<3> Dimension;
  typedef Dimension::Scalar Scalar;
  typedef Dimension::Vector Vector;
  typedef Dimension::Tensor Tensor;
  typedef Dimension::SymTensor SymTensor;

  typedef PhysicsSpace::Physics<Dimension>::TimeStepType TimeStepType;

  //! Constructor.
  //! \param G -- the gravitational constant.
  //! \param opening -- the opening ratio for approximating forces.
  //! \param softeningLength -- The Plummer softening length for the model.
  //! \param ftimestep -- safety factor to apply to pmom_i/F_i in setting time steps.
  OctTreeGravity(const double G,
                 const double opening,
                 const double softeningLength,
                 const double ftimestep);

  //! Destructor.
  virtual ~OctTreeGravity();

  //! This is the derivative method that all BodyForce classes must provide.
  virtual 
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBaseSpace::DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivs) const;

  //! Vote on the timestep.  This uses a velocity-limiting rule.
  virtual TimeStepType dt(const DataBaseSpace::DataBase<Dimension>& dataBase, 
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const;

  //! Initializations on problem start up.
  virtual void initializeProblemStartup(DataBaseSpace::DataBase<Dimension>& db);

  //! Initialize before we start a derivative evaluation.
  virtual void initialize(const Scalar time,
                          const Scalar dt,
                          const DataBaseSpace::DataBase<Dimension>& dataBase,
                          State<Dimension>& state,
                          StateDerivatives<Dimension>& derivs);
                       
  //! Return the total energy contribution due to the gravitational potential.
  virtual Scalar extraEnergy() const;

  //! Return the gravitational potential created by the particle distribution.
  const FieldSpace::FieldList<Dimension, Scalar>& potential() const;

  //! The gravitational constant we're using.
  double G() const;

  //! The opening angle threshold when we shift to tree cell approximations.
  double opening() const;
  void opening(const double x);

  //! The current softening length.
  double softeningLength() const;
  void softeningLength(const double x);

  //! The current time step scaling factor.
  double ftimestep() const;
  void ftimestep(const double x);

  //! The lower left corner of the computational cube that was last used.
  Vector xmin() const;

  //! The upper right corner of the computational cube that was last used.
  Vector xmax() const;

private:
  // Data types we use to build the internal tree structure.
  typedef uint32_t LevelKey;
  typedef uint64_t CellKey;
  typedef std::pair<LevelKey, CellKey> TreeKey;
  typedef std::pair<size_t, size_t> NodeID;

  static unsigned num1dbits;             // The number of bits we quantize 1D coordinates to.  We have to fit three of these in 64 bits.
  static CellKey max1dKey, max1dKey1;    // The maximum number of cells this corresponds to in a direction.

  // Cell holds the properties of cells in the tree.
  struct Cell {
    Vector xcm;                     // center of mass
    double M;                       // total mass
    std::vector<CellKey> daughters; // Keys of any daughter cells on level+1
    std::vector<NodeID> members;    // Any Spheral nodes that terminate in this cell.

    // Convenience constructors for OctTreeGravity::addNodeToTree.
    Cell(): xcm(), M(0.0), daughters(), members() {}
    Cell(const Vector& xi, const double mi, const NodeID& nid):
      xcm(xi), M(mi), daughters(), members(std::vector<NodeID>(1, nid)) {}
    Cell(const Vector& xi, const double mi, const CellKey& daughter):
      xcm(xi), M(mi), daughters(std::vector<CellKey>(1, daughter)), members() {}
  };

  typedef boost::unordered_map<TreeKey, Cell> Tree;

  // Private data.
  double mG, mOpening, mSofteningLength, mftimestep, mBoxLength;
  Vector mXmin, mXmax;
  Tree mTree;

  // The time step control info filled in during evaluateDerivatives.
  mutable size_t mdt_fieldi, mdt_nodei;
  mutable double mdt_veli, mdt_acci;

  // The potential fields filled in during evaluateDerivates.
  mutable FieldSpace::FieldList<Dimension, Scalar> mPotential;
  mutable Scalar mExtraEnergy;
  
  // Default constructor -- disabled.
  OctTreeGravity();

  // Copy constructor -- disabled.
  OctTreeGravity(const OctTreeGravity&);

  // Assignment operator -- disabled.
  OctTreeGravity& operator=(const OctTreeGravity&);

  // Build the key for a cell.
  TreeKey buildTreeKey(const LevelKey ilevel, const Vector& xi) const;

  // Add a cell key to the daughters of a cell.
  void addDaughter(Cell& cell, const CellKey daughterKey) const;

  // Add a node to the internal tree.
  void addNodeToTree(const size_t nodeListi,
                     const size_t i,
                     const double mi,
                     const Vector& xi);
};

}
}

#include "OctTreeGravityInline.hh"

#endif
