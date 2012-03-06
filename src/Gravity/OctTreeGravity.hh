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
                 const double softeningLength,
                 const double opening,
                 const double ftimestep);

  //! Destructor.
  virtual ~OctTreeGravity();

  //! We augment the generic body force state.
  virtual void registerState(DataBaseSpace::DataBase<Dimension>& dataBase,
                             State<Dimension>& state);

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
                       
  //! This package opts out of building connectivity.
  virtual bool requireConnectivity() const { return false; }

  //! Return the total energy contribution due to the gravitational potential.
  virtual Scalar extraEnergy() const;

  //! Return the gravitational potential created by the particle distribution.
  const FieldSpace::FieldList<Dimension, Scalar>& potential() const;

  //! Return a dump of the tree structure as a string.
  std::string dumpTree() const;

  //! Return a string describing the overall statistics of the tree.
  std::string dumpTreeStatistics() const;

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

  //! The last computed maximum tree cell density.
  double maxCellDensity() const;

private:
  // Data types we use to build the internal tree structure.
  typedef uint32_t LevelKey;
  typedef uint64_t CellKey;

  static unsigned num1dbits;                   // The number of bits we quantize 1D coordinates to.  We have to fit three of these in 64 bits.
  static CellKey max1dKey;                     // The maximum number of cells this corresponds to in a direction.
  static CellKey xkeymask, ykeymask, zkeymask; // Bit masks we can use to extract the coordinate specific indices from a cell key.

  //----------------------------------------------------------------------------
  // Cell holds the properties of cells in the tree.
  //----------------------------------------------------------------------------
  struct Cell {
    double M, Mglobal;               // total mass (and global sum)
    Vector xcm;                      // center of mass
    double rcm2cc2;                  // square of the distance between center of mass and geometric center
    std::vector<CellKey> daughters;  // Keys of any daughter cells on level+1
    std::vector<Cell*> daughterPtrs; // Pointers to the daughter cells.
    std::vector<double> masses;      // Masses of the nodes that terminate in this cell.
    std::vector<Vector> positions;   // Positions of the nodes that terminate in this cell.

    // Convenience constructors for OctTreeGravity::addNodeToTree.
    Cell(): M(0.0), Mglobal(0.0), xcm(), rcm2cc2(0.0), daughters(), daughterPtrs(), masses(), positions() {}
    Cell(const double mi, const Vector& xi):
      M(mi), Mglobal(mi), xcm(xi), rcm2cc2(0.0), daughters(), daughterPtrs(), masses(1, mi), positions(1, xi) {}
    Cell(const double mi, const Vector& xi, const CellKey& daughter):
      M(mi), Mglobal(mi), xcm(xi), rcm2cc2(0.0), daughters(1, daughter), daughterPtrs(), masses(), positions() {}
  };

  // Define the types we use to build the tree.
  typedef boost::unordered_map<CellKey, Cell> TreeLevel;
  typedef std::vector<TreeLevel> Tree;

  // Private data.
  double mG, mSofteningLength2, mOpening2, mftimestep, mBoxLength, mMaxCellDensity;
  Vector mXmin, mXmax;
  Tree mTree;

  // The potential fields filled in during evaluateDerivates.
  mutable FieldSpace::FieldList<Dimension, Scalar> mPotential;
  mutable Scalar mExtraEnergy;
  
  // Default constructor -- disabled.
  OctTreeGravity();

  // Copy constructor -- disabled.
  OctTreeGravity(const OctTreeGravity&);

  // Assignment operator -- disabled.
  OctTreeGravity& operator=(const OctTreeGravity&);

  // Build a cell key based on a level and position.
  void buildCellKey(const LevelKey ilevel,
                    const Vector& xi,
                    CellKey& key,
                    CellKey& ix,
                    CellKey& iy,
                    CellKey& iz) const;

  // Extract the individual coordinate indices from a cell key.
  void extractCellIndices(const CellKey& key,
                          CellKey& ix,
                          CellKey& iy,
                          CellKey& iz) const;

  // Add a cell key to the daughters of a cell.
  void addDaughter(Cell& cell, const CellKey daughterKey) const;

  // Add a node to the internal tree.
  void addNodeToTree(const double mi,
                     const Vector& xi);

  // Construct all the daughterPtrs in a tree.
  void constructDaughterPtrs(Tree& tree) const;

  // Walk a tree and apply it's forces to a set of points.
  void applyTreeForces(const Tree& tree,
                       const FieldSpace::FieldList<Dimension, Scalar>& mass,
                       const FieldSpace::FieldList<Dimension, Vector>& position,
                       FieldSpace::FieldList<Dimension, Vector>& DxDt,
                       FieldSpace::FieldList<Dimension, Vector>& DvDt,
                       FieldSpace::FieldList<Dimension, Scalar>& potential) const;

  // Methods to help serializing/deserializing Trees to buffers of char.
  void serialize(const Tree& tree, std::vector<char>& buffer) const;
  void serialize(const Cell& cell, std::vector<char>& buffer) const;

  // Unpack a tree from a buffer.
  void deserialize(Tree& tree, std::vector<char>::const_iterator& bufItr, const std::vector<char>::const_iterator& endItr) const;
  void deserialize(Cell& cell, std::vector<char>::const_iterator& bufItr, const std::vector<char>::const_iterator& endItr) const;
};

}
}

#include "OctTreeGravityInline.hh"

#endif
