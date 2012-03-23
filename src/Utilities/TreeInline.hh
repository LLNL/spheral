//---------------------------------Spheral++----------------------------------//
// Tree implementation.
//----------------------------------------------------------------------------//
#include "Geometry/Dimension.hh"

namespace Spheral {

// Put our private utility functions into an anonymous namespace.
namespace {

//------------------------------------------------------------------------------
// Dimension specialized methods for constructing keys.
// We use the convention of interleaving bits here, leading to Z or Morton
// ordering of the keys.
//------------------------------------------------------------------------------
template<typename Dimension, typename CellValue, typename LeafPolicy> struct CellKeyComputer;

// 1D
template<typename CellValue, typename LeafPolicy> 
struct
CellKeyComputer<Spheral::Dim<1>, CellValue, LeafPolicy> {
  static void tokey(const typename Tree<Spheral::Dim<1>, CellValue, LeafPolicy>::LevelKey ilevel,
                    const          Spheral::Dim<1>::Vector& xi,
                    const typename Tree<Spheral::Dim<1>, CellValue, LeafPolicy>::CellKey& key,
                    const typename Tree<Spheral::Dim<1>, CellValue, LeafPolicy>::CoordHash& ix,
                    const typename Tree<Spheral::Dim<1>, CellValue, LeafPolicy>::CoordHash& iy,
                    const typename Tree<Spheral::Dim<1>, CellValue, LeafPolicy>::CoordHash& iz) {
    typedef typename Tree<Spheral::Dim<1>, CellValue, LeafPolicy>::CoordHash CoordHash;
    REQUIRE(xi.x() >= mXmin.x() and xi.x() <= mXmax.x());
    REQUIRE(ilevel < Tree<Spheral::Dim<1>, CellValue, LeafPolicy>::num1dbits());
    const CoordHash ncell = (1U << ilevel);
    const CoordHash maxcell = ncell - 1U;
    ix = std::min(maxcell, CoordHash((xi.x() - mXmin.x())/mBoxLength * ncell));
    key = ix;
  }
  static void fromkey(const typename Tree<Spheral::Dim<1>, CellValue, LeafPolicy>::CellKey& key,
                      typename Tree<Spheral::Dim<1>, CellValue, LeafPolicy>::CoordHash& ix,
                      typename Tree<Spheral::Dim<1>, CellValue, LeafPolicy>::CoordHash& iy,
                      typename Tree<Spheral::Dim<1>, CellValue, LeafPolicy>::CoordHash& iz) {
    ix = key.to_ulong();
    iy = 0U;
    iz = 0U;
  }
};

// 2D
template<typename CellValue, typename LeafPolicy> 
struct
CellKeyComputer<Spheral::Dim<2>, CellValue, LeafPolicy> {
  static void tokey(const typename Tree<Spheral::Dim<2>, CellValue, LeafPolicy>::LevelKey ilevel,
                    const          Spheral::Dim<2>::Vector& xi,
                    const typename Tree<Spheral::Dim<2>, CellValue, LeafPolicy>::CellKey& key,
                    const typename Tree<Spheral::Dim<2>, CellValue, LeafPolicy>::CoordHash& ix,
                    const typename Tree<Spheral::Dim<2>, CellValue, LeafPolicy>::CoordHash& iy,
                    const typename Tree<Spheral::Dim<2>, CellValue, LeafPolicy>::CoordHash& iz) {
    typedef typename Tree<Spheral::Dim<2>, CellValue, LeafPolicy>::CoordHash CoordHash;
    REQUIRE(xi.x() >= mXmin.x() and xi.x() <= mXmax.x());
    REQUIRE(xi.y() >= mXmin.y() and xi.y() <= mXmax.y());
    REQUIRE(ilevel < Tree<Spheral::Dim<2>, CellValue, LeafPolicy>::num1dbits());
    const CoordHash ncell = (1U << ilevel);
    const CoordHash maxcell = ncell - 1U;
    ix = std::min(maxcell, CoordHash((xi.x() - mXmin.x())/mBoxLength * ncell));
    iy = std::min(maxcell, CoordHash((xi.y() - mXmin.y())/mBoxLength * ncell));
    key = 0U;
    for (unsigned i = 0; i != Tree<Dim<2>, CellValue, LeafPolicy>::num1dbits(); ++i) {
      key[2*i    ] = (ix >> i) & 1U;
      key[2*i + 1] = (iy >> i) & 1U;
    }
  }
  static void fromkey(const typename Tree<Spheral::Dim<2>, CellValue, LeafPolicy>::CellKey& key,
                      typename Tree<Spheral::Dim<2>, CellValue, LeafPolicy>::CoordHash& ix,
                      typename Tree<Spheral::Dim<2>, CellValue, LeafPolicy>::CoordHash& iy,
                      typename Tree<Spheral::Dim<2>, CellValue, LeafPolicy>::CoordHash& iz) {
    ix = 0U; iy = 0U; iz = 0U;
    unsigned jx, jy;
    for (unsigned i = 0; i != Tree<Dim<2>, CellValue, LeafPolicy>::num1dbits(); ++i) {
      jx = 2*i;
      jy = jx + 1;
      ix += ((key >> jx) & 1U) << i;
      iy += ((key >> jy) & 1U) << i;
    }
  }
};

// 3D
template<typename CellValue, typename LeafPolicy> 
struct
CellKeyComputer<Spheral::Dim<3>, CellValue, LeafPolicy>{
  static void tokey(const typename Tree<Spheral::Dim<3>, CellValue, LeafPolicy>::LevelKey ilevel,
                    const          Spheral::Dim<3>::Vector& xi,
                    const typename Tree<Spheral::Dim<3>, CellValue, LeafPolicy>::CellKey& key,
                    const typename Tree<Spheral::Dim<3>, CellValue, LeafPolicy>::CoordHash& ix,
                    const typename Tree<Spheral::Dim<3>, CellValue, LeafPolicy>::CoordHash& iy,
                    const typename Tree<Spheral::Dim<3>, CellValue, LeafPolicy>::CoordHash& iz) {
    typedef typename Tree<Spheral::Dim<3>, CellValue, LeafPolicy>::CoordHash CoordHash;
    REQUIRE(xi.x() >= mXmin.x() and xi.x() <= mXmax.x());
    REQUIRE(xi.y() >= mXmin.y() and xi.y() <= mXmax.y());
    REQUIRE(ilevel < Tree<Spheral::Dim<3>, CellValue, LeafPolicy>::num1dbits());
    const CoordHash ncell = (1U << ilevel);
    const CoordHash maxcell = ncell - 1U;
    ix = std::min(maxcell, CoordHash((xi.x() - mXmin.x())/mBoxLength * ncell));
    iy = std::min(maxcell, CoordHash((xi.y() - mXmin.y())/mBoxLength * ncell));
    iz = std::min(maxcell, CoordHash((xi.z() - mXmin.z())/mBoxLength * ncell));
    key = 0U;
    for (unsigned i = 0; i != Tree<Dim<3>, CellValue, LeafPolicy>::num1dbits(); ++i) {
      key[3*i    ] = (ix >> i) & 1U;
      key[3*i + 1] = (iy >> i) & 1U;
      key[3*i + 2] = (iz >> i) & 1U;
    }
  }
  static void fromkey(const typename Tree<Spheral::Dim<3>, CellValue, LeafPolicy>::CellKey& key,
                      typename Tree<Spheral::Dim<3>, CellValue, LeafPolicy>::CoordHash& ix,
                      typename Tree<Spheral::Dim<3>, CellValue, LeafPolicy>::CoordHash& iy,
                      typename Tree<Spheral::Dim<3>, CellValue, LeafPolicy>::CoordHash& iz) {
    ix = 0U; iy = 0U; iz = 0U;
    unsigned jx, jy, jz;
    for (unsigned i = 0; i != Tree<Dim<3>, CellValue, LeafPolicy>::num1dbits(); ++i) {
      jx = 3*i;
      jy = jx + 1;
      jz = jy + 1;
      ix += ((key >> jx) & 1U) << i;
      iy += ((key >> jy) & 1U) << i;
      iz += ((key >> jz) & 1U) << i;
    }
  }
};

} // end of anonymous namespace

//------------------------------------------------------------------------------
// Build a cell key from coordinate indices.
//------------------------------------------------------------------------------
template<typename Dimension, typename CellValue, typename LeafPolicy>
inline
void
Tree<Dimension, CellValue, LeafPolicy>::
buildCellKey(const typename Tree<Dimension, CellValue, LeafPolicy>::LevelKey ilevel,
             const typename Dimension::Vector& xi,
             typename Tree<Dimension, CellValue, LeafPolicy>::CellKey& key,
             typename Tree<Dimension, CellValue, LeafPolicy>::CoordHash& ix,
             typename Tree<Dimension, CellValue, LeafPolicy>::CoordHash& iy,
             typename Tree<Dimension, CellValue, LeafPolicy>::CoordHash& iz) const {
  // Dispatch to the specialized functor.
  CellKeyComputer<Dimension, CellValue, LeafPolicy>::tokey(ilevel, xi, key, ix, iy, iz);
}

//------------------------------------------------------------------------------
// Extract the individual coordinate indices from a cell index.
//------------------------------------------------------------------------------
template<typename Dimension, typename CellValue, typename LeafPolicy>
inline
void
Tree<Dimension, CellValue, LeafPolicy>::
extractCellIndices(const typename Tree<Dimension, CellValue, LeafPolicy>::CellKey& key,
                   typename Tree<Dimension, CellValue, LeafPolicy>::CoordHash& ix,
                   typename Tree<Dimension, CellValue, LeafPolicy>::CoordHash& iy,
                   typename Tree<Dimension, CellValue, LeafPolicy>::CoordHash& iz) const {
  // Dispatch to the specialized functor.
  CellKeyComputer<Dimension, CellValue, LeafPolicy>::fromkey(key, ix, iy, iz);
}

//------------------------------------------------------------------------------
// Add a daughter to a cell if not present.
//------------------------------------------------------------------------------
template<typename Dimension, typename CellValue, typename LeafPolicy>
inline
void
Tree<Dimension, CellValue, LeafPolicy>::
addDaughter(typename Tree<Dimension, CellValue, LeafPolicy>::Cell& cell,
            const Tree<Dimension, CellValue, LeafPolicy>::CellKey& daughterKey) const {
  if (std::find(cell.daughters.begin(), cell.daughters.end(), daughterKey) == cell.daughters.end())
    cell.daughters.push_back(daughterKey);
  ENSURE(cell.daughters.size() <= (1U << Dimension::nDim));
}

//------------------------------------------------------------------------------
// Add a node to the internal Tree structure.
//------------------------------------------------------------------------------
template<typename Dimension, typename CellValue, typename LeafPolicy>
template<typename NodeID, typename CellValueFactory>
inline
void
Tree<Dimension, CellValue, LeafPolicy>::
addNodeToTree(const Tree::Vector& xi, 
              const NodeID& nodeID,
              CellValueFactory& factory) {
  const unsigned nb1d = this->num1dbits();
  mTree.reserve(nb1d);                    // This is necessary to avoid memory errors!

  LevelKey ilevel = 0;
  bool terminated = false;
  CellKey key, parentKey, otherKey, ix, iy, iz;
  TreeLevel::iterator itr;
  while (ilevel < nb1d and not terminated) {

    // Do we need to add another level to the tree?
    if (ilevel == mTree.size() and ilevel < nb1d - 1) mTree.push_back(TreeLevel());

    // Create the key for the cell containing this particle on this level.
    buildCellKey(ilevel, xi, key, ix, iy, iz);
    itr = mTree[ilevel].find(key);

    // Check if this node should terminate on this level or not.
    terminated = (ilevel == nb1d - 1 or
                  LeafPolicy::terminate(xi, nodeID, ilevel, *this));

    if (itr == mTree[ilevel].end()) {

      // This is a new cell.
      mTree[ilevel][key] = Cell(key, xi, factory.newCellValue(nodeID));

    } else {

      // This is an existing cell.  Augment it's value as needed.
      Cell& cell = itr->second;
      factory.augmentCellValue(xi, nodeID, ilevel, cell, *this);

      if (terminated) {

        // If we're terminating in this cell just add to the member data.
        cell.positions.push_back(xi);
        factory.terminateNodeInCell(xi, nodeID, ilevel, cell, *this);

      } else {

        // Fork a new descendant on the next level.
        const LevelKey ilevel1 = ilevel + 1;
        CHECK(ilevel1 < Tree::num1dbits);
        if (ilevel1 == mTree.size()) mTree.push_back(TreeLevel());
        buildCellKey(ilevel1, cell.xcm, otherKey, ix, iy, iz);
          mTree[ilevel1][otherKey] = Cell(cell.M, cell.xcm, otherKey);
          cell.daughters = std::vector<CellKey>(1, otherKey);


        if (cell.daughters.find(key)


      // Is this cell a single leaf already?
      if (cell.positions.size() > 0) {
        CHECK(cell.masses.size() == cell.positions.size());
        CHECK(cell.daughters.size() == 0);

        // Yep, so we need to split it unless we're at the maximum refinement.
        if (ilevel < Tree::num1dbits - 1) {
          CHECK(cell.masses.size() == 1);
          const LevelKey ilevel1 = ilevel + 1;
          CHECK(ilevel1 < Tree::num1dbits);
          if (ilevel1 == mTree.size()) mTree.push_back(TreeLevel());
          buildCellKey(ilevel1, cell.xcm, otherKey, ix, iy, iz);
          mTree[ilevel1][otherKey] = Cell(cell.M, cell.xcm, otherKey);
          cell.daughters = std::vector<CellKey>(1, otherKey);
          cell.masses = std::vector<double>();
          cell.positions = std::vector<Vector>();

        } else {
          // If we've maxed out the levels, then we just huck this node in 
          // the members of this cell.
          cell.masses.push_back(mi);
          cell.positions.push_back(xi);
        }
      }

      // Increment the cell moments.
      cell.xcm = (cell.M*cell.xcm + mi*xi)/(cell.M + mi);
      cell.M += mi;
      cell.Mglobal = cell.M;
    }

    // Link this cell as a daughter of its parent.
    if (ilevel > 0) {
      CHECK(mTree[ilevel - 1].find(parentKey) != mTree[ilevel - 1].end());
      addDaughter(mTree[ilevel - 1][parentKey], key);
    }

    parentKey = key;
    ++ilevel;
  }

}

//------------------------------------------------------------------------------
// Construct the daughter pointers in a tree.
//------------------------------------------------------------------------------
inline
void
Tree::
constructDaughterPtrs(Tree::Tree& tree) const {
  const unsigned nlevels = tree.size();
  for (unsigned ilevel = 0; ilevel < nlevels - 1; ++ilevel) {
    const unsigned ilevel1 = ilevel + 1;
    for (TreeLevel::iterator itr = tree[ilevel].begin();
         itr != tree[ilevel].end();
         ++itr) {
      Cell& cell = itr->second;
      cell.daughterPtrs = std::vector<Cell*>();
      for (std::vector<CellKey>::const_iterator ditr = cell.daughters.begin();
           ditr != cell.daughters.end();
           ++ditr) {
        cell.daughterPtrs.push_back(&(tree[ilevel1][*ditr]));
      }
      CHECK(cell.daughters.size() == cell.daughterPtrs.size());
    }
  }
}

}
