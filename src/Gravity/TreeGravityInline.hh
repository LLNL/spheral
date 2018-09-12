//---------------------------------Spheral++----------------------------------//
// TreeGravity -- An implementation of the tree n-body gravity solver.
// Based on the original 3D only OctTreeGravity.
//
// Created by JMO, 2013-06-12
//----------------------------------------------------------------------------//

namespace Spheral {

//------------------------------------------------------------------------------
// Build a cell key from coordinate indices.
//------------------------------------------------------------------------------
template<>
inline
void
TreeGravity<Dim<2> >::
buildCellKey(const TreeGravity<Dim<2> >::LevelKey ilevel,
             const TreeGravity<Dim<2> >::Vector& xi,
             TreeGravity<Dim<2> >::CellKey& key,
             TreeGravity<Dim<2> >::CellKey& ix,
             TreeGravity<Dim<2> >::CellKey& iy,
             TreeGravity<Dim<2> >::CellKey& iz) const {
  REQUIRE(xi.x() >= mXmin.x() and xi.x() <= mXmax.x());
  REQUIRE(xi.y() >= mXmin.y() and xi.y() <= mXmax.y());
  const CellKey ncell = (1U << ilevel);
  const CellKey maxcell = ncell - 1U;
  ix = std::min(maxcell, CellKey((xi.x() - mXmin.x())/mBoxLength * ncell));
  iy = std::min(maxcell, CellKey((xi.y() - mXmin.y())/mBoxLength * ncell));
  iz = 0U;
  key = ((std::max(CellKey(0), std::min(max1dKey, iy)) <<   num1dbits) +
         (std::max(CellKey(0), std::min(max1dKey, ix))));
}

template<>
inline
void
TreeGravity<Dim<3> >::
buildCellKey(const TreeGravity<Dim<3> >::LevelKey ilevel,
             const TreeGravity<Dim<3> >::Vector& xi,
             TreeGravity<Dim<3> >::CellKey& key,
             TreeGravity<Dim<3> >::CellKey& ix,
             TreeGravity<Dim<3> >::CellKey& iy,
             TreeGravity<Dim<3> >::CellKey& iz) const {
  REQUIRE(xi.x() >= mXmin.x() and xi.x() <= mXmax.x());
  REQUIRE(xi.y() >= mXmin.y() and xi.y() <= mXmax.y());
  REQUIRE(xi.z() >= mXmin.z() and xi.z() <= mXmax.z());
  const CellKey ncell = (1U << ilevel);
  const CellKey maxcell = ncell - 1U;
  ix = std::min(maxcell, CellKey((xi.x() - mXmin.x())/mBoxLength * ncell));
  iy = std::min(maxcell, CellKey((xi.y() - mXmin.y())/mBoxLength * ncell));
  iz = std::min(maxcell, CellKey((xi.z() - mXmin.z())/mBoxLength * ncell));
  key = ((std::max(CellKey(0), std::min(max1dKey, iz)) << 2*num1dbits) +
         (std::max(CellKey(0), std::min(max1dKey, iy)) <<   num1dbits) +
         (std::max(CellKey(0), std::min(max1dKey, ix))));
}

//------------------------------------------------------------------------------
// Extract the individual coordinate indices from a cell index.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
TreeGravity<Dimension>::
extractCellIndices(const TreeGravity<Dimension>::CellKey& key,
                   TreeGravity<Dimension>::CellKey& ix,
                   TreeGravity<Dimension>::CellKey& iy,
                   TreeGravity<Dimension>::CellKey& iz) const {
  ix = key & xkeymask;
  iy = (key & ykeymask) >> num1dbits;
  iz = (key & zkeymask) >> 2*num1dbits;
}

//------------------------------------------------------------------------------
// Add a daughter to a cell if not present.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
TreeGravity<Dimension>::
addDaughter(TreeGravity<Dimension>::Cell& cell,
            const TreeGravity<Dimension>::CellKey daughterKey) const {
  if (std::find(cell.daughters.begin(), cell.daughters.end(), daughterKey) == cell.daughters.end())
    cell.daughters.push_back(daughterKey);
  ENSURE(cell.daughters.size() <= Dimension::pownu(2));
}

//------------------------------------------------------------------------------
// Add a node to the internal Tree structure.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
TreeGravity<Dimension>::
addNodeToTree(const double mi,
              const TreeGravity<Dimension>::Vector& xi,
              const TreeGravity<Dimension>::Vector& vi) {
  mTree.reserve(num1dbits); // This is necessary to avoid memory errors!

  LevelKey ilevel = 0;
  bool terminated = false;
  CellKey key, parentKey, otherKey, ix, iy, iz;
  typename TreeLevel::iterator itr;
  while (ilevel < TreeGravity<Dimension>::num1dbits and not terminated) {

    // Do we need to add another level to the tree?
    if (ilevel == mTree.size()) mTree.push_back(TreeLevel());

    // Create the key for the cell containing this particle on this level.
    buildCellKey(ilevel, xi, key, ix, iy, iz);
    itr = mTree[ilevel].find(key);

    if (itr == mTree[ilevel].end()) {
      // If this is an unregistered cell, add it with this node as the sole leaf
      // and we're done.
      terminated = true;
      mTree[ilevel][key] = Cell(mi, xi, vi, key);

    } else {
      Cell& cell = itr->second;

      // Is this cell a single leaf already?
      if (cell.masses.size() > 0) {
        CHECK(cell.masses.size() == cell.positions.size() and
              cell.masses.size() == cell.velocities.size());
        CHECK(cell.daughters.size() == 0);

        // Yep, so we need to split it unless we're at the maximum refinement.
        if (ilevel < TreeGravity<Dimension>::num1dbits - 1) {
          CHECK(cell.masses.size() == 1);
          const LevelKey ilevel1 = ilevel + 1;
          CHECK(ilevel1 < TreeGravity<Dimension>::num1dbits);
          if (ilevel1 == mTree.size()) mTree.push_back(TreeLevel());
          buildCellKey(ilevel1, cell.xcm, otherKey, ix, iy, iz);
          mTree[ilevel1][otherKey] = Cell(cell.M, cell.xcm, cell.vcm, otherKey);
          cell.daughters = std::vector<CellKey>(1, otherKey);
          cell.masses = std::vector<double>();
          cell.positions = std::vector<Vector>();
          cell.velocities = std::vector<Vector>();

        } else {
          // If we've maxed out the levels, then we just huck this node in 
          // the members of this cell.
          cell.masses.push_back(mi);
          cell.positions.push_back(xi);
          cell.velocities.push_back(vi);
        }
      }

      // Increment the cell moments.
      cell.xcm = (cell.M*cell.xcm + mi*xi)/(cell.M + mi);
      cell.vcm = (cell.M*cell.vcm + mi*vi)/(cell.M + mi);
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
template<typename Dimension>
inline
void
TreeGravity<Dimension>::
constructDaughterPtrs(TreeGravity<Dimension>::Tree& tree) const {
  const unsigned nlevels = tree.size();
  const unsigned n = nlevels > 0 ? nlevels - 1 : nlevels;
  unsigned ilevel, ilevel1;
  for (ilevel = 0; ilevel != n; ++ilevel) {
    ilevel1 = ilevel + 1;
    for (typename TreeLevel::iterator itr = tree[ilevel].begin();
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
