//---------------------------------Spheral++----------------------------------//
// OctTreeGravity -- An implementation of the tree n-body gravity solver.
//
// Created by JMO, 2012-02-28
//----------------------------------------------------------------------------//

namespace Spheral {
namespace GravitySpace {

//------------------------------------------------------------------------------
// Build a cell key from coordinate indices.
//------------------------------------------------------------------------------
inline
OctTreeGravity::CellKey
OctTreeGravity::
buildCellKey(const OctTreeGravity::CellKey ix,
             const OctTreeGravity::CellKey iy,
             const OctTreeGravity::CellKey iz) const {
  return (std::max(CellKey(0), std::min(max1dKey, iz << 2*num1dbits)) +
          std::max(CellKey(0), std::min(max1dKey, iy <<   num1dbits)) +
          std::max(CellKey(0), std::min(max1dKey, ix)));
}

//------------------------------------------------------------------------------
// Extract the individual coordinate indices from a cell index.
//------------------------------------------------------------------------------
inline
void
OctTreeGravity::
extractCellIndices(const OctTreeGravity::CellKey& key,
                   OctTreeGravity::CellKey& ix,
                   OctTreeGravity::CellKey& iy,
                   OctTreeGravity::CellKey& iz) const {
  ix = key % (1U << num1dbits);
  iy = (key >> num1dbits) % num1dbits;
  iz = (key >> 2*num1dbits);
}

//------------------------------------------------------------------------------
// Build the key for a cell based on a position and level.
//------------------------------------------------------------------------------
inline
OctTreeGravity::TreeKey
OctTreeGravity::
buildTreeKey(const OctTreeGravity::LevelKey ilevel,
             const OctTreeGravity::Vector& xi) const {
  REQUIRE(xi.x() >= mXmin.x() and xi.x() <= mXmax.x());
  REQUIRE(xi.y() >= mXmin.y() and xi.y() <= mXmax.y());
  REQUIRE(xi.z() >= mXmin.z() and xi.z() <= mXmax.z());
  const CellKey ckey = buildCellKey(CellKey((xi.x() - mXmin.x())/mBoxLength * max1dKey1 + 0.5),
                                    CellKey((xi.y() - mXmin.y())/mBoxLength * max1dKey1 + 0.5),
                                    CellKey((xi.z() - mXmin.z())/mBoxLength * max1dKey1 + 0.5));
  return std::make_pair(ilevel, ckey);
}

//------------------------------------------------------------------------------
// Add a daughter to a cell if not present.
//------------------------------------------------------------------------------
inline
void
OctTreeGravity::
addDaughter(OctTreeGravity::Cell& cell,
            const OctTreeGravity::CellKey daughterKey) const {
  if (std::find(cell.daughters.begin(), cell.daughters.end(), daughterKey) == cell.daughters.end())
    cell.daughters.push_back(daughterKey);
  ENSURE(cell.daughters.size() <= 8);
}

//------------------------------------------------------------------------------
// Add a node to the internal Tree structure.
//------------------------------------------------------------------------------
inline
void
OctTreeGravity::
addNodeToTree(const size_t nodeListi,
              const size_t i,
              const double mi,
              const OctTreeGravity::Vector& xi) {
  LevelKey ilevel = 0;
  const NodeID nodeID = std::make_pair(nodeListi, i);
  bool terminated = false;
  TreeKey key, parentKey;
  Tree::iterator itr;
  while (ilevel < OctTreeGravity::num1dbits and not terminated) {
    key = buildTreeKey(ilevel, xi);
    itr = mTree.find(key);

    if (itr == mTree.end()) {
      // If this is an unregistered cell, add it with this node as the sole leaf
      // and we're done.
      terminated = true;
      mTree[key] = Cell(xi, mi, nodeID);

    } else {
      Cell& cell = itr->second;

      // Is this cell a single leaf already?
      if (cell.members.size() > 0) {
        CHECK(cell.daughters.size() == 0);

        // Yep, so we need to split it unless we're at the maximum refinement.
        if (ilevel < OctTreeGravity::num1dbits - 1) {
          CHECK(cell.members.size() == 1);
          const LevelKey ilevel1 = ilevel + 1;
          CHECK(ilevel1 < OctTreeGravity::num1dbits);
          const TreeKey otherKey = buildTreeKey(ilevel1, cell.xcm);
          mTree[otherKey] = Cell(cell.xcm, cell.M, cell.members[0]);
          cell.daughters = std::vector<CellKey>(1, otherKey.second);
          cell.members = std::vector<NodeID>();

        } else {
          // If we've maxed out the levels, then we just huck this node in 
          // the members of this cell.
          cell.members.push_back(nodeID);
        }
      }

      // Increment the cell moments.
      cell.xcm = (cell.M*cell.xcm + mi*xi)/(cell.M + mi);
      cell.M += mi;
    }

    // Link this cell as a daughter of its parent.
    if (ilevel > 0) {
      CHECK(mTree.find(parentKey) != mTree.end());
      addDaughter(mTree[parentKey], key.second);
    }

    parentKey = key;
    ++ilevel;
  }

}

}
}
