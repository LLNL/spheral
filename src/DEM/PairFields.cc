//---------------------------------Spheral++----------------------------------//
// PairFields -- The DEM package for Spheral++.
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"

#include "Hydro/HydroFieldNames.hh"
#include "Hydro/PositionPolicy.hh"

#include "Physics/Physics.hh"

#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/DataBase.hh"

#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"

#include "Boundary/Boundary.hh"

#include "Neighbor/ConnectivityMap.hh"

#include "Utilities/timingUtilities.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/globalBoundingVolumes.hh"
#include "Utilities/Timer.hh"
#include "Utilities/globalNodeIDs.hh"
#include "Utilities/registerWithRedistribution.hh"

#include "DEM/PairFields.hh"
#include "DEM/DEMDimension.hh"
#include "DEM/DEMFieldNames.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

#include <limits.h>
#include <float.h>
#include <algorithm>
#include <fstream>
#include <map>
#include <vector>

using std::vector;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;


namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
template<typename Dimension>
PairFields<Dimension>::
PairFields(DataBase<Dimension>& dataBase):
  mUniqueIndices(FieldStorageType::CopyFields),
  mNeighborIndices(FieldStorageType::CopyFields),
  mShearDisplacement(FieldStorageType::CopyFields),
  mDDtShearDisplacement(FieldStorageType::CopyFields),
  mEquilibriumOverlap(FieldStorageType::CopyFields),
  mRestart(registerWithRestart(*this)),
  mRedistribute(registerWithRedistribution(*this,
                                           &PairFields<Dimension>::initializeBeforeRedistribution,
                                           &PairFields<Dimension>::finalizeAfterRedistribution)){
    mUniqueIndices = dataBase.newDEMFieldList(int(0),  DEMFieldNames::uniqueIndices);
    mUniqueIndices = globalNodeIDs<Dimension>(dataBase.nodeListBegin(),dataBase.nodeListEnd());
    mNeighborIndices = dataBase.newDEMFieldList(std::vector<int>(), DEMFieldNames::neighborIndices);
    mShearDisplacement = dataBase.newDEMFieldList(std::vector<Vector>(), DEMFieldNames::shearDisplacement);
    mDDtShearDisplacement = dataBase.newDEMFieldList(std::vector<Vector>(), IncrementFieldList<Dimension, Scalar>::prefix() + DEMFieldNames::shearDisplacement);
    mEquilibriumOverlap = dataBase.newDEMFieldList(std::vector<Scalar>(), DEMFieldNames::equilibriumOverlap);

}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
PairFields<Dimension>::
~PairFields() {
}



//------------------------------------------------------------------------------
// REdistribution methods
//------------------------------------------------------------------------------
template<typename Dimension>
void
PairFields<Dimension>::
finalizeAfterRedistribution() {
  std::cout << "SHABINGO!"<< std::endl;
}

template<typename Dimension>
void
PairFields<Dimension>::
initializeBeforeRedistribution(){
  std::cout << "HERE WE ARE BABY!"<< std::endl;
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PairFields<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  file.write(mUniqueIndices, pathName + "/uniqueIndices");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PairFields<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  file.read(mUniqueIndices, pathName + "/uniqueIndices");
}

}
