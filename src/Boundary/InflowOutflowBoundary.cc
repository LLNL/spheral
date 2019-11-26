//---------------------------------Spheral++----------------------------------//
// InflowOutflowBoundary -- creates inflow ghost images, which become internal nodes
// as they cross the specified boundary plane.
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "Boundary.hh"
#include "findNodesTouchingThroughPlanes.hh"
#include "mapPositionThroughPlanes.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Field/FieldBase.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Utilities/allReduce.hh"
#include "Utilities/planarReflectingOperator.hh"
#include "Utilities/DBC.hh"

#include "InflowOutflowBoundary.hh"

#include "boost/lexical_cast.hpp"
#include "boost/algorithm/string.hpp"

using std::vector;
using std::map;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

namespace {

//------------------------------------------------------------------------------
// Store the Field values of the given DataType for the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
void
storeFieldValues(const NodeList<Dimension>& nodeList,
                 const std::vector<int>& nodeIDs,
                 std::map<std::string, std::vector<DataType> >& values) {

  // Iterate over all the Fields defined on the NodeList.
  for (auto fieldItr = nodeList.registeredFieldsBegin();
       fieldItr != nodeList.registeredFieldsEnd();
       ++fieldItr) {

    // Determine if this Field is the type we're looking for.
    if (typeid(**fieldItr) == typeid(Field<Dimension, DataType>)) {
      const auto& field = (const Field<Dimension, DataType>&) **fieldItr;

      // Build a vector of the values of this field on the requested nodes.
      std::vector<DataType> vals;
      vals.reserve(nodeIDs.size());
      for (typename std::vector<int>::const_iterator nodeItr = nodeIDs.begin();
           nodeItr != nodeIDs.end();
           ++nodeItr) {
        CHECK2(*nodeItr >= 0 && *nodeItr < field.numElements(), *nodeItr << " " << field.numElements());
        vals.push_back(field(*nodeItr));
      }
      CHECK(vals.size() == nodeIDs.size());

      // Put this data into the result.
      // cerr << " --> ";
      // for (typename std::map<std::string, std::vector<DataType> >::const_iterator itr = values.begin();
      //      itr != values.end();
      //      ++itr) cerr << itr->first << " ";
      // cerr << endl;
      const auto key = StateBase<Dimension>::key(**fieldItr);
      CHECK2(values.find(key) == values.end(), key);
      values[key] = vals;
      // cerr << "Stored " << vals.size() << " values for " << key << endl;
    }
  }
}
  
//------------------------------------------------------------------------------
// Set the Field values in the given field using the given map.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
void
resetValues(Field<Dimension, DataType>& field,
            const std::vector<int>& nodeIDs,
            const std::map<std::string, std::vector<DataType> >& values,
            const bool dieOnMissingField) {

  const auto& nodeList = field.nodeList();

  // Find this Field in the set of stored values.
  const auto key = StateBase<Dimension>::key(field);
  auto itr = values.find(key);
  VERIFY2(itr != values.end() or not dieOnMissingField,
          "InflowOutflowBoundary error: " << key << " not found in stored field values.");

  // Now set the values.
  if (itr != values.end()) {
    const auto& newValues = itr->second;
    CHECK2(newValues.size() == nodeIDs.size(), key << " " << newValues.size() << " " << nodeIDs.size());
    for (int i = 0; i < nodeIDs.size(); ++i) {
      CHECK2(nodeIDs[i] >= 0 &&
             nodeIDs[i] < nodeList.numNodes(), nodeIDs[i] << " " << nodeList.numNodes());
      field(nodeIDs[i]) = newValues[i];
      // cerr << "Setting " << field.name() << "(" << nodeIDs[i] << ")" << endl;
    }
  }
}

//------------------------------------------------------------------------------
// Clear the values in the storage
//------------------------------------------------------------------------------
template<typename DataType>
void
clearValues(std::map<std::string, std::vector<DataType> >& values) {
  for (auto& pairvals: values) pairvals.second.clear();
}

//------------------------------------------------------------------------------
// Copy Field values of the given DataType from control to target IDs.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
void
copyFieldValues(NodeList<Dimension>& nodeList,
                const std::map<std::string, std::vector<DataType>>& storage,
                const std::vector<int>& controlIDs,
                const std::vector<int>& targetIDs) {
  REQUIRE(controlIDs.size() == targetIDs.size());
  const auto n = controlIDs.size();

  // Iterate over all the Fields defined on the NodeList.
  for (auto fieldItr = nodeList.registeredFieldsBegin();
       fieldItr != nodeList.registeredFieldsEnd();
       ++fieldItr) {

    // Determine if this Field is the type we're looking for.
    if (typeid(**fieldItr) == typeid(Field<Dimension, DataType>)) {
      auto& field = (Field<Dimension, DataType>&) **fieldItr;
      if ((*fieldItr)->name() != HydroFieldNames::position) {
        const auto key = StateBase<Dimension>::key(field);
        const auto itr = storage.find(key);
        if (itr != storage.end()) {
          // CHECK2(itr != storage.end(), "Key not found: " << key);
          const auto& vals = itr->second;

          // Copy the requested entries.
          for (auto k = 0; k < n; ++k) field[targetIDs[k]] = vals[controlIDs[k]];
        }
      }
    }
  }
}
  
}

//------------------------------------------------------------------------------
// Construct with the given NodeList and plane.
//------------------------------------------------------------------------------
template<typename Dimension>
InflowOutflowBoundary<Dimension>::
InflowOutflowBoundary(DataBase<Dimension>& dataBase,
                      const GeomPlane<Dimension>& plane,
                      const bool empty):
  Boundary<Dimension>(),
  Physics<Dimension>(),
  mDataBase(dataBase),
  mPlane(plane),
  mBoundaryCount(dataBase.numNodeLists()),
  mDT(1e100),
  mActive(false),
  mEmpty(empty),
  mNumInflowNodes(),
  mXmin(),
  mIntValues(),
  mScalarValues(),
  mVectorValues(),
  mTensorValues(),
  mSymTensorValues(),
  mThirdRankTensorValues(),
  mFourthRankTensorValues(),
  mFifthRankTensorValues(),
  mFacetedVolumeValues(),
  mVectorScalarValues(),
  mVectorVectorValues(),
  mRestart(registerWithRestart(*this)) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
InflowOutflowBoundary<Dimension>::~InflowOutflowBoundary() {
}

//------------------------------------------------------------------------------
// setGhostNodes
//------------------------------------------------------------------------------
template<typename Dimension>
void
InflowOutflowBoundary<Dimension>::
setGhostNodes(NodeList<Dimension>& nodeList) {
  this->addNodeList(nodeList);

  if (mActive) {
    auto& boundaryNodes = this->accessBoundaryNodes(nodeList);
    auto& cNodes = boundaryNodes.controlNodes;
    auto& gNodes = boundaryNodes.ghostNodes;
    const auto currentNumGhostNodes = nodeList.numGhostNodes();
    const auto firstNewGhostNode = nodeList.numNodes();
    // cerr << "Allocating new ghost nodes " << firstNewGhostNode << " -- " << (firstNewGhostNode + mNumInflowNodes[nodeList]) << endl;
    
    // Use the planar boundary to find the set of points that interact with
    // the entrance plane.  We make these the control nodes.
    cNodes = findNodesTouchingThroughPlanes(nodeList, mPlane, mPlane);

    // Create the ghost nodes.
    nodeList.numGhostNodes(currentNumGhostNodes + mNumInflowNodes[nodeList.name()]);
    gNodes = vector<int>(mNumInflowNodes[nodeList.name()]);
    for (auto i = 0; i < mNumInflowNodes[nodeList.name()]; ++i) gNodes[i] = firstNewGhostNode + i;
    this->updateGhostNodes(nodeList);
  }
}

//------------------------------------------------------------------------------
// updateGhostNodes
// For this boundary this just means moving the ghost points
//------------------------------------------------------------------------------
template<typename Dimension>
void
InflowOutflowBoundary<Dimension>::
updateGhostNodes(NodeList<Dimension>& nodeList) {
  if (mActive) {
    auto& pos = nodeList.positions();
    this->applyGhostBoundary(pos);
    this->applyGhostBoundary(nodeList.Hfield());

    const auto& boundaryNodes = this->accessBoundaryNodes(nodeList);
    const auto& cNodes = boundaryNodes.controlNodes;
    const auto& gNodes = boundaryNodes.ghostNodes;
    const auto& nhat = mPlane.normal();

    // Find how close the control nodes are to the entrance plane.
    Scalar xmin = 1e100;
    for (const auto i: cNodes) {
      const auto xd = mPlane.signedDistance(pos[i]);
      xmin = std::min(xmin, xd);
    }
    xmin = allReduce(xmin, MPI_MIN, Communicator::communicator());
    // CHECK(xmin >= 0.0);

    // Offset the current ghost points appropriately.
    const auto delta = (xmin < 1e100 ?
                        xmin - mXmin[nodeList.name()] :
                        0.0)*nhat;
    // cerr << " ************> " << xmin << " " << mXmin[nodeList.name()] << " " << nhat << " " << delta << endl;
    for (const auto i: gNodes) pos[i] += delta;

    // for (const auto i: gNodes) {
    //   cerr << " --> " << i << " " << pos(i) << " " << nodeList.Hfield()(i) << " : " << (pos(i) - pos(0)) << endl;
    // }
  }
}

//------------------------------------------------------------------------------
// Apply the boundary condition to fields of different DataTypes.
//------------------------------------------------------------------------------
// Specialization for int fields.
template<typename Dimension>
void
InflowOutflowBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, int>& field) const {
  if (mActive) {
    resetValues(field, this->ghostNodes(field.nodeList()), mIntValues, false);
  }
}

// Specialization for scalar fields.
template<typename Dimension>
void
InflowOutflowBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Scalar>& field) const {
  if (mActive) {
    resetValues(field, this->ghostNodes(field.nodeList()), mScalarValues, false);
  }
}

// Specialization for Vector fields.
template<typename Dimension>
void
InflowOutflowBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Vector>& field) const {
  if (mActive) {
    resetValues(field, this->ghostNodes(field.nodeList()), mVectorValues, false);
  }
}

// Specialization for Tensor fields.
template<typename Dimension>
void
InflowOutflowBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Tensor>& field) const {
  if (mActive) {
    resetValues(field, this->ghostNodes(field.nodeList()), mTensorValues, false);
  }
}

// Specialization for symmetric tensors.
template<typename Dimension>
void
InflowOutflowBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::SymTensor>& field) const {
  if (mActive) {
    resetValues(field, this->ghostNodes(field.nodeList()), mSymTensorValues, false);
  }
}

// Specialization for third rank tensors.
template<typename Dimension>
void
InflowOutflowBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::ThirdRankTensor>& field) const {
  if (mActive) {
    resetValues(field, this->ghostNodes(field.nodeList()), mThirdRankTensorValues, false);
  }
}


// Specialization for fourth rank tensors.
template<typename Dimension>
void
InflowOutflowBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::FourthRankTensor>& field) const {
  if (mActive) {
    resetValues(field, this->ghostNodes(field.nodeList()), mFourthRankTensorValues, false);
  }
}

// Specialization for fifth rank tensors.
template<typename Dimension>
void
InflowOutflowBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::FifthRankTensor>& field) const {
  if (mActive) {
    resetValues(field, this->ghostNodes(field.nodeList()), mFifthRankTensorValues, false);
  }
}

// Specialization for FacetedVolume.
template<typename Dimension>
void
InflowOutflowBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::FacetedVolume>& field) const {
  if (mActive) {
    resetValues(field, this->ghostNodes(field.nodeList()), mFacetedVolumeValues, false);
  }
}

// Specialization for vector<Scalar> fields.
template<typename Dimension>
void
InflowOutflowBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, std::vector<typename Dimension::Scalar>>& field) const {
  if (mActive) {
    resetValues(field, this->ghostNodes(field.nodeList()), mVectorScalarValues, false);
  }
}

// Specialization for vector<Vector> fields.
template<typename Dimension>
void
InflowOutflowBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, std::vector<typename Dimension::Vector>>& field) const {
  if (mActive) {
    resetValues(field, this->ghostNodes(field.nodeList()), mVectorVectorValues, false);
  }
}

//------------------------------------------------------------------------------
// Provide the setViolationNodes for a NodeList.
// In this case there shouldn't be any violation nodes -- if internal nodes
// are fighting upstream and through the entrance plane, there are bigger
// problems.
//------------------------------------------------------------------------------
template<typename Dimension>
void
InflowOutflowBoundary<Dimension>::
setViolationNodes(NodeList<Dimension>& nodeList) {
  this->addNodeList(nodeList);
}

//------------------------------------------------------------------------------
// Provide the updateViolationNodes for a NodeList, correcting positions and 
// H's.
//------------------------------------------------------------------------------
template<typename Dimension>
void
InflowOutflowBoundary<Dimension>::
updateViolationNodes(NodeList<Dimension>& nodeList) {
}

//------------------------------------------------------------------------------
// Apply the boundary condition to fields of different DataTypes.
//------------------------------------------------------------------------------
// Specialization for int fields.
template<typename Dimension>
void
InflowOutflowBoundary<Dimension>::
enforceBoundary(Field<Dimension, int>& field) const {
}

// Specialization for scalar fields.
template<typename Dimension>
void
InflowOutflowBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Scalar>& field) const {
}


// Specialization for Vector fields.
template<typename Dimension>
void
InflowOutflowBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Vector>& field) const {
}

// Specialization for Tensor fields.
template<typename Dimension>
void
InflowOutflowBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Tensor>& field) const {
}

// Specialization for symmetric tensors.
template<typename Dimension>
void
InflowOutflowBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::SymTensor>& field) const {
}

// Specialization for third rank tensors.
template<typename Dimension>
void
InflowOutflowBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::ThirdRankTensor>& field) const {
}

// Specialization for fourth rank tensors.
template<typename Dimension>
void
InflowOutflowBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::FourthRankTensor>& field) const {
}

// Specialization for fifth rank tensors.
template<typename Dimension>
void
InflowOutflowBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::FifthRankTensor>& field) const {
}

// Specialization for FacetedVolume.
template<typename Dimension>
void
InflowOutflowBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::FacetedVolume>& field) const {
}

//------------------------------------------------------------------------------
// On problem startup we take our snapshot of the state of the points that
// see/interact with the boundary plane.
//------------------------------------------------------------------------------
template<typename Dimension>
void
InflowOutflowBoundary<Dimension>::initializeProblemStartup() {

  // Clear any existing data.
  mIntValues.clear();
  mScalarValues.clear();
  mVectorValues.clear();
  mTensorValues.clear();
  mSymTensorValues.clear();
  mThirdRankTensorValues.clear();
  mFourthRankTensorValues.clear();
  mFifthRankTensorValues.clear();
  mFacetedVolumeValues.clear();
  mVectorScalarValues.clear();
  mVectorVectorValues.clear();

  // Check all NodeLists.
  for (auto itr = mDataBase.nodeListBegin(); itr < mDataBase.nodeListEnd(); ++itr) {
    auto& nodeList = **itr;
    // cerr << "--------------------------------------------------------------------------------" << endl
    //      << nodeList.name() << endl;

    // Use a planar boundary to figure out what sort of nodes are in range of the plane.
    // We use those to create a stencil of the in/outflow conditions.
    const auto& nhat = mPlane.normal();
    auto nodeIDs = findNodesTouchingThroughPlanes(nodeList, mPlane, mPlane, 1.0);
    if (mEmpty) nodeIDs.clear();

    // Remove any ghost nodes from other BCs.
    const auto firstGhostNode = nodeList.firstGhostNode();
    nodeIDs.erase(std::remove_if(nodeIDs.begin(), nodeIDs.end(), [&](const int& x) { return x >= firstGhostNode; }), nodeIDs.end());

    // cerr << "Node IDs: ";
    // std::copy(nodeIDs.begin(), nodeIDs.end(), std::ostream_iterator<int>(std::cerr, " "));
    // cerr << endl;

    // Now take a snapshot of the Fields.
    storeFieldValues(nodeList, nodeIDs, mIntValues);
    storeFieldValues(nodeList, nodeIDs, mScalarValues);
    storeFieldValues(nodeList, nodeIDs, mVectorValues);
    storeFieldValues(nodeList, nodeIDs, mTensorValues);
    storeFieldValues(nodeList, nodeIDs, mSymTensorValues);
    storeFieldValues(nodeList, nodeIDs, mThirdRankTensorValues);
    storeFieldValues(nodeList, nodeIDs, mFourthRankTensorValues);
    storeFieldValues(nodeList, nodeIDs, mFifthRankTensorValues);
    storeFieldValues(nodeList, nodeIDs, mFacetedVolumeValues);
    storeFieldValues(nodeList, nodeIDs, mVectorScalarValues);
    storeFieldValues(nodeList, nodeIDs, mVectorVectorValues);

    // Map the snapshot positions to put them outside the boundary.
    auto& pos = nodeList.positions();
    const auto poskey = StateBase<Dimension>::key(pos);
    CHECK(mVectorValues.find(poskey) != mVectorValues.end());
    auto& posvals = mVectorValues[poskey];
    const auto ni = nodeIDs.size();
    const GeomPlane<Dimension> exitPlane(mPlane.point(), -nhat);
    for (auto k = 0; k < ni; ++k) {
      const auto i = nodeIDs[k];
      posvals[k] = mapPositionThroughPlanes(pos[i], mPlane, mPlane);
      // cerr << "  Ghost position: " << i << " @ " << posvals[k] << endl;
    }

    // Determine the in/outflow velocity.
    const auto& vel = nodeList.velocity();
    Scalar vinflow = 0.0;
    for (const auto i: nodeIDs) {
      // CHECK(std::abs(vel[i].dot(nhat)/vel[i].magnitude() - 1.0) < 1.0e-5);
      vinflow += vel[i].dot(nhat);
    }
    vinflow = (allReduce(vinflow, MPI_SUM, Communicator::communicator())/
               std::max(1.0e-30, allReduce(double(nodeIDs.size()), MPI_SUM, Communicator::communicator())));  // Negative implies outflow
    // cerr << "Computed inflow velocity: " << vinflow << endl;

    // Figure out a timestep limit such that we don't move more than the ghost
    // node thickness.
    Scalar xmin = 1e100, xmax = -1e100;
    for (const auto i: nodeIDs) {
      const auto xd = mPlane.signedDistance(pos[i]);
      xmin = std::min(xmin, xd);
      xmax = std::max(xmax, xd);
    }
    xmin = allReduce(xmin, MPI_MIN, Communicator::communicator());
    xmax = allReduce(xmax, MPI_MAX, Communicator::communicator());
    mXmin[nodeList.name()] = xmin;
    mDT = std::min(mDT, std::abs(xmax - xmin)/std::max(1e-30, std::abs(vinflow)));   // Protect from negative outflow velocity
    // cerr << "Timestep constraint: " << mDT << endl;

    mNumInflowNodes[nodeList.name()] = nodeIDs.size();
    CHECK(mNumInflowNodes.size() == mDataBase.numNodeLists());

    // Turn the BC on.
    mActive = true;
  }
}

//------------------------------------------------------------------------------
// Physics::evaluateDerivatives
//------------------------------------------------------------------------------
template<typename Dimension>
void
InflowOutflowBoundary<Dimension>::evaluateDerivatives(const Scalar time,
                                                      const Scalar dt,
                                                      const DataBase<Dimension>& dataBase,
                                                      const State<Dimension>& state,
                                                      StateDerivatives<Dimension>& derivatives) const {
}

//------------------------------------------------------------------------------
// Physics::dt
//------------------------------------------------------------------------------
template<typename Dimension>
typename InflowOutflowBoundary<Dimension>::TimeStepType
InflowOutflowBoundary<Dimension>::dt(const DataBase<Dimension>& dataBase, 
                                     const State<Dimension>& state,
                                     const StateDerivatives<Dimension>& derivs,
                                     const Scalar currentTime) const {
  return TimeStepType(mDT, "InflowOutflowBoundary velocity constraint");
}

//------------------------------------------------------------------------------
// Physics::registerState
//------------------------------------------------------------------------------
template<typename Dimension>
void
InflowOutflowBoundary<Dimension>::registerState(DataBase<Dimension>& dataBase,
                                                State<Dimension>& state) {
}

//------------------------------------------------------------------------------
// Physics::registerDerivatives
//------------------------------------------------------------------------------
template<typename Dimension>
void
InflowOutflowBoundary<Dimension>::registerDerivatives(DataBase<Dimension>& dataBase,
                                                      StateDerivatives<Dimension>& derivs) {
}

//------------------------------------------------------------------------------
// Physics::finalize
// At the end of a step, any ghost points that have wandered inside the entrance
// plane become internal points.
//------------------------------------------------------------------------------
template<typename Dimension>
void
InflowOutflowBoundary<Dimension>::finalize(const Scalar time, 
                                           const Scalar dt,
                                           DataBase<Dimension>& dataBase, 
                                           State<Dimension>& state,
                                           StateDerivatives<Dimension>& derivatives) {

  // First check every NodeList for any inflow or outflow nodes.
  bool altered = false;
  for (auto itr = dataBase.nodeListBegin(); itr < dataBase.nodeListEnd(); ++itr) {
    auto& nodeList = **itr;
    bool nodeListAltered = false;
    // cerr << "--------------------------------------------------------------------------------" << endl
    //      << nodeList.name() << endl;

    // Find any ghost points that are inside the entrance plane now.  These are inflow.
    const auto& gNodes = this->ghostNodes(nodeList);
    auto& pos = nodeList.positions();
    vector<int> insideNodes;
    for (auto i: gNodes) {
      if (mPlane.compare(pos[i]) == -1) insideNodes.push_back(i - gNodes[0]);
    }
    const auto numNew = insideNodes.size();
    if (numNew > 0) {
      nodeListAltered = true;

      // cerr << "Promoting to internal: ";
      // for (auto i: insideNodes) cerr << " : " << gNodes[0] + i << " " << pos[gNodes[0] + i];
      // cerr << endl;

      // Allocate new internal nodes for those we're promoting.
      const auto firstID = nodeList.numInternalNodes();
      nodeList.numInternalNodes(firstID + numNew);

      // Update the positions.
      vector<int> newNodes(numNew);
      for (auto k = 0; k < numNew; ++k) {
        newNodes[k] = firstID + k;
        pos[firstID + k] = pos[gNodes[0] + insideNodes[k] + numNew];
      }

      // Copy all field values from ghosts to the new internal nodes.
      copyFieldValues<Dimension, int>            (nodeList, mIntValues, insideNodes, newNodes);
      copyFieldValues<Dimension, Scalar>         (nodeList, mScalarValues, insideNodes, newNodes);
      copyFieldValues<Dimension, Vector>         (nodeList, mVectorValues, insideNodes, newNodes);
      copyFieldValues<Dimension, Tensor>         (nodeList, mTensorValues, insideNodes, newNodes);
      copyFieldValues<Dimension, SymTensor>      (nodeList, mSymTensorValues, insideNodes, newNodes);
      copyFieldValues<Dimension, ThirdRankTensor>(nodeList, mThirdRankTensorValues, insideNodes, newNodes);
      copyFieldValues<Dimension, FourthRankTensor>(nodeList, mFourthRankTensorValues, insideNodes, newNodes);
      copyFieldValues<Dimension, FifthRankTensor>(nodeList, mFifthRankTensorValues, insideNodes, newNodes);
      copyFieldValues<Dimension, FacetedVolume>  (nodeList, mFacetedVolumeValues, insideNodes, newNodes);
      copyFieldValues<Dimension, vector<Scalar>> (nodeList, mVectorScalarValues, insideNodes, newNodes);
      copyFieldValues<Dimension, vector<Vector>> (nodeList, mVectorVectorValues, insideNodes, newNodes);

      // for (auto k = 0; k < numNew; ++k) {
      //   cerr << " assigning position " << newNodes[k] << " @ " << pos[newNodes[k]] << " vel=" << nodeList.velocity()[newNodes[k]] << endl;
      // }
    }

    // Look for any internal points that have exited through the plane.
    vector<int> outsideNodes;
    const auto ni = nodeList.numInternalNodes();
    for (auto i = 0; i < ni; ++i) {
      if (mPlane.compare(pos[i]) == 1) outsideNodes.push_back(i);
    }

    if (not outsideNodes.empty()) {
      // cerr << "Deleting outflow nodes: ";
      // for (auto i: outsideNodes) cerr << " : " << i << " " << pos[i];
      // cerr << endl;

      // Those nodes get deleted.
      nodeListAltered = true;
      nodeList.deleteNodes(outsideNodes);
    }

    if (nodeListAltered) {
      altered = true;
      nodeList.neighbor().updateNodes();
    }
  }
  altered = (allReduce((altered ? 1 : 0), MPI_MAX, Communicator::communicator()) == 1);

  // If any NodeLists were altered, recompute the boundary conditions.
  if (altered) {
    // Remove any old ghost node information from the NodeLists.
    for (auto nodeListItr = dataBase.fluidNodeListBegin();
         nodeListItr != dataBase.fluidNodeListEnd(); 
         ++nodeListItr) {
      (*nodeListItr)->numGhostNodes(0);
      (*nodeListItr)->neighbor().updateNodes();
    }
    for (auto boundaryItr = this->boundaryBegin(); 
         boundaryItr < this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->reset(dataBase);

    // Iterate over the boundaries and set their ghost node info.
    for (auto boundaryItr = this->boundaryBegin(); 
         boundaryItr < this->boundaryEnd();
         ++boundaryItr) {
      (*boundaryItr)->setAllGhostNodes(dataBase);
      (*boundaryItr)->finalizeGhostBoundary();
      for (auto nodeListItr = dataBase.fluidNodeListBegin();
           nodeListItr < dataBase.fluidNodeListEnd(); 
           ++nodeListItr) {
        (*nodeListItr)->neighbor().updateNodes();
      }
    }
  }

}

//------------------------------------------------------------------------------
// Return the keys for the Fields we have stored.
//------------------------------------------------------------------------------
template<typename Dimension>
std::vector<std::string>
InflowOutflowBoundary<Dimension>::storedKeys() const {
  vector<string> result;
  for (const auto& pairs: mIntValues)             result.push_back(pairs.first);
  for (const auto& pairs: mScalarValues)          result.push_back(pairs.first);
  for (const auto& pairs: mVectorValues)          result.push_back(pairs.first);
  for (const auto& pairs: mTensorValues)          result.push_back(pairs.first);
  for (const auto& pairs: mSymTensorValues)       result.push_back(pairs.first);
  for (const auto& pairs: mThirdRankTensorValues) result.push_back(pairs.first);
  for (const auto& pairs: mFourthRankTensorValues) result.push_back(pairs.first);
  for (const auto& pairs: mFifthRankTensorValues) result.push_back(pairs.first);
  for (const auto& pairs: mFacetedVolumeValues)   result.push_back(pairs.first);
  for (const auto& pairs: mVectorScalarValues)    result.push_back(pairs.first);
  for (const auto& pairs: mVectorVectorValues)    result.push_back(pairs.first);
  return result;
}

//------------------------------------------------------------------------------
// Clear the stored values.
//------------------------------------------------------------------------------
template<typename Dimension>
void
InflowOutflowBoundary<Dimension>::clearStoredValues() {
  for (auto& stuff: mNumInflowNodes) stuff.second = 0;
  clearValues(mIntValues);
  clearValues(mScalarValues);
  clearValues(mVectorValues);
  clearValues(mTensorValues);
  clearValues(mSymTensorValues);
  clearValues(mThirdRankTensorValues);
  clearValues(mFacetedVolumeValues);
  clearValues(mVectorScalarValues);
  clearValues(mVectorVectorValues);
}

//------------------------------------------------------------------------------
// Return a unique label for restart.
//------------------------------------------------------------------------------
template<typename Dimension>
std::string
InflowOutflowBoundary<Dimension>::label() const {
  return "InflowOutflowBoundary" + boost::lexical_cast<std::string>(mBoundaryCount);
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
InflowOutflowBoundary<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  file.write(mActive, pathName + "/active");
  file.write(mBoundaryCount, pathName + "/boundaryCount");

  vector<std::string> keys;
  for (const auto& p: mIntValues) {
    keys.push_back(p.first);
    file.write(p.second, pathName + "/IntValues/" + p.first);
  }
  file.write(keys, pathName + "/IntValues/keys");

  keys.clear();
  for (const auto& p: mScalarValues) {
    keys.push_back(p.first);
    file.write(p.second, pathName + "/ScalarValues/" + p.first);
  }
  file.write(keys, pathName + "/ScalarValues/keys");

  keys.clear();
  for (const auto& p: mVectorValues) {
    keys.push_back(p.first);
    file.write(p.second, pathName + "/VectorValues/" + p.first);
  }
  file.write(keys, pathName + "/VectorValues/keys");

  keys.clear();
  for (const auto& p: mTensorValues) {
    keys.push_back(p.first);
    file.write(p.second, pathName + "/TensorValues/" + p.first);
  }
  file.write(keys, pathName + "/TensorValues/keys");

  keys.clear();
  for (const auto& p: mSymTensorValues) {
    keys.push_back(p.first);
    file.write(p.second, pathName + "/SymTensorValues/" + p.first);
  }
  file.write(keys, pathName + "/SymTensorValues/keys");

  keys.clear();
  for (const auto& p: mThirdRankTensorValues) {
    keys.push_back(p.first);
    file.write(p.second, pathName + "/ThirdRankTensorValues/" + p.first);
  }
  file.write(keys, pathName + "/ThirdRankTensorValues/keys");

  keys.clear();
  for (const auto& p: mFourthRankTensorValues) {
    keys.push_back(p.first);
    file.write(p.second, pathName + "/FourthRankTensorValues/" + p.first);
  }
  file.write(keys, pathName + "/FourthRankTensorValues/keys");

  keys.clear();
  for (const auto& p: mFifthRankTensorValues) {
    keys.push_back(p.first);
    file.write(p.second, pathName + "/FifthRankTensorValues/" + p.first);
  }
  file.write(keys, pathName + "/FifthRankTensorValues/keys");
  
  keys.clear();
  for (const auto& p: mFacetedVolumeValues) {
    keys.push_back(p.first);
    file.write(p.second, pathName + "/FacetedVolumeValues/" + p.first);
  }
  file.write(keys, pathName + "/FacetedVolumeValues/keys");

  keys.clear();
  for (const auto& p: mVectorScalarValues) {
    keys.push_back(p.first);
    file.write(p.second, pathName + "/VectorScalarValues/" + p.first);
  }
  file.write(keys, pathName + "/VectorScalarValues/keys");

  keys.clear();
  for (const auto& p: mVectorVectorValues) {
    keys.push_back(p.first);
    file.write(p.second, pathName + "/VectorVectorValues/" + p.first);
  }
  file.write(keys, pathName + "/VectorVectorValues/keys");
}

//------------------------------------------------------------------------------
// Read the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
InflowOutflowBoundary<Dimension>::
restoreState(const FileIO& file, const string& pathName)  {
  file.read(mActive, pathName + "/active");
  file.read(mBoundaryCount, pathName + "/boundaryCount");

  vector<std::string> keys;
  file.read(keys, pathName + "/IntValues/keys");
  mIntValues.clear();
  for (const auto key: keys) {
    mIntValues[key] = std::vector<int>();
    file.read(mIntValues[key], pathName + "/IntValues/" + key);
  }

  keys.clear();
  file.read(keys, pathName + "/ScalarValues/keys");
  mScalarValues.clear();
  for (const auto key: keys) {
    mScalarValues[key] = std::vector<Scalar>();
    file.read(mScalarValues[key], pathName + "/ScalarValues/" + key);
  }

  keys.clear();
  file.read(keys, pathName + "/VectorValues/keys");
  mVectorValues.clear();
  for (const auto key: keys) {
    mVectorValues[key] = std::vector<Vector>();
    file.read(mVectorValues[key], pathName + "/VectorValues/" + key);
  }

  keys.clear();
  file.read(keys, pathName + "/TensorValues/keys");
  mTensorValues.clear();
  for (const auto key: keys) {
    mTensorValues[key] = std::vector<Tensor>();
    file.read(mTensorValues[key], pathName + "/TensorValues/" + key);
  }

  keys.clear();
  file.read(keys, pathName + "/SymTensorValues/keys");
  mSymTensorValues.clear();
  for (const auto key: keys) {
    mSymTensorValues[key] = std::vector<SymTensor>();
    file.read(mSymTensorValues[key], pathName + "/SymTensorValues/" + key);
  }

  keys.clear();
  file.read(keys, pathName + "/ThirdRankTensorValues/keys");
  mThirdRankTensorValues.clear();
  for (const auto key: keys) {
    mThirdRankTensorValues[key] = std::vector<ThirdRankTensor>();
    file.read(mThirdRankTensorValues[key], pathName + "/ThirdRankTensorValues/" + key);
  }

  keys.clear();
  file.read(keys, pathName + "/FourthRankTensorValues/keys");
  mFourthRankTensorValues.clear();
  for (const auto key: keys) {
    mFourthRankTensorValues[key] = std::vector<FourthRankTensor>();
    file.read(mFourthRankTensorValues[key], pathName + "/FourthRankTensorValues/" + key);
  }

  keys.clear();
  file.read(keys, pathName + "/FifthRankTensorValues/keys");
  mFifthRankTensorValues.clear();
  for (const auto key: keys) {
    mFifthRankTensorValues[key] = std::vector<FifthRankTensor>();
    file.read(mFifthRankTensorValues[key], pathName + "/FifthRankTensorValues/" + key);
  }
  
  keys.clear();
  file.read(keys, pathName + "/FacetedVolumeValues/keys");
  mFacetedVolumeValues.clear();
  for (const auto key: keys) {
    mFacetedVolumeValues[key] = std::vector<FacetedVolume>();
    file.read(mFacetedVolumeValues[key], pathName + "/FacetedVolumeValues/" + key);
  }

  keys.clear();
  file.read(keys, pathName + "/VectorScalarValues/keys");
  mVectorScalarValues.clear();
  for (const auto key: keys) {
    mVectorScalarValues[key] = std::vector<std::vector<Scalar> >();
    file.read(mVectorScalarValues[key], pathName + "/VectorScalarValues/" + key);
  }

  keys.clear();
  file.read(keys, pathName + "/VectorVectorValues/keys");
  mVectorVectorValues.clear();
  for (const auto key: keys) {
    mVectorVectorValues[key] = std::vector<std::vector<Vector> >();
    file.read(mVectorVectorValues[key], pathName + "/VectorVectorValues/" + key);
  }
}

}
