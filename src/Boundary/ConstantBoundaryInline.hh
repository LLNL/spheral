#include <typeinfo>

namespace Spheral {
namespace BoundarySpace {

//------------------------------------------------------------------------------
// Return the number of nodes this Boundary is going to create.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
ConstantBoundary<Dimension>::
numConstantNodes() const {
  return mNumConstantNodes;
}

//------------------------------------------------------------------------------
// Return the NodeList this Boundary is defined on.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const NodeSpace::NodeList<Dimension>&
ConstantBoundary<Dimension>::
nodeList() const {
  CHECK(mNodeListPtr != 0);
  return *mNodeListPtr;
}

//------------------------------------------------------------------------------
// Store the Field values of the given DataType for the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
inline
void
ConstantBoundary<Dimension>::
storeFieldValues(const NodeSpace::NodeList<Dimension>& nodeList,
                 const std::vector<int>& nodeIDs,
                 std::map<const FieldSpace::FieldBase<Dimension>*, std::vector<DataType> >& values) const {

  // Iterate over all the Fields defined on the NodeList.
  for (typename NodeSpace::NodeList<Dimension>::const_FieldBaseIterator
         fieldItr = nodeList.registeredFieldsBegin();
       fieldItr != nodeList.registeredFieldsEnd();
       ++fieldItr) {

    // Determine if this Field is the type we're looking for.
    if (typeid(**fieldItr) == typeid(FieldSpace::Field<Dimension, DataType>)) {
      const FieldSpace::Field<Dimension, DataType>& field = (const FieldSpace::Field<Dimension, DataType>&) **fieldItr;

      // Build a vector of the values of this field on the requested nodes.
      std::vector<DataType> vals;
      vals.reserve(nodeIDs.size());
      for (typename std::vector<int>::const_iterator nodeItr = nodeIDs.begin();
           nodeItr != nodeIDs.end();
           ++nodeItr) {
        CHECK(*nodeItr >= 0 && *nodeItr < field.numElements());
        vals.push_back(field(*nodeItr));
      }
      CHECK(vals.size() == mNumConstantNodes);

      // Put this data into the result.
      CHECK(values.find(*fieldItr) == values.end());
      values[*fieldItr] = vals;
    }
  }
}
  
//------------------------------------------------------------------------------
// Set the ghost values in the given field using the given map.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
inline
void
ConstantBoundary<Dimension>::
setGhostValues(FieldSpace::Field<Dimension, DataType>& field,
               const std::map<const FieldSpace::FieldBase<Dimension>*, 
                              std::vector<DataType> >& values) const {

  const NodeSpace::NodeList<Dimension>& nodeList = field.nodeList();
  if (&nodeList == mNodeListPtr) {

    // Get the set of ghost nodes.
    const std::vector<int>& ghostNodes = this->ghostNodes(nodeList);
    CHECK(ghostNodes.size() == mNumConstantNodes);

    // Find this Field in the set of stored values.
    const FieldSpace::FieldBase<Dimension>* fieldBasePtr = (const FieldSpace::FieldBase<Dimension>*) &field;
    typename std::map<const FieldSpace::FieldBase<Dimension>*, std::vector<DataType> >::const_iterator
      itr = values.find(fieldBasePtr);
//     CHECK(itr != values.end());
//     if (itr == values.end()) {
//       std::string message = "ConstantBoundary::setGhostValues: Attempt to set boundary values for unknown field.";
//       cerr << message << endl;
//       throw message;
//     }

    // Now set the ghost values.
    if (itr != values.end()) {
      const std::vector<DataType>& ghostValues = itr->second;;
      CHECK(ghostValues.size() == mNumConstantNodes);
      for (int i = 0; i < mNumConstantNodes; ++i) {
        CHECK(ghostNodes[i] >= nodeList.firstGhostNode() &&
              ghostNodes[i] < nodeList.numNodes());
        field(ghostNodes[i]) = ghostValues[i];
      }
    }

  }
}

}
}
