#include "Geometry/GeomPlane.hh"
#include "NodeList/NodeListRegistrar.hh"
#include "Utilities/DataTypeTraits.hh"
#include "Field/FieldList.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Write a FieldList of arbitrary DataType.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FileIO::write(const FieldList<Dimension, DataType>& fieldList,
              const std::string pathName) {

  const std::string divider = "|";

  // Is the FieldList responsible for it's own memory?  If so, we have to 
  // provide additional information so it can properly restore itself.
  if (fieldList.storageType() == FieldStorageType::CopyFields) {
    if (fieldList.numFields() > 0) {
      std::stringstream names;
      for (typename FieldList<Dimension, DataType>::const_iterator fieldItr = fieldList.begin();
           fieldItr != fieldList.end();
           ++fieldItr) {
        names << (**fieldItr).nodeList().name() << divider;
      }
      names << std::ends;
      write(names.str(), pathName + "/NodeListNames");
    } else {
      write("", pathName + "/NodeListNames");
    }
  }

  // Loop over each Field, and write each one using the descendent method.
  for (typename FieldList<Dimension, DataType>::const_iterator fieldItr = fieldList.begin();
       fieldItr != fieldList.end();
       ++fieldItr) {

    // Build a unique path name for the Field.
    std::stringstream varPath;
    int elementID = std::distance(fieldList.begin(), fieldItr);
    varPath << pathName << "/Field" << elementID;

    write(**fieldItr, varPath.str());
  }
}

//------------------------------------------------------------------------------
// Read a FieldList of arbitrary DataType.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FileIO::read(FieldList<Dimension, DataType>& fieldList,
             const std::string pathName) const {

  const std::string divider = "|";

  // Is the FieldList responsible for it's own memory?  If so, we have to 
  // first make sure it has memory for each of the NodeLists it's defined against.
  if (fieldList.storageType() == FieldStorageType::CopyFields) {
    // We need the NodeListRegistrar.
    const NodeListRegistrar<Dimension>& registrar = NodeListRegistrar<Dimension>::instance();
    const size_t numNodeLists = registrar.numNodeLists();
    const std::vector<std::string> registeredNames = registrar.registeredNames();

    // Read the set of NodeLists this FieldList is associated with.
    std::string names;
    read(names, pathName + "/NodeListNames");
    while (names.find(divider) != std::string::npos) {

      // Read the name of the next NodeList.
      const size_t len = names.find(divider);
      CHECK(len < names.size());
      const std::string name = names.substr(0, len);
      
      // Find the NodeList this is.
      const size_t nodeListi = std::distance(registeredNames.begin(), 
                                             find(registeredNames.begin(), registeredNames.end(), name));
      VERIFY(nodeListi < numNodeLists);
      const NodeList<Dimension>& nodeList = **(registrar.begin() + nodeListi);

      // If necessary, insert a field into this FieldList for this NodeList.
      if (not fieldList.haveNodeList(nodeList)) {
        fieldList.appendNewField("Unnamed Field", nodeList, DataTypeTraits<DataType>::zero());
      }

      // Remove this name from the list of potentials.
      names = names.substr(len + 1, names.size());
    }
  }

  // Loop over each Field, and read each one using the descendent method.
  for (typename FieldList<Dimension, DataType>::iterator fieldItr = fieldList.begin();
       fieldItr != fieldList.end();
       ++fieldItr) {

    // Build the path name for the individual Field.
    std::stringstream varPath;
    int elementID = std::distance(fieldList.begin(), fieldItr);
    varPath << pathName << "/Field" << elementID;

    read(**fieldItr, varPath.str());
  }
}

//------------------------------------------------------------------------------
// Write a Field of std::vector<DataType>.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FileIO::write(const Field<Dimension, std::vector<DataType> >& field,
              const std::string pathName) {

  // Build an array with the number of elements per node, and count the total number of elements.
  std::vector<int> numElementsPerNode;
  size_t totalNumElements = size_t(0);
  for (typename Field<Dimension, std::vector<DataType> >::const_iterator itr = field.internalBegin();
       itr != field.internalEnd();
       ++itr) {
    const int ni = (*itr).size();
    totalNumElements += ni;
    numElementsPerNode.push_back(ni);
  }
  CHECK(numElementsPerNode.size() == field.nodeList().numInternalNodes());

  // Serialize the elements into a flat array.
  std::vector<DataType> elements(totalNumElements);
  size_t offset = size_t(0);
  for (typename Field<Dimension, std::vector<DataType> >::const_iterator itr = field.internalBegin();
       itr != field.internalEnd();
       ++itr) {
    std::copy(itr->begin(), itr->end(), elements.begin() + offset);
    offset += itr->size();
  }
  CHECK(offset == totalNumElements);

  // Now we can use the available methods for writing vector<>'s to store
  // the data.
  write(numElementsPerNode, pathName + "/numElementsPerNode");
  write(elements, pathName + "/elements");
}

//------------------------------------------------------------------------------
// Read a Field of std::vector<DataType>.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FileIO::read(Field<Dimension, std::vector<DataType> >& field,
             const std::string pathName) const {

  // Read the serialized data back in.
  std::vector<int> numElementsPerNode;
  std::vector<DataType> elements;
  read(numElementsPerNode, pathName + "/numElementsPerNode");
  read(elements, pathName + "/elements");
  CHECK(numElementsPerNode.size() == field.nodeList().numInternalNodes());

  // Consistency check...
  BEGIN_CONTRACT_SCOPE
  {
    int numElements = 0;
    for (std::vector<int>::const_iterator itr = numElementsPerNode.begin();
         itr != numElementsPerNode.end();
         ++itr) numElements += *itr;
    CHECK(numElements == elements.size());
  }
  END_CONTRACT_SCOPE

  // Fill the Field back in.
  typename std::vector<DataType>::const_iterator elementItr = elements.begin();
  for (int i = 0; i != field.nodeList().numInternalNodes(); ++i) {
    field(i) = std::vector<DataType>();
    field(i).reserve(numElementsPerNode[i]);
    for (int j = 0; j != numElementsPerNode[i]; ++j) {
      CHECK(elementItr < elements.end());
      field(i).push_back(*elementItr);
      ++elementItr;
    }
  }
  CHECK(elementItr == elements.end());
}

//------------------------------------------------------------------------------
// Write a std::vector<DataType>.
//------------------------------------------------------------------------------
template<typename DataType>
inline
void
FileIO::write(const std::vector<DataType>& x, const std::string pathName) {
  const int numElements = x.size();
  write(numElements, pathName + "/numElements");
  for (int i = 0; i != numElements; ++i) {
    std::stringstream elementPathName;
    elementPathName << pathName << "/" << i << std::ends;
    write(x[i], elementPathName.str());
  }
}

//------------------------------------------------------------------------------
// Read a std::vector<DataType>.
//------------------------------------------------------------------------------
template<typename DataType>
inline
void
FileIO::read(std::vector<DataType>& x, const std::string pathName) const {
  int numElements;
  read(numElements, pathName + "/numElements");
  x.resize(numElements);
  CHECK(x.size() == numElements);
  for (int i = 0; i != numElements; ++i) {
    std::stringstream elementPathName;
    elementPathName << pathName << "/" << i << std::ends;
    read(x[i], elementPathName.str());
  }
}

}
