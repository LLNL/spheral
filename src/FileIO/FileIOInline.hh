#include "Geometry/GeomPlane.hh"
#include "NodeList/NodeListRegistrar.hh"
#include "Utilities/DataTypeTraits.hh"
#include "Field/FieldList.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Write Field as binary blob
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FileIO::writeFieldBlob(const Field<Dimension, DataType>& value,
                       const std::string pathName) {
  const auto n = value.nodeList().numInternalNodes();
  std::vector<int> inds(n);
  for (auto i = 0u; i < n; ++i) inds[i] = i;
  const auto blob = value.packValues(inds);
  const std::string sblob(blob.begin(), blob.end());
  this->write(sblob, pathName);
}

//------------------------------------------------------------------------------
// Read Field as binary blob
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FileIO::readFieldBlob(Field<Dimension, DataType>& value,
                      const std::string pathName) const {
  const auto n = value.nodeList().numInternalNodes();
  std::vector<int> inds(n);
  for (auto i = 0u; i < n; ++i) inds[i] = i;
  std::string sblob;
  this->read(sblob, pathName);
  std::vector<char> blob(sblob.begin(), sblob.end());
  value.unpackValues(inds, blob);
}

//------------------------------------------------------------------------------
// Write a FieldList of arbitrary DataType (private)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FileIO::writeFieldList(const FieldList<Dimension, DataType>& fieldList,
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
// Read a FieldList of arbitrary DataType. (private)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FileIO::readFieldList(FieldList<Dimension, DataType>& fieldList,
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
    // std::cerr << "  AFTER READING Field: " << (**fieldItr).size() << std::endl;
  }
}

//------------------------------------------------------------------------------
// Write a Field of std::vector<DataType>. (private)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FileIO::writeFieldVector(const Field<Dimension, std::vector<DataType> >& field,
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
// Read a Field of std::vector<DataType>. (private)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FileIO::readFieldVector(Field<Dimension, std::vector<DataType> >& field,
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
    CHECK(numElements == (int)elements.size());
  }
  END_CONTRACT_SCOPE

  // Fill the Field back in.
  typename std::vector<DataType>::const_iterator elementItr = elements.begin();
  for (auto i = 0u; i != field.nodeList().numInternalNodes(); ++i) {
    field(i) = std::vector<DataType>();
    field(i).reserve(numElementsPerNode[i]);
    for (auto j = 0; j != numElementsPerNode[i]; ++j) {
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
// Generic method
template<typename DataType>
inline
void
FileIO::write(const std::vector<DataType>& x, const std::string pathName) {
  std::vector<char> buf;
  packElement(x, buf);
  std::string bufstr(buf.begin(), buf.end());
  this->write(bufstr, pathName);
}

// Specialize for some types that can be treated as arrays of doubles more efficiently/portably
//..............................................................................
template<typename Value>
inline
void
FileIO::writeVector(const std::vector<Value>& x, const std::string pathName) {
  const auto n = x.size();
  const auto ne = Value::numElements;
  std::vector<double> buf(n*ne);
  for (auto i = 0u; i < n; ++i) std::copy(x[i].begin(), x[i].end(), &buf[i*ne]);
  this->write(buf, pathName);
}
template<> inline void FileIO::write<Dim<1>::Vector>         (const std::vector<Dim<1>::Vector>& x,          const std::string pathName) { this->writeVector(x, pathName); }
template<> inline void FileIO::write<Dim<1>::Tensor>         (const std::vector<Dim<1>::Tensor>& x,          const std::string pathName) { this->writeVector(x, pathName); }
template<> inline void FileIO::write<Dim<1>::SymTensor>      (const std::vector<Dim<1>::SymTensor>& x,       const std::string pathName) { this->writeVector(x, pathName); }
template<> inline void FileIO::write<Dim<1>::ThirdRankTensor>(const std::vector<Dim<1>::ThirdRankTensor>& x, const std::string pathName) { this->writeVector(x, pathName); }
template<> inline void FileIO::write<Dim<2>::Vector>         (const std::vector<Dim<2>::Vector>& x,          const std::string pathName) { this->writeVector(x, pathName); }
template<> inline void FileIO::write<Dim<2>::Tensor>         (const std::vector<Dim<2>::Tensor>& x,          const std::string pathName) { this->writeVector(x, pathName); }
template<> inline void FileIO::write<Dim<2>::SymTensor>      (const std::vector<Dim<2>::SymTensor>& x,       const std::string pathName) { this->writeVector(x, pathName); }
template<> inline void FileIO::write<Dim<2>::ThirdRankTensor>(const std::vector<Dim<2>::ThirdRankTensor>& x, const std::string pathName) { this->writeVector(x, pathName); }
template<> inline void FileIO::write<Dim<3>::Vector>         (const std::vector<Dim<3>::Vector>& x,          const std::string pathName) { this->writeVector(x, pathName); }
template<> inline void FileIO::write<Dim<3>::Tensor>         (const std::vector<Dim<3>::Tensor>& x,          const std::string pathName) { this->writeVector(x, pathName); }
template<> inline void FileIO::write<Dim<3>::SymTensor>      (const std::vector<Dim<3>::SymTensor>& x,       const std::string pathName) { this->writeVector(x, pathName); }
template<> inline void FileIO::write<Dim<3>::ThirdRankTensor>(const std::vector<Dim<3>::ThirdRankTensor>& x, const std::string pathName) { this->writeVector(x, pathName); }

//------------------------------------------------------------------------------
// Read a std::vector<DataType>.
//------------------------------------------------------------------------------
// Generic method
template<typename DataType>
inline
void
FileIO::read(std::vector<DataType>& x, const std::string pathName) const {
  std::string bufstr;
  this->read(bufstr, pathName);
  const std::vector<char> buf(bufstr.begin(), bufstr.end());
  auto itr = buf.begin();
  unpackElement(x, itr, buf.end());
  ENSURE(itr == buf.end());
}

// Specialize for some types that can be treated as arrays of doubles more efficiently/portably
//..............................................................................
template<typename Value>
inline
void
FileIO::readVector(std::vector<Value>& x, const std::string pathName) const {
  std::vector<double> buf;
  this->read(buf, pathName);
  const auto ne = Value::numElements;
  const auto n = buf.size()/ne;
  x.resize(n);
  for (auto i = 0u; i < n; ++i) std::copy(&buf[i*ne], &buf[i*ne] + ne, x[i].begin());
}
template<> inline void FileIO::read<Dim<1>::Vector>          (std::vector<Dim<1>::Vector>& x,          const std::string pathName) const { this->readVector(x, pathName); }
template<> inline void FileIO::read<Dim<1>::Tensor>          (std::vector<Dim<1>::Tensor>& x,          const std::string pathName) const { this->readVector(x, pathName); }
template<> inline void FileIO::read<Dim<1>::SymTensor>       (std::vector<Dim<1>::SymTensor>& x,       const std::string pathName) const { this->readVector(x, pathName); }
template<> inline void FileIO::read<Dim<1>::ThirdRankTensor> (std::vector<Dim<1>::ThirdRankTensor>& x, const std::string pathName) const { this->readVector(x, pathName); }
template<> inline void FileIO::read<Dim<2>::Vector>          (std::vector<Dim<2>::Vector>& x,          const std::string pathName) const { this->readVector(x, pathName); }
template<> inline void FileIO::read<Dim<2>::Tensor>          (std::vector<Dim<2>::Tensor>& x,          const std::string pathName) const { this->readVector(x, pathName); }
template<> inline void FileIO::read<Dim<2>::SymTensor>       (std::vector<Dim<2>::SymTensor>& x,       const std::string pathName) const { this->readVector(x, pathName); }
template<> inline void FileIO::read<Dim<2>::ThirdRankTensor> (std::vector<Dim<2>::ThirdRankTensor>& x, const std::string pathName) const { this->readVector(x, pathName); }
template<> inline void FileIO::read<Dim<3>::Vector>          (std::vector<Dim<3>::Vector>& x,          const std::string pathName) const { this->readVector(x, pathName); }
template<> inline void FileIO::read<Dim<3>::Tensor>          (std::vector<Dim<3>::Tensor>& x,          const std::string pathName) const { this->readVector(x, pathName); }
template<> inline void FileIO::read<Dim<3>::SymTensor>       (std::vector<Dim<3>::SymTensor>& x,       const std::string pathName) const { this->readVector(x, pathName); }
template<> inline void FileIO::read<Dim<3>::ThirdRankTensor> (std::vector<Dim<3>::ThirdRankTensor>& x, const std::string pathName) const { this->readVector(x, pathName); }

}
