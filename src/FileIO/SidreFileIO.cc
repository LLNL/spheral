//---------------------------------Spheral++----------------------------------//
// SidreFileIO -- Provide the interface to sidre file objects.
//
// Created by Mikhail Zakharchanka, 11/4/2021
//----------------------------------------------------------------------------//
#include "SidreFileIO.hh"
#include "axom/core/Path.hpp"

namespace Spheral
{

//------------------------------------------------------------------------------
// Specialize to read/write a std::string/std::vector type
//------------------------------------------------------------------------------
template <typename DataType>
void sidreWriteVec(axom::sidre::Group* grp, const DataType& value, const std::string& path)
{
  axom::sidre::DataTypeId dtype = DataTypeTraits<DataType>::axomTypeID();
  axom::sidre::Buffer* buff = grp->getDataStore()->createBuffer()
                                           ->allocate(dtype, value.size())
                                           ->copyBytesIntoBuffer((void*)value.data(), sizeof(typename DataType::value_type) * (value.size()));
  grp->createView(path, dtype, value.size(), buff);
}

template <typename DataType>
void sidreReadVec(axom::sidre::Group* grp, DataType& value, const std::string& path)
{
  int size = grp->getView(path)->getNumElements();
  value.resize(size);
  typename DataType::value_type* data = grp->getView(path)->getArray();
  value.assign(data, data + size);
}

//------------------------------------------------------------------------------
// Specialize to read/write a Geometric type
//------------------------------------------------------------------------------
template <typename DataType>
void sidreWriteGeom(axom::sidre::Group* grp, const DataType& value, const std::string& path)
{
  int size = int(DataType::numElements);
  axom::sidre::Buffer* buff = grp->getDataStore()->createBuffer()
                                          ->allocate(axom::sidre::DOUBLE_ID, size)
                                          ->copyBytesIntoBuffer((void*)value.begin(), sizeof(double) * (size));
  grp->createView(path, axom::sidre::DOUBLE_ID, size, buff);
}

template <typename DataType>
void sidreReadGeom(axom::sidre::Group* grp, DataType& value, const std::string& path)
{
  double* data = grp->getView(path)->getArray();
  for (int i = 0; i < int(DataType::numElements); ++i)
    value[i] = data[i];
}
//------------------------------------------------------------------------------
// Specialize to read/write a Field<scalar>
//------------------------------------------------------------------------------
template <typename Dimension, typename DataType,
          typename std::enable_if<std::is_arithmetic<DataType>::value,
                                  DataType>::type* = nullptr>
void sidreWriteField(axom::sidre::Group* grp,
                     const Spheral::Field<Dimension, DataType>& field,
                     const std::string& path)
{
  int size = field.numInternalElements();
  axom::sidre::DataTypeId dtype = field.getAxomTypeID();
  axom::sidre::Buffer* buff = grp->getDataStore()->createBuffer()
                                          ->allocate(dtype, size)
                                          ->copyBytesIntoBuffer((void*)&field[0], sizeof(DataType) * (size));
  grp->createView(path, dtype, size, buff);
}

template <typename Dimension, typename DataType,
          typename std::enable_if<std::is_arithmetic<DataType>::value,
                                  DataType>::type* = nullptr>
void sidreReadField(axom::sidre::Group* grp,
                     Spheral::Field<Dimension, DataType>& field,
                     const std::string& path)
{
  DataType* data = grp->getView(path)->getArray();
  for (int i = 0; i < grp->getView(path)->getNumElements(); ++i)
    field[i] = data[i];
}

//------------------------------------------------------------------------------
// Specialize to read/write a Field<Vector, Tensor, SymTensor>
//------------------------------------------------------------------------------
template <typename Dimension, typename DataType,
          typename std::enable_if<!std::is_arithmetic<DataType>::value,
                                  DataType>::type* = nullptr>
void sidreWriteField(axom::sidre::Group* grp,
                     const Spheral::Field<Dimension, DataType>& field,
                     const std::string& path)
{
  int size = field.numInternalElements();
  axom::sidre::DataTypeId dtype = field.getAxomTypeID();
  std::vector<double> fieldData(size * DataType::numElements);

  for (int i = 0; i < size; ++i)
    std::copy(field(i).begin(), field(i).end(), &fieldData[i * DataType::numElements]);

  axom::sidre::Buffer* buff = grp->getDataStore()->createBuffer()
                                          ->allocate(dtype, size * DataType::numElements)
                                          ->copyBytesIntoBuffer((void*)fieldData.data(), sizeof(double) * (size * DataType::numElements));
  grp->createView(path, dtype, size * DataType::numElements, buff);
}

template <typename Dimension, typename DataType,
          typename std::enable_if<!std::is_arithmetic<DataType>::value,
                                  DataType>::type* = nullptr>
void sidreReadField(axom::sidre::Group* grp,
                     Spheral::Field<Dimension, DataType>& field,
                     const std::string& path)
{
  double* data = grp->getView(path)->getArray();

  for (int i = 0; i < (grp->getView(path)->getNumElements() / (DataType::numElements)); ++i)
    for (int j = 0; j < int(DataType::numElements); ++j)
      field[i][j] = data[(i * (DataType::numElements)) + j];
}

//------------------------------------------------------------------------------
// Empty constructor.
//------------------------------------------------------------------------------
SidreFileIO::SidreFileIO():
    FileIO()
  {
    // mDataStorePtr = std::make_unique<axom::sidre::DataStore>();
  }

//------------------------------------------------------------------------------
// Construct and open the given file.
//------------------------------------------------------------------------------
SidreFileIO::SidreFileIO(const std::string fileName, AccessType access):
  FileIO(fileName, access)
{
  // mDataStorePtr = std::make_unique<axom::sidre::DataStore>();
  open(fileName, access);
  ENSURE(mFileOpen && mDataStorePtr != 0);
}

//------------------------------------------------------------------------------
// Construct and open the given file and set number of restart files
//------------------------------------------------------------------------------
SidreFileIO::SidreFileIO(const std::string fileName, AccessType access, int numFiles):
  FileIO(fileName, access)
{
  // mDataStorePtr = std::make_unique<axom::sidre::DataStore>();
  open(fileName, access);
  ENSURE(mFileOpen && mDataStorePtr != 0);
  if (numFiles > 0 && numFiles <= Process::getTotalNumberOfProcesses())
    numRestartFiles = numFiles;
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
SidreFileIO::~SidreFileIO()
{
  close();
}

//------------------------------------------------------------------------------
// Open a SiloFile file with the specified access.
//------------------------------------------------------------------------------
void SidreFileIO::open(const std::string fileName, AccessType access)
{
  VERIFY2(mDataStorePtr == 0 and mFileOpen == false,
          "ERROR: attempt to reopen SidreFileIO object.");

  mDataStorePtr = std::make_unique<axom::sidre::DataStore>();
  baseGroup = mDataStorePtr->getRoot();
  mFileName = fileName;

  if (access == AccessType::Read)
  {
#ifdef USE_MPI
    axom::sidre::IOManager reader(Communicator::communicator());
    reader.read(baseGroup, fileName + ".root");
#else
    baseGroup->load(fileName);
#endif // USE_MPI
  }

  VERIFY2(mDataStorePtr != 0, "SidreFileIO ERROR: unable to open " << fileName);
  mFileOpen = true;
}

//------------------------------------------------------------------------------
// Close the current file.
//------------------------------------------------------------------------------
void SidreFileIO::close()
{
  if (mDataStorePtr != 0)
  {
#ifdef USE_MPI
    axom::sidre::IOManager writer(Communicator::communicator());
    writer.write(baseGroup, numRestartFiles, mFileName, "sidre_hdf5");
#else
    baseGroup->save(mFileName);
#endif // USE_MPI
    mDataStorePtr.reset();
  }
  mFileOpen = false;
}

// void setGroup(axom::sidre::Group* group)
// {
//   baseGroup = group;
// }

//------------------------------------------------------------------------------
// Check if the specified path is in the file.
//------------------------------------------------------------------------------
bool SidreFileIO::pathExists(const std::string pathName) const
{
  return baseGroup->hasView(pathName);
}

//------------------------------------------------------------------------------
// Write an unsigned to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const unsigned& value, const std::string pathName)
{
  baseGroup->createViewScalar(pathName, value);
}

//------------------------------------------------------------------------------
// Write a size_t to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const size_t& value, const std::string pathName)
{
  baseGroup->createViewScalar(pathName, value);
}

//------------------------------------------------------------------------------
// Write an int to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const int& value, const std::string pathName)
{
  baseGroup->createViewScalar(pathName, value);
}

// ------------------------------------------------------------------------------
// Write a bool to the file.
// ------------------------------------------------------------------------------
void SidreFileIO::write(const bool& value, const std::string pathName)
{
  baseGroup->createViewScalar(pathName, (int)value);
}

//------------------------------------------------------------------------------
// Write a double to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const double& value, const std::string pathName)
{
  baseGroup->createViewScalar(pathName, value);
}

//------------------------------------------------------------------------------
// Write a std::string to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const std::string& value, const std::string pathName)
{
  // baseGroup->createViewString(pathName, value);
  sidreWriteVec(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Write a vector<int> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const std::vector<int>& value, const std::string pathName)
{
  sidreWriteVec(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Write a vector<double> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const std::vector<double>& value, const std::string pathName)
{
  sidreWriteVec(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Write a vector<std::string> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const std::vector<std::string>& value, const std::string pathName)
{
  axom::Path myPath = axom::Path(pathName);
  axom::sidre::Group* wholeField = baseGroup->createGroup(myPath.dirName());

  for (u_int i = 0; i < value.size(); ++i)
    wholeField->createViewString(myPath.baseName() + std::to_string(i), value[i]);
}

//------------------------------------------------------------------------------
// Write a Dim<1>::Vector to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<1>::Vector& value, const std::string pathName)
{
  sidreWriteGeom(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<1>::Tensor to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<1>::Tensor& value, const std::string pathName)
{
  sidreWriteGeom(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<1>::SymTensor to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<1>::SymTensor& value, const std::string pathName)
{
  sidreWriteGeom(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<1>::ThirdRankTensor to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<1>::ThirdRankTensor& value, const std::string pathName)
{
  sidreWriteGeom(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<2>::Vector to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<2>::Vector& value, const std::string pathName)
{
  sidreWriteGeom(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<2>::Tensor to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<2>::Tensor& value, const std::string pathName)
{
  sidreWriteGeom(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<2>::SymTensor to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<2>::SymTensor& value, const std::string pathName)
{
  sidreWriteGeom(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<2>::ThirdRankTensor to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<2>::ThirdRankTensor& value, const std::string pathName)
{
  sidreWriteGeom(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<3>::Vector to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<3>::Vector& value, const std::string pathName)
{
  sidreWriteGeom(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<3>::Tensor to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<3>::Tensor& value, const std::string pathName)
{
  sidreWriteGeom(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<3>::SymTensor to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<3>::SymTensor& value, const std::string pathName)
{
  sidreWriteGeom(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<3>::ThirdRankTensor to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<3>::ThirdRankTensor& value, const std::string pathName)
{
  sidreWriteGeom(baseGroup, value, pathName);
}


// ------------------------------------------------------------------------------
// Read an unsigned from the file.
// ------------------------------------------------------------------------------
void SidreFileIO::read(unsigned& value, const std::string pathName) const
{
  value = baseGroup->getView(pathName)->getScalar();
}

//------------------------------------------------------------------------------
// Read a size_t from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(size_t& value, const std::string pathName) const
{
  value = baseGroup->getView(pathName)->getScalar();
}

//------------------------------------------------------------------------------
// Read an int from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(int& value, const std::string pathName) const
{
  value = baseGroup->getView(pathName)->getScalar();
}

//------------------------------------------------------------------------------
// Read a bool from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(bool& value, const std::string pathName) const
{
  value = (int)baseGroup->getView(pathName)->getScalar();
}

//------------------------------------------------------------------------------
// Read a double from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(double& value, const std::string pathName) const
{
  value = baseGroup->getView(pathName)->getScalar();
}

//------------------------------------------------------------------------------
// Read a std::string from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(std::string& value, const std::string pathName) const
{
  // value = baseGroup->getView(pathName)->getString();
  sidreReadVec(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a vector<int> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(std::vector<int>& value, const std::string pathName) const
{
  sidreReadVec(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a vector<double> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(std::vector<double>& value, const std::string pathName) const
{
  sidreReadVec(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a vector<std::string> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(vector<std::string>& value, const std::string pathName) const
{
  axom::Path myPath = axom::Path(pathName);
  axom::sidre::Group* group = baseGroup->getGroup(myPath.dirName());
  int size = group->getNumViews();
  value.resize(size);

  for (int i = 0; i < size; ++i)
    value[i] = group->getView(myPath.baseName() + std::to_string(i))->getString();
}

//------------------------------------------------------------------------------
// Read a Dim<1>::Vector from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Dim<1>::Vector& value, const std::string pathName) const
{
  sidreReadGeom(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<1>::Tensor from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Dim<1>::Tensor& value, const std::string pathName) const
{
  sidreReadGeom(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<1>::SymTensor from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Dim<1>::SymTensor& value, const std::string pathName) const
{
  sidreReadGeom(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<1>::ThirdRankTensor from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Dim<1>::ThirdRankTensor& value, const std::string pathName) const
{
  sidreReadGeom(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<2>::Vector from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Dim<2>::Vector& value, const std::string pathName) const
{
  sidreReadGeom(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<2>::Tensor from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Dim<2>::Tensor& value, const std::string pathName) const
{
  sidreReadGeom(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<2>::SymTensor from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Dim<2>::SymTensor& value, const std::string pathName) const
{
  sidreReadGeom(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<2>::ThirdRankTensor from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Dim<2>::ThirdRankTensor& value, const std::string pathName) const
{
  sidreReadGeom(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<3>::Vector from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Dim<3>::Vector& value, const std::string pathName) const
{
  sidreReadGeom(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<3>::Tensor from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Dim<3>::Tensor& value, const std::string pathName) const
{
  sidreReadGeom(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<3>::SymTensor from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Dim<3>::SymTensor& value, const std::string pathName) const
{
  sidreReadGeom(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<3>::ThirdRankTensor from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Dim<3>::ThirdRankTensor& value, const std::string pathName) const
{
  sidreReadGeom(baseGroup, value, pathName);
}

#ifdef SPHERAL1D
//------------------------------------------------------------------------------
// Write a Field<Dim<1>, Dim<1>::Scalar> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<1>, Dim<1>::Scalar>& value, const std::string pathName)
{
  sidreWriteField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, Dim<1>::Vector> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<1>, Dim<1>::Vector>& value, const std::string pathName)
{
  sidreWriteField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, Dim<1>::Tensor> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<1>, Dim<1>::Tensor>& value, const std::string pathName)
{
  sidreWriteField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, Dim<1>::SymTensor> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<1>, Dim<1>::SymTensor>& value, const std::string pathName)
{
  sidreWriteField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, Dim<1>::ThirdRankTensor> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<1>, Dim<1>::ThirdRankTensor>& value, const std::string pathName)
{
  sidreWriteField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, int> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<1>, int>& value, const std::string pathName)
{
  sidreWriteField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, unsigned> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<1>, unsigned>& value, const std::string pathName)
{
  sidreWriteField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, Dim<1>::Scalar> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<1>, Dim<1>::Scalar>& value, const std::string pathName) const
{
  sidreReadField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, Dim<1>::Vector> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<1>, Dim<1>::Vector>& value, const std::string pathName) const
{
  sidreReadField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, Dim<1>::Tensor> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<1>, Dim<1>::Tensor>& value, const std::string pathName) const
{
  sidreReadField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, Dim<1>::SymTensor> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<1>, Dim<1>::SymTensor>& value, const std::string pathName) const
{
  sidreReadField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, Dim<1>::ThirdRankTensor> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<1>, Dim<1>::ThirdRankTensor>& value, const std::string pathName) const
{
  sidreReadField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, int> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<1>, int>& value, const std::string pathName) const
{
  sidreReadField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, unsigned> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<1>, unsigned>& value, const std::string pathName) const
{
  sidreReadField(baseGroup, value, pathName);
}
#endif

#ifdef SPHERAL2D
//------------------------------------------------------------------------------
// Write a Field<Dim<2>, Dim<2>::Scalar> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<2>, Dim<2>::Scalar>& value, const std::string pathName)
{
  sidreWriteField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, Dim<2>::Vector> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<2>, Dim<2>::Vector>& value, const std::string pathName)
{
  sidreWriteField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, Dim<2>::Tensor> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<2>, Dim<2>::Tensor>& value, const std::string pathName)
{
  sidreWriteField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, Dim<2>::SymTensor> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<2>, Dim<2>::SymTensor>& value, const std::string pathName)
{
  sidreWriteField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, Dim<2>::ThirdRankTensor> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<2>, Dim<2>::ThirdRankTensor>& value, const std::string pathName)
{
  sidreWriteField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, int> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<2>, int>& value, const std::string pathName)
{
  sidreWriteField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, unsigned> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<2>, unsigned>& value, const std::string pathName)
{
  sidreWriteField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, Dim<2>::Scalar> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<2>, Dim<2>::Scalar>& value, const std::string pathName) const
{
  sidreReadField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, Dim<2>::Vector> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<2>, Dim<2>::Vector>& value, const std::string pathName) const
{
  sidreReadField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, Dim<2>::Tensor> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<2>, Dim<2>::Tensor>& value, const std::string pathName) const
{
  sidreReadField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, Dim<2>::SymTensor> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<2>, Dim<2>::SymTensor>& value, const std::string pathName) const
{
  sidreReadField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, Dim<2>::ThirdRankTensor> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<2>, Dim<2>::ThirdRankTensor>& value, const std::string pathName) const
{
  sidreReadField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, int> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<2>, int>& value, const std::string pathName) const
{
  sidreReadField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, unsigned> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<2>, unsigned>& value, const std::string pathName) const
{
  sidreReadField(baseGroup, value, pathName);
}
#endif

#ifdef SPHERAL3D
//------------------------------------------------------------------------------
// Write a Field<Dim<3>, Dim<3>::Scalar> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<3>, Dim<3>::Scalar>& value, const std::string pathName)
{
  sidreWriteField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, Dim<3>::Vector> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<3>, Dim<3>::Vector>& value, const std::string pathName)
{
  sidreWriteField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, Dim<3>::Tensor> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<3>, Dim<3>::Tensor>& value, const std::string pathName)
{
  sidreWriteField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, Dim<3>::SymTensor> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<3>, Dim<3>::SymTensor>& value, const std::string pathName)
{
  sidreWriteField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, Dim<3>::ThirdRankTensor> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<3>, Dim<3>::ThirdRankTensor>& value, const std::string pathName)
{
  sidreWriteField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, int> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<3>, int>& value, const std::string pathName)
{
  sidreWriteField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, unsigned> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<3>, unsigned>& value, const std::string pathName)
{
  sidreWriteField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, Dim<3>::Scalar> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<3>, Dim<3>::Scalar>& value, const std::string pathName) const
{
  sidreReadField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, Dim<3>::Vector> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<3>, Dim<3>::Vector>& value, const std::string pathName) const
{
  sidreReadField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, Dim<3>::Tensor> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<3>, Dim<3>::Tensor>& value, const std::string pathName) const
{
  sidreReadField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, Dim<3>::SymTensor> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<3>, Dim<3>::SymTensor>& value, const std::string pathName) const
{
  sidreReadField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, Dim<3>::ThirdRankTensor> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<3>, Dim<3>::ThirdRankTensor>& value, const std::string pathName) const
{
  sidreReadField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, int> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<3>, int>& value, const std::string pathName) const
{
  sidreReadField(baseGroup, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, unsigned> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<3>, unsigned>& value, const std::string pathName) const
{
  sidreReadField(baseGroup, value, pathName);
}
#endif

}
