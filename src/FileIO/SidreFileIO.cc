//---------------------------------Spheral++----------------------------------//
// SidreFileIO -- Provide the interface to sidre file objects.
//
// Created by Mikhail Zakharchanka, 11/4/2021
//----------------------------------------------------------------------------//
#include "SidreFileIO.hh"
#include "Utilities/SidreDataCollection.hh"
#include "axom/core/Path.hpp"
#include <iostream>

namespace Spheral
{

//------------------------------------------------------------------------------
// Empty constructor.
//------------------------------------------------------------------------------
SidreFileIO::SidreFileIO():
  FileIO(),
  mDataStorePtr(0) {}

//------------------------------------------------------------------------------
// Construct and open the given file.
//------------------------------------------------------------------------------
SidreFileIO::SidreFileIO(const std::string fileName, AccessType access):
  FileIO(fileName, access),
  mDataStorePtr(0)
{
  open(fileName, access);
  //ENSURE(mFileOpen && mDataStorePtr != 0);
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

  mDataStorePtr = std::make_shared<axom::sidre::DataStore>();

  // Need this member var because save() needs to know the name too.
  mFileName = fileName;

  if (access == AccessType::Read) 
    mDataStorePtr->getRoot()->load(fileName);
  
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
    mDataStorePtr->getRoot()->save(mFileName, "sidre_hdf5");
    mDataStorePtr = 0;
  }
  mFileOpen = false;
}

//------------------------------------------------------------------------------
// Check if the specified path is in the file.
//------------------------------------------------------------------------------
bool SidreFileIO::pathExists(const std::string pathName) const
{
  return mFileOpen;
}

// 

//------------------------------------------------------------------------------
// Write an unsigned to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const unsigned& value, const std::string pathName)
{
  mDataStorePtr->getRoot()->createViewScalar(pathName, value);
}

//------------------------------------------------------------------------------
// Write a size_t to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const size_t& value, const std::string pathName)
{
  mDataStorePtr->getRoot()->createViewScalar(pathName, value);
}

//------------------------------------------------------------------------------
// Write an int to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const int& value, const std::string pathName)
{
  mDataStorePtr->getRoot()->createViewScalar(pathName, value);
}

// ------------------------------------------------------------------------------
// Write a bool to the file.
// ------------------------------------------------------------------------------
void SidreFileIO::write(const bool& value, const std::string pathName)
{
  mDataStorePtr->getRoot()->createViewScalar(pathName, (int)value);
}

//------------------------------------------------------------------------------
// Write a double to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const double& value, const std::string pathName)
{
  mDataStorePtr->getRoot()->createViewScalar(pathName, value);
}

//------------------------------------------------------------------------------
// Write a std::string to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const std::string& value, const std::string pathName) 
{
  mDataStorePtr->getRoot()->createViewString(pathName, value);
}

//------------------------------------------------------------------------------
// Write a vector<int> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const std::vector<int>& value, const std::string pathName)
{
  // std::cout << "This is writing a vector of ints and value = ";
  // if (value.size() == 0)
  //   std::cout << "value is empty.\n";
  // for (auto it = value.begin(); it != value.end(); ++it)
  //   std::cout << *it << " ";
  // std::cout << "\n";

  mDataStorePtr->getRoot()->createView(pathName, axom::sidre::INT_ID, value.size(), (void*)(&(*value.begin())));
  // mDataStorePtr->getRoot()->print();
  // std::cout << "This is how many views root has:" << mDataStorePtr->getRoot()->getNumViews() << "\n";
  // if (mDataStorePtr->getRoot()->getView(pathName)->isExternal())
  //   std::cout << "This is reading the correct external vector.\n";
  // std::cout << "This view has size of: " << mDataStorePtr->getRoot()->getView(pathName)->getNumElements() << "\n";
  // if (mDataStorePtr->getRoot()->hasChildView(pathName))
  //   std::cout << "Root does have the view " << pathName << "\n";

  // #if defined(AXOM_USE_HDF5)
  // std::cout << "Axom is using HDF5 as expected.\n";
  // #endif
}

//------------------------------------------------------------------------------
// Write a vector<double> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const std::vector<double>& value, const std::string pathName)
{
  mDataStorePtr->getRoot()->createView(pathName, axom::sidre::DOUBLE_ID, value.size(), (void*)(&(*value.begin())));
}

//------------------------------------------------------------------------------
// Write a vector<std::string> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const std::vector<std::string>& value, const std::string pathName)
{
  axom::Path myPath = axom::Path(pathName);
  std::pair<std::string, std::string> dirSplit = myPath.split();
  std::cout << "This is the split path, string 1: " << dirSplit.first << ". and string 2: " << dirSplit.second << "\n";
  //create a view with the pathname and then multiple views each containing one string
  axom::sidre::Group* group = mDataStorePtr->getRoot()->createGroup(pathName);
  // for (int i = 0; i < value.size(); ++i)
  //   group->createView(pathName + std::to_string(i), axom::sidre::INT8_ID, value[i].size(), (void*)(&(*value[i].begin())));

  // std::string totalString;
  // for (int i = 0; i < value.size(); ++i)
  //   totalString += value[i];

  // mDataStorePtr->getRoot()->createView(pathName, axom::sidre::INT8_ID, totalString.size(), (void*)(&(*totalString.begin())));
  mDataStorePtr->getRoot()->print();
  std::cout << std::endl;
}

//------------------------------------------------------------------------------
// Write a Dim<1>::Vector to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<1>::Vector& value, const std::string pathName)
{
  int size = int(Dim<1>::Vector::numElements);
  mDataStorePtr->getRoot()->createView(pathName, axom::sidre::DOUBLE_ID, size, (void*)(&(*value.begin())));
}

//------------------------------------------------------------------------------
// Write a Dim<1>::Tensor to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<1>::Tensor& value, const std::string pathName)
{
  int size = int(Dim<1>::Tensor::numElements);
  mDataStorePtr->getRoot()->createView(pathName, axom::sidre::DOUBLE_ID, size, (void*)(&(*value.begin())));
}

//------------------------------------------------------------------------------
// Write a Dim<1>::SymTensor to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<1>::SymTensor& value, const std::string pathName)
{
  int size = int(Dim<1>::SymTensor::numElements);
  mDataStorePtr->getRoot()->createView(pathName, axom::sidre::DOUBLE_ID, size, (void*)(&(*value.begin())));
}

//------------------------------------------------------------------------------
// Write a Dim<1>::ThirdRankTensor to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<1>::ThirdRankTensor& value, const std::string pathName)
{
  int size = int(Dim<1>::ThirdRankTensor::numElements);
  mDataStorePtr->getRoot()->createView(pathName, axom::sidre::DOUBLE_ID, size, (void*)(&(*value.begin())));
}

//------------------------------------------------------------------------------
// Write a Dim<2>::Vector to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<2>::Vector& value, const std::string pathName)
{
  int size = int(Dim<2>::Vector::numElements);
  mDataStorePtr->getRoot()->createView(pathName, axom::sidre::DOUBLE_ID, size, (void*)(&(*value.begin())));
}

//------------------------------------------------------------------------------
// Write a Dim<2>::Tensor to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<2>::Tensor& value, const std::string pathName)
{
  int size = int(Dim<2>::Tensor::numElements);
  mDataStorePtr->getRoot()->createView(pathName, axom::sidre::DOUBLE_ID, size, (void*)(&(*value.begin())));
}

//------------------------------------------------------------------------------
// Write a Dim<2>::SymTensor to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<2>::SymTensor& value, const std::string pathName)
{
  int size = int(Dim<2>::SymTensor::numElements);
  mDataStorePtr->getRoot()->createView(pathName, axom::sidre::DOUBLE_ID, size, (void*)(&(*value.begin())));
}

//------------------------------------------------------------------------------
// Write a Dim<2>::ThirdRankTensor to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<2>::ThirdRankTensor& value, const std::string pathName)
{
  int size = int(Dim<2>::ThirdRankTensor::numElements);
  mDataStorePtr->getRoot()->createView(pathName, axom::sidre::DOUBLE_ID, size, (void*)(&(*value.begin())));
}

//------------------------------------------------------------------------------
// Write a Dim<3>::Vector to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<3>::Vector& value, const std::string pathName)
{
  int size = int(Dim<3>::Vector::numElements);
  mDataStorePtr->getRoot()->createView(pathName, axom::sidre::DOUBLE_ID, size, (void*)(&(*value.begin())));
}

//------------------------------------------------------------------------------
// Write a Dim<3>::Tensor to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<3>::Tensor& value, const std::string pathName)
{
  int size = int(Dim<3>::Tensor::numElements);
  mDataStorePtr->getRoot()->createView(pathName, axom::sidre::DOUBLE_ID, size, (void*)(&(*value.begin())));
}

//------------------------------------------------------------------------------
// Write a Dim<3>::SymTensor to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<3>::SymTensor& value, const std::string pathName)
{
  int size = int(Dim<3>::SymTensor::numElements);
  mDataStorePtr->getRoot()->createView(pathName, axom::sidre::DOUBLE_ID, size, (void*)(&(*value.begin())));
}

//------------------------------------------------------------------------------
// Write a Dim<3>::ThirdRankTensor to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<3>::ThirdRankTensor& value, const std::string pathName)
{
  int size = int(Dim<3>::ThirdRankTensor::numElements);
  mDataStorePtr->getRoot()->createView(pathName, axom::sidre::DOUBLE_ID, size, (void*)(&(*value.begin())));
}


// ------------------------------------------------------------------------------
// Read an unsigned from the file.
// ------------------------------------------------------------------------------
void SidreFileIO::read(unsigned& value, const std::string pathName) const
{
  value = mDataStorePtr->getRoot()->getView(pathName)->getScalar();
}

//------------------------------------------------------------------------------
// Read a size_t from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(size_t& value, const std::string pathName) const
{
  value = mDataStorePtr->getRoot()->getView(pathName)->getScalar();
}

//------------------------------------------------------------------------------
// Read an int from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(int& value, const std::string pathName) const
{
  // mDataStorePtr->getRoot()->getView(pathName)->print();
  value = mDataStorePtr->getRoot()->getView(pathName)->getScalar();
}

//------------------------------------------------------------------------------
// Read a bool from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(bool& value, const std::string pathName) const
{
  value = (int)mDataStorePtr->getRoot()->getView(pathName)->getScalar();
}

//------------------------------------------------------------------------------
// Read a double from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(double& value, const std::string pathName) const
{
  value = mDataStorePtr->getRoot()->getView(pathName)->getScalar();
}

//------------------------------------------------------------------------------
// Read a std::string from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(std::string& value, const std::string pathName) const
{
  value = mDataStorePtr->getRoot()->getView(pathName)->getString();
}

//------------------------------------------------------------------------------
// Read a vector<int> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(std::vector<int>& value, const std::string pathName) const
{
  // std::cout << "---------------------------------------------------------------------------\n";

  // mDataStorePtr->print();
  // if (mDataStorePtr->getRoot()->getView(pathName)->isExternal())
  //   std::cout << "This is reading the correct external vector.\n";
  int size = mDataStorePtr->getRoot()->getView(pathName)->getNumElements();
  value.resize(size);
  // std::cout << value.size() << "\n";
  mDataStorePtr->getRoot()->getView(pathName)->setExternalDataPtr(static_cast<void*>(&value[0]));
  mDataStorePtr->getRoot()->loadExternalData(mFileName);
  // mDataStorePtr->getRoot()->print();
  // std::cout << std::endl;

  // for (int i = 0; i < size; ++i)
  //   std::cout << value[i] << " ";
  // std::cout << std::endl;
  // mDataStorePtr->getRoot()->getView(pathName)->print();
  // value = (std::vector<int>)mDataStorePtr->getRoot()->getView(pathName)->getData();
}

//------------------------------------------------------------------------------
// Read a vector<double> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(std::vector<double>& value, const std::string pathName) const
{
  int size = mDataStorePtr->getRoot()->getView(pathName)->getNumElements();
  value.resize(size);
  mDataStorePtr->getRoot()->getView(pathName)->setExternalDataPtr(static_cast<void*>(&value[0]));
  mDataStorePtr->getRoot()->loadExternalData(mFileName);
  // mDataStorePtr->getRoot()->getView(pathName)->print();
  // value = mDataStorePtr->getRoot()->getView(pathName)->getData();
}

//------------------------------------------------------------------------------
// Read a vector<std::string> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(vector<std::string>& value, const std::string pathName) const
{
  axom::sidre::Group* group = mDataStorePtr->getRoot()->getGroup(pathName);
  int stringCount = group->getNumViews();
  std::cout << "This is the amount of views in the group of strings: " << stringCount << std::endl;
  value.resize(stringCount);
  // int size = mDataStorePtr->getRoot()->getView(pathName)->getNumElements();
  // value.resize(size);
  // std::cout << value.size() << "\n";
  // mDataStorePtr->getRoot()->getView(pathName)->print();
}

//------------------------------------------------------------------------------
// Read a Dim<1>::Vector from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Dim<1>::Vector& value, const std::string pathName) const
{
  mDataStorePtr->getRoot()->getView(pathName)->setExternalDataPtr(static_cast<void*>(&value[0]));
  mDataStorePtr->getRoot()->loadExternalData(mFileName);
}

//------------------------------------------------------------------------------
// Read a Dim<1>::Tensor from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Dim<1>::Tensor& value, const std::string pathName) const
{
  mDataStorePtr->getRoot()->getView(pathName)->setExternalDataPtr(static_cast<void*>(&value[0]));
  mDataStorePtr->getRoot()->loadExternalData(mFileName);
}

//------------------------------------------------------------------------------
// Read a Dim<1>::SymTensor from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Dim<1>::SymTensor& value, const std::string pathName) const
{
  mDataStorePtr->getRoot()->getView(pathName)->setExternalDataPtr(static_cast<void*>(&value[0]));
  mDataStorePtr->getRoot()->loadExternalData(mFileName);
}

//------------------------------------------------------------------------------
// Read a Dim<1>::ThirdRankTensor from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Dim<1>::ThirdRankTensor& value, const std::string pathName) const
{
  mDataStorePtr->getRoot()->getView(pathName)->setExternalDataPtr(static_cast<void*>(&value[0]));
  mDataStorePtr->getRoot()->loadExternalData(mFileName);
}

//------------------------------------------------------------------------------
// Read a Dim<2>::Vector from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Dim<2>::Vector& value, const std::string pathName) const
{
  mDataStorePtr->getRoot()->getView(pathName)->setExternalDataPtr(static_cast<void*>(&value[0]));
  mDataStorePtr->getRoot()->loadExternalData(mFileName);
}

//------------------------------------------------------------------------------
// Read a Dim<2>::Tensor from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Dim<2>::Tensor& value, const std::string pathName) const
{
  mDataStorePtr->getRoot()->getView(pathName)->setExternalDataPtr(static_cast<void*>(&value[0]));
  mDataStorePtr->getRoot()->loadExternalData(mFileName);
}

//------------------------------------------------------------------------------
// Read a Dim<2>::SymTensor from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Dim<2>::SymTensor& value, const std::string pathName) const
{
  mDataStorePtr->getRoot()->getView(pathName)->setExternalDataPtr(static_cast<void*>(&value[0]));
  mDataStorePtr->getRoot()->loadExternalData(mFileName);
}

//------------------------------------------------------------------------------
// Read a Dim<2>::ThirdRankTensor from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Dim<2>::ThirdRankTensor& value, const std::string pathName) const
{
  mDataStorePtr->getRoot()->getView(pathName)->setExternalDataPtr(static_cast<void*>(&value[0]));
  mDataStorePtr->getRoot()->loadExternalData(mFileName);
}

//------------------------------------------------------------------------------
// Read a Dim<3>::Vector from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Dim<3>::Vector& value, const std::string pathName) const
{
  mDataStorePtr->getRoot()->getView(pathName)->setExternalDataPtr(static_cast<void*>(&value[0]));
  mDataStorePtr->getRoot()->loadExternalData(mFileName);
}

//------------------------------------------------------------------------------
// Read a Dim<3>::Tensor from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Dim<3>::Tensor& value, const std::string pathName) const
{
  mDataStorePtr->getRoot()->getView(pathName)->setExternalDataPtr(static_cast<void*>(&value[0]));
  mDataStorePtr->getRoot()->loadExternalData(mFileName);
}

//------------------------------------------------------------------------------
// Read a Dim<3>::SymTensor from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Dim<3>::SymTensor& value, const std::string pathName) const
{
  mDataStorePtr->getRoot()->getView(pathName)->setExternalDataPtr(static_cast<void*>(&value[0]));
  mDataStorePtr->getRoot()->loadExternalData(mFileName);
}

//------------------------------------------------------------------------------
// Read a Dim<3>::ThirdRankTensor from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Dim<3>::ThirdRankTensor& value, const std::string pathName) const
{
  mDataStorePtr->getRoot()->getView(pathName)->setExternalDataPtr(static_cast<void*>(&value[0]));
  mDataStorePtr->getRoot()->loadExternalData(mFileName);
}

#ifdef SPHERAL1D
//------------------------------------------------------------------------------
// Write a Field<Dim<1>, Dim<1>::Scalar> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<1>, Dim<1>::Scalar>& value, const std::string pathName)
{
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, Dim<1>::Vector> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<1>, Dim<1>::Vector>& value, const std::string pathName)
{
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, Dim<1>::Tensor> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<1>, Dim<1>::Tensor>& value, const std::string pathName)
{
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, Dim<1>::SymTensor> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<1>, Dim<1>::SymTensor>& value, const std::string pathName)
{
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, Dim<1>::ThirdRankTensor> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<1>, Dim<1>::ThirdRankTensor>& value, const std::string pathName)
{
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, int> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<1>, int>& value, const std::string pathName)
{
  axom::sidre::DataTypeId dtype = axom::sidre::INT_ID;
  axom::sidre::Group* wholeField = mDataStorePtr->getRoot()->createGroup(pathName);
  axom::IndexType num_elements = 1;

  for (u_int i = 0; i < value.size(); ++i)
  {
    const int *data = &(value[i]);
    wholeField->createView(pathName + std::to_string(i), dtype, num_elements, (void*)data);
  }
  std::cout << "This is a Field<Dim<1>, int> ===================================\n";
  mDataStorePtr->print();
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, unsigned> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<1>, unsigned>& value, const std::string pathName)
{
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, Dim<1>::Scalar> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<1>, Dim<1>::Scalar>& value, const std::string pathName) const
{
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, Dim<1>::Vector> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<1>, Dim<1>::Vector>& value, const std::string pathName) const
{
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, Dim<1>::Tensor> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<1>, Dim<1>::Tensor>& value, const std::string pathName) const
{
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, Dim<1>::SymTensor> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<1>, Dim<1>::SymTensor>& value, const std::string pathName) const
{
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, Dim<1>::ThirdRankTensor> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<1>, Dim<1>::ThirdRankTensor>& value, const std::string pathName) const
{
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, int> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<1>, int>& value, const std::string pathName) const
{
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, unsigned> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<1>, unsigned>& value, const std::string pathName) const
{
}
#endif

#ifdef SPHERAL2D
//------------------------------------------------------------------------------
// Write a Field<Dim<2>, Dim<2>::Scalar> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<2>, Dim<2>::Scalar>& value, const std::string pathName)
{
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, Dim<2>::Vector> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<2>, Dim<2>::Vector>& value, const std::string pathName)
{
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, Dim<2>::Tensor> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<2>, Dim<2>::Tensor>& value, const std::string pathName)
{
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, Dim<2>::SymTensor> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<2>, Dim<2>::SymTensor>& value, const std::string pathName)
{
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, Dim<2>::ThirdRankTensor> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<2>, Dim<2>::ThirdRankTensor>& value, const std::string pathName)
{
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, int> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<2>, int>& value, const std::string pathName)
{
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, unsigned> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<2>, unsigned>& value, const std::string pathName)
{
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, Dim<2>::Scalar> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<2>, Dim<2>::Scalar>& value, const std::string pathName) const
{
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, Dim<2>::Vector> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<2>, Dim<2>::Vector>& value, const std::string pathName) const
{
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, Dim<2>::Tensor> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<2>, Dim<2>::Tensor>& value, const std::string pathName) const
{
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, Dim<2>::SymTensor> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<2>, Dim<2>::SymTensor>& value, const std::string pathName) const
{

}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, Dim<2>::ThirdRankTensor> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<2>, Dim<2>::ThirdRankTensor>& value, const std::string pathName) const
{
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, int> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<2>, int>& value, const std::string pathName) const
{
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, unsigned> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<2>, unsigned>& value, const std::string pathName) const
{
}
#endif

#ifdef SPHERAL3D
//------------------------------------------------------------------------------
// Write a Field<Dim<3>, Dim<3>::Scalar> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<3>, Dim<3>::Scalar>& value, const std::string pathName)
{
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, Dim<3>::Vector> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<3>, Dim<3>::Vector>& value, const std::string pathName)
{
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, Dim<3>::Tensor> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<3>, Dim<3>::Tensor>& value, const std::string pathName)
{
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, Dim<3>::SymTensor> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<3>, Dim<3>::SymTensor>& value, const std::string pathName)
{
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, Dim<3>::ThirdRankTensor> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<3>, Dim<3>::ThirdRankTensor>& value, const std::string pathName)
{
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, int> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<3>, int>& value, const std::string pathName)
{
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, unsigned> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Field<Dim<3>, unsigned>& value, const std::string pathName)
{
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, Dim<3>::Scalar> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<3>, Dim<3>::Scalar>& value, const std::string pathName) const
{
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, Dim<3>::Vector> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<3>, Dim<3>::Vector>& value, const std::string pathName) const
{
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, Dim<3>::Tensor> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<3>, Dim<3>::Tensor>& value, const std::string pathName) const
{
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, Dim<3>::SymTensor> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<3>, Dim<3>::SymTensor>& value, const std::string pathName) const
{
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, Dim<3>::ThirdRankTensor> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<3>, Dim<3>::ThirdRankTensor>& value, const std::string pathName) const
{
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, int> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<3>, int>& value, const std::string pathName) const
{
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, unsigned> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(Field<Dim<3>, unsigned>& value, const std::string pathName) const
{
}
#endif


}