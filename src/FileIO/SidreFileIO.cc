//---------------------------------Spheral++----------------------------------//
// SidreFileIO -- Provide the interface to sidre file objects.
//
// Created by Mikhail Zakharchanka, 11/4/2021
//----------------------------------------------------------------------------//
#include "SidreFileIO.hh"

namespace Spheral
{

//------------------------------------------------------------------------------
// Empty constructor.
//------------------------------------------------------------------------------
SidreFileIO::SidreFileIO():
  FileIO(),
  mFilePtr(0) {}

//------------------------------------------------------------------------------
// Construct and open the given file.
//------------------------------------------------------------------------------
SidreFileIO::SidreFileIO(const std::string fileName, AccessType access):
  FileIO(fileName, access),
  mFilePtr(0)
{
  open(fileName, access);
  ENSURE(mFileOpen && mFilePtr != 0);
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
  VERIFY2(mFilePtr == 0 and mFileOpen == false,
          "ERROR: attempt to reopen SidreFileIO object.");

  // string fullFileName = fileName;
  // if (fullFileName.find(".silo") == string::npos)
  //   fullFileName += ".silo";


  // if (access == AccessType::Read) 
  //   mFilePtr = DBOpen(fullFileName.c_str(), DB_HDF5, DB_READ);
  // else
  //   mFilePtr = DBCreate(fullFileName.c_str(), DB_CLOBBER, DB_LOCAL, "Spheral++ restart file.", DB_HDF5);
  
  VERIFY2(mFilePtr != 0, "SidreFileIO ERROR: unable to open " << fileName);
  mFileOpen = true;
}

//------------------------------------------------------------------------------
// Close the current file.
//------------------------------------------------------------------------------
void SidreFileIO::close()
{
  if (mFilePtr != 0)
  {
    // VERIFY2(DBClose(mFilePtr) == 0,
    //         "SidreFileIO ERROR: unable to close file.");
    mFilePtr = 0;
  }
  mFileOpen = false;
}

//------------------------------------------------------------------------------
// Check if the specified path is in the file.
//------------------------------------------------------------------------------
bool SidreFileIO::pathExists(const std::string pathName) const
{
  return false;
}

//------------------------------------------------------------------------------
// Write an unsigned to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const unsigned& value, const std::string pathName)
{
  mFilePtr->getRoot()->createViewScalar(pathName, value);
}

// //------------------------------------------------------------------------------
// // Write a size_t to the file.
// //------------------------------------------------------------------------------
// void SidreFileIO::write(const size_t& value, const std::string pathName)
// {
//   mFilePtr->getRoot()->createViewScalar(pathName, value);
// }

//------------------------------------------------------------------------------
// Write an int to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const int& value, const std::string pathName)
{
  mFilePtr->getRoot()->createViewScalar(pathName, value);
}

//------------------------------------------------------------------------------
// Write a bool to the file.
//------------------------------------------------------------------------------
// void SidreFileIO::write(const bool& value, const std::string pathName)
// {
//   mFilePtr->getRoot()->createViewScalar(pathName, value);
// }

// //------------------------------------------------------------------------------
// // Write a double to the file.
// //------------------------------------------------------------------------------
// void SidreFileIO::write(const double& value, const string pathName)
// {
//   mFilePtr->getRoot()->createViewScalar(pathName, value);
// }

// //------------------------------------------------------------------------------
// // Write a string to the file.
// //------------------------------------------------------------------------------
// void SidreFileIO::write(const string& value, const string pathName) 
// {
//   mFilePtr->getRoot()->createViewString(pathName, value);
// }

// //------------------------------------------------------------------------------
// // Write a vector<int> to the file.
// //------------------------------------------------------------------------------
// void SidreFileIO::write(const std::vector<int>& value, const string pathName)
// {
//   mFilePtr->getRoot()->createView(pathName, int, value.size(), (void*)&value);
// }

// //------------------------------------------------------------------------------
// // Write a vector<double> to the file.
// //------------------------------------------------------------------------------
// void SidreFileIO::write(const std::vector<double>& value, const string pathName)
// {
//   mFilePtr->getRoot()->createView(pathName, double, value.size(), (void*)&value);
// }

// //------------------------------------------------------------------------------
// // Write a vector<string> to the file.
// //------------------------------------------------------------------------------
// void SidreFileIO::write(const std::vector<string>& value, const string pathName)
// {
  
// }


// ------------------------------------------------------------------------------
// Read an unsigned from the file.
// ------------------------------------------------------------------------------
void SidreFileIO::read(unsigned& value, const std::string pathName) const
{
  value = mFilePtr->getRoot()->getView(pathName)->getScalar();
}

// //------------------------------------------------------------------------------
// // Read a size_t from the file.
// //------------------------------------------------------------------------------
// void SidreFileIO::read(size_t& value, const string pathName) const
// {
//   value = mFilePtr->getRoot()->getView(pathName)->getScalar();
// }

//------------------------------------------------------------------------------
// Read an int to the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(int& value, const std::string pathName) const
{
  value = mFilePtr->getRoot()->getView(pathName)->getScalar();
}

// //------------------------------------------------------------------------------
// // Read a bool from the file.
// //------------------------------------------------------------------------------
// void SidreFileIO::read(bool& value, const std::string pathName) const
// {
//   value = mFilePtr->getRoot()->getView(pathName)->getScalar();
// }

// //------------------------------------------------------------------------------
// // Read a double from the file.
// //------------------------------------------------------------------------------
// void SidreFileIO::read(double& value, const string pathName) const
// {
//   value = mFilePtr->getRoot()->getView(pathName)->getScalar();
// }

// //------------------------------------------------------------------------------
// // Read a string from the file.
// //------------------------------------------------------------------------------
// void SidreFileIO::read(string& value, const string pathName) const
// {
//   value = mFilePtr->getRoot()->getView(pathName)->getString();

// //------------------------------------------------------------------------------
// // Read a vector<int> from the file.
// //------------------------------------------------------------------------------
// void SidreFileIO::read(std::vector<int>& value, const string pathName) const
// {
//   value = mFilePtr->getRoot()->getView(pathName)->getData();
// }

// //------------------------------------------------------------------------------
// // Read a vector<double> from the file.
// //------------------------------------------------------------------------------
// void SidreFileIO::read(std::vector<double>& value, const string pathName) const
// {
//   value = mFilePtr->getRoot()->getView(pathName)->getData();
// }

// //------------------------------------------------------------------------------
// // Read a vector<string> from the file.
// //------------------------------------------------------------------------------
// void SidreFileIO::read(vector<string>& value, const string pathName) const
// {
  
// }

}