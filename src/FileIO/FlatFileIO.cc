//---------------------------------Spheral++----------------------------------//
// FlatFileIO -- Provide the interface to FlatFile file objects.
//
// Created by JMO, Fri Apr 13 01:19:02 PDT 2001
//----------------------------------------------------------------------------//
#include "FlatFileIO.hh"
#include "Field/Field.hh"

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
using std::ios;

namespace Spheral {

//------------------------------------------------------------------------------
// Empty constructor.
//------------------------------------------------------------------------------
FlatFileIO::FlatFileIO():
  FileIO(),
  mPrecision(20),
  mFilePtr(0) {
}

//------------------------------------------------------------------------------
// Construct and open the given file.
//------------------------------------------------------------------------------
FlatFileIO::
FlatFileIO(const string fileName, AccessType access, FlatFileFormat format):
  FileIO(fileName, access),
  mPrecision(20),
  mFilePtr(0),
  mFileFormat(format) {
  open(fileName, access);
  ENSURE(mFileOpen && mFilePtr != 0);
}

// Do the same thing with an int in the place of the enum, to help Pyffle.
FlatFileIO::
FlatFileIO(const string fileName, int access, int format):
  FileIO(fileName, AccessType(access)),
  mPrecision(20),
  mFilePtr(0),
  mFileFormat(FlatFileFormat(format)) {
  open(fileName, AccessType(access));
  ENSURE(mFileOpen && mFilePtr != 0);
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
FlatFileIO::~FlatFileIO() {
  close();
}

//------------------------------------------------------------------------------
// Open a FlatFile file with the specified access.
//------------------------------------------------------------------------------
void
FlatFileIO::open(const string fileName, AccessType access) {

  // If we currently have a file open and attached to this object, close it!
  close();
  CHECK(mFilePtr == 0 && mFileOpen == false);

  // Build the file opening mode.
  ios::openmode mode;
  switch(access) {
  case AccessType::Create:
    mode = ios::out; //ios::trunc;
    break;
  case AccessType::Read:
    mode = ios::in;
    break;
  case AccessType::Write:
    mode = ios::out;
    break;
  case AccessType::ReadWrite:
    mode = ios::in | ios::out;
    break;
  default:
    VERIFY2(false, "Unhandled case in switch!");
  }

  if (mFileFormat == FlatFileFormat::binary) mode = mode | ios::binary;

  // Open a file stream and attach it to this objects pointer.
  mFilePtr = new std::fstream(fileName.c_str(), mode);
  mFileOpen = mFilePtr->is_open();
  ENSURE(mFilePtr != 0 && mFileOpen == true);

  // Set the precision we will be using to write to the file.
  (*mFilePtr).precision(mPrecision);
}

//------------------------------------------------------------------------------
// Close the current file.
//------------------------------------------------------------------------------
void
FlatFileIO::close() {

  if (mFilePtr != 0) {
    mFilePtr->flush();
    mFilePtr->close();
    delete mFilePtr;
    mFilePtr = 0;
  }
  mFileOpen = false;
}

//------------------------------------------------------------------------------
// Write an unsigned to the file.
//------------------------------------------------------------------------------
void
FlatFileIO::write(const unsigned value, const string pathName) {
  writeGenericType(value, pathName);
}

//------------------------------------------------------------------------------
// Write an int to the file.
//------------------------------------------------------------------------------
void
FlatFileIO::write(const int value, const string pathName) {
  writeGenericType(value, pathName);
}

//------------------------------------------------------------------------------
// Write a bool to the file.
//------------------------------------------------------------------------------
void
FlatFileIO::write(const bool value, const string pathName) {
  writeGenericType(value, pathName);
}

//------------------------------------------------------------------------------
// Write a Scalar to the file.
//------------------------------------------------------------------------------
void
FlatFileIO::
write(const double value, const string pathName) {
  writeGenericType(value, pathName);
}

//------------------------------------------------------------------------------
// Write a Vector to the file.
//------------------------------------------------------------------------------
void
FlatFileIO::
write(const Dim<1>::Vector& value, const string pathName) {
  writeGenericType(value, pathName);
}

void
FlatFileIO::
write(const Dim<2>::Vector& value, const string pathName) {
  writeGenericType(value, pathName);
}

void
FlatFileIO::
write(const Dim<3>::Vector& value, const string pathName) {
  writeGenericType(value, pathName);
}

//------------------------------------------------------------------------------
// Write a Tensor to the file.
//------------------------------------------------------------------------------
void
FlatFileIO::
write(const Dim<1>::Tensor& value, const string pathName) {
  writeGenericType(value, pathName);
}

void
FlatFileIO::
write(const Dim<2>::Tensor& value, const string pathName) {
  writeGenericType(value, pathName);
}

void
FlatFileIO::
write(const Dim<3>::Tensor& value, const string pathName) {
  writeGenericType(value, pathName);
}

//------------------------------------------------------------------------------
// Write a SymTensor to the file.
//------------------------------------------------------------------------------
void
FlatFileIO::
write(const Dim<1>::SymTensor& value, const string pathName) {
  writeGenericType(value, pathName);
}

void
FlatFileIO::
write(const Dim<2>::SymTensor& value, const string pathName) {
  writeGenericType(value, pathName);
}

void
FlatFileIO::
write(const Dim<3>::SymTensor& value, const string pathName) {
  writeGenericType(value, pathName);
}

//------------------------------------------------------------------------------
// Write a ThirdRankTensor to the file.
//------------------------------------------------------------------------------
void
FlatFileIO::
write(const Dim<1>::ThirdRankTensor& value, const string pathName) {
  writeGenericType(value, pathName);
}

void
FlatFileIO::
write(const Dim<2>::ThirdRankTensor& value, const string pathName) {
  writeGenericType(value, pathName);
}

void
FlatFileIO::
write(const Dim<3>::ThirdRankTensor& value, const string pathName) {
  writeGenericType(value, pathName);
}

//------------------------------------------------------------------------------
// Write a string to the file.
//------------------------------------------------------------------------------
void
FlatFileIO::
write(const string value, const string pathName) {
  writeGenericType(value, pathName);
}

//------------------------------------------------------------------------------
// Read an unsigned from the file.
//------------------------------------------------------------------------------
void
FlatFileIO::read(unsigned& value, const string pathName) const {
  readGenericType(value, pathName);
}

//------------------------------------------------------------------------------
// Read an int from the file.
//------------------------------------------------------------------------------
void
FlatFileIO::read(int& value, const string pathName) const {
  readGenericType(value, pathName);
}

//------------------------------------------------------------------------------
// Read a bool from the file.
//------------------------------------------------------------------------------
void
FlatFileIO::read(bool& value, const string pathName) const {
  readGenericType(value, pathName);
}

//------------------------------------------------------------------------------
// Read a Scalar from the file.
//------------------------------------------------------------------------------
void
FlatFileIO::read(double& value, const string pathName) const {
  readGenericType(value, pathName);
}

//------------------------------------------------------------------------------
// Read a Vector from the file.
//------------------------------------------------------------------------------
void
FlatFileIO::read(Dim<1>::Vector& value, const string pathName) const {
  readGenericType(value, pathName);
}

void
FlatFileIO::read(Dim<2>::Vector& value, const string pathName) const {
  readGenericType(value, pathName);
}

void
FlatFileIO::read(Dim<3>::Vector& value, const string pathName) const {
  readGenericType(value, pathName);
}

//------------------------------------------------------------------------------
// Read a Tensor from the file.
//------------------------------------------------------------------------------
void
FlatFileIO::read(Dim<1>::Tensor& value, const string pathName) const {
  readGenericType(value, pathName);
}

void
FlatFileIO::read(Dim<2>::Tensor& value, const string pathName) const {
  readGenericType(value, pathName);
}

void
FlatFileIO::read(Dim<3>::Tensor& value, const string pathName) const {
  readGenericType(value, pathName);
}

//------------------------------------------------------------------------------
// Read a SymTensor from the file.
//------------------------------------------------------------------------------
void
FlatFileIO::read(Dim<1>::SymTensor& value, const string pathName) const {
  readGenericType(value, pathName);
}

void
FlatFileIO::read(Dim<2>::SymTensor& value, const string pathName) const {
  readGenericType(value, pathName);
}

void
FlatFileIO::read(Dim<3>::SymTensor& value, const string pathName) const {
  readGenericType(value, pathName);
}

//------------------------------------------------------------------------------
// Read a ThirdRankTensor from the file.
//------------------------------------------------------------------------------
void
FlatFileIO::read(Dim<1>::ThirdRankTensor& value, const string pathName) const {
  readGenericType(value, pathName);
}

void
FlatFileIO::read(Dim<2>::ThirdRankTensor& value, const string pathName) const {
  readGenericType(value, pathName);
}

void
FlatFileIO::read(Dim<3>::ThirdRankTensor& value, const string pathName) const {
  readGenericType(value, pathName);
}

//------------------------------------------------------------------------------
// Read a string from the file.
//------------------------------------------------------------------------------
void
FlatFileIO::read(string& value, const string pathName) const {
  REQUIRE(readyToRead());

  // Set the current position to the beginning of the file.
//   mFilePtr->seekg(ios::beg);
  beginningOfFile();

  // Find the requested path, and read the value.
  findPathName(pathName);
  getline(*mFilePtr, value);

  // We pick up the leading space in the string, so get rid of that.
  if (value.size() > 0) value.erase(0,1);
}

//------------------------------------------------------------------------------
// Write a vector<int>.
//------------------------------------------------------------------------------
void
FlatFileIO::
write(const vector<int>& value, const string pathName) {
  writeGenericVector(value, pathName);
}

//------------------------------------------------------------------------------
// Write a vector<bool>.
//------------------------------------------------------------------------------
// void
// FlatFileIO::
// write(const vector<bool>& value, const string pathName) {

//   // Convert to a vector of int and store that.
//   vector<int> tmp(value.size());
//   for (int i = 0; i < value.size(); ++i) {
//     if (value[i]) {
//       tmp[i] = 1;
//     } else {
//       tmp[i] = 0;
//     }
//   }
//   writeGenericVector(tmp, pathName);
// }

//------------------------------------------------------------------------------
// Write a vector<Scalar>.
//------------------------------------------------------------------------------
void
FlatFileIO::
write(const vector<double>& value, const string pathName) {
  writeGenericVector(value, pathName);
}

//------------------------------------------------------------------------------
// Write a vector<Vector>.
//------------------------------------------------------------------------------
void
FlatFileIO::
write(const vector<Dim<1>::Vector>& value, const string pathName) {
  writeGenericVector(value, pathName);
}

void
FlatFileIO::
write(const vector<Dim<2>::Vector>& value, const string pathName) {
  writeGenericVector(value, pathName);
}

void
FlatFileIO::
write(const vector<Dim<3>::Vector>& value, const string pathName) {
  writeGenericVector(value, pathName);
}

//------------------------------------------------------------------------------
// Write a vector<Tensor>.
//------------------------------------------------------------------------------
void
FlatFileIO::
write(const vector<Dim<1>::Tensor>& value, const string pathName) {
  writeGenericVector(value, pathName);
}

void
FlatFileIO::
write(const vector<Dim<2>::Tensor>& value, const string pathName) {
  writeGenericVector(value, pathName);
}

void
FlatFileIO::
write(const vector<Dim<3>::Tensor>& value, const string pathName) {
  writeGenericVector(value, pathName);
}

//------------------------------------------------------------------------------
// Write a vector<SymTensor>.
//------------------------------------------------------------------------------
void
FlatFileIO::
write(const vector<Dim<1>::SymTensor>& value, const string pathName) {
  writeGenericVector(value, pathName);
}

void
FlatFileIO::
write(const vector<Dim<2>::SymTensor>& value, const string pathName) {
  writeGenericVector(value, pathName);
}

void
FlatFileIO::
write(const vector<Dim<3>::SymTensor>& value, const string pathName) {
  writeGenericVector(value, pathName);
}

//------------------------------------------------------------------------------
// Write a vector<ThirdRankTensor>.
//------------------------------------------------------------------------------
void
FlatFileIO::
write(const vector<Dim<1>::ThirdRankTensor>& value, const string pathName) {
  writeGenericVector(value, pathName);
}

void
FlatFileIO::
write(const vector<Dim<2>::ThirdRankTensor>& value, const string pathName) {
  writeGenericVector(value, pathName);
}

void
FlatFileIO::
write(const vector<Dim<3>::ThirdRankTensor>& value, const string pathName) {
  writeGenericVector(value, pathName);
}

//------------------------------------------------------------------------------
// Write a vector<string>.
//------------------------------------------------------------------------------
void
FlatFileIO::
write(const vector<std::string>& value, const string pathName) {
  writeGenericVector(value, pathName);
}

//------------------------------------------------------------------------------
// Read a vector<int>.
//------------------------------------------------------------------------------
void
FlatFileIO::
read(vector<int>& value, const string pathName) const {
  readGenericVector(value, pathName);
}

//------------------------------------------------------------------------------
// Read a vector<bool>.
//------------------------------------------------------------------------------
// void
// FlatFileIO::
// read(vector<bool>& value, const string pathName) const {

//   // Read to a temporary vector<int>, and convert it.
//   vector<int> tmp(value.size());
//   readGenericVector(tmp, pathName);
//   value.resize(tmp.size());
//   for (int i = 0; i < value.size(); ++i) {
//     if (tmp[i] == 1) {
//       value[i] = true;
//     } else {
//       value[i] = false;
//     }
//   }
// }

//------------------------------------------------------------------------------
// Read a vector<Scalar>.
//------------------------------------------------------------------------------
void
FlatFileIO::
read(vector<double>& value, const string pathName) const {
  readGenericVector(value, pathName);
}

//------------------------------------------------------------------------------
// Read a vector<Vector>.
//------------------------------------------------------------------------------
void
FlatFileIO::
read(vector<Dim<1>::Vector>& value, const string pathName) const {
  readGenericVector(value, pathName);
}

void
FlatFileIO::
read(vector<Dim<2>::Vector>& value, const string pathName) const {
  readGenericVector(value, pathName);
}

void
FlatFileIO::
read(vector<Dim<3>::Vector>& value, const string pathName) const {
  readGenericVector(value, pathName);
}

//------------------------------------------------------------------------------
// Read a vector<Tensor>.
//------------------------------------------------------------------------------
void
FlatFileIO::
read(vector<Dim<1>::Tensor>& value, const string pathName) const {
  readGenericVector(value, pathName);
}

void
FlatFileIO::
read(vector<Dim<2>::Tensor>& value, const string pathName) const {
  readGenericVector(value, pathName);
}

void
FlatFileIO::
read(vector<Dim<3>::Tensor>& value, const string pathName) const {
  readGenericVector(value, pathName);
}

//------------------------------------------------------------------------------
// Read a vector<SymTensor>.
//------------------------------------------------------------------------------
void
FlatFileIO::
read(vector<Dim<1>::SymTensor>& value, const string pathName) const {
  readGenericVector(value, pathName);
}

void
FlatFileIO::
read(vector<Dim<2>::SymTensor>& value, const string pathName) const {
  readGenericVector(value, pathName);
}

void
FlatFileIO::
read(vector<Dim<3>::SymTensor>& value, const string pathName) const {
  readGenericVector(value, pathName);
}

//------------------------------------------------------------------------------
// Read a vector<ThirdRankTensor>.
//------------------------------------------------------------------------------
void
FlatFileIO::
read(vector<Dim<1>::ThirdRankTensor>& value, const string pathName) const {
  readGenericVector(value, pathName);
}

void
FlatFileIO::
read(vector<Dim<2>::ThirdRankTensor>& value, const string pathName) const {
  readGenericVector(value, pathName);
}

void
FlatFileIO::
read(vector<Dim<3>::ThirdRankTensor>& value, const string pathName) const {
  readGenericVector(value, pathName);
}

//------------------------------------------------------------------------------
// Read a vector<string>.
//------------------------------------------------------------------------------
void
FlatFileIO::
read(vector<std::string>& value, const string pathName) const {
  readGenericVector(value, pathName);
}

//------------------------------------------------------------------------------
// 1d fields (write).
//------------------------------------------------------------------------------
#ifdef SPHERAL1D
// Scalar
void
FlatFileIO::
write(const Field<Dim<1>, Dim<1>::Scalar>& value, const string pathName) {
  writeGenericType(value, pathName);
}

//Vector
void
FlatFileIO::
write(const Field<Dim<1>, Dim<1>::Vector>& value, const string pathName) {
  writeGenericType(value, pathName);
}

// Tensor
void
FlatFileIO::
write(const Field<Dim<1>, Dim<1>::Tensor>& value, const string pathName) {
  writeGenericType(value, pathName);
}

// SymTensor
void
FlatFileIO::
write(const Field<Dim<1>, Dim<1>::SymTensor>& value, const string pathName) {
  writeGenericType(value, pathName);
}

// ThirdRankTensor
void
FlatFileIO::
write(const Field<Dim<1>, Dim<1>::ThirdRankTensor>& value, const string pathName) {
  writeGenericType(value, pathName);
}

// int
void
FlatFileIO::
write(const Field<Dim<1>, int>& value, const string pathName) {
  writeGenericType(value, pathName);
}

//------------------------------------------------------------------------------
// 1d fields (read).
//------------------------------------------------------------------------------
// Scalar
void
FlatFileIO::
read(Field<Dim<1>, Dim<1>::Scalar>& value, const string pathName) const {
  readGenericType(value, pathName);
}

// Vector
void
FlatFileIO::
read(Field<Dim<1>, Dim<1>::Vector>& value, const string pathName) const {
  readGenericType(value, pathName);
}

// Tensor
void
FlatFileIO::
read(Field<Dim<1>, Dim<1>::Tensor>& value, const string pathName) const {
  readGenericType(value, pathName);
}

// SymTensor
void
FlatFileIO::
read(Field<Dim<1>, Dim<1>::SymTensor>& value, const string pathName) const {
  readGenericType(value, pathName);
}

// ThirdRankTensor
void
FlatFileIO::
read(Field<Dim<1>, Dim<1>::ThirdRankTensor>& value, const string pathName) const {
  readGenericType(value, pathName);
}

// int
void
FlatFileIO::
read(Field<Dim<1>, int>& value, const string pathName) const {
  readGenericType(value, pathName);
}
#endif

//------------------------------------------------------------------------------
// 2d fields (write).
//------------------------------------------------------------------------------
#ifdef SPHERAL2D
// Scalar
void
FlatFileIO::
write(const Field<Dim<2>, Dim<2>::Scalar>& value, const string pathName) {
  writeGenericType(value, pathName);
}

//Vector
void
FlatFileIO::
write(const Field<Dim<2>, Dim<2>::Vector>& value, const string pathName) {
  writeGenericType(value, pathName);
}

// Tensor
void
FlatFileIO::
write(const Field<Dim<2>, Dim<2>::Tensor>& value, const string pathName) {
  writeGenericType(value, pathName);
}

// SymTensor
void
FlatFileIO::
write(const Field<Dim<2>, Dim<2>::SymTensor>& value, const string pathName) {
  writeGenericType(value, pathName);
}

// ThirdRankTensor
void
FlatFileIO::
write(const Field<Dim<2>, Dim<2>::ThirdRankTensor>& value, const string pathName) {
  writeGenericType(value, pathName);
}

// int
void
FlatFileIO::
write(const Field<Dim<2>, int>& value, const string pathName) {
  writeGenericType(value, pathName);
}

//------------------------------------------------------------------------------
// 2d fields (read).
//------------------------------------------------------------------------------
// Scalar
void
FlatFileIO::
read(Field<Dim<2>, Dim<2>::Scalar>& value, const string pathName) const {
  readGenericType(value, pathName);
}

// Vector
void
FlatFileIO::
read(Field<Dim<2>, Dim<2>::Vector>& value, const string pathName) const {
  readGenericType(value, pathName);
}

// Tensor
void
FlatFileIO::
read(Field<Dim<2>, Dim<2>::Tensor>& value, const string pathName) const {
  readGenericType(value, pathName);
}

// SymTensor
void
FlatFileIO::
read(Field<Dim<2>, Dim<2>::SymTensor>& value, const string pathName) const {
  readGenericType(value, pathName);
}

// ThirdRankTensor
void
FlatFileIO::
read(Field<Dim<2>, Dim<2>::ThirdRankTensor>& value, const string pathName) const {
  readGenericType(value, pathName);
}

// int
void
FlatFileIO::
read(Field<Dim<2>, int>& value, const string pathName) const {
  readGenericType(value, pathName);
}
#endif

//------------------------------------------------------------------------------
// 3d fields (write).
//------------------------------------------------------------------------------
#ifdef SPHERAL3D
// Scalar
void
FlatFileIO::
write(const Field<Dim<3>, Dim<3>::Scalar>& value, const string pathName) {
  writeGenericType(value, pathName);
}

//Vector
void
FlatFileIO::
write(const Field<Dim<3>, Dim<3>::Vector>& value, const string pathName) {
  writeGenericType(value, pathName);
}

// Tensor
void
FlatFileIO::
write(const Field<Dim<3>, Dim<3>::Tensor>& value, const string pathName) {
  writeGenericType(value, pathName);
}

// SymTensor
void
FlatFileIO::
write(const Field<Dim<3>, Dim<3>::SymTensor>& value, const string pathName) {
  writeGenericType(value, pathName);
}

// ThirdRankTensor
void
FlatFileIO::
write(const Field<Dim<3>, Dim<3>::ThirdRankTensor>& value, const string pathName) {
  writeGenericType(value, pathName);
}

// int
void
FlatFileIO::
write(const Field<Dim<3>, int>& value, const string pathName) {
  writeGenericType(value, pathName);
}

//------------------------------------------------------------------------------
// 3d fields (read).
//------------------------------------------------------------------------------
// Scalar
void
FlatFileIO::
read(Field<Dim<3>, Dim<3>::Scalar>& value, const string pathName) const {
  readGenericType(value, pathName);
}

// Vector
void
FlatFileIO::
read(Field<Dim<3>, Dim<3>::Vector>& value, const string pathName) const {
  readGenericType(value, pathName);
}

// Tensor
void
FlatFileIO::
read(Field<Dim<3>, Dim<3>::Tensor>& value, const string pathName) const {
  readGenericType(value, pathName);
}

// SymTensor
void
FlatFileIO::
read(Field<Dim<3>, Dim<3>::SymTensor>& value, const string pathName) const {
  readGenericType(value, pathName);
}

// ThirdRankTensor
void
FlatFileIO::
read(Field<Dim<3>, Dim<3>::ThirdRankTensor>& value, const string pathName) const {
  readGenericType(value, pathName);
}

// int
void
FlatFileIO::
read(Field<Dim<3>, int>& value, const string pathName) const {
  readGenericType(value, pathName);
}
#endif

//------------------------------------------------------------------------------
// Access the precision we're using to write to the file.
//------------------------------------------------------------------------------
int
FlatFileIO::
precision() const {
  return mPrecision;
}

void
FlatFileIO::
setPrecision(int precision) {
  REQUIRE(precision > 0);
  mPrecision = precision;
}

//------------------------------------------------------------------------------
// Generic function to write an arbitrary DataType to a FlatFile file.
//------------------------------------------------------------------------------
template<typename DataType>
void
FlatFileIO::
writeGenericType(const DataType& value,
		 const string pathName) {
  REQUIRE(readyToWrite());
//   mFilePtr->seekp(ios::end);
  *mFilePtr << pathName << " " << value << endl;
}

//------------------------------------------------------------------------------
// Generic function to read an arbitrary DataType from a FlatFile file.
//------------------------------------------------------------------------------
template<typename DataType>
void
FlatFileIO::
readGenericType(DataType& value,
		const string pathName) const {


  REQUIRE(readyToRead());

  // Set the current position to the beginning of the file.
//   mFilePtr->seekg(ios::beg);
  beginningOfFile();

  // Find the requested path, and read the value.
  findPathName(pathName);
  *mFilePtr >> value;

}

//------------------------------------------------------------------------------
// Generic function to write a vector<DataType> to a FlatFile file.
//------------------------------------------------------------------------------
template<typename DataType>
void
FlatFileIO::
writeGenericVector(const vector<DataType>& value,
                   const string pathName) {
  REQUIRE(readyToWrite());
//   mFilePtr->seekp(ios::end);
  *mFilePtr << pathName << " " << value.size();
  for (typename vector<DataType>::const_iterator itr = value.begin();
       itr < value.end();
       ++itr) *mFilePtr << " " << *itr;
  *mFilePtr << endl;
}

//------------------------------------------------------------------------------
// Generic function to read a vector<DataType> from a FlatFile file.
//------------------------------------------------------------------------------
template<typename DataType>
void
FlatFileIO::
readGenericVector(vector<DataType>& value,
                  const string pathName) const {


  REQUIRE(readyToRead());

  // Set the current position to the beginning of the file.
//   mFilePtr->seekg(ios::beg);
  beginningOfFile();

  // Find the requested path, and read the value.
  findPathName(pathName);

  // Read the size of the stored vector data.
  int size;
  *mFilePtr >> size;
  value.resize(size);

  // Now read in the vector values.
  for (typename vector<DataType>::iterator itr = value.begin();
       itr < value.end();
       ++itr) *mFilePtr >> *itr;

}

//------------------------------------------------------------------------------
// A function to determine if the current file is ready to write to.
//------------------------------------------------------------------------------
bool
FlatFileIO::readyToWrite() const {
  return (mFilePtr != 0 &&
	  (access() == AccessType::Write || access() == AccessType::ReadWrite || access() == AccessType::Create));
}

//------------------------------------------------------------------------------
// A function to determine if the current file is ready to read from.
//------------------------------------------------------------------------------
bool
FlatFileIO::readyToRead() const {
  return (mFilePtr != 0 && 
	  (access() == AccessType::Read || access() == AccessType::ReadWrite || access() == AccessType::Create));
}

//------------------------------------------------------------------------------
// Scan the file until we find the given path name.
//------------------------------------------------------------------------------
void
FlatFileIO::findPathName(const string pathName) const {
  REQUIRE(readyToRead());


  // Set the pointer to the beginning of the file stream.
//   mFilePtr->seekg(ios::beg);
  beginningOfFile();

  // Scan the file line by line until we find the line containing 
  // the path, or we hit the end.
  string currentPath;
  while (!mFilePtr->eof() && currentPath != pathName) {
    currentPath = "";
    char thpt = '0';
    int i = 0;
    while (!mFilePtr->eof() && thpt != '\n' && i < pathName.size()) {
      mFilePtr->get(thpt);
      currentPath += thpt;
      ++i;
    }
//     cerr << pathName << " " << currentPath << " " << (pathName == currentPath) << endl;

    // If this line doesn't contain the path, skip to the next line.
    if (currentPath != pathName) {
      while (!mFilePtr->eof() && thpt != '\n') mFilePtr->get(thpt);
    }
  }

  VERIFY2(currentPath == pathName,
          "FlatFileIO::findPathName ERROR: couldn't find path " << pathName);
}

//------------------------------------------------------------------------------
// Position the current pointer at the beginning of the file.
//------------------------------------------------------------------------------
void
FlatFileIO::beginningOfFile() const {

//   mFilePtr->seekg(0);
//   ENSURE(mFilePtr->tellg() == 0);

//   close();
//   open(fileName(), access());

  delete mFilePtr;
  mFilePtr = new std::fstream(fileName().c_str(), ios::in);
}

}
