//---------------------------------Spheral++----------------------------------//
// SiloFileIO -- Provide the interface to silo file objects.
//
// Created by JMO, Sat Feb  7 23:06:03 PST 2015
//----------------------------------------------------------------------------//
#include "SiloFileIO.hh"
#include "Field/Field.hh"

#include "boost/algorithm/string/replace.hpp"

#include <algorithm>
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

namespace {
//------------------------------------------------------------------------------
// Mangle a path name for silo consumption.
//------------------------------------------------------------------------------
string SILO_mangle(const string& x) {
  string result = x;
  boost::replace_all(result, " ", "_");
  boost::replace_all(result, ".", "_p_");
  boost::replace_all(result, ",", "_c_");
  boost::replace_all(result, "|", "_P_");
  return result;
}

}

//------------------------------------------------------------------------------
// Empty constructor.
//------------------------------------------------------------------------------
SiloFileIO::SiloFileIO():
  FileIO(),
  mFilePtr(0) {
}

//------------------------------------------------------------------------------
// Construct and open the given file.
//------------------------------------------------------------------------------
SiloFileIO::
SiloFileIO(const string fileName, AccessType access):
  FileIO(fileName, access),
  mFilePtr(0) {
  open(fileName, access);
  ENSURE(mFileOpen && mFilePtr != 0);
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
SiloFileIO::~SiloFileIO() {
  close();
}

//------------------------------------------------------------------------------
// Open a SiloFile file with the specified access.
//------------------------------------------------------------------------------
void
SiloFileIO::open(const string fileName, AccessType access) {
  VERIFY2(mFilePtr == 0 and mFileOpen == false,
          "ERROR: attempt to reopen SiloFileIO object.");

  string fullFileName = fileName;
  if (fullFileName.find(".silo") == string::npos) {
    fullFileName += ".silo";
  }

  if (access == AccessType::Read) {
    mFilePtr = DBOpen(fullFileName.c_str(), DB_HDF5, DB_READ);
  } else {
    mFilePtr = DBCreate(fullFileName.c_str(), DB_CLOBBER, DB_LOCAL, "Spheral++ restart file.", DB_HDF5);
  }
  VERIFY2(mFilePtr != 0, "SiloFileIO ERROR: unable to open " << fullFileName);
  mFileOpen = true;
}

//------------------------------------------------------------------------------
// Close the current file.
//------------------------------------------------------------------------------
void
SiloFileIO::close() {
  if (mFilePtr != 0) {
    VERIFY2(DBClose(mFilePtr) == 0,
            "SiloFileIO ERROR: unable to close file.");
    mFilePtr = 0;
  }
  mFileOpen = false;
}

//------------------------------------------------------------------------------
// Write an unsigned to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const unsigned& value, const string pathName) {
  const string varname = this->setDir(pathName);
  int dims[1] = {1};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &value, dims, 1, DB_INT) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write an int to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const int& value, const string pathName) {
  const string varname = this->setDir(pathName);
  int dims[1] = {1};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &value, dims, 1, DB_INT) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a bool to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const bool& value, const string pathName) {
  const string varname = this->setDir(pathName);
  int dims[1] = {1};
  int ivalue = value ? 1 : 0;
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &ivalue, dims, 1, DB_INT) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a double to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const double& value, const string pathName) {
  const string varname = this->setDir(pathName);
  int dims[1] = {1};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &value, dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a string to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const string& value, const string pathName) {
  const char* cvalue = value.c_str();
  int dims[1] = {int(strlen(cvalue))};
  this->write(dims[0], pathName + "/size");
  if (dims[0] > 0) {
    const string varname = this->setDir(pathName + "/value");
    VERIFY2(DBWrite(mFilePtr, varname.c_str(), cvalue, dims, 1, DB_CHAR) == 0,
            "SiloFileIO ERROR: unable to write variable " << pathName);
  }
}

//------------------------------------------------------------------------------
// Write a Dim<1>::Vector to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<1>::Vector& value, const string pathName) {
  const string varname = this->setDir(pathName);
  int dims[1] = {int(Dim<1>::Vector::numElements)};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<1>::Tensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<1>::Tensor& value, const string pathName) {
  const string varname = this->setDir(pathName);
  int dims[1] = {int(Dim<1>::Tensor::numElements)};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<1>::SymTensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<1>::SymTensor& value, const string pathName) {
  const string varname = this->setDir(pathName);
  int dims[1] = {int(Dim<1>::SymTensor::numElements)};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<1>::ThirdRankTensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<1>::ThirdRankTensor& value, const string pathName) {
  const string varname = this->setDir(pathName);
  int dims[1] = {int(Dim<1>::ThirdRankTensor::numElements)};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<2>::Vector to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<2>::Vector& value, const string pathName) {
  const string varname = this->setDir(pathName);
  int dims[1] = {int(Dim<2>::Vector::numElements)};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<2>::Tensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<2>::Tensor& value, const string pathName) {
  const string varname = this->setDir(pathName);
  int dims[1] = {int(Dim<2>::Tensor::numElements)};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<2>::SymTensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<2>::SymTensor& value, const string pathName) {
  const string varname = this->setDir(pathName);
  int dims[1] = {int(Dim<2>::SymTensor::numElements)};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<2>::ThirdRankTensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<2>::ThirdRankTensor& value, const string pathName) {
  const string varname = this->setDir(pathName);
  int dims[1] = {int(Dim<2>::ThirdRankTensor::numElements)};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<3>::Vector to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<3>::Vector& value, const string pathName) {
  const string varname = this->setDir(pathName);
  int dims[1] = {int(Dim<3>::Vector::numElements)};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<3>::Tensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<3>::Tensor& value, const string pathName) {
  const string varname = this->setDir(pathName);
  int dims[1] = {int(Dim<3>::Tensor::numElements)};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<3>::SymTensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<3>::SymTensor& value, const string pathName) {
  const string varname = this->setDir(pathName);
  int dims[1] = {int(Dim<3>::SymTensor::numElements)};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<3>::ThirdRankTensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<3>::ThirdRankTensor& value, const string pathName) {
  const string varname = this->setDir(pathName);
  int dims[1] = {int(Dim<3>::ThirdRankTensor::numElements)};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Read an unsigned from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(unsigned& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  // VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &value) == 0,
  //         "SiloFileIO ERROR: unable to read variable " << pathName);
  value = *static_cast<unsigned*>(DBGetVar(mFilePtr, varname.c_str()));
}

//------------------------------------------------------------------------------
// Read an int from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(int& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  // CHECK2(DBReadVar(mFilePtr, varname.c_str(), &value) == 0,
  //         "SiloFileIO ERROR: unable to read variable " << pathName);
  value = *static_cast<int*>(DBGetVar(mFilePtr, varname.c_str()));
}

//------------------------------------------------------------------------------
// Read a bool from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(bool& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  const int ivalue = *static_cast<int*>(DBGetVar(mFilePtr, varname.c_str()));
  // VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &ivalue) == 0,
  //         "SiloFileIO ERROR: unable to read variable " << pathName);
  value = (ivalue == 1 ? true : false);
}

//------------------------------------------------------------------------------
// Read a double from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(double& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  value = *static_cast<double*>(DBGetVar(mFilePtr, varname.c_str()));
  // VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &value) == 0,
  //         "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a string from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(string& value, const string pathName) const {
  int valsize;
  this->read(valsize, pathName + "/size");
  if (valsize == 0) {
    value = "";
  } else {
    const string varname = this->setDir(pathName + "/value");
    // char* cvalue = (char*) DBGetVar(mFilePtr, varname.c_str());
    // VERIFY2(cvalue != NULL,
    //         "SiloFileIO ERROR: unable to read variable " << pathName);
    // value = string(cvalue);
    // CHECK2(value.size() == valsize, value << " " << value.size() << " " << valsize);
    // free(cvalue);
    char cvalue[valsize + 1];  // Do we need to allow space for trailing null?
    VERIFY2(DBReadVar(mFilePtr, varname.c_str(), cvalue) == 0,
            "SiloFileIO ERROR: failed to read variable " << pathName);
    value = string(&cvalue[0], &cvalue[valsize]);
  }
}

//------------------------------------------------------------------------------
// Read a Dim<1>::Vector from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<1>::Vector& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<1>::Tensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<1>::Tensor& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<1>::SymTensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<1>::SymTensor& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<1>::ThirdRankTensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<1>::ThirdRankTensor& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<2>::Vector from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<2>::Vector& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<2>::Tensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<2>::Tensor& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<2>::SymTensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<2>::SymTensor& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<2>::ThirdRankTensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<2>::ThirdRankTensor& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<3>::Vector from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<3>::Vector& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<3>::Tensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<3>::Tensor& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<3>::SymTensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<3>::SymTensor& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<3>::ThirdRankTensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<3>::ThirdRankTensor& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

#ifdef SPHERAL1D
//------------------------------------------------------------------------------
// Write a Field<Dim<1>, Dim<1>::Scalar> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<1>, Dim<1>::Scalar>& value, const string pathName) {
  this->writeField(value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, Dim<1>::Vector> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<1>, Dim<1>::Vector>& value, const string pathName) {
  this->writeField(value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, Dim<1>::Tensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<1>, Dim<1>::Tensor>& value, const string pathName) {
  this->writeField(value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, Dim<1>::SymTensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<1>, Dim<1>::SymTensor>& value, const string pathName) {
  this->writeField(value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, Dim<1>::ThirdRankTensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<1>, Dim<1>::ThirdRankTensor>& value, const string pathName) {
  this->writeField(value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, int> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<1>, int>& value, const string pathName) {
  this->writeField(value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, Dim<1>::Scalar> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<1>, Dim<1>::Scalar>& value, const string pathName) const {
  this->readField(value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, Dim<1>::Vector> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<1>, Dim<1>::Vector>& value, const string pathName) const {
  this->readField(value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, Dim<1>::Tensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<1>, Dim<1>::Tensor>& value, const string pathName) const {
  this->readField(value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, Dim<1>::SymTensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<1>, Dim<1>::SymTensor>& value, const string pathName) const {
  this->readField(value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, Dim<1>::ThirdRankTensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<1>, Dim<1>::ThirdRankTensor>& value, const string pathName) const {
  this->readField(value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, int> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<1>, int>& value, const string pathName) const {
  this->readField(value, pathName);
}
#endif

#ifdef SPHERAL2D
//------------------------------------------------------------------------------
// Write a Field<Dim<2>, Dim<2>::Scalar> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<2>, Dim<2>::Scalar>& value, const string pathName) {
  this->writeField(value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, Dim<2>::Vector> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<2>, Dim<2>::Vector>& value, const string pathName) {
  this->writeField(value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, Dim<2>::Tensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<2>, Dim<2>::Tensor>& value, const string pathName) {
  this->writeField(value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, Dim<2>::SymTensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<2>, Dim<2>::SymTensor>& value, const string pathName) {
  this->writeField(value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, Dim<2>::ThirdRankTensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<2>, Dim<2>::ThirdRankTensor>& value, const string pathName) {
  this->writeField(value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, int> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<2>, int>& value, const string pathName) {
  this->writeField(value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, Dim<2>::Scalar> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<2>, Dim<2>::Scalar>& value, const string pathName) const {
  this->readField(value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, Dim<2>::Vector> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<2>, Dim<2>::Vector>& value, const string pathName) const {
  this->readField(value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, Dim<2>::Tensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<2>, Dim<2>::Tensor>& value, const string pathName) const {
  this->readField(value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, Dim<2>::SymTensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<2>, Dim<2>::SymTensor>& value, const string pathName) const {
  this->readField(value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, Dim<2>::ThirdRankTensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<2>, Dim<2>::ThirdRankTensor>& value, const string pathName) const {
  this->readField(value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, int> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<2>, int>& value, const string pathName) const {
  this->readField(value, pathName);
}
#endif

#ifdef SPHERAL3D
//------------------------------------------------------------------------------
// Write a Field<Dim<3>, Dim<3>::Scalar> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<3>, Dim<3>::Scalar>& value, const string pathName) {
  this->writeField(value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, Dim<3>::Vector> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<3>, Dim<3>::Vector>& value, const string pathName) {
  this->writeField(value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, Dim<3>::Tensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<3>, Dim<3>::Tensor>& value, const string pathName) {
  this->writeField(value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, Dim<3>::SymTensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<3>, Dim<3>::SymTensor>& value, const string pathName) {
  this->writeField(value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, Dim<3>::ThirdRankTensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<3>, Dim<3>::ThirdRankTensor>& value, const string pathName) {
  this->writeField(value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, int> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<3>, int>& value, const string pathName) {
  this->writeField(value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, Dim<3>::Scalar> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<3>, Dim<3>::Scalar>& value, const string pathName) const {
  this->readField(value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, Dim<3>::Vector> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<3>, Dim<3>::Vector>& value, const string pathName) const {
  this->readField(value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, Dim<3>::Tensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<3>, Dim<3>::Tensor>& value, const string pathName) const {
  this->readField(value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, Dim<3>::SymTensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<3>, Dim<3>::SymTensor>& value, const string pathName) const {
  this->readField(value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, Dim<3>::ThirdRankTensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<3>, Dim<3>::ThirdRankTensor>& value, const string pathName) const {
  this->readField(value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, int> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<3>, int>& value, const string pathName) const {
  this->readField(value, pathName);
}
#endif

//------------------------------------------------------------------------------
// Set the directory in the current silo file to the path portion of a pathName,
// and return the variable name section.
//------------------------------------------------------------------------------
// Non-const version creates directory.
string 
SiloFileIO::setDir(const string& ipathName) {
  REQUIRE(mFileOpen and mFilePtr != 0);

  // Mangle the path into something acceptable to silo.
  string pathName = SILO_mangle(ipathName);

  // We always start from the top.
  VERIFY2(DBSetDir(mFilePtr, "/") == 0,
          "SiloFileIO ERROR: unable to change to path /");

  // Split up the path into directories and the var name.
  vector<string> components = this->splitPathComponents(pathName);
  CHECK(components.size() > 0);

  // Set the path and return just the variable name.
  for (unsigned i = 0; i != components.size() - 1; ++i) {
    const string dirName = components[i];

    // Check if this directory already exists or not.
    DBtoc* toc = DBGetToc(mFilePtr);
    bool exists = false;
    size_t j = 0;
    while (not exists and j != toc->ndir) {
      exists = (dirName == string(toc->dir_names[j]));
      ++j;
    }
    if (not exists) {
      VERIFY2(DBMkDir(mFilePtr, dirName.c_str()) == 0,
              "SiloFileIO ERROR: unable to create path " << dirName);
    }
    VERIFY2(DBSetDir(mFilePtr, dirName.c_str()) == 0,
            "SiloFileIO ERROR: unable to change to path " << dirName);
  }

  return components.back();
}

// const version only changes directory.
string 
SiloFileIO::setDir(const string& ipathName) const {
  REQUIRE(mFileOpen and mFilePtr != 0);

  // Mangle the path into something acceptable to silo.
  string pathName = SILO_mangle(ipathName);

  // We always start from the top.
  VERIFY2(DBSetDir(mFilePtr, "/") == 0,
          "SiloFileIO ERROR: unable to change to path /");

  // If there is no absolute path specified, we just go with the cwd.
  const size_t i = pathName.find_last_of("/");
  if (i >= pathName.size()) {
    return pathName;
  }

  // Otherwise set the path and return just the variable name.
  const string dirName = pathName.substr(0,i);
  const string varName = pathName.substr(i+1);
  VERIFY2(DBSetDir(mFilePtr, dirName.c_str()) == 0,
          "SiloFileIO ERROR: unable to change to path " << dirName);
  return varName;
}

//------------------------------------------------------------------------------
// Worker method to write a generic Field
//------------------------------------------------------------------------------
template<typename FieldType>
void
SiloFileIO::
writeField(const FieldType& value, 
           const string pathName) {
  const unsigned n = value.numInternalElements();
  std::vector<int> packIndices(n);
  for (auto i = 0; i < n; ++i) packIndices[i] = i;
  const auto buffer = packFieldValues(value, packIndices);
  this->write(value.name(), pathName + "/name");
  this->write(n, pathName + "/size");
  this->write(std::string(buffer.begin(), buffer.end()), pathName + "/values");
}
  
template<typename FieldType>
void
SiloFileIO::
readField(FieldType& value, 
          const string pathName) const {
  unsigned n;
  this->read(n, pathName + "/size");
  CHECK(value.numInternalElements() == n);
  std::string name;
  this->read(name, pathName + "/name");
  std::vector<int> packIndices(n);
  for (auto i = 0; i < n; ++i) packIndices[i] = i;
  std::string stuff;
  this->read(stuff, pathName + "/values");
  unpackFieldValues(value, packIndices, std::vector<char>(stuff.begin(), stuff.end()));
  value.name(name);
}
  
}
