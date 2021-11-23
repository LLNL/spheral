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

  // std::string fullFileName = fileName;
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

//------------------------------------------------------------------------------
// Write a size_t to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const size_t& value, const std::string pathName)
{
  mFilePtr->getRoot()->createViewScalar(pathName, value);
}

//------------------------------------------------------------------------------
// Write an int to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const int& value, const std::string pathName)
{
  mFilePtr->getRoot()->createViewScalar(pathName, value);
}

// ------------------------------------------------------------------------------
// Write a bool to the file.
// ------------------------------------------------------------------------------
void SidreFileIO::write(const bool& value, const std::string pathName)
{
  mFilePtr->getRoot()->createViewScalar(pathName, value);
}

//------------------------------------------------------------------------------
// Write a double to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const double& value, const std::string pathName)
{
  mFilePtr->getRoot()->createViewScalar(pathName, value);
}

//------------------------------------------------------------------------------
// Write a std::string to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const std::string& value, const std::string pathName) 
{
  mFilePtr->getRoot()->createViewString(pathName, value);
}

//------------------------------------------------------------------------------
// Write a vector<int> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const std::vector<int>& value, const std::string pathName)
{
  // mFilePtr->getRoot()->createView(pathName, int, value.size(), (void*)&value);
}

//------------------------------------------------------------------------------
// Write a vector<double> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const std::vector<double>& value, const std::string pathName)
{
  // mFilePtr->getRoot()->createView(pathName, double, value.size(), (void*)&value);
}

//------------------------------------------------------------------------------
// Write a vector<std::string> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const std::vector<std::string>& value, const std::string pathName)
{
  
}

//------------------------------------------------------------------------------
// Write a Dim<1>::Vector to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<1>::Vector& value, const std::string pathName)
{

}

//------------------------------------------------------------------------------
// Write a Dim<1>::Tensor to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<1>::Tensor& value, const std::string pathName)
{

}

//------------------------------------------------------------------------------
// Write a Dim<1>::SymTensor to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<1>::SymTensor& value, const std::string pathName)
{

}

//------------------------------------------------------------------------------
// Write a Dim<1>::ThirdRankTensor to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<1>::ThirdRankTensor& value, const std::string pathName)
{

}

//------------------------------------------------------------------------------
// Write a Dim<2>::Vector to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<2>::Vector& value, const std::string pathName)
{

}

//------------------------------------------------------------------------------
// Write a Dim<2>::Tensor to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<2>::Tensor& value, const std::string pathName)
{

}

//------------------------------------------------------------------------------
// Write a Dim<2>::SymTensor to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<2>::SymTensor& value, const std::string pathName)
{

}

//------------------------------------------------------------------------------
// Write a Dim<2>::ThirdRankTensor to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<2>::ThirdRankTensor& value, const std::string pathName)
{

}

//------------------------------------------------------------------------------
// Write a Dim<3>::Vector to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<3>::Vector& value, const std::string pathName)
{
}

//------------------------------------------------------------------------------
// Write a Dim<3>::Tensor to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<3>::Tensor& value, const std::string pathName)
{
}

//------------------------------------------------------------------------------
// Write a Dim<3>::SymTensor to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<3>::SymTensor& value, const std::string pathName)
{
}

//------------------------------------------------------------------------------
// Write a Dim<3>::ThirdRankTensor to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const Dim<3>::ThirdRankTensor& value, const std::string pathName)
{

}


// ------------------------------------------------------------------------------
// Read an unsigned from the file.
// ------------------------------------------------------------------------------
void SidreFileIO::read(unsigned& value, const std::string pathName) const
{
  value = mFilePtr->getRoot()->getView(pathName)->getScalar();
}

//------------------------------------------------------------------------------
// Read a size_t from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(size_t& value, const std::string pathName) const
{
  value = mFilePtr->getRoot()->getView(pathName)->getScalar();
}

//------------------------------------------------------------------------------
// Read an int to the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(int& value, const std::string pathName) const
{
  value = mFilePtr->getRoot()->getView(pathName)->getScalar();
}

//------------------------------------------------------------------------------
// Read a bool from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(bool& value, const std::string pathName) const
{
  // value = mFilePtr->getRoot()->getView(pathName)->getScalar();
}

//------------------------------------------------------------------------------
// Read a double from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(double& value, const std::string pathName) const
{
  value = mFilePtr->getRoot()->getView(pathName)->getScalar();
}

//------------------------------------------------------------------------------
// Read a std::string from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(std::string& value, const std::string pathName) const
{
  // value = mFilePtr->getRoot()->getView(pathName)->getString();
}

//------------------------------------------------------------------------------
// Read a vector<int> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(std::vector<int>& value, const std::string pathName) const
{
  // value = mFilePtr->getRoot()->getView(pathName)->getData();
}

//------------------------------------------------------------------------------
// Read a vector<double> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(std::vector<double>& value, const std::string pathName) const
{
  // value = mFilePtr->getRoot()->getView(pathName)->getData();
}

//------------------------------------------------------------------------------
// Read a vector<std::string> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(vector<std::string>& value, const std::string pathName) const
{
  
}

//------------------------------------------------------------------------------
// Read a Dim<1>::Vector from the file.
//------------------------------------------------------------------------------
void
SidreFileIO::read(Dim<1>::Vector& value, const std::string pathName) const {

}

//------------------------------------------------------------------------------
// Read a Dim<1>::Tensor from the file.
//------------------------------------------------------------------------------
void
SidreFileIO::read(Dim<1>::Tensor& value, const std::string pathName) const {

}

//------------------------------------------------------------------------------
// Read a Dim<1>::SymTensor from the file.
//------------------------------------------------------------------------------
void
SidreFileIO::read(Dim<1>::SymTensor& value, const std::string pathName) const {

}

//------------------------------------------------------------------------------
// Read a Dim<1>::ThirdRankTensor from the file.
//------------------------------------------------------------------------------
void
SidreFileIO::read(Dim<1>::ThirdRankTensor& value, const std::string pathName) const {

}

//------------------------------------------------------------------------------
// Read a Dim<2>::Vector from the file.
//------------------------------------------------------------------------------
void
SidreFileIO::read(Dim<2>::Vector& value, const std::string pathName) const {

}

//------------------------------------------------------------------------------
// Read a Dim<2>::Tensor from the file.
//------------------------------------------------------------------------------
void
SidreFileIO::read(Dim<2>::Tensor& value, const std::string pathName) const {

}

//------------------------------------------------------------------------------
// Read a Dim<2>::SymTensor from the file.
//------------------------------------------------------------------------------
void
SidreFileIO::read(Dim<2>::SymTensor& value, const std::string pathName) const {

}

//------------------------------------------------------------------------------
// Read a Dim<2>::ThirdRankTensor from the file.
//------------------------------------------------------------------------------
void
SidreFileIO::read(Dim<2>::ThirdRankTensor& value, const std::string pathName) const {

}

//------------------------------------------------------------------------------
// Read a Dim<3>::Vector from the file.
//------------------------------------------------------------------------------
void
SidreFileIO::read(Dim<3>::Vector& value, const std::string pathName) const {

}

//------------------------------------------------------------------------------
// Read a Dim<3>::Tensor from the file.
//------------------------------------------------------------------------------
void
SidreFileIO::read(Dim<3>::Tensor& value, const std::string pathName) const {

}

//------------------------------------------------------------------------------
// Read a Dim<3>::SymTensor from the file.
//------------------------------------------------------------------------------
void
SidreFileIO::read(Dim<3>::SymTensor& value, const std::string pathName) const {

}

//------------------------------------------------------------------------------
// Read a Dim<3>::ThirdRankTensor from the file.
//------------------------------------------------------------------------------
void
SidreFileIO::read(Dim<3>::ThirdRankTensor& value, const std::string pathName) const {

}

#ifdef SPHERAL1D
//------------------------------------------------------------------------------
// Write a Field<Dim<1>, Dim<1>::Scalar> to the file.
//------------------------------------------------------------------------------
void
SidreFileIO::write(const Field<Dim<1>, Dim<1>::Scalar>& value, const std::string pathName) {
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, Dim<1>::Vector> to the file.
//------------------------------------------------------------------------------
void
SidreFileIO::write(const Field<Dim<1>, Dim<1>::Vector>& value, const std::string pathName) {
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, Dim<1>::Tensor> to the file.
//------------------------------------------------------------------------------
void
SidreFileIO::write(const Field<Dim<1>, Dim<1>::Tensor>& value, const std::string pathName) {
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, Dim<1>::SymTensor> to the file.
//------------------------------------------------------------------------------
void
SidreFileIO::write(const Field<Dim<1>, Dim<1>::SymTensor>& value, const std::string pathName) {
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, Dim<1>::ThirdRankTensor> to the file.
//------------------------------------------------------------------------------
void
SidreFileIO::write(const Field<Dim<1>, Dim<1>::ThirdRankTensor>& value, const std::string pathName) {
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, int> to the file.
//------------------------------------------------------------------------------
void
SidreFileIO::write(const Field<Dim<1>, int>& value, const std::string pathName) {
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, unsigned> to the file.
//------------------------------------------------------------------------------
void
SidreFileIO::write(const Field<Dim<1>, unsigned>& value, const std::string pathName) {
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, Dim<1>::Scalar> from the file.
//------------------------------------------------------------------------------
void
SidreFileIO::read(Field<Dim<1>, Dim<1>::Scalar>& value, const std::string pathName) const {
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, Dim<1>::Vector> from the file.
//------------------------------------------------------------------------------
void
SidreFileIO::read(Field<Dim<1>, Dim<1>::Vector>& value, const std::string pathName) const {
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, Dim<1>::Tensor> from the file.
//------------------------------------------------------------------------------
void
SidreFileIO::read(Field<Dim<1>, Dim<1>::Tensor>& value, const std::string pathName) const {
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, Dim<1>::SymTensor> from the file.
//------------------------------------------------------------------------------
void
SidreFileIO::read(Field<Dim<1>, Dim<1>::SymTensor>& value, const std::string pathName) const {
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, Dim<1>::ThirdRankTensor> from the file.
//------------------------------------------------------------------------------
void
SidreFileIO::read(Field<Dim<1>, Dim<1>::ThirdRankTensor>& value, const std::string pathName) const {
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, int> from the file.
//------------------------------------------------------------------------------
void
SidreFileIO::read(Field<Dim<1>, int>& value, const std::string pathName) const {
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, unsigned> from the file.
//------------------------------------------------------------------------------
void
SidreFileIO::read(Field<Dim<1>, unsigned>& value, const std::string pathName) const {
}
#endif

#ifdef SPHERAL2D
//------------------------------------------------------------------------------
// Write a Field<Dim<2>, Dim<2>::Scalar> to the file.
//------------------------------------------------------------------------------
void
SidreFileIO::write(const Field<Dim<2>, Dim<2>::Scalar>& value, const std::string pathName) {
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, Dim<2>::Vector> to the file.
//------------------------------------------------------------------------------
void
SidreFileIO::write(const Field<Dim<2>, Dim<2>::Vector>& value, const std::string pathName) {
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, Dim<2>::Tensor> to the file.
//------------------------------------------------------------------------------
void
SidreFileIO::write(const Field<Dim<2>, Dim<2>::Tensor>& value, const std::string pathName) {
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, Dim<2>::SymTensor> to the file.
//------------------------------------------------------------------------------
void
SidreFileIO::write(const Field<Dim<2>, Dim<2>::SymTensor>& value, const std::string pathName) {
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, Dim<2>::ThirdRankTensor> to the file.
//------------------------------------------------------------------------------
void
SidreFileIO::write(const Field<Dim<2>, Dim<2>::ThirdRankTensor>& value, const std::string pathName) {
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, int> to the file.
//------------------------------------------------------------------------------
void
SidreFileIO::write(const Field<Dim<2>, int>& value, const std::string pathName) {
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, unsigned> to the file.
//------------------------------------------------------------------------------
void
SidreFileIO::write(const Field<Dim<2>, unsigned>& value, const std::string pathName) {
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, Dim<2>::Scalar> from the file.
//------------------------------------------------------------------------------
void
SidreFileIO::read(Field<Dim<2>, Dim<2>::Scalar>& value, const std::string pathName) const {
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, Dim<2>::Vector> from the file.
//------------------------------------------------------------------------------
void
SidreFileIO::read(Field<Dim<2>, Dim<2>::Vector>& value, const std::string pathName) const {
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, Dim<2>::Tensor> from the file.
//------------------------------------------------------------------------------
void
SidreFileIO::read(Field<Dim<2>, Dim<2>::Tensor>& value, const std::string pathName) const {
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, Dim<2>::SymTensor> from the file.
//------------------------------------------------------------------------------
void
SidreFileIO::read(Field<Dim<2>, Dim<2>::SymTensor>& value, const std::string pathName) const {

}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, Dim<2>::ThirdRankTensor> from the file.
//------------------------------------------------------------------------------
void
SidreFileIO::read(Field<Dim<2>, Dim<2>::ThirdRankTensor>& value, const std::string pathName) const {
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, int> from the file.
//------------------------------------------------------------------------------
void
SidreFileIO::read(Field<Dim<2>, int>& value, const std::string pathName) const {
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, unsigned> from the file.
//------------------------------------------------------------------------------
void
SidreFileIO::read(Field<Dim<2>, unsigned>& value, const std::string pathName) const {
}
#endif

#ifdef SPHERAL3D
//------------------------------------------------------------------------------
// Write a Field<Dim<3>, Dim<3>::Scalar> to the file.
//------------------------------------------------------------------------------
void
SidreFileIO::write(const Field<Dim<3>, Dim<3>::Scalar>& value, const std::string pathName) {
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, Dim<3>::Vector> to the file.
//------------------------------------------------------------------------------
void
SidreFileIO::write(const Field<Dim<3>, Dim<3>::Vector>& value, const std::string pathName) {
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, Dim<3>::Tensor> to the file.
//------------------------------------------------------------------------------
void
SidreFileIO::write(const Field<Dim<3>, Dim<3>::Tensor>& value, const std::string pathName) {
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, Dim<3>::SymTensor> to the file.
//------------------------------------------------------------------------------
void
SidreFileIO::write(const Field<Dim<3>, Dim<3>::SymTensor>& value, const std::string pathName) {
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, Dim<3>::ThirdRankTensor> to the file.
//------------------------------------------------------------------------------
void
SidreFileIO::write(const Field<Dim<3>, Dim<3>::ThirdRankTensor>& value, const std::string pathName) {
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, int> to the file.
//------------------------------------------------------------------------------
void
SidreFileIO::write(const Field<Dim<3>, int>& value, const std::string pathName) {
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, unsigned> to the file.
//------------------------------------------------------------------------------
void
SidreFileIO::write(const Field<Dim<3>, unsigned>& value, const std::string pathName) {
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, Dim<3>::Scalar> from the file.
//------------------------------------------------------------------------------
void
SidreFileIO::read(Field<Dim<3>, Dim<3>::Scalar>& value, const std::string pathName) const {
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, Dim<3>::Vector> from the file.
//------------------------------------------------------------------------------
void
SidreFileIO::read(Field<Dim<3>, Dim<3>::Vector>& value, const std::string pathName) const {
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, Dim<3>::Tensor> from the file.
//------------------------------------------------------------------------------
void
SidreFileIO::read(Field<Dim<3>, Dim<3>::Tensor>& value, const std::string pathName) const {
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, Dim<3>::SymTensor> from the file.
//------------------------------------------------------------------------------
void
SidreFileIO::read(Field<Dim<3>, Dim<3>::SymTensor>& value, const std::string pathName) const {
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, Dim<3>::ThirdRankTensor> from the file.
//------------------------------------------------------------------------------
void
SidreFileIO::read(Field<Dim<3>, Dim<3>::ThirdRankTensor>& value, const std::string pathName) const {
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, int> from the file.
//------------------------------------------------------------------------------
void
SidreFileIO::read(Field<Dim<3>, int>& value, const std::string pathName) const {
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, unsigned> from the file.
//------------------------------------------------------------------------------
void
SidreFileIO::read(Field<Dim<3>, unsigned>& value, const std::string pathName) const {
}
#endif


}