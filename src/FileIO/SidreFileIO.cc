//---------------------------------Spheral++----------------------------------//
// SidreFileIO -- Provide the interface to sidre file objects.
//
// Created by Mikhail Zakharchanka, 11/4/2021
//----------------------------------------------------------------------------//
#include "SidreFileIO.hh"

namespace Spheral
{

//------------------------------------------------------------------------------
// Write an unsigned to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const unsigned& value, const string pathName)
{
  m_datastore_ptr->getRoot()->createViewScalar(pathName, value);
}

//------------------------------------------------------------------------------
// Write a size_t to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const size_t& value, const string pathName)
{
  m_datastore_ptr->getRoot()->createViewScalar(pathName, value);
}

//------------------------------------------------------------------------------
// Write an int to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const int& value, const string pathName)
{
  m_datastore_ptr->getRoot()->createViewScalar(pathName, value);
}

//------------------------------------------------------------------------------
// Write a bool to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const bool& value, const string pathName)
{
  m_datastore_ptr->getRoot()->createViewScalar(pathName, value);
}

//------------------------------------------------------------------------------
// Write a double to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const double& value, const string pathName)
{
  m_datastore_ptr->getRoot()->createViewScalar(pathName, value);
}

//------------------------------------------------------------------------------
// Write a string to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const string& value, const string pathName) 
{
  m_datastore_ptr->getRoot()->createViewString(pathName, value);
}

//------------------------------------------------------------------------------
// Write a vector<int> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const std::vector<int>& value, const string pathName)
{
  m_datastore_ptr->getRoot()->createView(pathName, int, value.size(), (void*)&value);
}

//------------------------------------------------------------------------------
// Write a vector<double> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const std::vector<double>& value, const string pathName)
{
  m_datastore_ptr->getRoot()->createView(pathName, double, value.size(), (void*)&value);
}

//------------------------------------------------------------------------------
// Write a vector<string> to the file.
//------------------------------------------------------------------------------
void SidreFileIO::write(const std::vector<string>& value, const string pathName)
{
  
}


//------------------------------------------------------------------------------
// Read an unsigned from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(unsigned& value, const string pathName) const
{
  value = m_datastore_ptr->getRoot()->getView(pathName)->getScalar();
}

//------------------------------------------------------------------------------
// Read a size_t from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(size_t& value, const string pathName) const
{
  value = m_datastore_ptr->getRoot()->getView(pathName)->getScalar();
}

//------------------------------------------------------------------------------
// Read an int to the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(int& value, const string pathName) const
{
  value = m_datastore_ptr->getRoot()->getView(pathName)->getScalar();
}

//------------------------------------------------------------------------------
// Read a bool from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(bool& value, const string pathName) const
{
  value = m_datastore_ptr->getRoot()->getView(pathName)->getScalar();
}

//------------------------------------------------------------------------------
// Read a double from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(double& value, const string pathName) const
{
  value = m_datastore_ptr->getRoot()->getView(pathName)->getScalar();
}

//------------------------------------------------------------------------------
// Read a string from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(string& value, const string pathName) const
{
  value = m_datastore_ptr->getRoot()->getView(pathName)->getString();

//------------------------------------------------------------------------------
// Read a vector<int> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(std::vector<int>& value, const string pathName) const
{
  value = m_datastore_ptr->getRoot()->getView(pathName)->getData();
}

//------------------------------------------------------------------------------
// Read a vector<double> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(std::vector<double>& value, const string pathName) const
{
  value = m_datastore_ptr->getRoot()->getView(pathName)->getData();
}

//------------------------------------------------------------------------------
// Read a vector<string> from the file.
//------------------------------------------------------------------------------
void SidreFileIO::read(vector<string>& value, const string pathName) const
{
  
}

}