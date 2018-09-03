//---------------------------------Spheral++----------------------------------//
// DbFileIO -- Provide the database interface for I/O
//
// Created by Jeffrey Johnson and Charlie Crabb, Mon Aug 27 11:01:58 PDT 2001
//----------------------------------------------------------------------------//

#include "DbFileIO.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

// #ifdef GNUCXX
// #include <strstream>
// #define STRINGSTREAM ostrstream
// #else
#include <sstream>
#define STRINGSTREAM stringstream
// #endif

namespace Spheral {


#ifdef USE_POSTGRES
//------------------------------------------------------------------------------
// Static database pointer needs to be defined.
//------------------------------------------------------------------------------
template<typename Dimension>
PgDatabase*
DbFileIO<Dimension>::sDatabasePtr = 0;
#endif

//------------------------------------------------------------------------------
// Empty Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
DbFileIO<Dimension>::DbFileIO():
FileIO<Dimension>() {

  // The connection will be made the first time it is accessed, per the 
  // mConnection() method.
}

//------------------------------------------------------------------------------
// Construct with the given file name and access type.
//------------------------------------------------------------------------------
template<typename Dimension>
DbFileIO<Dimension>::DbFileIO(const string filename) :
FileIO<Dimension>(filename, ReadWrite) {

  // The connection will be made the first time it is accessed, per the 
  // mConnection() method.
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
DbFileIO<Dimension>::~DbFileIO() {

  // Close the database connection.
  close();
}

//------------------------------------------------------------------------------
// Open the database given a "filename" and an access type.
// NOTE: The access typ e is ignored. :-P
//------------------------------------------------------------------------------
template<typename Dimension>
void
DbFileIO<Dimension>::open(const string fileName, AccessType access) {
#ifdef USE_POSTGRES
  mFileName = fileName;
  mAccess = access;
  mFileOpen = ! mConnection().ConnectionBad();
#endif
}

//------------------------------------------------------------------------------
// Close the database.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DbFileIO<Dimension>::close() {
#ifdef USE_POSTGRES
  mFileOpen = false;
  if (sDatabasePtr)
  {
    delete sDatabasePtr;
  } // end if
#endif
}

#ifdef USE_POSTGRES
//------------------------------------------------------------------------------
// Return a valid database connection.
//------------------------------------------------------------------------------
template <typename Dimension>
PgDatabase&
DbFileIO<Dimension>::mConnection() const {
  if (!sDatabasePtr)
  {
    // Make a connection string with which to connect to the database.
    std::strstream connectStr;
    connectStr //<< "host = localhost "
      << "dbname = spheral" << '\0';
    sDatabasePtr = new PgDatabase(
        const_cast<const char*>(connectStr.str()));

    // Check the connection.
    if (sDatabasePtr->ConnectionBad())
    {
      cerr << "Error: " << sDatabasePtr->ErrorMessage() << endl;
    } // end if
  } // end if 

  return *sDatabasePtr;
}

//------------------------------------------------------------------------------
// Return a string containing query results for a given strstream containing 
// a query.
//------------------------------------------------------------------------------
template <typename Dimension>
const char*
DbFileIO<Dimension>::mSingleResult(strstream& query) const {
  // Make sure that the strstream is null-terminated.
  query << '\0';
  // Get results from the query.
  int status = 
    mConnection().ExecTuplesOk(const_cast<const char*>(query.str()));

  if (!status)
  {
    cerr << "Error: " << mConnection().ErrorMessage() << endl;
  } // end if

  // Load the resulting value into value.
  if (mConnection().Tuples())
  {
    return mConnection().GetValue(0, 0);
  } // end if
  else
  {
    return "";
  } // end else
} // end mSingleResult

//------------------------------------------------------------------------------
// Send a query to the database.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DbFileIO<Dimension>::mSendQuery(strstream& query) const {
  // Make sure that the strstream is null-terminated.
  query << '\0';
  int status = 
    mConnection().ExecCommandOk(const_cast<const char*>(query.str()));
  // HACK: We really shouldn't have to check the length of the error message,
  //       but Postgres is being stupid and setting status to nonzero even
  //       though it seems that no error is occurring. :-/
  if (!status && strlen(mConnection().ErrorMessage()))
  {
    cerr << "Error: " << mConnection().ErrorMessage() << endl;
  } // end if
} // end mSendQuery
#endif

//------------------------------------------------------------------------------
// Write an int to the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DbFileIO<Dimension>::write(const int& value, const string pathName) {
#ifdef USE_POSTGRES
  // Form a command to write our value to the database.
  strstream query;
  query << "SELECT write_integer('" << mFileName << "', '" << pathName << "', " << value << ");" << '\0';

  // Send the query to the database.
  mSendQuery(query);
#endif
}

//------------------------------------------------------------------------------
// Write a bool to the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DbFileIO<Dimension>::write(const bool& value, const string pathName) {
#ifdef USE_POSTGRES
  // Form a command to write our value to the database.
  strstream query;
  query << "SELECT write_boolean('" << mFileName << "', '" << pathName << "', " << value << ");" << '\0';

  // Send the query to the database.
  mSendQuery(query);
#endif
}

//------------------------------------------------------------------------------
// Write a Scalar to the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DbFileIO<Dimension>::
write(const typename Dimension::Scalar& value, const string pathName) {
#ifdef USE_POSTGRES
  // Form a command to write our value to the database.
  strstream query;
  query << "SELECT write_scalar('" << mFileName << "', '" << pathName << "', " << value << ");" << '\0';

  // Send the query to the database.
  mSendQuery(query);
#endif
}

//------------------------------------------------------------------------------
// Write a Vector to the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DbFileIO<Dimension>::
write(const typename Dimension::Vector& value, const string pathName) {
#ifdef USE_POSTGRES
  // Construct a string representation of the vector as follows:
  // (v1, v2, ..., vN) for an N-dimensional vector v.
  strstream valueStr;
  valueStr << '(';
  for (int i = 0; i < Dimension::nDim; ++i)
  {
    if (i > 0)
    {
      valueStr << ',';
    } // end if
    valueStr << value(i);
  } // end for
  valueStr << ')' << '\0';

  // Form a command to write our value to the database.
  strstream query;
  query << "SELECT write_vector" << Dimension::nDim << "d('" << mFileName << "', '" << pathName << "', '" << valueStr.str()
    << "');" << '\0';

  // Send the query to the database.
  mSendQuery(query);
#endif
}

//------------------------------------------------------------------------------
// Write a Tensor to the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DbFileIO<Dimension>::
write(const typename Dimension::Tensor& value, const string pathName) {
#ifdef USE_POSTGRES
  // Construct a string representation of the tensor as follows:
  // (t11, t12, ..., t1N, t21, ..., t2N, ..., tNN) 
  // for a general N-dimensional tensor t.
  strstream valueStr;
  valueStr << '(';
  for (int i = 0; i < Dimension::nDim; ++i)
  {
    for (int j = 0; j < Dimension::nDim; ++j)
    {
      if ((i > 0) || (j > 0))
      {
        valueStr << ',';
      } // end if
      valueStr << value(i,j);
    } // end for
  } // end for
  valueStr << ')' << '\0';

  // Form a command to write our value to the database.
  strstream query;
  query << "SELECT write_tensor" << Dimension::nDim << "d('" << mFileName << "', '" << pathName << "', '" << valueStr.str()
    << "');" << '\0';

  // Send the query to the database.
  mSendQuery(query);
#endif
}

//------------------------------------------------------------------------------
// Write a SymTensor to the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DbFileIO<Dimension>::
write(const typename Dimension::SymTensor& value, const string pathName) {
#ifdef USE_POSTGRES
  // Construct a string representation of the tensor as in the example below:
  // In 3 dimensions:
  // (t11, t12, t13, t22, t23, t33)
  strstream valueStr;
  valueStr << '(';
  for (int i = 0; i < Dimension::nDim; ++i)
  {
    for (int j = 0; j <= i; ++j)
    {
      if ((i > 0) || (j > 0))
      {
        valueStr << ',';
      } // end if
      valueStr << value(j,i);
    } // end for
  } // end for
  valueStr << ')' << '\0';

  // Form a command to write our value to the database.
  strstream query;
  query << "SELECT write_symtensor" << Dimension::nDim << "d('" << mFileName << "', '" << pathName << "', '" << valueStr.str()
    << "');" << '\0';

  // Send the query to the database.
  mSendQuery(query);
#endif
}

//------------------------------------------------------------------------------
// Write a string to the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DbFileIO<Dimension>::
write(const string& value, const string pathName) {
#ifdef USE_POSTGRES
  // Form a command to write our value to the database.
  strstream query;
  query << "SELECT write_string('" << mFileName << "', '" << pathName << "', '" << value << "');" << '\0';

  // Send the query to the database.
  mSendQuery(query);
#endif
}

//------------------------------------------------------------------------------
// Read an int from the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DbFileIO<Dimension>::read(int& value, const string pathName) const {
#ifdef USE_POSTGRES
  // Set up the query.
  strstream query;
  query << "SELECT read_integer('" << mFileName << "', '"
    << pathName << "');" << '\0';

  // Get results from the query.
  value = atoi(mSingleResult(query));
#endif
}

//------------------------------------------------------------------------------
// Read a bool from the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DbFileIO<Dimension>::read(bool& value, const string pathName) const {
#ifdef USE_POSTGRES
  // Set up the query.
  strstream query;
  query << "SELECT read_boolean('" << mFileName << "', '"
    << pathName << "');" << '\0';

  // Get results from the query.
  value = static_cast<bool>(atoi(mSingleResult(query)));
#endif
}

//------------------------------------------------------------------------------
// Read a Scalar from the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DbFileIO<Dimension>::read(typename Dimension::Scalar& value, const string pathName) const {
#ifdef USE_POSTGRES
  // Set up the query.
  strstream query;
  query << "SELECT read_scalar('" << mFileName << "', '"
    << pathName << "');" << '\0';

  // Get results from the query.
  value = atof(mSingleResult(query));
#endif
}

//------------------------------------------------------------------------------
// Read a Vector from the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DbFileIO<Dimension>::read(typename Dimension::Vector& value, const string pathName) const {
#ifdef USE_POSTGRES
  // Set up the query.
  strstream query;
  query << "SELECT read_vector" << Dimension::nDim << "d('" << mFileName << "', '" << pathName << "');" << '\0';

  // Get the result string from the query and parse it into value.
  string valueString = mSingleResult(query);
  mParseVector(value, valueString);
#endif
}

//------------------------------------------------------------------------------
// Read a Tensor from the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DbFileIO<Dimension>::read(typename Dimension::Tensor& value, const string pathName) const {
#ifdef USE_POSTGRES
  // Set up the query.
  strstream query;
  query << "SELECT read_tensor" << Dimension::nDim << "d('" << mFileName << "', '" << pathName << "');" << '\0';

  // Get the result string from the query and parse it into value.
  string valueString = mSingleResult(query);
  mParseTensor(value, valueString);
#endif
}

//------------------------------------------------------------------------------
// Read a SymTensor from the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DbFileIO<Dimension>::read(typename Dimension::SymTensor& value, const string pathName) const {
#ifdef USE_POSTGRES
  // Set up the query.
  strstream query;
  query << "SELECT read_symtensor" << Dimension::nDim << "d('" << mFileName << "', '" << pathName << "');" << '\0';

  // Get the result string from the query and parse it into value.
  string valueString = mSingleResult(query);
  mParseSymTensor(value, valueString);
#endif
}

//------------------------------------------------------------------------------
// Read a string from the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DbFileIO<Dimension>::read(string& value, const string pathName) const {
#ifdef USE_POSTGRES
  // Set up the query.
  strstream query;
  query << "SELECT read_string('" << mFileName << "', '" << pathName << "');" << '\0';

  // Get the result string from the query.
  value = mSingleResult(query);
#endif
}

//------------------------------------------------------------------------------
// Write a Scalar Field to the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DbFileIO<Dimension>::
write(const Field<Dimension, typename Dimension::Scalar>& value, const string pathName) {
#ifdef USE_POSTGRES
  // Delete any existing information about this field.
  strstream query;
  query << "SELECT delete_scalarfield('" << mFileName << "', '"
    << pathName << "');" << '\0';
  mSendQuery(query);

  // Insert the new information about this field.
  strstream query2;
  query2 << "INSERT INTO scalarfields (tag_id, name_id, elements) "
    << "     VALUES (get_tag_id('" << mFileName 
    << "'), get_name_id('" << pathName << "'), '{";
  for (int i = 0; i < value.numElements(); ++i)
  {
    if (i != 0)
    {
      query2 << ",";
    } // end if
    query2 << value(i);
  } // end for
  query2 << "}');" << '\0';
  mSendQuery(query2);
#endif
}

//------------------------------------------------------------------------------
// Write a Vector Field to the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DbFileIO<Dimension>::
write(const Field<Dimension, typename Dimension::Vector>& value, const string pathName) {
#ifdef USE_POSTGRES
  // Delete any existing information about this field.
  strstream query;
  query << "SELECT delete_vector" << Dimension::nDim
    << "dfield('" << mFileName << "', '"
    << pathName << "');" << '\0';
  mSendQuery(query);

  // Insert the new information about this field.
  // The rather strange syntax for the elements of the field is:
  // '{"(v11, v12, ..., v1N)", ..., "(vN1, vN2, ..., VMN)"}'
  strstream query2;
  query2 << "INSERT INTO vector" << Dimension::nDim
    << "dfields (tag_id, name_id, elements) "
    << "     VALUES (get_tag_id('" << mFileName 
    << "'), get_name_id('" << pathName << "'), '{";
  for (int i = 0; i < value.numElements(); ++i)
  {
    if (i != 0)
    {
      query2 << ",";
    } // end if
    query2 << "\"(";
    for (int j = 0; j < Dimension::nDim; ++j)
    {
      if (j != 0)
        query2 << ',';
      query2 << (value(i))(j);
    } // end for
    query2 << ")\"";
  } // end for
  query2 << "}');" << '\0';
  mSendQuery(query2);
#endif
}

//------------------------------------------------------------------------------
// Write a Tensor Field to the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DbFileIO<Dimension>::
write(const Field<Dimension, typename Dimension::Tensor>& value, const string pathName) {
#ifdef USE_POSTGRES
  // Delete any existing information about this field.
  strstream query;
  query << "SELECT delete_tensor" << Dimension::nDim
    << "dfield('" << mFileName << "', '"
    << pathName << "');" << '\0';
  mSendQuery(query);

  // Insert the new information about this field.
  // The syntax for inserting tensor fields is similar to that of vector 
  // fields.  See above.
  strstream query2;
  query2 << "INSERT INTO tensor" << Dimension::nDim
    << "dfields(tag_id, name_id, elements) "
    << "     VALUES (get_tag_id('" << mFileName 
    << "'), get_name_id('" << pathName << "'), '{";
  for (int i = 0; i < value.numElements(); ++i)
  {
    if (i != 0)
    {
      query2 << ",";
    } // end if
    query2 << "\"(";
    for (int j = 0; j < Dimension::nDim; ++j)
    {
      if (j != 0)
        query2 << ',';
      for (int k = 0; k < Dimension::nDim; ++k)
      {
        if (k != 0)
          query2 << ',';
        query2 << (value(i))(j, k);
      } // end for
    } // end for
    query2 << ")\"";
  } // end for
  query2 << "}');" << '\0';
  mSendQuery(query2);
#endif
}

//------------------------------------------------------------------------------
// Write a SymTensor Field to the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DbFileIO<Dimension>::
write(const Field<Dimension, typename Dimension::SymTensor>& value, const string pathName) {
#ifdef USE_POSTGRES
  // Delete any existing information about this field.
  strstream query;
  query << "SELECT delete_symtensor" << Dimension::nDim
    << "dfield('" << mFileName << "', '"
    << pathName << "');" << '\0';
  mSendQuery(query);

  // Insert the new information about this field.
  strstream query2;
  query2 << "INSERT INTO symtensor" << Dimension::nDim
    << "dfields(tag_id, name_id, elements) "
    << "     VALUES (get_tag_id('" << mFileName 
    << "'), get_name_id('" << pathName << "'), '{";
  for (int i = 0; i < value.numElements(); ++i)
  {
    if (i != 0)
    {
      query2 << ',';
    } // end if
    query2 << "\"(";
    for (int j = 0; j < Dimension::nDim; ++j)
    {
      if (j != 0)
        query2 << ',';
      for (int k = 0; k <= j; ++k)
      {
        if (k != 0)
          query2 << ',';
        query2 << (value(i))(j, k);
      } // end for
    } // end for
    query2 << ")\"";
  } // end for
  query2 << "}');" << '\0';
  mSendQuery(query2);
#endif
}

//------------------------------------------------------------------------------
// Read a Scalar Field from the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DbFileIO<Dimension>::
read(Field<Dimension, typename Dimension::Scalar>& value, const string pathName) const {
#ifdef USE_POSTGRES
  // Get the elements of the given field from the database.
  strstream query;
  query << "SELECT elements "
    << "  FROM scalarfields "
    << " WHERE tag_id = get_tag_id('" << mFileName
    << "') AND name_id = get_name_id('" << pathName
    << "');" << '\0';

  // Get the value string and parse it into the given field.
  string valueString = mSingleResult(query);
  // The format of the field is "'{f1,f2,...,fN}'".
  size_t index = 0, firstDelim = 0, secondDelim = 0;
  int i = 0;
  while (index < valueString.length())
  {
    // Find a delimiter in the string and set firstDelim equal to its 
    // offset.
    firstDelim = valueString.find('"', secondDelim + 1);
    index = firstDelim;
    if (index > valueString.length())
      break;

    // Find the following delimiter in the string and set secondDelim equal 
    // to its offset.
    secondDelim = valueString.find('"', firstDelim + 1);

    // Form a substring bounded by but not including our delimiters.
    string subString = valueString.substr(firstDelim + 1, 
        (secondDelim - firstDelim) - 1);

    // Scan the number contained within the substring into the ith  
    // component of value.
    sscanf(subString.c_str(), "%lf", &(value(i)));

    // Increment the vector component index i.
    ++i;
  } // end while
#endif
}

//------------------------------------------------------------------------------
// Read a Vector Field from the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DbFileIO<Dimension>::
read(Field<Dimension, typename Dimension::Vector>& value, const string pathName) const {
#ifdef USE_POSTGRES
  // Get the elements of the given field from the database.
  strstream query;
  query << "SELECT elements "
    << "  FROM vector" << Dimension::nDim << "dfields "
    << " WHERE tag_id = get_tag_id('" << mFileName
    << "') AND name_id = get_name_id('" << pathName
    << "');" << '\0';
  string valueString = mSingleResult(query);
  // The field string is "{'vector1','vector2',...,'vectorN'}"
  size_t index = 0, firstDelim = 0, secondDelim = 0;
  int i = 0;
  while (index < valueString.length())
  {
    // Find a delimiter in the string and set firstDelim equal to its 
    // offset.
    firstDelim = valueString.find('"', secondDelim + 1);
    index = firstDelim;
    if (index > valueString.length())
      break;

    // Find the following delimiter in the string and set secondDelim equal 
    // to its offset.
    secondDelim = valueString.find('"', firstDelim + 1);

    // Form a substring bounded by but not including our delimiters.
    string subString = valueString.substr(firstDelim + 1, 
        (secondDelim - firstDelim) - 1);

    // Parse the vector contained within the substring.
    mParseVector(value(i), subString);

    // Increment the vector component index i.
    i++;
  } // end while
#endif
}

//------------------------------------------------------------------------------
// Read a Tensor Field from the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DbFileIO<Dimension>::
read(Field<Dimension, typename Dimension::Tensor>& value, const string pathName) const {
#ifdef USE_POSTGRES
  // Get the elements of the given field from the database.
  strstream query;
  query << "SELECT elements "
    << "  FROM tensor" << Dimension::nDim << "dfields "
    << " WHERE tag_id = get_tag_id('" << mFileName
    << "') AND name_id = get_name_id('" << pathName
    << "');" << '\0';
  string valueString = mSingleResult(query);
  // The field string is "{'tensor1','tensor2',...,'tensorN'}"
  size_t index = 0, firstDelim = 0, secondDelim = 0;
  int i = 0;
  while (index < valueString.length())
  {
    // Find a delimiter in the string and set firstDelim equal to its 
    // offset.
    firstDelim = valueString.find('"', secondDelim + 1);
    index = firstDelim;
    if (index > valueString.length())
      break;

    // Find the following delimiter in the string and set secondDelim equal 
    // to its offset.
    secondDelim = valueString.find('"', firstDelim + 1);

    // Form a substring bounded by but not including our delimiters.
    string subString = valueString.substr(firstDelim + 1, 
        (secondDelim - firstDelim) - 1);

    // Parse the tensor contained within the substring.
    mParseTensor(value(i), subString);

    // Increment the vector component index i.
    i++;
  } // end while
#endif
}

//------------------------------------------------------------------------------
// Read a SymTensor Field from the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DbFileIO<Dimension>::
read(Field<Dimension, typename Dimension::SymTensor>& value, const string pathName) const {
#ifdef USE_POSTGRES
  // Get the elements of the given field from the database.
  strstream query;
  query << "SELECT elements "
    << "  FROM symtensor" << Dimension::nDim << "dfields "
    << " WHERE tag_id = get_tag_id('" << mFileName
    << "') AND name_id = get_name_id('" << pathName
    << "');" << '\0';
  string valueString = mSingleResult(query);
  // The field string is "{'tensor1','tensor2',...,'tensorN'}"
  size_t index = 0, firstDelim = 0, secondDelim = 0;
  int i = 0;
  while (index < valueString.length())
  {
    // Find a delimiter in the string and set firstDelim equal to its 
    // offset.
    firstDelim = valueString.find('"', secondDelim + 1);
    index = firstDelim;
    if (index > valueString.length())
      break;

    // Find the following delimiter in the string and set secondDelim equal 
    // to its offset.
    secondDelim = valueString.find('"', firstDelim + 1);

    // Form a substring bounded by but not including our delimiters.
    string subString = valueString.substr(firstDelim + 1, 
        (secondDelim - firstDelim) - 1);

    // Parse the tensor contained within the substring.
    mParseSymTensor(value(i), subString);

    // Increment the vector component index i.
    i++;
  } // end while
#endif
}

#ifdef USE_POSTGRES
//------------------------------------------------------------------------------
// Parse a Vector from a given value string.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DbFileIO<Dimension>::
mParseVector(typename Dimension::Vector& value,
             const string valueString) const {
  // The format of the vector string is "(v0,v1,...,vn-1)"
  switch(Dimension::nDim)
  {
    case 1:
      sscanf(valueString.c_str(), "(%lf)", &(value(0)));
      break;
    case 2:
      sscanf(valueString.c_str(), "(%lf, %lf)", &(value(0)), &(value(1)));
      break;
    case 3:
      sscanf(valueString.c_str(), "(%lf, %lf, %lf)", &(value(0)), &(value(1)), &(value(2)));
      break;
    default:
      break;
  } // end switch
} // end mParseVector

//------------------------------------------------------------------------------
// Parse a Tensor from a given value string.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DbFileIO<Dimension>::
mParseTensor(typename Dimension::Tensor& value,
	     const string valueString) const {
  // The format of the tensor string is "(t00,t01,...,tn-1tn-1)"
  switch(Dimension::nDim)
  {
    case 1:
      sscanf(valueString.c_str(), "(%lf)", &value(0,0));
      break;
    case 2:
      sscanf(valueString.c_str(), "(%lf,%lf,%lf,%lf)", 
          &value(0,0), &value(0,1), &value(1,0), &value(1,1));
      break;
    case 3:
      sscanf(valueString.c_str(), "(%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf)",
          &value(0,0), &value(0,1), &value(0,2), 
          &value(1,0), &value(1,1), &value(1,2),
          &value(2,0), &value(2,1), &value(2,2));
      break;
    default:
      break;
  } // end switch
} // end mParseTensor

//------------------------------------------------------------------------------
// Parse a SymTensor from a given value string.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DbFileIO<Dimension>::
mParseSymTensor(typename Dimension::SymTensor& value,
		const string valueString) const {
  // The format of the tensor string is "(t00,t01,...,tn-1tn-1)"
  switch(Dimension::nDim)
  {
    case 1:
      sscanf(valueString.c_str(), "(%lf)", &value(0,0));
      break;
    case 2:
      sscanf(valueString.c_str(), "(%lf,%lf,%lf)", 
          &value(0,0), &value(1,0), &value(1,1));
      break;
    case 3:
      sscanf(valueString.c_str(), "(%lf,%lf,%lf,%lf,%lf,%lf)",
          &value(0,0), &value(1,0), &value(1,1),
          &value(2,0), &value(2,1), &value(2,2));
      break;
    default:
      break;
  } // end switch
} // end mParseSymTensor

#endif // #ifdef USE_POSTGRES

} // end namespace Spheral

//------------------------------------------------------------------------------
// Explicit instantiation
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
template class Spheral::DbFileIO< Dim<1> >;
template class Spheral::DbFileIO< Dim<2> >;
template class Spheral::DbFileIO< Dim<3> >;


