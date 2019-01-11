//---------------------------------Spheral++----------------------------------//
// DbFileIO -- Provide the interface to a Postgres database.
//
// Created by Jeffrey Johnson and Charlie Crabb , Mon Aug 27 14:04:00 PDT 2001
// Notice that the conditional compile guards are placed so that this class 
// will exist without Postgres, but will not do anything.
//----------------------------------------------------------------------------//
#ifndef DbFileIO_HH
#define DbFileIO_HH

//#include <strstream>
#include <sstream>
#include "FileIO.hh"

#if USE_POSTGRES
// This brings in the PostgreSQL C++ support.  You will need to have built 
// libpq++ within your Postgres distribution.
#include <libpq++.h>
#endif

namespace Spheral {

template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class DbFileIO: public FileIO<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors.
  DbFileIO();
  DbFileIO(const string filename);

  // Destructor.
  virtual ~DbFileIO();

  // All File objects must provide methods to open and close the files.
  virtual void open(const string fileName, AccessType access); 
  virtual void close();

  //******************************************************************************
  // Methods all FileIO descendent classes must provide.
  //******************************************************************************
  // All FileIO objects had better be able to read and write the primitive 
  // DataTypes.
  virtual void write(const int& value, const string pathName); 
  virtual void write(const bool& value, const string pathName); 
  virtual void write(const Scalar& value, const string pathName); 
  virtual void write(const Vector& value, const string pathName); 
  virtual void write(const Tensor& value, const string pathName); 
  virtual void write(const SymTensor& value, const string pathName); 
  virtual void write(const string& value, const string pathName); 

  virtual void read(int& value, const string pathName) const; 
  virtual void read(bool& value, const string pathName) const; 
  virtual void read(Scalar& value, const string pathName) const; 
  virtual void read(Vector& value, const string pathName) const; 
  virtual void read(Tensor& value, const string pathName) const; 
  virtual void read(SymTensor& value, const string pathName) const;
  virtual void read(string& value, const string pathName) const;

  // Require that all FileIO objects provide methods to read and write
  // Fields of specific DataTypes.
  virtual void write(const Field<Dimension, Scalar>& field, const string pathName);
  virtual void write(const Field<Dimension, Vector>& field, const string pathName);
  virtual void write(const Field<Dimension, Tensor>& field, const string pathName); 
  virtual void write(const Field<Dimension, SymTensor>& field, const string pathName);

  virtual void read(Field<Dimension, Scalar>& field, const string pathName) const; 
  virtual void read(Field<Dimension, Vector>& field, const string pathName) const; 
  virtual void read(Field<Dimension, Tensor>& field, const string pathName) const; 
  virtual void read(Field<Dimension, SymTensor>& field, const string pathName) const; 
  //******************************************************************************

#ifdef USE_POSTGRES

protected:
  //--------------------------- Protected Interface ---------------------------//

  // Method to provide connection given the proper info.
  PgDatabase& mConnection() const;

  // Method to return a scalar result from a query.
  // Use this for reading values from the database.
  const char* mSingleResult(std::strstream& query) const;

  // Method to send a simple query to the database.
  // Use this for modifying the database.
  void mSendQuery(std::strstream& query) const;

  // Parse a vector from a vector string.
  void mParseVector(typename Dimension::Vector& value,
                    const string valueString) const;

  // Parse a tensor from a tensor string.
  void mParseTensor(typename Dimension::Tensor& value,
                    const string valueString) const;

  // Parse a vector from a vector string.
  void mParseSymTensor(typename Dimension::SymTensor& value,
                       const string valueString) const;

private:
  //--------------------------- Private Interface ---------------------------//
  // Don't allow assignment.
  DbFileIO& operator=(const DbFileIO& rhs);

  // Database connection.  Access this with mConnection() above.
  static PgDatabase* sDatabasePtr;
#endif // #ifdef USE_POSTGRES

};

}

#endif

