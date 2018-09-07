//---------------------------------------------------------------------------
//
// AggregateField -- a collection of fields from several NodeLists.
// You can use an AggregateField just as you would a Field; the interface 
// is almost identical.
//
// $Id: AggregateField.hh 613 2003-07-31 22:44:32Z mikeowen $
// 
//---------------------------------------------------------------------------

#ifndef AGGREGATEFIELD_HH
#define AGGREGATEFIELD_HH

#include "FieldSet.hh"

namespace Spheral {

template <typename Dimension, typename DataType>
class AggregateField: public FieldBase<AggregateField<Dimension, DataType> > {
  public:

  //! Default constructor.
  AggregateField();

  //! Constructor with a set of fields.
  explicit AggregateField(FieldSet& fields);

  //! Destructor.
  virtual ~AggregateField();

  private:

  // No copy constructor.
  AggregateField(const AggregateField&);
  
  // No assignment operator.
  AggregateField& operator=(const AggregateField&);

}; // end class AggregateField
  
} // end namespace Spheral

#endif
