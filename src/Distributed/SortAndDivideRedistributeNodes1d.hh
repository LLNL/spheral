//---------------------------------Spheral++----------------------------------//
// SortAndDivideRedistributeNodes1d -- 1-D implementation of the sort and 
// divide algorithm for domain decomposition.
//
// Created by JMO, Thu Dec 2 11:01:53 2004
//----------------------------------------------------------------------------//
#ifndef Spheral_SortAndDivideRedistributeNodes1d_hh
#define Spheral_SortAndDivideRedistributeNodes1d_hh

#include "SortAndDivideRedistributeNodes.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {

template<typename Dimension> class DataBase;
template<typename Dimension> class NodeList;
template<typename Dimension> class Boundary;
template<typename Dimension, typename DataType> class FieldList;

class SortAndDivideRedistributeNodes1d: public SortAndDivideRedistributeNodes<Dim<1> > {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef Dim<1> Dimension;
  typedef Dim<1>::Scalar Scalar;
  typedef Dim<1>::Vector Vector;
  typedef Dim<1>::Tensor Tensor;
  typedef Dim<1>::SymTensor SymTensor;

  // Constructors
  SortAndDivideRedistributeNodes1d(const double HExtent);

  // Destructor
  virtual ~SortAndDivideRedistributeNodes1d();

  // Given a Spheral++ data base of NodeLists, repartition it among the processors.
  virtual void redistributeNodes(DataBase<Dim<1> >& dataBase,
                                 std::vector<Boundary<Dim<1> >*> boundaries = std::vector<Boundary<Dim<1> >*>());


private:
  //--------------------------- Private Interface ---------------------------//
  // No default constructor, copy, or assignment operations.
  SortAndDivideRedistributeNodes1d();
  SortAndDivideRedistributeNodes1d(const SortAndDivideRedistributeNodes1d&);
  SortAndDivideRedistributeNodes1d& operator=(const SortAndDivideRedistributeNodes1d&);

};

}

#endif
