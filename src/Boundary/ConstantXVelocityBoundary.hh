//---------------------------------Spheral++----------------------------------//
// ConstantXVelocityBoundary -- A boundary condition to enforce a constant 
// x velocity component on a given set of nodes.
//
// This boundary is very specialized -- it explicitly works on only one 
// NodeList.
//
// Created by JMO, Tue Oct 12 09:57:21 PDT 2004
//----------------------------------------------------------------------------//

#ifndef ConstantXVelocityBoundary_HH
#define ConstantXVelocityBoundary_HH

#include "ConstantVelocityBoundary.hh"

namespace Spheral {
  namespace NodeSpace {
    template<typename Dimension> class NodeList;
  }
  namespace FieldSpace {
    template<typename Dimension, typename DataType> class Field;
    template<typename Dimension, typename DataType> class FieldList;
  }
  namespace DataBaseSpace {
    template<typename Dimension> class DataBase;
  }
  namespace FileIOSpace {
    class FileIO;
  }
}

namespace Spheral {
namespace BoundarySpace {

template<typename Dimension>
class ConstantXVelocityBoundary: 
    public ConstantVelocityBoundary<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  // Constructors and destructors.
  ConstantXVelocityBoundary(const NodeSpace::NodeList<Dimension>& nodeList,
                            const std::vector<int>& nodeIndicies);
  virtual ~ConstantXVelocityBoundary();

  //**********************************************************************
  // Override the vector enforceBoundary method.
  // Use the given NodeList's neighbor object to select the ghost nodes.
  virtual void enforceBoundary(FieldSpace::Field<Dimension, Vector>& field) const;

  // Restart methods.
  virtual std::string label() const { return "ConstantXVelocityBoundary"; }
};

}
}

#else

// Forward declaration.
namespace Spheral {
  namespace BoundarySpace {
    template<typename Dimension> class ConstantXVelocityBoundary;
  }
}

#endif
