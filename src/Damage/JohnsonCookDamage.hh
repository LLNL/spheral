//---------------------------------Spheral++----------------------------------//
// JohnsonCookDamage -- an implementation of a Johnson-Cook damage law.
//
// Created by JMO, Mon Jul  9 08:21:23 PDT 2018
//----------------------------------------------------------------------------//
#ifndef __Spheral_JohnsonCookDamage_hh__
#define __Spheral_JohnsonCookDamage_hh__

#include <vector>

#include "Physics/Physics.hh"
#include "DataOutput/registerWithRestart.hh"

// Forward declarations.
namespace Spheral {
  template<typename Dimension> class State;
  template<typename Dimension> class StateDerivatives;
  namespace NodeSpace {
    template<typename Dimension> class SolidNodeList;
  }
  namespace DataBaseSpace {
    template<typename Dimension> class DataBase;
  }
  namespace FieldSpace {
    template<typename Dimension, typename DataType> class Field;
    template<typename Dimension, typename DataType> class FieldList;
  }
  namespace KernelSpace {
    template<typename Dimension> class TableKernel;
  }
  namespace FileIOSpace {
    class FileIO;
  }
}

namespace Spheral {
namespace PhysicsSpace {

template<typename Dimension>
class JohnsonCookDamage: public Physics<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs.
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;
  typedef FieldSpace::Field<Dimension, std::vector<double> > FlawStorageType;

  // Constructors, destructor.
  JohnsonCookDamage(NodeSpace::SolidNodeList<Dimension>& nodeList,
                    const unsigned seed,
                    const bool domainIndependent,
                    const double D1,
                    const double aD1,
                    const double bD1,
                    const double eps0D1,
                    const double D2,
                    const double aD2,
                    const double bD2,
                    const double eps0D2,
                    const double D3,
                    const double D4,
                    const double D5,
                    const double epsilondot0,
                    const double Tcrit,
                    const double sigmamax,
                    const double efailmin);
  virtual ~JohnsonCookDamage();

  // Attributes.
  const NodeSpace::SolidNodeList<Dimension>& nodeList() const;
  const FieldSpace::Field<Dimension, Scalar>& failureStrain() const;
  const FieldSpace::Field<Dimension, Scalar>& D1() const;
  const FieldSpace::Field<Dimension, Scalar>& D2() const;
  double D3() const;
  double D4() const;
  double D5() const;
  double epsilondot0() const;
  double Tcrit() const;
  double sigmamax() const;
  double efailmin() const;

  // Register our state and updates.
  virtual void registerState(DataBaseSpace::DataBase<Dimension>& dataBase,
                             State<Dimension>& state) override;

  //**************************************************************************
  // Restart methods.
  virtual std::string label() const { return "JohnsonCookDamage"; }
  virtual void dumpState(FileIOSpace::FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIOSpace::FileIO& file, const std::string& pathName);
  //**************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  NodeSpace::SolidNodeList<Dimension>& mNodeList;
  FieldSpace::Field<Dimension, Scalar> mD1, mD2, mFailureStrain;
  double mD3, mD4, mD5, mepsilondot0, mTcrit, msigmamax, mefailmin;

  // The restart registration.
  DataOutput::RestartRegistrationType mRestart;

  // No default constructor, copying or assignment.
  JohnsonCookDamage();
  JohnsonCookDamage(const JohnsonCookDamage&);
  JohnsonCookDamage& operator=(const JohnsonCookDamage&);
};

}
}

#include "JohnsonCookDamageInline.hh"

#else

// Forward declaration.
namespace Spheral {
  namespace PhysicsSpace {
    template<typename Dimension> class JohnsonCookDamage;
  }
}

#endif

