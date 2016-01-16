//---------------------------------Spheral++----------------------------------//
// An implementation of the Cullen Dehnen Viscosity, with Hopkins alterations
// Computes a correction for the scaling.
// References:
// Cullen, L., & Dehnen, W. 2010, MNRAS, 408, 669
// Hopkins arXiv:1409.7395
//----------------------------------------------------------------------------//
#ifndef __CullenDehnenViscosity__
#define __CullenDehnenViscosity__

#include "ArtificialViscosity.hh"
#include "Physics/Physics.hh"

namespace Spheral {
  template<typename Dimension> class State;
  template<typename Dimension> class StateDerivatives;
  namespace NodeSpace {
    template<typename Dimension> class SmoothingScaleBase;
  }
  namespace ArtificialViscositySpace {
    template<typename Dimension> class ArtificialViscosity;
  }
  namespace KernelSpace {
    template<typename Dimension> class TableKernel;
  }
  namespace DataBaseSpace {
    template<typename Dimension> class DataBase;
  }
  namespace FieldSpace {
    template<typename Dimension, typename DataType> class Field;
    template<typename Dimension, typename DataType> class FieldList;
  }
  namespace FileIOSpace {
    class FileIO;
  }
}

namespace Spheral {
namespace ArtificialViscositySpace {
    
template<typename Dimension>
class CullenDehnenViscosity: public PhysicsSpace::Physics<Dimension>{
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename PhysicsSpace::Physics<Dimension>::TimeStepType TimeStepType;
        
  typedef typename PhysicsSpace::Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;
    
  // Constructors & Destructors
  CullenDehnenViscosity(ArtificialViscosity<Dimension>& q,
                        const KernelSpace::TableKernel<Dimension>& W,
                        const Scalar alphMax,
                        const Scalar alphMin,
                        const Scalar betaC,
                        const Scalar betaD,
                        const Scalar betaE,
                        const Scalar fKern,
                        const bool boolHopkins,
                        const bool reproducingKernelGradient);
  virtual ~CullenDehnenViscosity();
    
  //............................................................................
  // Physics methods.
  // Increment the derivatives.
  virtual
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBaseSpace::DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const;
  virtual
  void finalizeDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBaseSpace::DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const;
  virtual
  void finalize(const Scalar time,
                const Scalar dt,
                DataBaseSpace::DataBase<Dimension>& dataBase,
                State<Dimension>& state,
                StateDerivatives<Dimension>& derivs);
  virtual   
  void applyGhostBoundaries(State<Dimension>& state,
                            StateDerivatives<Dimension>& derivs);
  virtual
  void enforceBoundaries(State<Dimension>& state,
                         StateDerivatives<Dimension>& derivs);

    
  // Vote on a time step.
  virtual TimeStepType dt(const DataBaseSpace::DataBase<Dimension>& dataBase,
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const;
    
  // Register our state.
  virtual void registerState(DataBaseSpace::DataBase<Dimension>& dataBase,
                             State<Dimension>& state);
    
  // Register the derivatives/change fields for updating state.
  virtual void registerDerivatives(DataBaseSpace::DataBase<Dimension>& dataBase,
                                   StateDerivatives<Dimension>& derivs);
    
  // Do any required one-time initializations on problem start up.
  virtual void initializeProblemStartup(DataBaseSpace::DataBase<Dimension>& dataBase);
  //............................................................................
    
  // Restart methods.
  virtual std::string label() const { return "CullenDehnenViscosity"; }
  virtual void dumpState(FileIOSpace::FileIO& file, std::string pathName) const;
  virtual void restoreState(const FileIOSpace::FileIO& file, std::string pathName);

  Scalar alphMax() const;
  Scalar alphMin() const;
  Scalar betaE() const;
  Scalar betaD() const;
  Scalar betaC() const;
  Scalar fKern() const;
  bool boolHopkins() const;
  bool reproducingKernelGradient() const;
    
  void alphMax(Scalar val);
  void alphMin(Scalar val);
  void betaE(Scalar val);
  void betaD(Scalar val);
  void betaC(Scalar val);
  void fKern(Scalar val);
  void boolHopkins(bool val);
  void reproducingKernelGradient(bool val);

  // Access the stored interpolation kernels.
  const KernelSpace::TableKernel<Dimension>& kernel() const;
  const FieldSpace::FieldList<Dimension, Vector>&    PrevDvDt() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    PrevDivV() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    CullAlpha() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    PrevDivV2() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    CullAlpha2() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    DalphaDt() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    alphaLocal() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    alpha0() const;
    
private:
  //--------------------------- Private Interface ---------------------------//
  // The restart registration.
  DataOutput::RestartRegistrationType mRestart;
        
        
  CullenDehnenViscosity();
  CullenDehnenViscosity(const CullenDehnenViscosity&);
  CullenDehnenViscosity& operator=(const CullenDehnenViscosity&) const;

  Scalar malphMax, malphMin, mbetaE, mbetaD, mbetaC, mfKern;
  bool mboolHopkins;//Use Hopkins Reformulation
  bool mReproducingKernelGradient;  // Use reproducing kernels to estimate gradients.
  ArtificialViscosity<Dimension>& myq;
  const KernelSpace::TableKernel<Dimension>& mKernel;
  FieldSpace::FieldList<Dimension, Vector>    mPrevDvDt;//Will enroll as state fields
  FieldSpace::FieldList<Dimension, Scalar>    mPrevDivV;
  FieldSpace::FieldList<Dimension, Scalar>    mCullAlpha;
  FieldSpace::FieldList<Dimension, Scalar>    mPrevDivV2;//Will enroll as derivative fields so that we can store the previous time step value.
  FieldSpace::FieldList<Dimension, Scalar>    mCullAlpha2;

  FieldSpace::FieldList<Dimension, Scalar>    mDalphaDt;     // Time derivative of alpha
  FieldSpace::FieldList<Dimension, Scalar>    mAlphaLocal;   // Alpha local to be filled in derivatives
  FieldSpace::FieldList<Dimension, Scalar>    mAlpha0;       // The Hopkins form actually evolves alpha0
};
    
}
}

#else

namespace Spheral {
  namespace ArtificialViscositySpace {
    // Forward declaration.
    template<typename Dimension> class CullenDehnenViscosity;
  }
}

#endif
