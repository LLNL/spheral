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
  template<typename Dimension> class SmoothingScaleBase;
  template<typename Dimension> class ArtificialViscosity;
  template<typename Dimension> class TableKernel;
  template<typename Dimension> class DataBase;
  template<typename Dimension, typename DataType> class Field;
  template<typename Dimension, typename DataType> class FieldList;
  class FileIO;
}

namespace Spheral {
    
template<typename Dimension>
class CullenDehnenViscosity: public Physics<Dimension>{
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Physics<Dimension>::TimeStepType TimeStepType;
        
  typedef typename Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;
    
  // Constructors & Destructors
  CullenDehnenViscosity(ArtificialViscosity<Dimension>& q,
                        const TableKernel<Dimension>& W,
                        const Scalar alphMax,
                        const Scalar alphMin,
                        const Scalar betaC,
                        const Scalar betaD,
                        const Scalar betaE,
                        const Scalar fKern,
                        const bool boolHopkins);
  virtual ~CullenDehnenViscosity();
    
  //............................................................................
  // Physics methods.
  // Increment the derivatives.
  virtual
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const;
  virtual
  void finalizeDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const;
  virtual
  void finalize(const Scalar time,
                const Scalar dt,
                DataBase<Dimension>& dataBase,
                State<Dimension>& state,
                StateDerivatives<Dimension>& derivs);
  virtual   
  void applyGhostBoundaries(State<Dimension>& state,
                            StateDerivatives<Dimension>& derivs);
  virtual
  void enforceBoundaries(State<Dimension>& state,
                         StateDerivatives<Dimension>& derivs);

    
  // Vote on a time step.
  virtual TimeStepType dt(const DataBase<Dimension>& dataBase,
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const;
    
  // Register our state.
  virtual void registerState(DataBase<Dimension>& dataBase,
                             State<Dimension>& state);
    
  // Register the derivatives/change fields for updating state.
  virtual void registerDerivatives(DataBase<Dimension>& dataBase,
                                   StateDerivatives<Dimension>& derivs);
    
  // Do any required one-time initializations on problem start up.
  virtual void initializeProblemStartup(DataBase<Dimension>& dataBase);
  //............................................................................
    
  // Restart methods.
  virtual std::string label() const { return "CullenDehnenViscosity"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);

  Scalar alphMax() const;
  Scalar alphMin() const;
  Scalar betaE() const;
  Scalar betaD() const;
  Scalar betaC() const;
  Scalar fKern() const;
  bool boolHopkins() const;
    
  void alphMax(Scalar val);
  void alphMin(Scalar val);
  void betaE(Scalar val);
  void betaD(Scalar val);
  void betaC(Scalar val);
  void fKern(Scalar val);
  void boolHopkins(bool val);

  // Access the stored interpolation kernels.
  const TableKernel<Dimension>& kernel() const;
  const FieldList<Dimension, Vector>&    PrevDvDt() const;
  const FieldList<Dimension, Scalar>&    PrevDivV() const;
  const FieldList<Dimension, Scalar>&    CullAlpha() const;
  const FieldList<Dimension, Scalar>&    PrevDivV2() const;
  const FieldList<Dimension, Scalar>&    CullAlpha2() const;
  const FieldList<Dimension, Scalar>&    DalphaDt() const;
  const FieldList<Dimension, Scalar>&    alphaLocal() const;
  const FieldList<Dimension, Scalar>&    alpha0() const;
  const FieldList<Dimension, Scalar>&    R() const;
  const FieldList<Dimension, Scalar>&    vsig() const; 
private:
  //--------------------------- Private Interface ---------------------------//

  CullenDehnenViscosity();
  CullenDehnenViscosity(const CullenDehnenViscosity&);
  CullenDehnenViscosity& operator=(const CullenDehnenViscosity&) const;

  FieldList<Dimension, Vector>    mPrevDvDt;//Will enroll as state fields
  FieldList<Dimension, Scalar>    mPrevDivV;
  FieldList<Dimension, Scalar>    mCullAlpha;
  FieldList<Dimension, Scalar>    mPrevDivV2;//Will enroll as derivative fields so that we can store the previous time step value.
  FieldList<Dimension, Scalar>    mCullAlpha2;

  FieldList<Dimension, Scalar>    mDalphaDt;     // Time derivative of alpha
  FieldList<Dimension, Scalar>    mAlphaLocal;   // Alpha local to be filled in derivatives
  FieldList<Dimension, Scalar>    mAlpha0;       // The Hopkins form actually evolves alpha0
  
  FieldList<Dimension, Scalar>    mR;    //Will enroll as derivs fields
  FieldList<Dimension, Scalar>    mVsig; //Will enroll as derivs fields

  Scalar malphMax, malphMin, mbetaC, mbetaD, mbetaE, mfKern;
  bool mboolHopkins;//Use Hopkins Reformulation
  ArtificialViscosity<Dimension>& myq;
  const TableKernel<Dimension>& mKernel;

  // The restart registration.
  RestartRegistrationType mRestart;
};
    
}

#include "CullenDehnenViscosityInline.hh"

#else

namespace Spheral {
  // Forward declaration.
  template<typename Dimension> class CullenDehnenViscosity;
}

#endif
