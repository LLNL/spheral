//---------------------------------Spheral++----------------------------------//
// A simple form of the reducing artificial viscosity from Morris & Monaghan.
// Computes a correction scaling for the typical Monaghan & Gingold visc.
// References:
//   Monaghan, J. J, & Gingold, R. A. 1983, J. Comput. Phys., 52, 374
//   Monaghan, J. J. 1992, ARA&A, 30, 543
//   Morris, J. P., & Monaghan, J. J. 1997, J. Comput. Phys., 136, 41
//
// Created by CDR, Aug 8th, 2014
//----------------------------------------------------------------------------//
#ifndef __MorrisMonaghanReducingViscosity__
#define __MorrisMonaghanReducingViscosity__

#include "ArtificialViscosity.hh"
#include "Physics/Physics.hh"

namespace Spheral {
namespace ArtificialViscositySpace {
    
template<typename Dimension>
    class MorrisMonaghanReducingViscosity: public PhysicsSpace::Physics<Dimension>{
public:
    //--------------------------- Public Interface ---------------------------//
    typedef typename Dimension::Scalar Scalar;
    typedef typename Dimension::Vector Vector;
    typedef typename Dimension::Tensor Tensor;
    typedef typename Dimension::SymTensor SymTensor;
    typedef typename PhysicsSpace::Physics<Dimension>::TimeStepType TimeStepType;
    
    // Constructors & Destructors
    MorrisMonaghanReducingViscosity(ArtificialViscosity<Dimension>& q,
                                    const Scalar nh,
                                    const Scalar aMin,
                                    const Scalar aMax);
    virtual ~MorrisMonaghanReducingViscosity();
    
    //............................................................................
    // Physics methods.
    // Increment the derivatives.
    virtual
    void evaluateDerivatives(const Scalar time,
                             const Scalar dt,
                             const DataBaseSpace::DataBase<Dimension>& dataBase,
                             const State<Dimension>& state,
                             StateDerivatives<Dimension>& derivatives) const;
    
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
    virtual std::string label() const { return "MorrisMonaghanReducingViscosity"; }
        
        
    const FieldSpace::FieldList<Dimension, Scalar>& DrvAlphaDt() const;
    Scalar nh() const;
    Scalar aMin() const;
    Scalar aMax() const;
    
    void aMin(Scalar val);
    void aMax(Scalar val);
    void nh(Scalar val);
    
private:
    //--------------------------- Private Interface ---------------------------//
    
    MorrisMonaghanReducingViscosity();
    MorrisMonaghanReducingViscosity(const MorrisMonaghanReducingViscosity&);
    MorrisMonaghanReducingViscosity& operator=(const MorrisMonaghanReducingViscosity&) const;
    
    Scalar mnh,maMin,maMax;
    FieldSpace::FieldList<Dimension, Scalar> mDrvAlphaDt;
    ArtificialViscosity<Dimension>& myq;
};
    
}
}

#else

namespace Spheral {
    namespace ArtificialViscositySpace {
        // Forward declaration.
        template<typename Dimension> class MorrisMonaghanReducingViscosity;
    }
}

#endif
