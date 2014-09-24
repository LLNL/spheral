//---------------------------------Spheral++----------------------------------//
// ArtificialConduction -- Artificial smoothing of energy discontinuities
//
//
// Created by CDR, 9/24/2014
//----------------------------------------------------------------------------//

#ifndef ArtificialConduction_HH
#define ArtificialConduction_HH

#include "Physics.hh"

namespace Spheral {
    template<typename Dimension> class State;
    template<typename Dimension> class StateDerivatives;
    
    namespace PhysicsSpace {
        template<typename Dimension>
        class ArtificialConduction: public Physics<Dimension> {
        public:
            //--------------------------- Public Interface ---------------------------//
            typedef typename Dimension::Scalar Scalar;
            typedef typename Dimension::Vector Vector;
            typedef typename Dimension::Tensor Tensor;
            typedef typename Dimension::SymTensor SymTensor;
            
            // Constructors
            ArtificialConduction(const bool gradPMode);
            
            // Destructor
            virtual ~ArtificialConduction();
            
            // Provide default methods for creating and registering an energy derivative.
            virtual void registerState(DataBaseSpace::DataBase<Dimension>& dataBase,
                                       State<Dimension>& state);
            virtual void registerDerivatives(DataBaseSpace::DataBase<Dimension>& dataBase,
                                             StateDerivatives<Dimension>& derivs);
            
            
            // Do any required one-time initializations on problem start up.
            virtual void initializeProblemStartup(DataBaseSpace::DataBase<Dimension>& dataBase);
            
            // Access to the derivative field(s).
            const FieldSpace::FieldList<Dimension, Vector>& DepsDt() const;
            
            // Access to vsig and gradPMode
            Scalar vsig() const;
            bool gradPMode() const;
            
            void gradPMode(bool val);
            
        private:
            //--------------------------- Private Interface ---------------------------//
            // Our derivative field(s).
            FieldSpace::FieldList<Dimension, Vector> mDepsDt;
            
            const bool mGradPMode;
            Scalar mVsig;

        };
    }
}

#else

// Fwd Declaration
namespace Spheral{
    namespace PhysicsSpace {
        template<typename Dimension> class ArtificialConduction;
    }
}

#endif
