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
            ArtificialConduction(const Scalar alphaArCond);
            
            // Destructor
            virtual ~ArtificialConduction();
            
            // Do any required one-time initializations on problem start up.
            virtual void initializeProblemStartup(DataBaseSpace::DataBase<Dimension>& dataBase);
            
            // Provide default methods for registering and iterating derivatives.
            virtual void registerDerivatives(DataBaseSpace::DataBase<Dimension>& dataBase,
                                             StateDerivatives<Dimension>& derivs);
            virtual
            void evaluateDerivatives(const Scalar time,
                                     const Scalar dt,
                                     const DataBaseSpace::DataBase<Dimension>& dataBase,
                                     const State<Dimension>& state,
                                     StateDerivatives<Dimension>& derivatives) const;
            
            
            
            // Accessor Fns
            
        private:
            //--------------------------- Private Interface ---------------------------//
            // Our derivative field(s).
            FieldSpace::FieldList<Dimension, Scalar> mGradP;
            Scalar mAlphaArCond;

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
