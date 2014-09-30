//---------------------------------Spheral++----------------------------------//
// ArtificialConductionPolicy -- Override the default energy policy in the
// presence of artificial conduction.
//
// Created by CDR, 9/30/2014
//----------------------------------------------------------------------------//

#ifndef __ArtificialConductionPolicy_hh__
#define __ArtificialConductionPolicy_hh__

#include <string>
#include "Database/FieldListUpdatePolicyBase.hh"

namespace Spheral {
    
    // Forward declarations.
    template<typename Dimension> class State;
    template<typename Dimension> class StateDerivatives;
    namespace NodeSpace {
        template<typename Dimension> class FluidNodeList;
    }
    namespace FieldSpace {
        template<typename Dimension, typename DataType> class Field;
    }
    
    template<typename Dimension>
    class ArtificialConductionPolicy: public FieldListUpdatePolicyBase<Dimension, ValueType> {
    public:
        //--------------------------- Public Interface ---------------------------//
        // Useful typedefs
        typedef typename Dimension::Scalar Scalar;
        typedef typename Dimension::SymTensor SymTensor;
        typedef typename UpdatePolicyBase<Dimension>::KeyType KeyType;
        typedef typename State<Dimension>::PolicyPointer PolicyPointer;
        
        // Constructors, destructor.
        ArtificialConductionPolicy(PolicyPointer& energyPolicy);
        virtual ~ArtificialConductionPolicy();
        
        // Overload the methods describing how to update Fields.
        virtual void update(const KeyType& key,
                            State<Dimension>& state,
                            StateDerivatives<Dimension>& derivs,
                            const double multiplier,
                            const double t,
                            const double dt);
        
        // Equivalence.
        virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const;
        
    private:
        //--------------------------- Private Interface ---------------------------//
        ArtificialConductionPolicy(const ArtificialConductionPolicy& rhs);
        ArtificialConductionPolicy& operator=(const ArtificialConductionPolicy& rhs);
    };
    
}

#else

// Forward declaration.
namespace Spheral {
    template<typename Dimension> class ArtificialConductionPolicy;
}

#endif