//---------------------------------Spheral++----------------------------------//
// ArtificialConductionPolicy -- Override the default energy policy in the
// presence of artificial conduction.
//----------------------------------------------------------------------------//
#include "ArtificialConductionPolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "NodeList/FluidNodeList.hh"
#include "Utilities/DBC.hh"

namespace Spheral {
    
    using NodeSpace::NodeList;
    using NodeSpace::FluidNodeList;
    using FieldSpace::FieldList;
    
    //------------------------------------------------------------------------------
    // Constructor.
    //------------------------------------------------------------------------------
    template<typename Dimension, typename Value>
    ArtificialConductionPolicy<Dimension>::
    ArtificialConductionPolicy(State<Dimension>::PolicyPointer& energyPolicy):
    FieldListUpdatePolicyBase<Dimension>(),
    mEnergyPolicy(energyPolicy){
        
    }
    
    //------------------------------------------------------------------------------
    // Destructor.
    //------------------------------------------------------------------------------
    template<typename Dimension, typename Value>
    ArtificialConductionPolicy<Dimension>::
    ~ArtificialConductionPolicy() {
    }
    
    //------------------------------------------------------------------------------
    // Update the field.
    //------------------------------------------------------------------------------
    template<typename Dimension, typename Value>
    void
    ArtificialConductionPolicy<Dimension>::
    update(const KeyType& key,
           State<Dimension>& state,
           StateDerivatives<Dimension>& derivs,
           const double multiplier,
           const double t,
           const double dt) {
        KeyType fieldKey, nodeListKey;
        StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
        REQUIRE(fieldKey == HydroFieldNames::specificThermalEnergy and
                nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
        
        
        
        FieldSpace::FieldList<Dimension, Value> f = state.fields(fieldKey, Value());
        const FieldSpace::FieldList<Dimension, Value> df = derivs.fields(incrementKey, Value());
        CHECK(f.size() == df.size());
        
        // Have the base class update the energy.
        mEnergyPolicy<Dimension>::update(key, state, derivs, multiplier, t, dt);
        
        // Get the artycond depsdt from derivs

        
        // Loop over the internal values of the field.
        const unsigned numNodeLists = f.size();
        for (unsigned k = 0; k != numNodeLists; ++k) {
            const unsigned n = f[k]->numInternalElements();
            for (unsigned i = 0; i != n; ++i) {
                f(k, i) += df(k,i) * dt;
            }
        }

    }
    
    //------------------------------------------------------------------------------
    // Equivalence operator.
    //------------------------------------------------------------------------------
    template<typename Dimension, typename Value>
    bool
    ArtificialConductionPolicy<Dimension>::
    operator==(const UpdatePolicyBase<Dimension>& rhs) const {
        
        // We're only equal if the other guy is also an increment operator.
        const ArtificialConductionPolicy<Dimension>* rhsPtr = dynamic_cast<const ArtificialConductionPolicy<Dimension>*>(&rhs);
        if (rhsPtr == 0) {
            return false;
        } else {
            return true;
        }
    }
    
}

