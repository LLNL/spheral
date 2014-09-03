//
//  HelmholtzEquationOfState.h
//  
//
//  Created by Raskin, Cody Dantes on 8/28/14.
//
//

#ifndef ____HelmholtzEquationOfState_hh__
#define ____HelmholtzEquationOfState_hh__

#include "EquationOfState.hh"

namespace Spheral {
    namespace Material {
        
        template<typename Dimension>
        class HelmholtzEquationOfState: public EquationOfState<Dimension> {
        
        public:
            //--------------------------- Public Interface ---------------------------//
            typedef typename Dimension::Scalar Scalar;
            typedef typename Dimension::Vector Vector;
            typedef typename Dimension::Tensor Tensor;
            typedef typename Dimension::SymTensor SymTensor;
            
            // Constructors, destructors.
            HelmholtzEquationOfState(const PhysicalConstants& constants,
                                     const double minimumPressure,
                                     const double maximumPressure,
                                     const double minimumTemperature,
                                     const double maximumTemperature,
                                     const MaterialPressureMinType minPressureType,
                                     const Scalar abar0,
                                     const Scalar zbar0);
            ~HelmholtzEquationOfState();

            
            // We require any equation of state to define the following properties.
            virtual void setPressure(FieldSpace::Field<Dimension, Scalar>& Pressure,
                                     const FieldSpace::Field<Dimension, Scalar>& massDensity,
                                     const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy) const;
            
            virtual void setTemperature(FieldSpace::Field<Dimension, Scalar>& temperature,
                                        const FieldSpace::Field<Dimension, Scalar>& massDensity,
                                        const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy) const;
            
            virtual void setSpecificThermalEnergy(FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy,
                                                  const FieldSpace::Field<Dimension, Scalar>& massDensity,
                                                  const FieldSpace::Field<Dimension, Scalar>& temperature) const;
            
            virtual void setSpecificHeat(FieldSpace::Field<Dimension, Scalar>& specificHeat,
                                         const FieldSpace::Field<Dimension, Scalar>& massDensity,
                                         const FieldSpace::Field<Dimension, Scalar>& temperature) const;
            
            virtual void setSoundSpeed(FieldSpace::Field<Dimension, Scalar>& soundSpeed,
                                       const FieldSpace::Field<Dimension, Scalar>& massDensity,
                                       const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy) const;
            
            virtual void setGammaField(FieldSpace::Field<Dimension, Scalar>& gamma,
                                       const FieldSpace::Field<Dimension, Scalar>& massDensity,
                                       const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy) const;
            
            virtual void setBulkModulus(FieldSpace::Field<Dimension, Scalar>& bulkModulus,
                                        const FieldSpace::Field<Dimension, Scalar>& massDensity,
                                        const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy) const;
            
            // Some of the following methods are disabled
            virtual Scalar pressure(const Scalar massDensity,
                                    const Scalar specificThermalEnergy) const { VERIFY2(false, "HelmholtzEquationOfState does not support individual state calls."); }
            
            virtual Scalar temperature(const Scalar massDensity,
                                       const Scalar specificThermalEnergy) const { VERIFY2(false, "HelmholtzEquationOfState does not support individual state calls."); }
            
            virtual Scalar specificThermalEnergy(const Scalar massDensity,
                                                 const Scalar temperature) const { VERIFY2(false, "HelmholtzEquationOfState does not support individual state calls."); }
            
            virtual Scalar specificHeat(const Scalar massDensity,
                                        const Scalar temperature) const { VERIFY2(false, "HelmholtzEquationOfState does not support individual state calls."); }
            
            virtual Scalar soundSpeed(const Scalar massDensity,
                                      const Scalar specificThermalEnergy) const { VERIFY2(false, "HelmholtzEquationOfState does not support individual state calls."); }
            
            // Get the effective gamma (ratio of specific heats) for this eos.
            virtual Scalar gamma(const Scalar massDensity,
                                 const Scalar specificThermalEnergy) const { VERIFY2(false, "HelmholtzEquationOfState does not support individual state calls."); }
            
            // Get the bulk modulus.
            virtual Scalar bulkModulus(const Scalar massDensity,
                                       const Scalar specificThermalEnergy) const { VERIFY2(false, "HelmholtzEquationOfState does not support individual state calls."); }

            
            const FieldSpace::FieldList<Dimension, Scalar>& abar() const;
            const FieldSpace::FieldList<Dimension, Scalar>& zbar() const;
            
            
        private:
            //--------------------------- Private Interface ---------------------------//
            
            const FieldSpace::FieldList<Dimension, Scalar>* mabarPtr;
            const FieldSpace::FieldList<Dimension, Scalar>* mzbarPtr;
            
            FieldSpace::FieldList<Dimension, Scalar> myAbar;
            FieldSpace::FieldList<Dimension, Scalar> myZbar;
            FieldSpace::FieldList<Dimension, Scalar> mySpecificThermalEnergy;
            FieldSpace::FieldList<Dimension, Scalar> myMassDensity;
            FieldSpace::FieldList<Dimension, Scalar> myTemperature;
            FieldSpace::FieldList<Dimension, Scalar> myPressure;
            FieldSpace::FieldList<Dimension, Scalar> mySoundSpeed;
            FieldSpace::FieldList<Dimension, Scalar> myGamma;
            
            Scalar mabar0, mzbar0, mTmin, mTmax, mPmin, mPmax;
            bool needUpdate;

        };
    }
}

#else
        
namespace Spheral {
    namespace Material {
        // Forward declaration.
        template<typename Dimension> class HelmholtzEquationOfState;
    }
}

#endif
