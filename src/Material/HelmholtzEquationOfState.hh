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
                                     const MaterialPressureMinType minPressureType);
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
            
            // We also want the equivalent functions for individual calculations.
            virtual Scalar pressure(const Scalar massDensity,
                                    const Scalar specificThermalEnergy) const;
            
            virtual Scalar temperature(const Scalar massDensity,
                                       const Scalar specificThermalEnergy) const;
            
            virtual Scalar specificThermalEnergy(const Scalar massDensity,
                                                 const Scalar temperature) const;
            
            virtual Scalar specificHeat(const Scalar massDensity,
                                        const Scalar temperature) const;
            
            virtual Scalar soundSpeed(const Scalar massDensity,
                                      const Scalar specificThermalEnergy) const;
            
            virtual Scalar gamma(const Scalar massDensity,
                                 const Scalar specificThermalEnergy) const;
            
            virtual Scalar bulkModulus(const Scalar massDensity,
                                       const Scalar specificThermalEnergy) const;

            
            const FieldSpace::FieldList<Dimension, Scalar>& abar() const;
            const FieldSpace::FieldList<Dimension, Scalar>& zbar() const;
            
            
        private:
            //--------------------------- Private Interface ---------------------------//
            
            FieldSpace::FieldList<Dimension, Scalar>& mabar;
            FieldSpace::FieldList<Dimension, Scalar>& mzbar;
            

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
