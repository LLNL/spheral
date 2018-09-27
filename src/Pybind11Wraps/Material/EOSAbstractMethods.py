#-------------------------------------------------------------------------------
# Helper to inject common virtual methods for equations of state
#-------------------------------------------------------------------------------
from PYB11Generator import *
import inspect

@PYB11ignore
class EOSAbstractMethods:

    @PYB11const
    def setPressure(self,
                    Pressure = "ScalarField&",
                    massDensity = "const ScalarField&",
                    specificThermalEnergy = "const ScalarField&"):
        return "void"

    @PYB11const
    def setTemperature(self,
                       temperature = "ScalarField&",
                       massDensity = "const ScalarField&",
                       specificThermalEnergy = "const ScalarField&"):
        return "void"

    @PYB11const
    def setSpecificThermalEnergy(self,
                                 specificThermalEnergy = "ScalarField&",
                                 massDensity = "const ScalarField&",
                                 temperature = "const ScalarField&"):
        return "void"

    @PYB11const
    def setSpecificHeat(self,
                        specificHeat = "ScalarField&",
                        massDensity = "const ScalarField&",
                        temperature = "const ScalarField&"):
        return "void"

    @PYB11const
    def setSoundSpeed(self,
                      soundSpeed = "ScalarField&",
                      massDensity = "const ScalarField&",
                      specificThermalEnergy = "const ScalarField&"):
        return "void"

    @PYB11const
    def setGammaField(self,
                      gamma = "ScalarField&",
                      massDensity = "const ScalarField&",
                      specificThermalEnergy = "const ScalarField&"):
        return "void"

    @PYB11const
    def setBulkModulus(self,
                       bulkModulus = "ScalarField&",
                       massDensity = "const ScalarField&",
                       specificThermalEnergy = "const ScalarField&"):
        return "void"

    @PYB11const
    def setEntropy(self,
                   entropy = "ScalarField&",
                   massDensity = "const ScalarField&",
                   specificThermalEnergy = "const ScalarField&"):
        return "void"

    @PYB11const
    def valid(self):
        return "bool"
