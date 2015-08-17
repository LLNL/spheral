#-------------------------------------------------------------------------------
# MaterialEquationsOfState.
# This is a backwards compatibility module for the equations of state in Spheral.
# I used to encode the units as a template parameter, but we've switched to 
# passing in the units as a constructor argument.  The interfaces & names for 
# the equations of state provided here emulate the original interfaces.
#-------------------------------------------------------------------------------
from MaterialUnits import MKS, CGS, Cosmological, Solar
from SpheralModules.Spheral.Material import PhysicalConstants
from spheralDimensions import spheralDimensions

dims = spheralDimensions()
for dim in dims:
    exec("from SpheralModules.Spheral.Material import GammaLawGas%(dim)id, PolytropicEquationOfState%(dim)id, IsothermalEquationOfState%(dim)id" % {"dim" : dim})

EOSFactoryString = """
#-------------------------------------------------------------------------------
# GammaLawGas
#-------------------------------------------------------------------------------
class GammaLawGas%(units)s%(dim)id(GammaLawGas%(dim)id):
    def __init__(self, 
                 gamma,
                 mu,
                 minimumPressure = -1e200,
                 maximumPressure =  1e200):
        self._units = %(units)s()
        GammaLawGas%(dim)id.__init__(self,
                                     gamma,
                                     mu, 
                                     self._units,
                                     minimumPressure,
                                     maximumPressure)
        return

#-------------------------------------------------------------------------------
# PolytropicEquationOfState
#-------------------------------------------------------------------------------
class PolytropicEquationOfState%(units)s%(dim)id(PolytropicEquationOfState%(dim)id):
    def __init__(self, 
                 K,
                 index,
                 mu,
                 minimumPressure = -1e200,
                 maximumPressure =  1e200):
        self._units = %(units)s()
        PolytropicEquationOfState%(dim)id.__init__(self,
                                                   K,
                                                   index,
                                                   mu, 
                                                   self._units,
                                                   minimumPressure,
                                                   maximumPressure)
        return

#-------------------------------------------------------------------------------
# IsothermalEquationOfState
#-------------------------------------------------------------------------------
class IsothermalEquationOfState%(units)s%(dim)id(IsothermalEquationOfState%(dim)id):
    def __init__(self, 
                 K,
                 mu,
                 minimumPressure = -1e200,
                 maximumPressure =  1e200):
        self._units = %(units)s()
        IsothermalEquationOfState%(dim)id.__init__(self,
                                                   K,
                                                   mu, 
                                                   self._units,
                                                   minimumPressure,
                                                   maximumPressure)
        return

"""

#-------------------------------------------------------------------------------
# Create the different instantiations.
#-------------------------------------------------------------------------------
for dim in dims:
    for units in ("MKS", "CGS", "Cosmological", "Solar"):
        exec(EOSFactoryString % {"dim"   : dim,
                                 "units" : units})
