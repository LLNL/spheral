#-------------------------------------------------------------------------------
# SolidMaterialEquationsOfState.
# This is a backwards compatibility module for the equations of state in Spheral.
# I used to encode the units as a template parameter, but we've switched to 
# passing in the units as a constructor argument.  The interfaces & names for 
# the equations of state provided here emulate the original interfaces.
#-------------------------------------------------------------------------------
from SpheralCompiledPackages import PhysicalConstants
from MaterialUnits import MKS, CGS, Cosmological, Solar
from SolidMaterialUnits import CGuS

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

for dim in dims:
    exec("""
from SpheralCompiledPackages import (LinearPolynomialEquationOfState%(dim)sd,
                                     GruneisenEquationOfState%(dim)sd,
                                     MurnaghanEquationOfState%(dim)sd,
                                     TillotsonEquationOfState%(dim)sd,
                                     SteinbergGuinanStrength%(dim)sd)
""" % {"dim" : dim})

EOSFactoryString = """
#-------------------------------------------------------------------------------
# LinearPolynomialEquationOfState
#-------------------------------------------------------------------------------
class LinearPolynomialEquationOfState%(units)s%(dim)id(LinearPolynomialEquationOfState%(dim)id):
    def __init__(self, 
                 referenceDensity,
                 etamin,
                 etamax,
                 a0,
                 a1,
                 a2,
                 a3,
                 b0,
                 b1,
                 b2,
                 atomicWeight,
                 externalPressure = 0.0,
                 minimumPressure = -1e200,
                 maximumPressure =  1e200):
        self._units = %(units)s()
        LinearPolynomialEquationOfState%(dim)id.__init__(self,
                                                         referenceDensity,
                                                         etamin,
                                                         etamax,
                                                         a0,
                                                         a1,
                                                         a2,
                                                         a3,
                                                         b0,
                                                         b1,
                                                         b2,
                                                         atomicWeight,
                                                         self._units,
                                                         externalPressure,
                                                         minimumPressure,
                                                         maximumPressure)
        return

#-------------------------------------------------------------------------------
# GruneisenEquationOfState
#-------------------------------------------------------------------------------
class GruneisenEquationOfState%(units)s%(dim)id(GruneisenEquationOfState%(dim)id):
    def __init__(self, 
                 referenceDensity,
                 etamin,
                 etamax,
                 C0,
                 S1,
                 S2,
                 S3,
                 gamma0,
                 b,
                 atomicWeight,
                 externalPressure = 0.0,
                 minimumPressure = -1e200,
                 maximumPressure =  1e200):
        self._units = %(units)s()
        GruneisenEquationOfState%(dim)id.__init__(self,
                                                  referenceDensity,
                                                  etamin,
                                                  etamax,
                                                  C0,
                                                  S1,
                                                  S2,
                                                  S3,
                                                  gamma0,
                                                  b,
                                                  atomicWeight,
                                                  self._units,
                                                  externalPressure,
                                                  minimumPressure,
                                                  maximumPressure)
        return

#-------------------------------------------------------------------------------
# MurnaghanEquationOfState
#-------------------------------------------------------------------------------
class MurnaghanEquationOfState%(units)s%(dim)id(MurnaghanEquationOfState%(dim)id):
    def __init__(self, 
                 referenceDensity,
                 etamin,
                 etamax,
                 n,
                 K,
                 atomicWeight,
                 externalPressure = 0.0,
                 minimumPressure = -1e200,
                 maximumPressure =  1e200):
        self._units = %(units)s()
        MurnaghanEquationOfState%(dim)id.__init__(self,
                                                 referenceDensity,
                                                 etamin,
                                                 etamax,
                                                 n,
                                                 K,
                                                 atomicWeight,
                                                 self._units,
                                                 externalPressure,
                                                 minimumPressure,
                                                 maximumPressure)
        return

#-------------------------------------------------------------------------------
# TillotsonEquationOfState
#-------------------------------------------------------------------------------
class TillotsonEquationOfState%(units)s%(dim)id(TillotsonEquationOfState%(dim)id):
    def __init__(self, 
                 referenceDensity,
                 etamin,
                 etamax,
                 etamin_solid,
                 etamax_solid,
                 a,
                 b,
                 A,
                 B,
                 alpha,
                 beta,
                 eps0,
                 epsLiquid,
                 epsVapor,
                 atomicWeight,
                 externalPressure = 0.0,
                 minimumPressure = -1e200,
                 maximumPressure =  1e200):
        self._units = %(units)s()
        TillotsonEquationOfState%(dim)id.__init__(self,
                                                  referenceDensity,
                                                  etamin,
                                                  etamax,
                                                  etamin_solid,
                                                  etamax_solid,
                                                  a,
                                                  b,
                                                  A,
                                                  B,
                                                  alpha,
                                                  beta,
                                                  eps0,
                                                  epsLiquid,
                                                  epsVapor,
                                                  atomicWeight,
                                                  self._units,
                                                  externalPressure,
                                                  minimumPressure,
                                                  maximumPressure)
        return

#-------------------------------------------------------------------------------
# SteinbergGuinanStrength
#-------------------------------------------------------------------------------
class SteinbergGuinanStrength%(units)s%(dim)id(SteinbergGuinanStrength%(dim)id):
    def __init__(self, 
                 eos,
                 G0,
                 A,
                 B,
                 Y0,
                 Ymax,
                 Yp,
                 beta,
                 gamma0,
                 nhard,
                 coldEnergyFit,
                 meltEnergyFit):
        SteinbergGuinanStrength%(dim)id.__init__(self,
                                                 eos,
                                                 G0,
                                                 A,
                                                 B,
                                                 Y0,
                                                 Ymax,
                                                 Yp,
                                                 beta,
                                                 gamma0,
                                                 nhard,
                                                 coldEnergyFit,
                                                 meltEnergyFit)
        return

# #-------------------------------------------------------------------------------
# # SteinbergGuinanLundStrength
# #-------------------------------------------------------------------------------
# class SteinbergGuinanLundStrength%(units)s%(dim)id(SteinbergGuinanLundStrength%(dim)id):
#     def __init__(self, 
#                  eos,
#                  G0,
#                  A,
#                  B,
#                  Y0,
#                  Ymax,
#                  Yp,
#                  beta,
#                  gamma0,
#                  nhard,
#                  C1,
#                  C2,
#                  UK,
#                  YP,
#                  YTmax,
#                  coldEnergyFit,
#                  meltEnergyFit):
#         SteinbergGuinanStrength%(dim)id.__init__(self,
#                                                  eos,
#                                                  G0,
#                                                  A,
#                                                  B,
#                                                  Y0,
#                                                  Ymax,
#                                                  Yp,
#                                                  beta,
#                                                  gamma0,
#                                                  nhard,
#                                                  C1,
#                                                  C2,
#                                                  UK,
#                                                  YP,
#                                                  YTmax,
#                                                  coldEnergyFit,
#                                                  meltEnergyFit)
#         return

"""

#-------------------------------------------------------------------------------
# Create the different instantiations.
#-------------------------------------------------------------------------------
for dim in dims:
    for units in ("MKS", "CGS", "Cosmological", "Solar", "CGuS"):
        exec(EOSFactoryString % {"dim"   : dim,
                                 "units" : units})
