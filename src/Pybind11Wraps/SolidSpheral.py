#-------------------------------------------------------------------------------
# Master import file to import the Spheral packages and a standard set of
# helper extensions, including the optional solid material and strength
# extensions.
#-------------------------------------------------------------------------------
from Spheral import *
from SpheralModules.Spheral.SolidMaterial import *
from SolidNodeLists import *
from GradyKippTensorDamage import *
from SolidSPHHydros import *
from SolidSPHHydrosRZ import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

from spheralDimensions import spheralDimensions
if 2 in spheralDimensions():
    from SolidSPHHydrosRZ import *
    from SolidCRKSPHHydrosRZ import *

# ------------------------------------------------------------------------------
# Import the SolidMaterial python extensions.
# ------------------------------------------------------------------------------
from SolidMaterialUnits import *
from SolidMaterialEquationsOfState import *

# ------------------------------------------------------------------------------
# Import our shadow layers for augmenting C++ types.
# ------------------------------------------------------------------------------
for shadowedthing in ("TillotsonEquationOfState",
                      "ConstantStrength"):
    for dim in dims:
        exec("from Shadow%(thing)s import %(thing)s%(dim)sd" % {"thing" : shadowedthing,
                                                                "dim"   : dim})
