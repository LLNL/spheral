#-------------------------------------------------------------------------------
# A mock physics package to mess around with the CSPH corrections.
#-------------------------------------------------------------------------------
from Spheral1d import *

class CSPH_mod_package(Physics):

    def __init__(self):
        Physics.__init__(self)
        return

    def evaluateDerivatives(self, t, dt, db, state, derivs):
        return

    def dt(self, db, state, derivs, t):
        return pair_double_string(1e100, "No vote")

    def registerState(self, dt, state):
        return

    def registerDerivatives(self, db, derivs):
        return

    def label(self):
        return "CSPH_mod_package"

    def initialize(self, t, dt, db, state, derivs):

        # Grab the CSPH arrays.
        A0_fl = state.scalarFields(HydroFieldNames.A0_CSPH)
        A_fl = state.scalarFields(HydroFieldNames.A_CSPH)
        B_fl = state.vectorFields(HydroFieldNames.B_CSPH)

        A0 = A0_fl[0]
        A = A_fl[0]
        B = B_fl[0]

        print "A", A.internalValues()
        return

