#-------------------------------------------------------------------------------
# RKCoefficients
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11template("Dimension")
class RKCoefficients:
    "Carries the RK correction coefficients around"

    correctionOrder = PYB11readwrite(doc="The correction order of the coefficients")
    coeffs = PYB11readwrite(doc="The coefficients vector")
