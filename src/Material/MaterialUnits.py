from SpheralCompiledPackages import PhysicalConstants

#-------------------------------------------------------------------------------
# MKS units.
#-------------------------------------------------------------------------------
class MKS(PhysicalConstants):
    def __init__(self):
        PhysicalConstants.__init__(self, 
                                   1.0,   # Unit length (m)
                                   1.0,   # Unit mass (kg)
                                   1.0,   # Unit time (sec)
                                   1.0,   # Unit temp (kelvin)
                                   1.0)   # Unit charge (coulomb)
        return

#-------------------------------------------------------------------------------
# CGS units.
#-------------------------------------------------------------------------------
class CGS(PhysicalConstants):
    def __init__(self):
        PhysicalConstants.__init__(self, 
                                   0.01,  # Unit length (m)
                                   0.001, # Unit mass (kg)
                                   1.0,   # Unit time (sec)
                                   1.0,   # Unit temp (kelvin)
                                   1.0)   # Unit charge (coulomb)
        return

#-------------------------------------------------------------------------------
# Cosmological units (Mpc, Mmsun, Myr)
#-------------------------------------------------------------------------------
class Cosmological(PhysicalConstants):
    def __init__(self):
        PhysicalConstants.__init__(self, 
                                   3.08567757e22, # Unit length (m)
                                   1.9891e36,     # Unit mass (kg)
                                   3.155674e19,   # Unit time (sec)
                                   1.0,   # Unit temp (kelvin)
                                   1.0)   # Unit charge (coulomb)
        return

#-------------------------------------------------------------------------------
# Solar units. (AU, Msun, yr)
#-------------------------------------------------------------------------------
class Solar(PhysicalConstants):
    def __init__(self):
        PhysicalConstants.__init__(self, 
                                   149597870700.0,   # Unit length (m)
                                   1.98892e30,       # Unit mass (kg)
                                   365.25*3600*24,   # Unit time (sec)
                                   1.0,   # Unit temp (kelvin)
                                   1.0)   # Unit charge (coulomb)
        return

