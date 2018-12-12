from SpheralCompiledPackages import PhysicalConstants

#-------------------------------------------------------------------------------
# CGuS units.
#-------------------------------------------------------------------------------
class CGuS(PhysicalConstants):
    def __init__(self):
        PhysicalConstants.__init__(self, 
                                   0.01,   # Unit length (m)
                                   0.001,  # Unit mass (kg)
                                   1.0e-6) # Unit time (sec)
        return
