# This bit of code solves the equation for hydrostatic equilibrium to return
# a density profile function to be used in various NodeGenerators
# NOTE: the eos tuple should contain at least one pair of eos and the range
# to which it applies. ex. (tillotsoneos, [rmin,rmax], gammaeos, [rmax,rmax2])

from math import *

#-------------------------------------------------------------------------------
# 3-D solvers
#-------------------------------------------------------------------------------
class HydroStaticProfileConstantTemp3d():
    from Spheral3d import *
    def __init__(self,
                 rho0,
                 rMax,
                 temp,
                 eostup,
                 units,
                 y0=0,
                 nbins=1000):

        self.y0     = y0
        self.rMax   = rMax
        self.nbins  = nbins
        self.rho0   = rho0
        self.soln   = []
        
        eoscount    = len(eostup)/2
        
        r   = self.rMax
        rho = self.rho0
        dr  = self.rMax/self.nbins
        y   = self.y0
        dy  = 0
        
        nodes   = makeVoidNodeList("nodes", numInternal=1)
        ef      = ScalarField("eps", nodes)
        Kf      = ScalarField("mod", nodes)
        rhof    = ScalarField("rho", nodes)
        tempf   = ScalarField("temp", nodes)
        
        tempf[0] = temp
        
        while(r>0):
            # get the eos for this radius
            if(eoscount>1):
                for i in xrange(eoscount):
                    ermin = eostup[2*i+1][0]
                    ermax = eostup[2*i+1][1]
                    if(r<=ermax and r>=ermin):
                        eos = eostup[2*i]
                        break
            else:
                eos = eostup[0]
        
            rhof[0] = rho
            eos.setSpecificThermalEnergy(ef,rhof,tempf)
            e       = ef[0]
            eos.setBulkModulus(Kf,rhof,ef)
            K       = Kf[0]
            dy      = dr*(2.0/rho*y*y - 2.0/r*y - units.G/K*4.0*pi*pow(rho,3.0))
            self.soln.append([r,rho])
            y       = y + dy
            rho     = rho - y*dr
            r       = r - dr

        self.soln.sort()

    def __call__(self,r):
        rho = self.rho0
        for i in xrange(len(self.soln)):
            if(self.soln[i][0] > r):
                if(i>0):
                    f1  = self.soln[i][1]
                    f0  = self.soln[i-1][1]
                    r1  = self.soln[i][0]
                    r0  = self.soln[i-1][0]
                    rho = (f1-f0)*(r-r1)/(r1-r0)+f1
                else:
                    rho = self.soln[0][1]
                break
        return rho

#-------------------------------------------------------------------------------
# 2-D solvers
#-------------------------------------------------------------------------------
class HydroStaticProfileConstantTemp2d():
    from Spheral2d import *
    def __init__(self,
                 rho0,
                 rMax,
                 temp,
                 eostup,
                 units,
                 y0=0,
                 nbins=1000):
        
        self.y0     = y0
        self.rMax   = rMax
        self.nbins  = nbins
        self.rho0   = rho0
        eoscount    = len(eostup)/2
        self.soln   = []
        
        r   = self.rMax
        rho = self.rho0
        dr  = self.rMax/self.nbins
        y   = self.y0
        dy  = 0
        
        nodes   = makeVoidNodeList("nodes", numInternal=1)
        ef      = ScalarField("eps", nodes)
        Kf      = ScalarField("mod", nodes)
        rhof    = ScalarField("rho", nodes)
        tempf   = ScalarField("temp", nodes)
        
        tempf[0] = temp
        
        while(r>0):
            # get the eos for this radius
            if(eoscount>1):
                for i in xrange(eoscount):
                    ermin = eostup[2*i+1][0]
                    ermax = eostup[2*i+1][1]
                    if(r<=ermax and r>=ermin):
                        eos = eostup[2*i]
                        break
            else:
                eos = eostup[0]
            
            rhof[0] = rho
            eos.setSpecificThermalEnergy(ef,rhof,tempf)
            e       = ef[0]
            eos.setBulkModulus(Kf,rhof,ef)
            K       = Kf[0]
            dy      = dr*(2.0/rho*y*y - 1.0/r*y - units.G/K*2.0*pi*pow(rho,3.0))
            self.soln.append([r,rho])
            y       = y + dy
            rho     = rho - y*dr
            r       = r - dr
        
        self.soln.sort()
    
    def __call__(self,r):
        rho = self.rho0
        for i in xrange(len(self.soln)):
            if(self.soln[i][0] > r):
                if(i>0):
                    f1  = self.soln[i][1]
                    f0  = self.soln[i-1][1]
                    r1  = self.soln[i][0]
                    r0  = self.soln[i-1][0]
                    rho = (f1-f0)*(r-r1)/(r1-r0)+f1

                else:
                    rho = self.soln[0][1]
                break
        return rho


