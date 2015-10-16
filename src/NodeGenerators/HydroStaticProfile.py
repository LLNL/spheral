# This bit of code solves the equation for hydrostatic equilibrium to return
# a density profile function to be used in various NodeGenerators
# NOTE: the eos tuple should contain at least one pair of eos and the range
# to which it applies. ex. (tillotsoneos, [rmin,rmax], gammaeos, [rmax,rmax2])

from math import *

#-------------------------------------------------------------------------------
# 3-D solvers
#-------------------------------------------------------------------------------
class LaneEmdenSolver():
    def __init__(self,
                 rhoc,
                 rMax,
                 Tc,
                 gamma,
                 kappa,
                 eostup,
                 units,
                 nbins=1000):
        
        self.rho0 = rhoc
        n = 1.0/(gamma-1.0)
        G = units.G
        a = sqrt((n+1)*kappa*pow(rhoc,1/n-1)/(4.0*pi*G))
        
        eoscount    = len(eostup)/2
        sMax        = rMax / a
        
        self.soln   = []
        
        from SolidSpheral3d import makeVoidNodeList
        from SolidSpheral3d import ScalarField
        
        nodes   = makeVoidNodeList("nodes", numInternal=1)
        ef      = ScalarField("eps", nodes)
        Kf      = ScalarField("mod", nodes)
        rhof    = ScalarField("rho", nodes)
        tempf   = ScalarField("temp", nodes)
        
        storedE0 = 0
        self.e0  = 0
        
        th = 1.0
        y  = 0
        s  = 0
        ds = sMax / nbins
        
        while (s<sMax):
            r = a * s
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


            Temp    = th*Tc
            rho     = rhoc*pow(th,1.0/(gamma-1.0))
            
            tempf[0]= Temp
            rhof[0] = rho
            
            eos.setSpecificThermalEnergy(ef,rhof,tempf)
            e       = ef[0]
            
            if(storedE0==0):
                self.e0     = e
                storedE0    = 1
            
            self.soln.append([r,rho,e])
            
            # now advance to next values using current s
            
            dy      = ds * (-2.0*y/s - th**n)
            y       = y + dy
            th      = th + ds*y
            s       = s + ds

        self.soln.sort()

    def __call__(self,r):
        return self.density(r)

    def density(self,r):
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
                
                
                
    def energy(self,r):
        e   = self.e0
        for i in xrange(len(self.soln)):
            if(self.soln[i][0] > r):
                if(i>0):
                    e1  = self.soln[i][2]
                    e0  = self.soln[i-1][2]
                    r1  = self.soln[i][0]
                    r0  = self.soln[i-1][0]
                    e   = (e1-e0)*(r-r1)/(r1-r0)+e1
                else:
                    e   = self.soln[0][2]
                break
        return e


class WeppnerSolver():
    def __init__(self,
                 rhoc,
                 rMax,
                 Tc,
                 gamma,
                 kappa,
                 eostup,
                 units,
                 nbins=1000):
        nothing = 0
    def __call__(self,r):
        return 0

class HydroStaticProfileConstantTemp3d():
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
        
        from SolidSpheral3d import makeVoidNodeList
        from SolidSpheral3d import ScalarField

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
            
            print "dy, dr, rho, y, G, K = {0:3.3e} {1:3.3e} {2:3.3e} {3:3.3e} {4:3.3e} {5:3.3e}".format(dy,dr,rho,y,units.G,K)
            
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
        
        from SolidSpheral2d import makeVoidNodeList
        from SolidSpheral2d import ScalarField

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


