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
                for i in range(eoscount):
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
        for i in range(len(self.soln)):
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
        for i in range(len(self.soln)):
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

class EarthLikeProfileConstantTemp3d():
    # this version will first solve inward to get central density, then outward to fix the total mass
    def __init__(self,
                 rho0,      # Density at Radius
                 rMax,      # Radius
                 M0,        # Mass at Radius
                 temp,      # Temperature throughout
                 eostup,    # tuple that indicates how materials/eos change
                 units,
                 nbins=1000):
        
        self.soln = []
        self.rho0 = rho0
        
        from SolidSpheral3d import makeVoidNodeList
        from SolidSpheral3d import ScalarField
        
        eoscount    = len(eostup)/2
        
        nodes   = makeVoidNodeList("nodes", numInternal=1)
        ef      = ScalarField("eps", nodes)
        Kf      = ScalarField("mod", nodes)
        Pf      = ScalarField("pressure", nodes)
        rhof    = ScalarField("rho", nodes)
        tempf   = ScalarField("temp", nodes)
        
        # get the eos for this radius
        if(eoscount>1):
            eos = eostup[2*(eoscount-1)]
        else:
            eos = eostup[0]
        
        rhof[0] = rho0
        eos.setSpecificThermalEnergy(ef,rhof,tempf)
        e       = ef[0]
        eos.setBulkModulus(Kf,rhof,ef)
        K       = Kf[0]
        
        y0  = -M0*units.G/(rMax**2)*(rho0**2)/K
        
        r   = rMax
        rho = rho0
        dr  = rMax/nbins
        y   = y0
        dy  = 0
        
        tempf[0] = temp
        eosold = None
        #eosSwitch = False
        step = 0

        while(r>0):
            # get the eos for this radius
            if(eoscount>1):
                for i in range(eoscount):
                    ermin = eostup[2*i+1][0]
                    ermax = eostup[2*i+1][1]
                    if(r<=ermax and r>=ermin):
                        if eos is not None:
                            eosold = eos
                        eos = eostup[2*i]
                        break
            else:
                if eos is not None:
                    eosold = eos
                eos = eostup[0]

            if step>0:
                #print "Switching eos at r=%e" % r
                #print eosold,eos
                # compute a new rho based on pressure and the new eos
                # first get old rho -> pressure
                rhof[0] = rho
                eosold.setSpecificThermalEnergy(ef,rhof,tempf)
                eosold.setPressure(Pf,rhof,ef)
                P = Pf[0]
                #print "P=%e" % P
                # now root-find for new density based on this pressure
                tol = 0.0001
                d = 1.0
                iter = 0
                #print "old rho was %f" % rho
                rhoold = rho
                #rho = 10.0*rho
                while ((abs(d)>tol) and (iter < 1000)):
                    if (d>1.0 and iter ==1 ):
                        print("i think eos just changed?")
                        #rho = 10.0*rhoold
                    rhof[0] = rho
                    eos.setSpecificThermalEnergy(ef,rhof,tempf)
                    eos.setPressure(Pf,rhof,ef)
                    eos.setBulkModulus(Kf,rhof,ef)
                    Pn = Pf[0]
                    Kn = Kf[0]
                    d = (Pn-P)/Kn
                    #print "Pn=%e Kn=%e d=%e rho=%e e=%e" % (Pn,Kn,d,rho,ef[0])
                    rho *= (1.0-d)
                    iter += 1
                #print "new rho is %f after %d iterations, d was %f" % (rho,iter,d)
                    


            rhof[0] = rho
            eos.setSpecificThermalEnergy(ef,rhof,tempf)
            e       = ef[0]
            eos.setBulkModulus(Kf,rhof,ef)
            K       = Kf[0]
            
            #print "dy, dr, rho, y, r, K = {0:3.3e} {1:3.3e} {2:3.3e} {3:3.3e} {4:3.3e} {5:3.3e}".format(dy,dr,rho,y,r,K)

            dy      = dr*(3.0/rho*y*y - 2.0/r*y - units.G/K*4.0*pi*pow(rho,3.0))
            #self.soln.append([r,rho])
            r       = r - dr
            if (r>0):
                rho = rho - y*dr
                y   = y + dy
            step    += 1
        
        print("\n\n\n\nNow Forward...\n\n\n\n")
        # got central density, now solve outward until Mtot = M0
        self.soln.append([0,rho])
        Mt  = 0
        r   = dr
        step = 0
        eosold = None
        
        
        while(Mt<=M0):
            Mt = Mt + 4.0*pi*r*r*rho*dr
            # get the eos for this radius
            if(eoscount>1):
                for i in range(eoscount):
                    ermin = eostup[2*i+1][0]
                    ermax = eostup[2*i+1][1]
                    if(r<=ermax and r>=ermin):
                        if eos is not None:
                            eosold = eos
                        eos = eostup[2*i]
                        break
            else:
                if eos is not None:
                    eosold = eos
                eos = eostup[0]
            
            if step>0:
                #print "Switching eos at r=%e" % r
                #print eosold,eos
                # compute a new rho based on pressure and the new eos
                # first get old rho -> pressure
                rhof[0] = rho
                eosold.setSpecificThermalEnergy(ef,rhof,tempf)
                eosold.setPressure(Pf,rhof,ef)
                P = Pf[0]
                #print "P=%e" % P
                # now root-find for new density based on this pressure
                tol = 0.0001
                d = 1.0
                iter = 0
                #print "old rho was %f" % rho
                while ((abs(d)>tol) and (iter < 1000)):
                    if (abs(d) > 1 and iter ==1):
                        print("i think the eos just changed?")
                    rhof[0] = rho
                    eos.setSpecificThermalEnergy(ef,rhof,tempf)
                    eos.setPressure(Pf,rhof,ef)
                    eos.setBulkModulus(Kf,rhof,ef)
                    Pn = Pf[0]
                    #print "Pn=%e" % Pn
                    Kn = Kf[0]
                    d = (Pn-P)/Kn
                    rho *= (1.0-d)
                    iter += 1
            #print "new rho is %f after %d iterations, d was %f" % (rho,iter,d)
            
            rhof[0] = rho
            eos.setSpecificThermalEnergy(ef,rhof,tempf)
            e       = ef[0]
            eos.setBulkModulus(Kf,rhof,ef)
            K       = Kf[0]
            
            #print "dy, dr, rho, y, r, Mt, K = {0:3.3e} {1:3.3e} {2:3.3e} {3:3.3e} {4:3.3e} {5:3.3e} {6:3.3e}".format(dy,dr,rho,y,r,Mt,K)
            dy      = dr*(3.0/rho*y*y - 2.0/r*y - units.G/K*4.0*pi*pow(rho,3.0))
            #self.soln.append([r,rho])
            y       = y + dy
            rho     = rho + y*dr
            self.soln.append([r,rho])
            r       = r + dr
            step    += 1
        
        self.soln.sort()
        self.rMax = r - dr  # to call inside the script to reset rMax

        totM = 0.0
        for i in range(len(self.soln)-1):
            r1 = self.soln[i+1][0]
            r0 = self.soln[i][0]
            f1 = self.soln[i+1][1]*r1*r1
            f0 = self.soln[i][0]*r0*r0
            totM += 4.0*pi*0.5*(f1+f0)*(r1-r0)
        self.totalMass = totM

    def __call__(self,r):
        rho = self.rho0
        for i in range(len(self.soln)):
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

class HydroStaticProfileConstantTemp3d():
    # this version will first solve inward to get central density, then outward to fix the total mass
    def __init__(self,
                 rho0,      # Density at Radius
                 rMax,      # Radius
                 M0,        # Mass at Radius
                 temp,      # Temperature throughout
                 eostup,    # tuple that indicates how materials/eos change
                 units,
                 nbins=1000):
        
        self.soln = []
        self.rho0 = rho0
        
        from SolidSpheral3d import makeVoidNodeList
        from SolidSpheral3d import ScalarField
        
        eoscount    = len(eostup)/2
        
        nodes   = makeVoidNodeList("nodes", numInternal=1)
        ef      = ScalarField("eps", nodes)
        Kf      = ScalarField("mod", nodes)
        Pf      = ScalarField("pressure", nodes)
        rhof    = ScalarField("rho", nodes)
        tempf   = ScalarField("temp", nodes)
        
        # get the eos for this radius
        if(eoscount>1):
            eos = eostup[2*(eoscount-1)]
        else:
            eos = eostup[0]
        
        rhof[0] = rho0
        eos.setSpecificThermalEnergy(ef,rhof,tempf)
        e       = ef[0]
        eos.setBulkModulus(Kf,rhof,ef)
        K       = Kf[0]
        
        y0  = -M0*units.G/(rMax**2)*(rho0**2)/K
        
        r   = rMax
        rho = rho0
        dr  = rMax/nbins
        y   = y0
        dy  = 0
        
        tempf[0] = temp
        eosold = None
        eosSwitch = False
        
        while(r>0):
            # get the eos for this radius
            if(eoscount>1):
                for i in range(eoscount):
                    ermin = eostup[2*i+1][0]
                    ermax = eostup[2*i+1][1]
                    if(r<=ermax and r>=ermin):
                        if r!=rMax:
                            if ((eos != eosold) and (eosold is not None) and (eosSwitch == False)):
                                eosSwitch = True
                            else:
                                eosSwitch = False
                                eosold = eos
                        
                        eos = eostup[2*i]
                        break
            else:
                eos = eostup[0]
            
            if eosSwitch:
                print("Switching eos at r=%e" % r)
                #print eosold,eos
                # compute a new rho based on pressure and the new eos
                # first get old rho -> pressure
                rhof[0] = rho
                eosold.setSpecificThermalEnergy(ef,rhof,tempf)
                eosold.setPressure(Pf,rhof,ef)
                P = Pf[0]
                #print "P=%e" % P
                # now root-find for new density based on this pressure
                tol = 0.0001
                d = 1.0
                iter = 0
                print("old rho was %f" % rho)
                while ((abs(d)>tol) and (iter < 1000)):
                    rhof[0] = rho
                    eos.setSpecificThermalEnergy(ef,rhof,tempf)
                    eos.setPressure(Pf,rhof,ef)
                    eos.setBulkModulus(Kf,rhof,ef)
                    Pn = Pf[0]
                    #print "Pn=%e" % Pn
                    Kn = Kf[0]
                    d = (Pn-P)/Kn
                    rho *= (1.0-d)
                    iter += 1
                print("new rho is %f after %d iterations, d was %f" % (rho,iter,d))
            
            
            rhof[0] = rho
            eos.setSpecificThermalEnergy(ef,rhof,tempf)
            e       = ef[0]
            eos.setBulkModulus(Kf,rhof,ef)
            K       = Kf[0]
            
            #print "dy, dr, rho, y, r, K = {0:3.3e} {1:3.3e} {2:3.3e} {3:3.3e} {4:3.3e} {5:3.3e}".format(dy,dr,rho,y,r,K)
            
            dy      = dr*(2.0/rho*y*y - 2.0/r*y - units.G/K*4.0*pi*pow(rho,3.0))
            #self.soln.append([r,rho])
            y       = y + dy
            rho     = rho - y*dr
            r       = r - dr
        
        #print "Now Forward..."
        # got central density, now solve outward until Mtot = M0
        self.soln.append([0,rho])
        Mt  = 0
        r   = dr
        eosSwitch = False
        eosold = None
        
        
        while(Mt<=M0):
            Mt = Mt + 4.0*pi*r*r*rho*dr
            # get the eos for this radius
            if(eoscount>1):
                for i in range(eoscount):
                    ermin = eostup[2*i+1][0]
                    ermax = eostup[2*i+1][1]
                    if(r<=ermax and r>=ermin):
                        if r>0:
                            if ((eos != eosold) and (eosold is not None) and (eosSwitch == False)):
                                eosSwitch = True
                            else:
                                eosSwitch = False
                                eosold = eos
                        
                        eos = eostup[2*i]
                        break
            else:
                eos = eostup[0]
            
            if eosSwitch:
                print("Switching eos at r=%e" % r)
                #print eosold,eos
                # compute a new rho based on pressure and the new eos
                # first get old rho -> pressure
                rhof[0] = rho
                eosold.setSpecificThermalEnergy(ef,rhof,tempf)
                eosold.setPressure(Pf,rhof,ef)
                P = Pf[0]
                #print "P=%e" % P
                # now root-find for new density based on this pressure
                tol = 0.0001
                d = 1.0
                iter = 0
                print("old rho was %f" % rho)
                while ((abs(d)>tol) and (iter < 1000)):
                    rhof[0] = rho
                    eos.setSpecificThermalEnergy(ef,rhof,tempf)
                    eos.setPressure(Pf,rhof,ef)
                    eos.setBulkModulus(Kf,rhof,ef)
                    Pn = Pf[0]
                    #print "Pn=%e" % Pn
                    Kn = Kf[0]
                    d = (Pn-P)/Kn
                    rho *= (1.0-d)
                    iter += 1
                print("new rho is %f after %d iterations, d was %f" % (rho,iter,d))
            
            rhof[0] = rho
            eos.setSpecificThermalEnergy(ef,rhof,tempf)
            e       = ef[0]
            eos.setBulkModulus(Kf,rhof,ef)
            K       = Kf[0]
            
            #print "dy, dr, rho, y, r, Mt, K = {0:3.3e} {1:3.3e} {2:3.3e} {3:3.3e} {4:3.3e} {5:3.3e} {6:3.3e}".format(dy,dr,rho,y,r,Mt,K)
            dy      = dr*(2.0/rho*y*y - 2.0/r*y - units.G/K*4.0*pi*pow(rho,3.0))
            #self.soln.append([r,rho])
            y       = y + dy
            rho     = rho + y*dr
            self.soln.append([r,rho])
            r       = r + dr
        
        self.soln.sort()
        self.rMax = r - dr  # to call inside the script to reset rMax
    
    def __call__(self,r):
        rho = self.rho0
        for i in range(len(self.soln)):
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


class OldHydroStaticProfileConstantTemp3d():
    def __init__(self,
                 rho0,      # Density at Radius
                 rMax,      # Radius
                 M0,        # Mass at Radius
                 temp,      # Temperature throughout
                 eostup,    # tuple that indicates how materials/eos change
                 units,
                 nbins=1000):
        
        self.soln = []
        self.rho0 = rho0
        
        from SolidSpheral3d import makeVoidNodeList
        from SolidSpheral3d import ScalarField
        
        eoscount    = len(eostup)/2
        
        nodes   = makeVoidNodeList("nodes", numInternal=1)
        ef      = ScalarField("eps", nodes)
        Kf      = ScalarField("mod", nodes)
        rhof    = ScalarField("rho", nodes)
        tempf   = ScalarField("temp", nodes)
        
        # get the eos for this radius
        if(eoscount>1):
            eos = eostup[2*(eoscount-1)]
        else:
            eos = eostup[0]
        
        rhof[0] = rho0
        eos.setSpecificThermalEnergy(ef,rhof,tempf)
        e       = ef[0]
        eos.setBulkModulus(Kf,rhof,ef)
        K       = Kf[0]
        
        y0  = -M0*units.G/(rMax**2)*(rho0**2)/K
        
        r   = rMax
        rho = rho0
        dr  = rMax/nbins
        y   = y0
        dy  = 0
        
        tempf[0] = temp
        
        while(r>0):
            # get the eos for this radius
            if(eoscount>1):
                for i in range(eoscount):
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
            
            print("dy, dr, rho, y, G, K = {0:3.3e} {1:3.3e} {2:3.3e} {3:3.3e} {4:3.3e} {5:3.3e}".format(dy,dr,rho,y,units.G,K))
            
            dy      = dr*(2.0/rho*y*y - 2.0/r*y - units.G/K*4.0*pi*pow(rho,3.0))
            self.soln.append([r,rho])
            y       = y + dy
            rho     = rho - y*dr
            r       = r - dr
        
        self.soln.sort()
    
    def __call__(self,r):
        rho = self.rho0
        for i in range(len(self.soln)):
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
                 rho0,      # Density at Radius
                 rMax,      # Radius
                 M0,        # Mass at Radius
                 temp,      # Temperature throughout
                 eostup,    # tuple that indicates how materials/eos change
                 units,
                 nbins=1000):
        
        self.soln = []
        self.rho0 = rho0
        
        from SolidSpheral2d import makeVoidNodeList
        from SolidSpheral2d import ScalarField

        eoscount    = len(eostup)/2

        nodes   = makeVoidNodeList("nodes", numInternal=1)
        ef      = ScalarField("eps", nodes)
        Kf      = ScalarField("mod", nodes)
        rhof    = ScalarField("rho", nodes)
        tempf   = ScalarField("temp", nodes)

        # get the eos for this radius
        if(eoscount>1):
            eos = eostup[2*(eoscount-1)]
        else:
            eos = eostup[0]

        rhof[0] = rho0
        eos.setSpecificThermalEnergy(ef,rhof,tempf)
        e       = ef[0]
        eos.setBulkModulus(Kf,rhof,ef)
        K       = Kf[0]

        y0  = -M0*units.G/(rMax)*(rho0**2)/K
        
        r   = rMax
        rho = rho0
        dr  = rMax/nbins
        y   = y0
        dy  = 0
        
        tempf[0] = temp
        
        while(r>0):
            # get the eos for this radius
            if(eoscount>1):
                for i in range(eoscount):
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
            #self.soln.append([r,rho])
            y       = y + dy
            rho     = rho - y*dr
            r       = r - dr
    
        print("Now Forward...")
        # got central density, now solve outward until Mtot = M0
        self.soln.append([0,rho])
        Mt  = 0
        r   = dr
        
        while(Mt<=M0):
            Mt = Mt + 2.0*pi*r*rho*dr
            # get the eos for this radius
            if(eoscount>1):
                for i in range(eoscount):
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
            
            print("dy, dr, rho, y, r, Mt, K = {0:3.3e} {1:3.3e} {2:3.3e} {3:3.3e} {4:3.3e} {5:3.3e} {6:3.3e}".format(dy,dr,rho,y,r,Mt,K))
            dy      = dr*(2.0/rho*y*y - 2.0/r*y - units.G/K*4.0*pi*pow(rho,3.0))
            #self.soln.append([r,rho])
            y       = y + dy
            rho     = rho + y*dr
            self.soln.append([r,rho])
            r       = r + dr
        
        self.soln.sort()
        self.rMax = r - dr  # to call inside the script to reset rMax
    
        self.soln.sort()
    
    def __call__(self,r):
        rho = self.rho0
        for i in range(len(self.soln)):
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

class OldHydroStaticProfileConstantTemp2d():
    def __init__(self,
                 rho0,      # Density at Radius
                 rMax,      # Radius
                 M0,        # Mass at Radius
                 temp,      # Temperature throughout
                 eostup,    # tuple that indicates how materials/eos change
                 units,
                 nbins=1000):
        
        self.soln = []
        self.rho0 = rho0
        
        from SolidSpheral2d import makeVoidNodeList
        from SolidSpheral2d import ScalarField
        
        eoscount    = len(eostup)/2
        
        nodes   = makeVoidNodeList("nodes", numInternal=1)
        ef      = ScalarField("eps", nodes)
        Kf      = ScalarField("mod", nodes)
        rhof    = ScalarField("rho", nodes)
        tempf   = ScalarField("temp", nodes)
        
        # get the eos for this radius
        if(eoscount>1):
            eos = eostup[2*(eoscount-1)]
        else:
            eos = eostup[0]
        
        rhof[0] = rho0
        eos.setSpecificThermalEnergy(ef,rhof,tempf)
        e       = ef[0]
        eos.setBulkModulus(Kf,rhof,ef)
        K       = Kf[0]
        
        y0  = -M0*units.G/(rMax)*(rho0**2)/K
        
        r   = rMax
        rho = rho0
        dr  = rMax/nbins
        y   = y0
        dy  = 0
        
        tempf[0] = temp
        
        while(r>0):
            # get the eos for this radius
            if(eoscount>1):
                for i in range(eoscount):
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
        for i in range(len(self.soln)):
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

class EarthLikeProfileConstantTemp2d():
    # this version will first solve inward to get central density, then outward to fix the total mass
    def __init__(self,
                 rho0,      # Density at Radius
                 rMax,      # Radius
                 M0,        # Mass at Radius
                 temp,      # Temperature throughout
                 eostup,    # tuple that indicates how materials/eos change
                 units,
                 nbins=1000):
        
        self.soln = []
        self.rho0 = rho0
        
        from SolidSpheral2d import makeVoidNodeList
        from SolidSpheral2d import ScalarField
        
        eoscount    = len(eostup)/2
        
        nodes   = makeVoidNodeList("nodes", numInternal=1)
        ef      = ScalarField("eps", nodes)
        Kf      = ScalarField("mod", nodes)
        Pf      = ScalarField("pressure", nodes)
        rhof    = ScalarField("rho", nodes)
        tempf   = ScalarField("temp", nodes)
        
        # get the eos for this radius
        if(eoscount>1):
            eos = eostup[2*(eoscount-1)]
        else:
            eos = eostup[0]
        
        rhof[0] = rho0
        eos.setSpecificThermalEnergy(ef,rhof,tempf)
        e       = ef[0]
        eos.setBulkModulus(Kf,rhof,ef)
        K       = Kf[0]
        
        y0  = -M0*units.G/(rMax)*(rho0**2)/K
        
        r   = rMax
        rho = rho0
        dr  = rMax/nbins
        y   = y0
        dy  = 0
        
        tempf[0] = temp
        eosold = None
        #eosSwitch = False
        step = 0
        
        while(r>0):
            # get the eos for this radius
            if(eoscount>1):
                for i in range(eoscount):
                    ermin = eostup[2*i+1][0]
                    ermax = eostup[2*i+1][1]
                    if(r<=ermax and r>=ermin):
                        if eos is not None:
                            eosold = eos
                        eos = eostup[2*i]
                        break
            else:
                if eos is not None:
                    eosold = eos
                eos = eostup[0]
            
            if step>0:
                #print "Switching eos at r=%e" % r
                #print eosold,eos
                # compute a new rho based on pressure and the new eos
                # first get old rho -> pressure
                rhof[0] = rho
                eosold.setSpecificThermalEnergy(ef,rhof,tempf)
                eosold.setPressure(Pf,rhof,ef)
                P = Pf[0]
                #print "P=%e" % P
                # now root-find for new density based on this pressure
                tol = 0.0001
                d = 1.0
                iter = 0
                #print "old rho was %f" % rho
                while ((abs(d)>tol) and (iter < 1000)):
                    rhof[0] = rho
                    eos.setSpecificThermalEnergy(ef,rhof,tempf)
                    eos.setPressure(Pf,rhof,ef)
                    eos.setBulkModulus(Kf,rhof,ef)
                    Pn = Pf[0]
                    #print "Pn=%e" % Pn
                    Kn = Kf[0]
                    d = (Pn-P)/Kn
                    rho *= (1.0-d)
                    iter += 1
            #print "new rho is %f after %d iterations, d was %f" % (rho,iter,d)
            
            
            rhof[0] = rho
            eos.setSpecificThermalEnergy(ef,rhof,tempf)
            e       = ef[0]
            eos.setBulkModulus(Kf,rhof,ef)
            K       = Kf[0]
            
            #print "dy, dr, rho, y, r, K = {0:3.3e} {1:3.3e} {2:3.3e} {3:3.3e} {4:3.3e} {5:3.3e}".format(dy,dr,rho,y,r,K)
            
            dy      = dr*(2.0/rho*y*y - 1.0/r*y - units.G/K*2.0*pi*pow(rho,2.0))
            #self.soln.append([r,rho])
            y       = y + dy
            rho     = rho - y*dr
            r       = r - dr
            step    += 1
        
        #print "Now Forward..."
        # got central density, now solve outward until Mtot = M0
        self.soln.append([0,rho])
        Mt  = 0
        r   = dr
        step = 0
        eosold = None
        
        
        while(Mt<=M0):
            Mt = Mt + 4.0*pi*r*r*rho*dr
            # get the eos for this radius
            if(eoscount>1):
                for i in range(eoscount):
                    ermin = eostup[2*i+1][0]
                    ermax = eostup[2*i+1][1]
                    if(r<=ermax and r>=ermin):
                        if eos is not None:
                            eosold = eos
                        eos = eostup[2*i]
                        break
            else:
                if eos is not None:
                    eosold = eos
                eos = eostup[0]
            
            if step>0:
                #print "Switching eos at r=%e" % r
                #print eosold,eos
                # compute a new rho based on pressure and the new eos
                # first get old rho -> pressure
                rhof[0] = rho
                eosold.setSpecificThermalEnergy(ef,rhof,tempf)
                eosold.setPressure(Pf,rhof,ef)
                P = Pf[0]
                #print "P=%e" % P
                # now root-find for new density based on this pressure
                tol = 0.0001
                d = 1.0
                iter = 0
                #print "old rho was %f" % rho
                while ((abs(d)>tol) and (iter < 1000)):
                    rhof[0] = rho
                    eos.setSpecificThermalEnergy(ef,rhof,tempf)
                    eos.setPressure(Pf,rhof,ef)
                    eos.setBulkModulus(Kf,rhof,ef)
                    Pn = Pf[0]
                    #print "Pn=%e" % Pn
                    Kn = Kf[0]
                    d = (Pn-P)/Kn
                    rho *= (1.0-d)
                    iter += 1
            #print "new rho is %f after %d iterations, d was %f" % (rho,iter,d)
            
            rhof[0] = rho
            eos.setSpecificThermalEnergy(ef,rhof,tempf)
            e       = ef[0]
            eos.setBulkModulus(Kf,rhof,ef)
            K       = Kf[0]
            
            #print "dy, dr, rho, y, r, Mt, K = {0:3.3e} {1:3.3e} {2:3.3e} {3:3.3e} {4:3.3e} {5:3.3e} {6:3.3e}".format(dy,dr,rho,y,r,Mt,K)
            dy      = dr*(2.0/rho*y*y - 1.0/r*y - units.G/K*2.0*pi*pow(rho,2.0))
            #self.soln.append([r,rho])
            y       = y + dy
            rho     = rho + y*dr
            self.soln.append([r,rho])
            r       = r + dr
            step    += 1
        
        self.soln.sort()
        self.rMax = r - dr  # to call inside the script to reset rMax
        
        totM = 0.0
        for i in range(len(self.soln)-1):
            r1 = self.soln[i+1][0]
            r0 = self.soln[i][0]
            f1 = self.soln[i+1][1]*r1*r1
            f0 = self.soln[i][0]*r0*r0
            totM += 4.0*pi*0.5*(f1+f0)*(r1-r0)
        self.totalMass = totM
    
    def __call__(self,r):
        rho = self.rho0
        for i in range(len(self.soln)):
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



