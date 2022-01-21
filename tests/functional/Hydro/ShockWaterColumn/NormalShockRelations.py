class NormalShockRelations:
    def __init__(self,
                 M = 1.0,      # Mach number
                 gamma = 1.4): # ratio of specic heats
        self.M = M
        self.g = gamma

    def totalPressure(self): # total pressure in front of shock
        return (1+(self.g-1.0)/2.0*self.M**2) ** (self.g/(self.g-1.0)) 

    def pressureRatio(self): # pressure ratio across shock
        return (2.0*self.g*self.M**2 - self.g + 1.0)/(self.g+1.0)

    def densityRatio(self): # pressure ratio across shock
        return ((self.g+1.0)*self.M**2)/(2.0+(self.g-1.0)*self.M**2)

    def velocityRatio(self): # pressure ratio across shock
        return (2.0+(self.g-1.0)*self.M**2)/((self.g+1.0)*self.M**2)

    def totalPressureRatio(self):
        num1 = (self.g+1.0)*self.M**2
        num2 = (self.g+1.0)
        denom1 = (self.g-1.0)*self.M**2 +2.0
        denom2 = 2.0*self.g*self.M**2 - self.g + 1.0
        exp1 = self.g/(self.g-1.0)
        exp2 = 1.0/(self.g-1.0)
        return (num1/denom1)**exp1 * (num2/denom2)**exp2

    def StagnationPressure(self):
        return self.totalPressure()*self.totalPressureRatio()

    def StagnationPressure2(self):
        num1 = (self.g+1.0)**2 * self.M**2
        num2 = 1.0 - self.g + 2.0*self.g*self.M**2
        denom1 = 4.0*self.g*self.M**2 - 2.0*(self.g-1.0)
        denom2 = self.g+1.0
        exp1 = self.g/(self.g-1.0)
        exp2 = 1.0
        return (num1/denom1)**exp1 * (num2/denom2)**exp2



# if we don't want to do class initializations to get fun
def totalPressure(M,g): 
    return (1+(g-1.0)/2.0*M**2) ** (g/(g-1.0)) 

def pressureRatio(M,g): # pressure ratio across shock
    return (2.0*g*M**2 - g + 1.0)/(g+1.0)
    
def totalPressureRatio(M,g):
    num1 = (g+1.0)*M**2
    num2 = (g+1.0)
    denom1 = (g-1.0)*M**2 +2.0
    denom2 = 2.0*g*M**2 - g + 1.0
    exp1 = g/(g-1.0)
    exp2 = 1.0/(g-1.0)
    return (num1/denom1)**exp1 * (num2/denom2)**exp2

def StagnationPressure(M,g):
    return totalPressure(M,g)*totalPressureRatio(M,g)

def StagnationPressure2(M,g):
    num1 = (g+1.0)**2 * M**2
    num2 = 1.0 - g + 2.0*g*M**2
    denom1 = 4.0*g*M**2 - 2.0*(g-1.0)
    denom2 = g+1.0
    exp1 = g/(g-1.0)
    exp2 = 1.0
    return (num1/denom1)**exp1 * (num2/denom2)**exp2