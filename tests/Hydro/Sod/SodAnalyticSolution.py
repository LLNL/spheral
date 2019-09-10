from math import *

################################################################################
# Provide the fortran sign function.
################################################################################
def sign(x,y):
    sgn = y/(abs(y) + 1e-20)
    return x*sgn

################################################################################
# Find the roots of a function -- taken from Numerical Recipes.
################################################################################
def zbrent(func,x1,x2,tol):
    itmax = 100
    eps = 3.e-8
    a=x1
    b=x2
    fa=func(a)
    fb=func(b)
    if(fb*fa > 0.):
        print 'root must be bracketed for zbrent.'
        return 0.0
    fc=fb
    for iter in xrange(itmax):
        if(fb*fc >= 0.):
            c=a
            fc=fa
            d=b-a
            e=d
        if(abs(fc) > abs(fb)):
            a=b
            b=c
            c=a
            fa=fb
            fb=fc
            fc=fa
        tol1=2.*eps*abs(b)+0.5*tol
        xm=.5*(c-b)
        if(abs(xm) <= tol1 or fb == 0.):
            return b
        if (abs(e) >= tol1 and abs(fa) > abs(fb)):
            s=fb/fa
            if(a == c):
                p=2.*xm*s
                q=1.-s
            else:
                q=fa/fc
                r=fb/fc
                p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
                q=(q-1.)*(r-1.)*(s-1.)
            if(p > 0.):
                q=-q
            p=abs(p)
            if (2.*p < min(3.*xm*q-abs(tol1*q),abs(e*q))):
                e=d
                d=p/q
            else:
                d=xm
                e=d
        else:
            d=xm
            e=d
        a=b
        fa=fb
        if(abs(d) > tol1):
            b=b+d
        else:
            b = b + sign(tol1,xm)
        fb=func(b)
    
    print 'zbrent exceeding maximum iterations.'
    return b

################################################################################
# A class to provide the analytic solution for the Sod problem
################################################################################
class SodSolution:

    # Constructor.
    def __init__(self,
                 nPoints=101,
                 gamma=1.4,
                 rho1=1.0, P1=1.0, v1=0.0,
                 rho2=0.25, P2=0.1795, v2=0.0,
                 x0=-0.5, x1=0.0, x2=0.5,
                 h1=0.01, h2=0.01):
        self.gamma = gamma
        self.nPoints = nPoints
        self.rho1 = rho1
        self.P1 = P1
        self.v1 = v1
        self.rho2 = rho2
        self.P2 = P2
        self.v2 = v2
        self.x0 = x0
        self.x1 = x1
        self.x2 = x2
        self.h1 = h1
        self.h2 = h2
        self.dx = (x2 - x0)/(nPoints - 1)
        self.cs1 = sqrt(gamma*P1/rho1)
        self.cs2 = sqrt(gamma*P2/rho2)
        self.mu2 = (gamma - 1.0)/(gamma + 1.0)
        self.Pm = zbrent(self.func, P2, P1, 0.5*(P1+P2)*0.0001)
        self.rhom1 = rho1 * (self.Pm/P1)**(1.0/gamma)
        self.rhom2 = rho2 * (self.Pm + self.mu2*P2)/(P2 + self.mu2*self.Pm)
        self.vm = 2.*self.cs1/(gamma - 1.) * (1. - (self.Pm/P1)**(0.5*(gamma-1.0)/gamma))
        self.vt = self.cs1 - self.vm/(1. - self.mu2)
        thpt = 1. - rho2/self.rhom2
        self.vs = self.vm*thpt/(1.0e-50 + thpt*thpt)
        self.x = []
        for i in xrange(self.nPoints):
            self.x.append(x0 + i*self.dx)

        return

    def solution(self, time,
                 x = None):

        # Did the user specify the x coordinates?
        if not x is None:
            self.x = x
            self.nPoints = len(x)
        assert len(self.x) == self.nPoints

        # Make the return arrays.
        v = [0.0]*self.nPoints
        u = [0.0]*self.nPoints
        rho = [0.0]*self.nPoints
        P = [0.0]*self.nPoints
        h = [0.0]*self.nPoints

        i = 0
        for xx in self.x:
            x = xx - self.x1
            if (x <= -self.cs1*time):
                rho[i] = self.rho1
                P[i] = self.P1
                v[i] = self.v1
                h[i] = self.h1
            elif (x <= -self.vt*time):
                thpt = (-self.mu2*x/(self.cs1*time) + (1. - self.mu2))
                rho[i] = self.rho1*thpt**(2./(self.gamma - 1.))
                P[i] = self.P1*thpt**(2.0*self.gamma/(self.gamma - 1.))
                v[i] = (1. - self.mu2)*(x/time + self.cs1)
                h[i] = self.h1*self.rho1/rho[i]
            elif (x <= self.vm*time):
                rho[i] = self.rhom1
                P[i] = self.Pm
                v[i] = self.vm
                h[i] = self.h1*self.rho1/rho[i]
            elif (x <= self.vs*time):
                rho[i] = self.rhom2
                P[i] = self.Pm
                v[i] = self.vm
                h[i] = self.h2*self.rho2/rho[i]
            else:
                rho[i] = self.rho2
                P[i] = self.P2
                v[i] = 0.
                h[i] = self.h2
            u[i] = P[i]/((self.gamma - 1.0)*rho[i])
            i = i + 1

        return self.x, v, u, rho, P, h

    def func(self, Pm):
        return ((Pm/self.P2 - 1.)*sqrt((1. - self.mu2)/
                                       (self.gamma*(Pm/self.P2 + self.mu2))) - 
                2./(self.gamma - 1.)*self.cs1/self.cs2*
                (1. - (Pm/self.P1)**(0.5*(self.gamma - 1.)/self.gamma)))
