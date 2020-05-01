from math import sqrt

class Vector2d:
    def __init__(self,x,y):
        self.x = float(x)
        self.y = float(y)
        self.items = [x,y]
    def __getitem__(self,i):
        return self.items[i]
    def __sub__(self,v):
        return Vector2d(self.x-v.x,self.y-v.y)
    def __add__(self,v):
        return Vector2d(self.x+v.x,self.y+v.y)
    def __str__(self):
        return "<"+str(self.x)+","+str(self.y)+">"
    def __repr__(self):
        return str(self)
    def __mul__(self,f):
        return Vector2d(self.x*f,self.y*f)
    def __truediv__(self,f):
        return Vector2d(self.x/f,self.y/f)
    def magnitude(self):
        return sqrt(self.x**2+self.y**2)
    def norm(self):
        mag = self.magnitude()
        return Vector2d(self.x/mag,self.y/mag)
    def lerp(self,v2,t):
        return v2*t+self*(1.0-t)

def vectorfromtuple(tup):
    return Vector2d(tup[0],tup[1])

def interpVector(v1,v2,t):
    return v2*t+v1*(1.0-t)