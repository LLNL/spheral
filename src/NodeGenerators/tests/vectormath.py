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
    if len(tup) == 2:
        return Vector2d(tup[0],tup[1])
    elif len(tup) == 3:
        return Vector3d(tup[0],tup[1],tup[2])

def interpVector(v1,v2,t):
    return v2*t+v1*(1.0-t)

class Vector3d:
    def __init__(self,x,y=0,z=0):
        if isinstance(x,int) or isinstance(x,float):
            self.x = float(x)
            self.y = float(y)
            self.z = float(z)
        elif isinstance(x,tuple) or isinstance(x,list) or isinstance(x,vector3d):
            self.x = float(x[0])
            self.y = float(x[1])
            self.z = float(x[2])
        self.magnitude = sqrt(self.x*self.x+self.y*self.y+self.z*self.z)
    def __getitem__(self, key):
        if( key == 0):
            return self.x        
        elif( key == 1):            
            return self.y     
        elif( key == 2):            
            return self.z      
        else:            
            raise Exception("Invalid key to Point")
    def __setitem__(self, key, value):        
        if( key == 0):            
            self.x = value        
        elif( key == 1):            
            self.y = value  
        elif( key == 2):            
            self.z = value      
        else:            
            raise Exception("Invalid key to Point")
    def __add__(self,val):
        return Vector3d(self.x+val.x,self.y+val.y,self.z+val.z)
    def __sub__(self,val):
        return Vector3d(self.x-val.x,self.y-val.y,self.z-val.z)
    def __str__(self):
        return "<"+str(self.x)+","+str(self.y)+","+str(self.z)+">"
    def dot(self,val):
        return self.x*val.x+self.y*val.y+self.z*val.z
    def __mul__(self,val):
        return Vector3d(self.x*val,self.y*val,self.z*val)
    def __truediv__(self,val):
        return Vector3d(self.x/val,self.y/val,self.z/val)
    def cross(self,val):
        x = self.y*val.z - self.z*val.y
        y = self.z*val.x - self.x*val.z
        z = self.x*val.y - self.y*val.x
        return Vector3d(x,y,z)
    def normal(self):
        return Vector3d(self.x,self.y,self.z)/self.magnitude