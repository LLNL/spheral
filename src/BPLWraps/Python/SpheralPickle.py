#-------------------------------------------------------------------------------
# Helper to define methods for pickling/unpickling Spheral types.
# This can be important with pyMPI, since pyMPI relies on picking data types
# for communication!
#-------------------------------------------------------------------------------
import copy_reg

#-------------------------------------------------------------------------------
# Vector1d
#-------------------------------------------------------------------------------
from Geometry import GeomVector1d as Vector1d

def construct_Vector1d(x):
    return Vector1d(x)

def reduce_Vector1d(obj):
    return construct_Vector1d, (obj.x,)

copy_reg.pickle(type(Vector1d()), reduce_Vector1d, construct_Vector1d)

#-------------------------------------------------------------------------------
# Vector2d
#-------------------------------------------------------------------------------
from Geometry import GeomVector2d as Vector2d

def construct_Vector2d(x, y):
    return Vector2d(x, y)

def reduce_Vector2d(obj):
    return construct_Vector2d, (obj.x, obj.y)

copy_reg.pickle(type(Vector2d()), reduce_Vector2d, construct_Vector2d)

#-------------------------------------------------------------------------------
# Vector3d
#-------------------------------------------------------------------------------
from Geometry import GeomVector3d as Vector3d

def construct_Vector3d(x, y, z):
    return Vector3d(x, y, z)

def reduce_Vector3d(obj):
    return construct_Vector3d, (obj.x, obj.y, obj.z)

copy_reg.pickle(type(Vector3d()), reduce_Vector3d, construct_Vector3d)

#-------------------------------------------------------------------------------
# Tensor1d
#-------------------------------------------------------------------------------
from Geometry import GeomTensor1d as Tensor1d

def construct_Tensor1d(xx):
    return Tensor1d(xx)

def reduce_Tensor1d(obj):
    return construct_Tensor1d, (obj.xx,)

copy_reg.pickle(type(Tensor1d()), reduce_Tensor1d, construct_Tensor1d)

#-------------------------------------------------------------------------------
# Tensor2d
#-------------------------------------------------------------------------------
from Geometry import GeomTensor2d as Tensor2d

def construct_Tensor2d(xx, xy,
                       yx, yy):
    return Tensor2d(xx, xy,
                    yx, yy)

def reduce_Tensor2d(obj):
    return construct_Tensor2d, (obj.xx, obj.xy,
                                obj.yx, obj.yy)

copy_reg.pickle(type(Tensor2d()), reduce_Tensor2d, construct_Tensor2d)

#-------------------------------------------------------------------------------
# Tensor3d
#-------------------------------------------------------------------------------
from Geometry import GeomTensor3d as Tensor3d

def construct_Tensor3d(xx, xy, xz,
                       yx, yy, yz,
                       zx, zy, zz):
    return Tensor3d(xx, xy, xz,
                    yx, yy, yz,
                    zx, zy, zz)

def reduce_Tensor3d(obj):
    return construct_Tensor3d, (obj.xx, obj.xy, obj.xz,
                                obj.yx, obj.yy, obj.yz,
                                obj.zx, obj.zy, obj.zz)

copy_reg.pickle(type(Tensor3d()), reduce_Tensor3d, construct_Tensor3d)

#-------------------------------------------------------------------------------
# SymTensor1d
#-------------------------------------------------------------------------------
from Geometry import GeomSymTensor1d as SymTensor1d

def construct_SymTensor1d(xx):
    return SymTensor1d(xx)

def reduce_SymTensor1d(obj):
    return construct_SymTensor1d, (obj.xx,)

copy_reg.pickle(type(SymTensor1d()), reduce_SymTensor1d, construct_SymTensor1d)

#-------------------------------------------------------------------------------
# SymTensor2d
#-------------------------------------------------------------------------------
from Geometry import GeomSymTensor2d as SymTensor2d

def construct_SymTensor2d(xx, xy,
                          yx, yy):
    return SymTensor2d(xx, xy,
                       yx, yy)

def reduce_SymTensor2d(obj):
    return construct_SymTensor2d, (obj.xx, obj.xy,
                                   obj.yx, obj.yy)

copy_reg.pickle(type(SymTensor2d()), reduce_SymTensor2d, construct_SymTensor2d)

#-------------------------------------------------------------------------------
# SymTensor3d
#-------------------------------------------------------------------------------
from Geometry import GeomSymTensor3d as SymTensor3d

def construct_SymTensor3d(xx, xy, xz,
                          yx, yy, yz,
                          zx, zy, zz):
    return SymTensor3d(xx, xy, xz,
                       yx, yy, yz,
                       zx, zy, zz)

def reduce_SymTensor3d(obj):
    return construct_SymTensor3d, (obj.xx, obj.xy, obj.xz,
                                   obj.yx, obj.yy, obj.yz,
                                   obj.zx, obj.zy, obj.zz)

copy_reg.pickle(type(SymTensor3d()), reduce_SymTensor3d, construct_SymTensor3d)

