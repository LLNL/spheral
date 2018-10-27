#-------------------------------------------------------------------------------
# Helper to define methods for pickling/unpickling Spheral types.
# This can be important with pyMPI, since pyMPI relies on picking data types
# for communication!
#-------------------------------------------------------------------------------
import copy_reg, pickle
from SpheralCompiledPackages import *

#-------------------------------------------------------------------------------
# Vector1d
#-------------------------------------------------------------------------------
def reduce_Vector1d(obj):
    return Vector1d, (obj.x,)
copy_reg.pickle(Vector1d, reduce_Vector1d)

#-------------------------------------------------------------------------------
# Vector2d
#-------------------------------------------------------------------------------
def reduce_Vector2d(obj):
    return Vector2d, (obj.x, obj.y)
copy_reg.pickle(Vector2d, reduce_Vector2d)

#-------------------------------------------------------------------------------
# Vector3d
#-------------------------------------------------------------------------------
def reduce_Vector3d(obj):
    return Vector3d, (obj.x, obj.y, obj.z)
copy_reg.pickle(Vector3d, reduce_Vector3d)

#-------------------------------------------------------------------------------
# Tensor1d
#-------------------------------------------------------------------------------
def reduce_Tensor1d(obj):
    return Tensor1d, (obj.xx,)
copy_reg.pickle(Tensor1d, reduce_Tensor1d)

#-------------------------------------------------------------------------------
# Tensor2d
#-------------------------------------------------------------------------------
def reduce_Tensor2d(obj):
    return Tensor2d, (obj.xx, obj.xy,
                      obj.yx, obj.yy)
copy_reg.pickle(Tensor2d(), reduce_Tensor2d)

#-------------------------------------------------------------------------------
# Tensor3d
#-------------------------------------------------------------------------------
def reduce_Tensor3d(obj):
    return Tensor3d, (obj.xx, obj.xy, obj.xz,
                      obj.yx, obj.yy, obj.yz,
                      obj.zx, obj.zy, obj.zz)
copy_reg.pickle(Tensor3d, reduce_Tensor3d)

#-------------------------------------------------------------------------------
# SymTensor1d
#-------------------------------------------------------------------------------
def reduce_SymTensor1d(obj):
    return SymTensor1d, (obj.xx,)
copy_reg.pickle(SymTensor1d, reduce_SymTensor1d)

#-------------------------------------------------------------------------------
# SymTensor2d
#-------------------------------------------------------------------------------
def reduce_SymTensor2d(obj):
    return SymTensor2d, (obj.xx, obj.xy,
                         obj.yx, obj.yy)
copy_reg.pickle(SymTensor2d, reduce_SymTensor2d)

#-------------------------------------------------------------------------------
# SymTensor3d
#-------------------------------------------------------------------------------
def reduce_SymTensor3d(obj):
    return SymTensor3d, (obj.xx, obj.xy, obj.xz,
                         obj.yx, obj.yy, obj.yz,
                         obj.zx, obj.zy, obj.zz)
copy_reg.pickle(SymTensor3d, reduce_SymTensor3d)

#-------------------------------------------------------------------------------
# ThirdRankTensor1d
#-------------------------------------------------------------------------------
def construct_ThirdRankTensor1d(x00):
    result = ThirdRankTensor1d()
    for i in xrange(1):
        result[i] = eval("x%02i" % i)
    return result

def reduce_ThirdRankTensor1d(obj):
    return construct_ThirdRankTensor1d, tuple(obj[i] for i in xrange(1))

copy_reg.pickle(type(ThirdRankTensor1d()), reduce_ThirdRankTensor1d, construct_ThirdRankTensor1d)

#-------------------------------------------------------------------------------
# ThirdRankTensor2d
#-------------------------------------------------------------------------------
def construct_ThirdRankTensor2d(x00, x01, x02, x03, x04, x05, x06, x07):
    result = ThirdRankTensor2d()
    for i in xrange(8):
        result[i] = eval("x%02i" % i)
    return result

def reduce_ThirdRankTensor2d(obj):
    return construct_ThirdRankTensor2d, tuple(obj[i] for i in xrange(8))

copy_reg.pickle(type(ThirdRankTensor2d()), reduce_ThirdRankTensor2d, construct_ThirdRankTensor2d)

#-------------------------------------------------------------------------------
# ThirdRankTensor3d
#-------------------------------------------------------------------------------
def construct_ThirdRankTensor3d(x00, x01, x02, x03, x04, x05, x06, x07, x08, x09,
                                x10, x11, x12, x13, x14, x15, x16, x17, x18, x19,
                                x20, x21, x22, x23, x24, x25, x26):
    result = ThirdRankTensor3d()
    for i in xrange(27):
        result[i] = eval("x%02i" % i)
    return result

def reduce_ThirdRankTensor3d(obj):
    return construct_ThirdRankTensor3d, tuple(obj[i] for i in xrange(27))

copy_reg.pickle(type(ThirdRankTensor3d()), reduce_ThirdRankTensor3d, construct_ThirdRankTensor3d)

#-------------------------------------------------------------------------------
# Box1d
#-------------------------------------------------------------------------------
def construct_Box1d(encoded_string):
    return unpackElementBox1d(encoded_string)

def reduce_Box1d(obj):
    return construct_Box1d, (packElementBox1d(obj),)

copy_reg.pickle(type(Box1d()), reduce_Box1d, construct_Box1d)

#-------------------------------------------------------------------------------
# Polygon
#-------------------------------------------------------------------------------
def construct_Polygon(encoded_string):
    return unpackElementPolygon(encoded_string)

def reduce_Polygon(obj):
    return construct_Polygon, (packElementPolygon(obj),)

copy_reg.pickle(type(Polygon()), reduce_Polygon, construct_Polygon)

#-------------------------------------------------------------------------------
# Polyhedron
#-------------------------------------------------------------------------------
def construct_Polyhedron(encoded_string):
    return unpackElementPolyhedron(encoded_string)

def reduce_Polyhedron(obj):
    return construct_Polyhedron, (packElementPolyhedron(obj),)

copy_reg.pickle(type(Polyhedron()), reduce_Polyhedron, construct_Polyhedron)

#------------------------------------------------------------------------------
# std::vectors
#------------------------------------------------------------------------------
vector_template = """
def construct_vector_of_%(value_type)s(encoded_string):
    return string2vector_of_%(value_type)s(encoded_string)
def reduce_vector_of_%(value_type)s(obj):
    return construct_vector_of_%(value_type)s, (vector2string(obj, 20),)
copy_reg.constructor(construct_vector_of_%(value_type)s)
copy_reg.pickle(vector_of_%(value_type)s, reduce_vector_of_%(value_type)s, construct_vector_of_%(value_type)s)
"""

for t in ("int", "unsigned", "ULL", "double", "string"):
          # "Vector1d", "Vector2d", "Vector3d", 
          # "Tensor1d", "Tensor2d", "Tensor3d", 
          # "SymTensor1d", "SymTensor2d", "SymTensor3d", 
          # "ThirdRankTensor1d", "ThirdRankTensor2d", "ThirdRankTensor3d"):
    exec(vector_template % {"value_type" : t})

