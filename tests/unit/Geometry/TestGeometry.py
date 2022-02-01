# Python script to test Spheral++ Geometry package.

from Spheral import *

print "-------------------- Vector Tests --------------------"
print "Vector3d(): ", Vector3d()
#print "Vector3d(2.0): ", Vector3d(2.0)

r1 = Vector3d(1, 2, 3)
r2 = Vector3d(5, 10, 20)
print "r1 = ", r1, "\tr2 = ", r2
print "Elements of r1: ", r1.x, r1.y, r1.z

r2.x = 20
print "r2.x = 20", r2

r3 = Vector3d(r1)
print "r3 = ", r3
r3.Zero()
print "r3.Zero() = ", r3

print "r1.dot(r2) = ", r1.dot(r2)
print "r1.cross(r2) = ", r1.cross(r2)
print "r1.magnitude = ", r1.magnitude()
print "r1.magnitude2 = ", r1.magnitude2()
print "r1.minElement = ", r1.minElement()
print "r1.maxElement = ", r1.maxElement()
print "r1.sumElements = ", r1.sumElements()

print "-------------------- Tensor Tests --------------------"
print "Empty tensor constructor: ", Tensor3d()
## print "Tensor3d(5.0) = ", Tensor3d(5.0)

t1 = Tensor3d(1, 2, 3, 4, 5, 6, 7, 8, 9)
t2 = Tensor3d(1, 0, 0, 0, 1, 0, 0, 0, 2)
print "t1 = ", t1
print "t2 = ", t2

print "Elements of t1: ", t1.xx, t1.xy, t1.xz, t1.yx, t1.yy, t1.yz, t1.zx, t1.zy, t1.zz

t2.zz = 10
print "t2.zz = 10, t2 = ", t2

print "t1.getRow(1) = ", t1.getRow(1)
print "t1.getColumn(2) = ", t2.getColumn(2)
r0 = Vector3d(-1, -2, -3)
t3 = Tensor3d(t1)
t3.setRow(1, r0)
print "t1.setRow(1, (-1, -2, -3)) = ", t1
t3 = Tensor3d(t1)
t3.setColumn(1, Vector3d(-1, -2, -3))
print "t1.setColumn(1, (-1, -2, -3)), t1 = ", t1

t3 = Tensor3d(t1)
t3.Zero()
print "t3.Zero() = ", t3

print "t1.Symmetric = ", t1.Symmetric()
print "t1.SkewSymmetric = ", t1.SkewSymmetric()
print "t1.Transpose = ", t1.Transpose()
print "t1.Trace = ", t1.Trace()
print "t1.Determinant = ", t1.Determinant()

print "t1*t2 = ", t1*t2
print "t1*r1 = ", t1*r1
print "t1.doubledot(t2) = ", t1.doubledot(t2)

