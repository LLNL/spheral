t2d = [["xx", "xy"], ["yx", "yy"]]
t3d = [["xx", "xy", "xz"], ["yx", "yy", "yz"], ["zx", "zy", "zz"]]

print "matrix products:\n\n"
for i in xrange(2):
    for j in xrange(2):
        result = t2d[i][j] + " = "
        for k in xrange(2):
            result += "(this->m%s)*(rhs.%s()) + " % (t2d[i][k], t2d[k][j])
        print result

print "\n\n"
for i in xrange(3):
    for j in xrange(3):
        result = t3d[i][j] + " = "
        for k in xrange(3):
            result += "(this->m%s)*(rhs.%s()) + " % (t3d[i][k], t3d[k][j])
        print result

print "\n\ndouble dots:\n\n"

result = "dd2d = "
for i in xrange(2):
    for j in xrange(2):
        result += "(this->m%s)*(rhs.%s()) + " % (t2d[i][j], t2d[j][i])
print result

result = "dd3d = "
for i in xrange(3):
    for j in xrange(3):
        result += "(this->m%s)*(rhs.%s()) + " % (t3d[i][j], t3d[j][i])
print result

