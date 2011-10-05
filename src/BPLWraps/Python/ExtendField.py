# ------------------------------------------------------------------------------
# Define the base class for extending field types.
# ------------------------------------------------------------------------------
class ExtendFieldBase:
    RealFieldType = None
    RealVectorType = None

    def __init__(self):
        return

    # Calculate the actual index into this container for the given python
    # style index.
    def boundIndex(self, index):
        if index < 0:
            index += len(self)
        return max(0, min(len(self), index))

    # Return safe values for the min, max, and stride from a slice
    def processSlice(self, obj):
        if obj.start:
            minIndex = self.boundIndex(obj.start)
        else:
            minIndex = 0
        if obj.stop:
            maxIndex = self.boundIndex(obj.stop)
        else:
            maxIndex = len(self)
        if obj.step:
            stride = obj.step
        else:
            stride = 1
        return minIndex, maxIndex, stride

    # Overload the getitem method and take care of the special case of slices.
    def __getitem__(self, arg):
        if type(arg) != type(slice(1)):
            return self.RealFieldType.__getitem__(self, arg)

        else:
            minIndex, maxIndex, stride = self.processSlice(arg)
            sliceSize = (maxIndex - minIndex) / stride
            result = self.RealVectorType(sliceSize)
            j = 0
            for i in xrange(minIndex, maxIndex, stride):
                result[j] = self[i]
            return result

    # Overload the setitem method to take care of slices.
    def __setitem__(self, arg, val):
        if type(arg) != type(slice(1)):
            self.RealFieldType.__setitem__(self, arg, val)

        else:
            minIndex, maxIndex, stride = self.processSlice(arg)
            sliceSize = (maxIndex - minIndex) / stride
            assert len(val) == sliceSize
            j = 0
            for i in xrange(minIndex, maxIndex, stride):
                self.RealFieldType[i] = val[j]

# ------------------------------------------------------------------------------
# Define the base class for extending field types.
# ------------------------------------------------------------------------------
import CXXWraps

from Field import ScalarField1d as RealScalarField1d
class ScalarField1d(ExtendFieldBase, RealScalarField1d):
    RealFieldType = RealScalarField1d
    RealVectorType = CXXWraps.vector_of_double
  
    def __init__(self, nodeList):
        print 'calling ExtendFieldBase.__init__'
        ExtendFieldBase.__init__(self)
        print 'calling RealScalarField1d.__init__'
        RealScalarField1d.__init__(self, nodeList)
        print 'done'
        return

##    def __init__(self, nodeList, field):
##        print 'calling ExtendFieldBase.__init__'
##        ExtendFieldBase.__init__(self)
##        print 'calling RealScalarField1d.__init__'
##        RealScalarField1d.__init__(self, nodeList, field)
##        print 'done'
##        return
