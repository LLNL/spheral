from Spheral import *
from SpheralTestUtilities import *

################################################################################
fileName = "test.h5"
access = 0           # Create

file = HDF5IO2d(fileName, access)
output('file')

################################################################################
title('2-D HDF5 I/O Test: write atomic data types')

groupName = "/Root/AtomicTypes/"

fileInt = 42
fileScalar = 42.5
fileVector = Vector2d(1.3, 4.92)
fileTensor = Tensor2d(0, 1, 2, 3)
fileSymTensor = SymTensor2d(0, 1, 1, 2)
fileString = "Booga booga"

output("fileInt")
file.writeInt(fileInt, groupName + "fileInt")
output("fileScalar")
file.writeScalar(fileScalar, groupName + "fileScalar")
output("fileVector")
file.writeVector(fileVector, groupName + "fileVector")
output("fileTensor")
file.writeTensor(fileTensor, groupName + "fileTensor")
output("fileSymTensor")
file.writeSymTensor(fileSymTensor, groupName + "fileSymTensor")
output("fileString")
file.writeString(fileString, groupName + "fileString")

################################################################################
title('2-D HDF5 I/O Test: write Field data types')

groupName = "/Root/FieldTypes/"

n = 50
nodes = SphNodeList2d(n)

scalarField = ScalarField2d(nodes)
vectorField = VectorField2d(nodes)
tensorField = TensorField2d(nodes)
symTensorField = SymTensorField2d(nodes)

for i in xrange(n):
    scalarField[i] = i
    vectorField[i] = (i, 2*i)
    tensorField[i] = Tensor2d(i, 2*i, 3*i, 4*i)
    symTensorField[i] = SymTensor2d(i, 2*i, 2*i, 3*i)

output("scalarField")
file.writeScalarField(scalarField, groupName + "scalarField")
output("vectorField")
file.writeVectorField(vectorField, groupName + "vectorField")
output("tensorField")
file.writeTensorField(tensorField, groupName + "tensorField")
output("symTensorField")
file.writeSymTensorField(symTensorField, groupName + "symTensorField")

################################################################################
title('2-D HDF5 I/O Test: read atomic data types')

# Close the current File
file = None

# Open it up again in read only mode.
file = HDF5IO2d(fileName, 1)

groupName = "/Root/AtomicTypes/"

readInt = 0
file.readInt(readInt, groupName + "fileInt")
output("fileInt, readInt")
##assert(fileInt == readInt)

##fileInt = 42
##fileScalar = 42.5
##fileVector = Vector2d(1.3, 4.92)
##fileTensor = Tensor2d(0, 1, 2, 3)
##fileSymTensor = SymTensor2d(0, 1, 1, 2)
##fileString = "Booga booga"

##output("fileInt")
##file.writeInt(fileInt, groupName + "fileInt")
##output("fileScalar")
##file.writeScalar(fileScalar, groupName + "fileScalar")
##output("fileVector")
##file.writeVector(fileVector, groupName + "fileVector")
##output("fileTensor")
##file.writeTensor(fileTensor, groupName + "fileTensor")
##output("fileSymTensor")
##file.writeSymTensor(fileSymTensor, groupName + "fileSymTensor")
##output("fileString")
##file.writeString(fileString, groupName + "fileString")

################################################################################
